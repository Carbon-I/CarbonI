using CarbonI
using Unitful
using DiffResults
using DelimitedFiles
using ForwardDiff
using InstrumentOperator
using Interpolations
using Unitful, PrettyTables, Printf, Markdown

# Map Symbol -> Val for PrettyTables
_pt_backend(sym::Symbol) = Val(sym)

function show_required_precisions(σ_ch4_req, σ_co2_req, σ_co_req;
                                  uch4=u"ppb", uco2=u"ppm", uco=u"ppb",
                                  sigdigits::Int=3,
                                  render::Symbol = :markdown)

    ch4 = uconvert(uch4, σ_ch4_req)
    co2 = uconvert(uco2, σ_co2_req)
    co  = uconvert(uco,  σ_co_req)

    gases = ["CH₄", "CO₂", "CO"]
    vals  = [
        @sprintf("%.*g", sigdigits, ustrip(ch4)),
        @sprintf("%.*g", sigdigits, ustrip(co2)),
        @sprintf("%.*g", sigdigits, ustrip(co))
    ]
    units = [ string(unit(ch4)), string(unit(co2)), string(unit(co)) ]

    data = hcat(gases, vals, units)

    if render === :markdown
        io = IOBuffer()
        pretty_table(io, data;
                     header   = ["Gas", "Required single-measurement precision", "Unit"],
                     alignment = :l,
                     backend  = _pt_backend(:markdown))
        md = String(take!(io))
        display("text/markdown", md)   # <-- forces Markdown rendering in VS Code/Jupyter
    elseif render === :html
        io = IOBuffer()
        pretty_table(io, data;
                     header   = ["Gas", "Required single-measurement precision", "Unit"],
                     alignment = :l,
                     backend  = _pt_backend(:html))
        html = String(take!(io))
        display("text/html", html)     # <-- forces HTML rendering
    else
        # plain text fallback
        pretty_table(data;
                     header   = ["Gas", "Required single-measurement precision", "Unit"],
                     alignment = :l,
                     backend  = _pt_backend(:text))
    end
end

function md_expected_precisions(error_cbe, sqrtN; sigdigits=3)
    ch4 = error_cbe["ch4"] / sqrtN                 # ppb
    co2 = error_cbe["co2"] / 1000 / sqrtN          # ppm
    co  = error_cbe["co"]  / sqrtN                 # ppb

    row(gas, v, unit) = "| $(gas) | " * @sprintf("%.*g", sigdigits, v) * " | $(unit) |"

    md = """
    | Gas | Expected single-measurement precision | Unit |
    |-----|--------------------------------------:|------|
    $(row("CH₄", ch4, "ppb"))
    $(row("CO₂", co2, "ppm"))
    $(row("CO",  co,  "ppb"))
    """

    display("text/markdown", md)
end

using Printf, Markdown

"""
    req_vs_exp_table(rows; sigdigits=3, render=:markdown)

`rows` is a Vector of NamedTuples like:
    (gas="CH₄", unit="ppb", req=..., exp=...)

`render` can be :markdown, :html, or :text.
"""
function req_vs_exp_table(rows::Vector{<:NamedTuple};
                          sigdigits::Int = 3,
                          render::Symbol = :markdown)

    fmt(v) = @sprintf("%.*g", sigdigits, v)
    margin(req, exp) = 100 * (exp - req) / exp
    margin2(req, exp) = 100 * (req - exp) / req

    header = """
    | Gas | Required σ | Expected σ | Margin wrt to CBE | Margin wrt to REQ |  Unit |
    |-----|-----------:|-----------:|-------:|-------:|------|
    """

    body = IOBuffer()
    for r in rows
        m = margin(r.req, r.exp)
        m2 = margin2(r.req, r.exp)
        @printf(body, "| %s | %s | %s | %+0.1f%% | %+0.1f%% |%s |\n",
                r.gas, fmt(r.req), fmt(r.exp), m, m2, r.unit)
    end
    table = header * String(take!(body))

    if render === :markdown
        display("text/markdown", table)
    elseif render === :html
        # very simple HTML wrapper (optional)
        html = replace(table, "\n" => "<br>")
        display("text/html", html)
    else
        print(table)
    end
end

include(joinpath(dirname(pathof(CarbonI)), "forwardModel.jl"))

function calc_rel_error(specs, x, solarIrr, refl, sza, σ_matrix, profile, h, ins, Sₐ; return_F=false)

	# Run Forward model and Jacobian generation together:
	result = DiffResults.JacobianResult(zeros(length(specs.instrument_wl)),x);

	ForwardDiff.jacobian!(result, (x) -> 
	forward_model_x_(x, 
				sun=solarIrr, 
				instrument=specs.instrument_kernel,
				reflectance=refl, 
				sza=sza, 
				σ_matrix=σ_matrix, 
				profile=profile,
				wl=specs.modelling_wl), x);
	K = DiffResults.jacobian(result);
	F = DiffResults.value(result);  

	# Compute noise equivalent radiance:
	nesr = InstrumentOperator.noise_equivalent_radiance(ins,(specs.instrument_wl)u"nm", (F)u"mW/m^2/nm/sr");
	nesr_ = nesr./1u"mW/m^2/nm/sr"

	# Generate S\_epsilon matrix
	Se = Diagonal(nesr_.^2);

	# Compute the Gain Matrix:
	#G = inv(K'inv(Se)K + inv(Sₐ))K'inv(Se);

	# Posterior covariance matrix (Profile retrievals, 10 layers):
	Ŝ = inv(K'inv(Se)K + inv(Sₐ));

	# Apply column operator to get to total column uncertainties:
	error_ppb = Dict()
	for spec in ["co2","co213", "ch4", "h2o", "co", "hdo", "n2o", "c2h6"]
		error_ppb[spec] = sqrt(h[spec]' * Ŝ * h[spec])*1e9
	end
    if return_F
        return error_ppb, F
    else
        return error_ppb
    end
	
	
end


function setup_data(scenario, specs)

	# Load spectroscopies:
	co2, ch4, h2o, hdo, n2o, co, co2_iso2, c2h6 = CarbonI.loadXSModels();


	DS = Dataset(CarbonI.solar_file);
	wlSol = 1e3*DS["wl"][:]
	solar_irr = 1e3*DS["solar_irr"][:] # convert to mW/m2/nm
	close(DS)

	# optionally load a custom reflectance
	#DS = Dataset("data/reflectance_cube_all_1nm_from450nm.h5")
	#r = DS["reflectance_cube"][:]
	#close(DS)
	
	hitran_array = (co2, h2o, ch4, co, n2o, hdo, co2_iso2, c2h6);
	# Precompute the cross sections:
	σ_matrix_hr = CarbonI.compute_profile_crossSections(scenario.profile_hr, hitran_array , specs.modelling_wl);

	nL = length(scenario.profile_hr.T)

	vmr_co2 = zeros(nL) .+ 407e-6
	vmr_ch4 = zeros(nL) .+ 1.8e-6
	vmr_ch4[1:3] .= 1.4e-6
	vmr_h2o = scenario.profile_hr.vcd_h2o ./ scenario.profile_hr.vcd_dry
	vmr_co  = zeros(nL) .+ 100e-9
	vmr_n2o = zeros(nL) .+ 337e-9
	vmr_n2o[1:3] .= 100e-9
	vmr_hdo = vmr_h2o * 0.9
	vmr_c2h6 = zeros(nL) .+ 1.0e-9
	vmrs = [vmr_co2, vmr_h2o, vmr_ch4,vmr_co, vmr_n2o, vmr_hdo, vmr_co2, vmr_c2h6];

	sol  = CubicSplineInterpolation(range(wlSol[1],wlSol[end], length=length(wlSol)),solar_irr, extrapolation_bc=Interpolations.Flat());

	# Reduce dimensions, group layers together to get  layers of equal pressure difference:
	n_layers = 10;

	# Reduce profile and cross sections to fewer dimensions:
	profile, σ_matrix, indis, gasProfiles = CarbonI.reduce_profile(n_layers, scenario.profile_hr, σ_matrix_hr,vmrs);

	# Number of Legendre Polynomials for the spectrally resolved albedo term
	nLeg = 10

	# Define state vector elements for legendre polynomials
	xPoly = zeros(nLeg).+eps()

	# First element has to be unity here
	xPoly[1] = 1.0

	# Define full state vector (includes trace gas VMRs and Legendre Polynomial)
	x = [reduce(vcat,gasProfiles) ; xPoly ];

	# Define Covariance matrices and column kernel operator
	# Get prior covariance matrix:
	n_state = length(x);
	Sₐ = zeros(n_state,n_state);
	rel_error = 0.0001;
	rel_error = 0.02;
	# could try modifying this to something like 0.02 (TBD)
	# currently good for point surfaces and surface components

	# vcd_ratio = profile_caltech.vcd_dry ./ mean(profile_caltech.vcd_dry)

	# Fill the diagonal for the trace gases:
	for i=1:80
		Sₐ[i,i] = (rel_error*x[i])^2   
	end
	# CO2 at surface, 100% error
	Sₐ[10,10] = (200x[10])^2
	Sₐ[20,20] = (200x[20])^2
	Sₐ[30,30] = (200x[30])^2
	Sₐ[40,40] = (200x[40])^2
	Sₐ[50,50] = (200x[50])^2
	Sₐ[60,60] = (200x[60])^2
	Sₐ[70,70] = (200x[70])^2
	Sₐ[80,80] = (20000x[80])^2
	# Put in arbitrarily high numbers for the polynomial term, so these won't be constrained at all! 
	for i=81:n_state
		Sₐ[i,i] = 1e2;
	end
	ratio = profile.vcd_dry/sum(profile.vcd_dry);

	h = Dict()
	for spec in ["co2","co213", "ch4", "h2o", "co", "hdo", "n2o", "c2h6"]
		h[spec] = zeros(length(x));
	end

	h["co2"][1:10] .= ratio;
	h["h2o"][11:20] .= ratio;
	h["ch4"][21:30] .= ratio;
	h["co"][31:40] .= ratio;
	h["n2o"][41:50] .= ratio;
	h["hdo"][51:60] .= ratio;
	h["co213"][61:70] .= ratio;
	h["c2h6"][71:80] .= ratio;

	# Load tropical albedo
	clima_alb = readdlm(CarbonI.albedo_file,',', skipstart=1);

	# Generate albedo function generator
	soil = CubicSplineInterpolation(300:2400,clima_alb[:,2]/1.16, extrapolation_bc=Interpolations.Flat());

	# Create solar spectrum at Forward model resolution
	solarIrr = sol(specs.modelling_wl);


	return soil, x, solarIrr, σ_matrix, profile, h, Sₐ

end

