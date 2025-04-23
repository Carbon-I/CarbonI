using CarbonI

scenario = CarbonI.stressing_scenario()
grid = 1980:10:2400
spl = CubicSplineInterpolation(grid, scenario.surface_albedo(grid)*alb_scaling)

cM, wl_ci = CarbonI.create_carbonI_conv_matrix_cbe(wl)
ind = 1:70
wll = wl_ci[ind]
wll = 2030:0.1:2380

l = length(wll);
index = range(-1,1,l)

# Get Spectrum to fit:

#spec = spl(wll)
spec = scenario.surface_albedo(wll)

# Define Polynomial Terms
deg_range = 1:60
mods = zeros(l, length(deg_range));
fitSTD = zeros(length(deg_range));


for (i,deg) in enumerate(deg_range)
    K = zeros(l,deg)
    x = zeros(deg)

    for i=1:deg
        x.=0;
        x[i] = 1;
        pp = Legendre(x)
        K[:,i] = pp.(index)
    end

    # Fit
    a = K\spec;

    # recosntruct
    modeled = K*a
    mods[:,i] = modeled
    fitSTD[i] = sqrt(sum((spec-modeled).^2)/l)
end

Plots.plot(wll,spec./mods[:,2], label="degree 2")
Plots.plot!(wll,spec./mods[:,3], label="degree 3")
Plots.plot!(wll,spec./mods[:,5], label="degree 5")
Plots.plot!(wll,spec./mods[:,10], label="degree 10")
Plots.plot!(wll,spec./mods[:,15], label="degree 15")
Plots.plot!(wll,spec./mods[:,20], label="degree 20")
