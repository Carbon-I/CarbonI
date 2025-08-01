# This Cell is used to set up the environment for the Carbon-I Open Source project.
using Pkg
Pkg.activate("../.."); # Actives the environment at the root of the project
using CarbonI, vSmartMOM

# Add packages for plotting, etc.
using Printf, CairoMakie, DelimitedFiles, Statistics, Interpolations, NCDatasets, InstrumentOperator, LinearAlgebra, Unitful, UnitfulEquivalences, LaTeXStrings

include(joinpath("../../src/Plots", "CI_colorsNew.jl"));
include(joinpath(dirname(pathof(CarbonI)), "readSun.jl"))
include(joinpath(dirname(pathof(CarbonI)), "Requirements", "common.jl"))
set_theme!(theme_ggplot2())

using PrettyTables, Printf

"""
    show_pixels_table(N_max_req, N_eff_req, N_max_cbe, N_eff_cbe;
                      grid_label="12 km box", render=:markdown, pct_digits=1)

Show a table with Scenario, Max pixels, Cloud-free %, and Effective pixels
for Required vs CBE. `render` = :markdown (default), :html, or :text.
"""
function show_pixels_table(N_max_req, N_eff_req, N_max_cbe, N_eff_cbe;
                           grid_label::AbstractString="12 km box",
                           render::Symbol=:markdown, pct_digits::Int=1)

    req_max = round(Int, N_max_req);  req_eff = round(Int, N_eff_req)
    cbe_max = round(Int, N_max_cbe);  cbe_eff = round(Int, N_eff_cbe)

    f_req = req_max == 0 ? 0.0 : req_eff / req_max * 100
    f_cbe = cbe_max == 0 ? 0.0 : cbe_eff / cbe_max * 100
    f_req_str = @sprintf("%.*f%%", pct_digits, f_req)
    f_cbe_str = @sprintf("%.*f%%", pct_digits, f_cbe)

    # --- make a 2×4 table (not a vector of tuples) ---
    data = Any[
        "Required"  req_max  f_req_str  req_eff;
        "CBE"       cbe_max  f_cbe_str  cbe_eff
    ]

    header = ["Scenario", "Max pixels in $grid_label", "Cloud-free %", "Effective pixels"]
    aligns = [:l, :r, :r, :r]

    if render === :markdown || render === :html
        io = IOBuffer()
        backend = render === :markdown ? Val(:markdown) : Val(:html)
        pretty_table(io, data; header=header, alignment=aligns, backend=backend)
        display(render === :markdown ? "text/markdown" : "text/html", String(take!(io)))
    else
        pretty_table(data; header=header, alignment=aligns, backend=Val(:text))
    end
end

# Example call (using your variables):
# show_pixels_table(N_max_global_req, N_eff_global_req, N_max_global_cbe, N_eff_global_cbe;
#                   grid_label="12 km box", render=:markdown)

# using PrettyTables
"""
    show_max_pixels_table(N_max_req, N_max_cbe; grid_label="12 km box", render=:markdown)

Pretty-print a table of the maximum number of ground pixels for the Required vs CBE scenarios.

- `grid_label`: shown in the column header (e.g., "12 km box")
- `render`: :markdown (default), :html, or :text
"""
function show_max_pixels_table(N_max_req, N_max_cbe; grid_label="12 km box", render::Symbol=:markdown)
    # Ensure integer counts
    vals = [round(Int, N_max_req), round(Int, N_max_cbe)]

    # Build a 2D array so header length matches column count
    data = hcat(["Required", "CBE"], string.(vals))

    header = ["Scenario", "Max pixels in $(grid_label)"]

    if render === :markdown
        io = IOBuffer()
        pretty_table(io, data; header=header, alignment=:l, backend=Val(:markdown))
        display("text/markdown", String(take!(io)))
    elseif render === :html
        io = IOBuffer()
        pretty_table(io, data; header=header, alignment=:l, backend=Val(:html))
        display("text/html", String(take!(io)))
    else
        pretty_table(data; header=header, alignment=:l, backend=Val(:text))
    end
end

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