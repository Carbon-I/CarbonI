
using Unitful

Base.@kwdef mutable struct Scenario
    ### Constants; these are design characteristics
    # not expected to change
    lat::Float64
    lon::Float64

    ### Dependent on build outcome; each has a 
    # required value and a current best estimate (CBE)
end


function reference_scenario()
    ref = Scenario(
        lat = 0.0,
        lon = -62.0
    )
    return ref
end

function cbe_instrument()
    inst = Scenario(
        FPA_quantum_efficiency = 0.85,
        bench_efficiency = 0.72,
        readout_noise = 100.0,
        dark_current = 5e3u"1/s"
    )
    return inst
end
