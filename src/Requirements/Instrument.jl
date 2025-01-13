
using Unitful

Base.@kwdef mutable struct InstrumentSpecs
    ### Constants; these are design characteristics
    # not expected to change
    ET::Float64 = 44.0u"ms"
    SSI::Float64 = 1.4u"nm"
    Pitch::Float64 = 18.0u"um"
    Fnumber::Float64 = 2.2

    ### Dependent on build outcome; each has a 
    # required value and a current best estimate (CBE)
    FPA_quantum_efficiency::Float64 
    bench_efficiency::Float64
    readout_noise::Float64
    dark_current
end


function requirement_instrument()
    inst = InstrumentSpecs(
        FPA_quantum_efficiency = 0.80,
        bench_efficiency = 0.72,
        readout_noise = 120,
        dark_current = 10e3u"1/s"
    )
    return inst
end

function cbe_instrument()
    inst = InstrumentSpecs(
        FPA_quantum_efficiency = 0.85,
        bench_efficiency = 0.72,
        readout_noise = 100.0,
        dark_current = 5e3u"1/s"
    )
    return inst
end
