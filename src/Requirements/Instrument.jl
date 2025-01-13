
using Unitful


Base.@kwdef mutable struct InstrumentSpecs

    ### Constants; these are design characteristics
    # not expected to change
    ET            # Exposure time
	SSI           # Spectral resolution
	Pitch         # Pixel pitch
	Fnumber       # F-number

    ### Dependent on build outcome; each has a 
    # required value and a current best estimate (CBE)
    FPA_quantum_efficiency
    bench_efficiency
    readout_noise
    dark_current

end


# Provide an optional constuctor with the constant 
# values locked, for convenience
Base.@kwdef mutable struct InstrumentSpecs
    ET = 44.0u"ms",
    SSI = (2*0.7)u"nm",
    Pitch = 18.0u"μm",
    Fnumber = 2.2,
    FPA_quantum_efficiency = a,
    bench_efficiency = b,
    readout_noise = c,
    dark_current = d
end

function requirement_specs()
    inst = InstrumentSpecs(
        FPA_quantum_efficiency = 0.85,
        bench_efficiency = 0.72,
        readout_noise = 70.0,
        dark_current = 100.0u"1/s"
    )
    return inst
end

function cbe_specs()
    inst = InstrumentSpecs
    inst.bench_efficiency = 0.72
    inst.FPA_quantum_efficiency = 0.85
    return inst
end