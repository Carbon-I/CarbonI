function create_carbonI_conv_matrix(wl::StepRangeLen{FT}) where FT
    # Define a Fixed instrument:
    FWHM  = 2.2  # 
    SSI  = 0.7
    Δwl = wl.step.hi
    kern1 = CarbonI.box_kernel(2*SSI, Δwl)
    kern2 = CarbonI.gaussian_kernel(FWHM, Δwl)
    kernf = imfilter(kern1, kern2)
    
    # Hardcoded for Carbon-I
    lociBox = CarbonI.KernelInstrument(kernf, collect(2035:SSI:2380));
    # Generate convolution matrix:
    cM = CarbonI.generate_conv_matrix(lociBox,wl, Δwl)
    return cM, lociBox.ν_out
end

function create_carbonI_conv_matrix_cbe(wl::StepRangeLen{FT}) where FT
    # Define a Fixed instrument:
    FWHM  = 1.5  # 
    SSI  = 0.7
    Δwl = wl.step.hi
    kern1 = CarbonI.box_kernel(2*SSI, Δwl)
    kern2 = CarbonI.gaussian_kernel(FWHM, Δwl)
    kernf = imfilter(kern1, kern2)
    
    # Hardcoded for Carbon-I
    lociBox = CarbonI.KernelInstrument(kernf, collect(2035:SSI:2380));
    # Generate convolution matrix:
    cM = CarbonI.generate_conv_matrix(lociBox,wl, Δwl)
    return cM, lociBox.ν_out
end

function create_emit_conv_matrix(wl::StepRangeLen{FT}) where FT
    Δwl = wl.step.hi
    wl_emit = collect(2035:7.5:2380)
    emitBox = CarbonI.KernelInstrument(CarbonI.box_kernel(8.5, Δwl), wl_emit);
    cM = CarbonI.generate_conv_matrix(emitBox,wl, Δwl)
    return cM, emitBox.ν_out
end
