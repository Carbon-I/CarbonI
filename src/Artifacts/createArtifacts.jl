using Pkg.Artifacts

deploy = true

local_download_url = "/home/cfranken/code/gitHub/CarbonI/"
target_url         = "http://web.gps.caltech.edu/~cfranken/Carbon-I/"
filename           = "MERRA2_300.tavg3_3d_asm_Nv.20100610.nc4"
name               = "merra"

# Create Artifact
hash = create_artifact() do artifact_dir
            # Copy in weights
    cp(joinpath(local_download_url, filename), joinpath(artifact_dir, filename))
end

# Spit tarballs to be hosted out to local temporary directory:
if deploy
    tarball_hash = archive_artifact(hash, joinpath(local_download_url, "$(name).tar.gz"))

    # Calculate tarball url
    tarball_url = target_url*name*".tar.gz"
    # Bind this to an Artifacts.toml file
    @info("Binding $(name) in Artifacts.toml...")
    bind_artifact!(joinpath(@__DIR__, "Artifacts.toml"), name, hash; download_info=[(tarball_url, tarball_hash)], lazy=true, force=true)
end

# Copy file into www directory manually afterwards!
# cp merra.tar.gz ~/www/Carbon-I/

# Cross section database:
names = ["co2_model.jld2", "sco2_v52.jld2", "ch4_model.jld2", "h2o_model.jld2", "hdo_model.jld2", "n2o_model.jld2", "co_model.jld2", "co2_model_iso2.jld2", "c2h6_model.jld2"]

# This is the path to the Artifacts.toml we will manipulate
artifact_toml = joinpath(@__DIR__, "Artifacts.toml")
local_download_url = "/home/cfranken/code/gitHub/CarbonI/data_xs"
# Query the `Artifacts.toml` file for the hash bound to the name "xs"
# (returns `nothing` if no such binding exists)
# Create Artifact
xs_hash = create_artifact() do artifact_dir
    # Copy in weights
    for filename in names
        cp(joinpath(local_download_url, filename), joinpath(artifact_dir, filename))
    end
end


# Spit tarballs to be hosted out to local temporary directory:
if deploy
    tarball_hash = archive_artifact(xs_hash, joinpath(local_download_url, "xs.tar.gz"))

    # Calculate tarball url
    tarball_url = target_url*"xs.tar.gz"
    # Bind this to an Artifacts.toml file
    @info("Binding xs in Artifacts.toml...")
    bind_artifact!(joinpath(@__DIR__, "Artifacts.toml"), name, xs_hash; download_info=[(tarball_url, tarball_hash)], lazy=true, force=true)
end