################################################################################
# Spherical Earth constants (in kilometers):
################################################################################
const R_EARTH = 6378.137  # Approx. Earth radius

################################################################################
# Helper rotation matrices about x- and y- axes (intrinsic rotations):
################################################################################
using LinearAlgebra
using StaticArrays

""" Rotation matrix about the x-axis by angle θ. """
function rotX(θ)
    c = cos(θ)
    s = sin(θ)
    return [
        1.0   0.0   0.0
        0.0    c    -s
        0.0    s     c
    ]
end

""" Rotation matrix about the y-axis by angle θ. """
function rotY(θ)
    c = cos(θ)
    s = sin(θ)
    return [
         c   0.0    s
        0.0  1.0   0.0
        -s   0.0    c
    ]
end

################################################################################
# Convert from geodetic latitude/longitude (spherical approximation) to ECEF.
#   lat, lon in radians
#   radius = R_EARTH + satellite altitude (km) if for a satellite.
################################################################################
function latlon_to_ecef(lat::Float64, lon::Float64, radius::Float64)
    x = radius * cos(lat) * cos(lon)
    y = radius * cos(lat) * sin(lon)
    z = radius * sin(lat)
    return SVector{3}(x, y, z)
end

################################################################################
# Convert ECEF to geodetic latitude/longitude on a spherical Earth.
#   returns (lat, lon) in radians.
################################################################################
function ecef_to_latlon(pos::SVector{3,Float64})
    x, y, z = pos
    # Spherical Earth:
    r_xy = sqrt(x^2 + y^2)  # projection in equatorial plane
    lat  = atan(z, r_xy)
    lon  = atan(y, x)
    return (lat, lon)
end

################################################################################
# Build local East, North, Up (ENU) axes at the sub-satellite point (lat0, lon0).
#   lat0, lon0 in radians
#   Returns (e, n, u) as normalized 3D vectors in ECEF frame.
################################################################################
function build_enu(lat0::Float64, lon0::Float64)
    # Up vector (local zenith)
    u = SVector{3}( cos(lat0)*cos(lon0),
                    cos(lat0)*sin(lon0),
                    sin(lat0) )

    # East vector
    e = SVector{3}(-sin(lon0), cos(lon0), 0.0)

    # North vector
    n = cross(u, e)  # or cross(e, u), depending on sign convention
    return (e, n, u)
end

################################################################################
# Main function: given:
#   lat0, lon0: sub-satellite lat/lon in radians
#   h: satellite altitude above Earth (km)
#   pitch, roll: angles in radians
#
# Returns:
#   dLat, dLon - the difference (in degrees) from (lat0, lon0) to the target
#                intersection point on Earth’s surface
################################################################################
function relative_latlon_offset(
        lat0::Float64,
        lon0::Float64,
        h::Float64,
        pitch::Float64,
        roll::Float64;
        R::Float64 = R_EARTH
    )

    # 1) Satellite position in ECEF:
    r_satellite = R + h
    sat_ecef = latlon_to_ecef(lat0, lon0, r_satellite)

    # 2) Local ENU axes at sub-satellite point:
    ê, n̂, û = build_enu(lat0, lon0)

    # 3) Define the nadir direction in local coordinates.
    #    In our local ENU, the “down” (nadir) direction is -û.
    #    We'll store it as a 3×1 vector: v_local = [0, 0, -1] is simpler
    #    but we need to align with the actual ENU basis if applying rotation
    #    in that local coordinate frame.
    # 
    #    We'll do this:
    #      - start with the local-vector [0,0,-1],
    #      - apply roll about x, then pitch about y,
    #      - then map it into ECEF.
    #
    #    So the local tangent coordinate system is X = ê, Y = n̂, Z = û.
    #    v_local = (0, 0, -1).
    #
    v_local_nadir = SVector{3}(0.0, 0.0, -1.0)

    # 4) Build the total rotation in local coordinates:
    #    roll about x => Rx(roll)
    #    pitch about y => Ry(pitch)
    #    final local vector = Ry(pitch)*Rx(roll)*v_local_nadir
    #
    #    (You might swap the order of rotations if your pitch/roll definitions
    #     differ from the convention used here.)
    #
    v_local = rotY(pitch) * (rotX(roll) * v_local_nadir)

    # 5) Convert this local vector into ECEF by linear combination of (ê, n̂, û):
    #    v_ecef = v_local[1]*ê + v_local[2]*n̂ + v_local[3]*û
    #
    v_ecef = v_local[1]*ê + v_local[2]*n̂ + v_local[3]*û

    # 6) Solve for intersection of ray with Earth (sphere of radius R):
    #    We have parametric line:
    #      X(t) = sat_ecef + t * v_ecef
    #    We want norm(X(t)) = R
    #
    #    => ||sat_ecef + t v_ecef||^2 = R^2
    #    => let S = sat_ecef, V = v_ecef:
    #       (S + tV)·(S + tV) = R^2
    #       S·S + 2t(S·V) + t^2 (V·V) = R^2
    #    Solve for t and pick the smaller positive root (the intersection “below” the satellite).
    #
    S = sat_ecef
    V = v_ecef
    S_dot_S = dot(S,S)
    S_dot_V = dot(S,V)
    V_dot_V = dot(V,V)

    # Quadratic: a*t^2 + b*t + c = 0
    a = V_dot_V
    b = 2.0 * S_dot_V
    c = S_dot_S - R^2

    # Solve:
    discriminant = b^2 - 4*a*c
    if discriminant < 0
        error("No real intersection found (discriminant < 0). Check geometry.")
    end

    t1 = (-b + sqrt(discriminant)) / (2a)
    t2 = (-b - sqrt(discriminant)) / (2a)

    # We want the intersection *in front* of the satellite along v_ecef
    # Typically t will be negative if the direction is 'down' toward Earth,
    # but it depends how you define v_ecef.  Check which root is physically correct.
    #
    # Generally, one root will be between 0 and 1*(some range),
    # or might be negative. Evaluate both and pick whichever is valid and is closer to the satellite.
    #
    # Because v_local_nadir = (0,0,-1), the direction is indeed “down.”
    # But let's just pick the smaller positive solution or the larger negative solution
    # depending on sign conventions. In many setups, the correct intersection is
    # the *larger* negative t if v is downward. You can test & pick the one that
    # yields a smaller ||X - S|| in magnitude. We'll do a simple check here:
    #
    candidates = [t1, t2]
    # pick the one that yields the smaller distance from sat_ecef to the Earth:
    function distance_to_sat(t) 
        return norm(S + t*V - S) 
    end
    t_candidates = filter(t -> !isnan(t), candidates)
    if isempty(t_candidates)
        error("No valid intersection parameter found.")
    end
    t_best = t_candidates[argmin(distance_to_sat.(t_candidates))]

    # Intersection point on Earth:
    X_int = S + t_best * V

    # 7) Convert intersection ECEF -> lat/lon (in radians):
    (lat_int, lon_int) = ecef_to_latlon(X_int)

    # 8) Return the difference (in degrees) from the sub-satellite point
    dLat = (lat_int - lat0) * (180/pi)
    dLon = (lon_int - lon0) * (180/pi)

    return dLat, dLon
end


################################################################################
# Example usage
################################################################################

# Suppose the sub-satellite point is at latitude=10°, longitude=45°,
# altitude = 700 km, with pitch=2°, roll=-1° (all angles in degrees).
# Let's compute the relative lat/lon shift (in degrees).
#


function example()
    lat0_deg  = 10.0
    lon0_deg  = 45.0
    alt_km    = 668.0
    pitch_deg = 0.0
    roll_deg  = 30.0

    # Convert input angles to radians for internal trigonometry:
    deg2rad = π/180
    lat0  = lat0_deg * deg2rad
    lon0  = lon0_deg * deg2rad
    pitch = pitch_deg * deg2rad
    roll  = roll_deg * deg2rad

    ifov = 51.7/1000000

    dLat, dLon = relative_latlon_offset(lat0, lon0, alt_km, pitch, roll)
    dLatIFOV, dLonIFOV = relative_latlon_offset(lat0, lon0, alt_km, pitch+ifov, roll+ifov)
    println("Relative offset in latitude = $dLat degrees")
    println("Relative offset in longitude = $dLon degrees")

    ΔLat = (dLat - dLatIFOV)*110e3
    ΔLon = (dLon - dLonIFOV)*110e3
    println("Relative footprint size in latitude = $ΔLat m")
    println("Relative footprint size in longitude = $ΔLon m")
end

example()