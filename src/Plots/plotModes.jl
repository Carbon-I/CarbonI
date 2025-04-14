using GLMakie

# Define Earth's radius in kilometers
earth_radius = 6371.0

# Create a circle representing the Earth's curvature
θ = LinRange(0, 2π, 200)
earth_x = earth_radius .* cos.(θ)
earth_y = earth_radius .* sin.(θ)

# Define the target point on the Earth's surface (example: 30° from the positive x-axis)
target_angle = deg2rad(90)
target_x = earth_radius * cos(target_angle)
target_y = earth_radius * sin(target_angle)

# Define the satellite position above the Earth's surface (using 688 km altitude).
# Example: place the satellite along a 45° direction.
sat_angle = deg2rad(90)
sat_altitude = 688.0
sat_radius = earth_radius + sat_altitude
sat_x = sat_radius * cos(sat_angle)
sat_y = sat_radius * sin(sat_angle)

# The nadir point is directly below the satellite on the Earth's surface.
nadir_x = earth_radius * cos(sat_angle)
nadir_y = earth_radius * sin(sat_angle)

# Define a simple glint point: here, a slight offset from the target angle (offset by 5°)
glint_angle = deg2rad(95)
glint_x = earth_radius * cos(glint_angle)
glint_y = earth_radius * sin(glint_angle)

# Set up the figure and axis
fig = Figure(resolution = (800, 600))
ax = Axis(fig[1, 1],
    xlabel = "X (km)",
    ylabel = "Y (km)",
    title = "Satellite Geometry with Target, Glint & Nadir",
    aspect = DataAspect()
)

# Plot the Earth as a circle
lines!(ax, earth_x, earth_y, color = :black, linewidth=4, label = "Earth")

# Plot the target point on the Earth's surface
#scatter!(ax, [target_x], [target_y], color = :red, markersize = 12, label = "Target")

# Plot the glint point and draw dashed lines indicating the glint geometry
scatter!(ax, [glint_x], [glint_y], color = :orange, markersize = 12, label = "Glint")
#lines!(ax, [target_x, glint_x], [target_y, glint_y], color = :orange, linestyle = :dash)
lines!(ax, [sat_x, glint_x], [sat_y, glint_y], color = :orange,  linewidth=3)

# Plot the satellite position
scatter!(ax, [sat_x], [sat_y], color = :green, markersize = 12, label = "Satellite")

# Draw a dashed line from the satellite to its nadir on the Earth
lines!(ax, [sat_x, nadir_x], [sat_y, nadir_y], color = :green, linewidth=3, label = "Nadir")
xlims!(-1000, 1000)
ylims!(6150, 7300)
# Display legend and the final figure
axislegend(ax)
fig


using GLMakie

# Define Earth and satellite parameters (in kilometers)
earth_radius = 6371.0
sat_altitude = 688.0
sat_radius = earth_radius + sat_altitude

# Create a parametric mesh for the Earth (a sphere)
u = LinRange(0, 2π, 50)
v = LinRange(0, π, 25)
x = [earth_radius*sin(φ)*cos(θ) for φ in v, θ in u]
y = [earth_radius*sin(φ)*sin(θ) for φ in v, θ in u]
z = [earth_radius*cos(φ) for φ in v, θ in u]

# Set up the figure and a 3D axis
fig = Figure(resolution = (800, 600))
ax = LScene(fig[1, 1])
cc = Makie.Camera3D(scene.scene, projectiontype = Makie.Perspective)

# Plot the Earth as a semi-transparent blue surface
surface!(ax, x, y, z, color = :blue, transparency = true, alpha = 0.3)

# Define key points on the Earth's surface (assumed here on the equatorial plane, z = 0)
# Target: located at 30° from the x-axis
target_angle_deg = 30.0
target_angle = deg2rad(target_angle_deg)
target = Point3f(earth_radius*cos(target_angle), earth_radius*sin(target_angle), 0)

# Glint: an example offset (5° offset) from the target location
glint_angle_deg = target_angle_deg + 5.0
glint_angle = deg2rad(glint_angle_deg)
glint = Point3f(earth_radius*cos(glint_angle), earth_radius*sin(glint_angle), 0)

# Satellite: positioned at a 45° angle (in the equatorial plane) at 688 km altitude
sat_angle_deg = 45.0
sat_angle = deg2rad(sat_angle_deg)
satellite = Point3f(sat_radius*cos(sat_angle), sat_radius*sin(sat_angle), 0)

# Nadir: the point on the Earth's surface directly below the satellite
nadir = Point3f(earth_radius*cos(sat_angle), earth_radius*sin(sat_angle), 0)

# Markers for the points
scatter!(ax, [target], markersize = 15, color = :red, label = "Target")
scatter!(ax, [glint], markersize = 15, color = :orange, label = "Glint")
scatter!(ax, [satellite], markersize = 15, color = :green, label = "Satellite")
scatter!(ax, [nadir], markersize = 15, color = :purple, label = "Nadir")

# Draw lines illustrating the viewing geometry
# Line from target to glint and from satellite to glint (glint geometry)
lines!(ax, [target, glint], linestyle = :dash, color = :orange)
lines!(ax, [satellite, glint], linestyle = :dash, color = :orange)
# Line from satellite to its nadir (nadir viewing geometry)
lines!(ax, [satellite, nadir], linestyle = :dash, color = :green)

# Set the camera view to a slanted angle for an intuitive 3D perspective.
# The vector below sets an arbitrary camera position that you can adjust as needed.
center!(scene.scene)

# Display the figure
fig