using NCDatasets, CairoMakie, Colors

f = Dataset("data/Sentinel2_GreekFire_30m.nc")

cc = to_color((:red,0.4))
greek300 = f["data"][:]
carbonI = greek300[1:3500,end:-1:1,1]
r   = greek300[1:3500,end:-1:1,2]
g  = greek300[1:3500,end:-1:1,3] 
b   = greek300[1:3500,end:-1:1,4]
o2 = greek300[1:3500,end:-1:1,5]
wco2 = greek300[1:3500,end:-1:1,6]

zoomX = 1900:2300
zoomY = 1450:1850
# Define the coordinates of the rectangle
x_rect = zoomX[1]  # x-coordinate of the lower-left corner
y_rect = zoomY[1] # y-coordinate of the lower-left corner
width_rect = zoomX[end]-zoomX[1]  # Width of the rectangle
height_rect = zoomY[end]-zoomY[1]  # Height of the rectangle
# Calculate the corners of the rectangle
rect_corners = [
    Point(x_rect, y_rect),  # Bottom-left
    Point(x_rect + width_rect, y_rect),  # Bottom-right
    Point(x_rect + width_rect, y_rect + height_rect),  # Top-right
    Point(x_rect, y_rect + height_rect),  # Top-left
    Point(x_rect, y_rect)  # Closing the loop back to bottom-left
]

R = replace!(r, NaN => 0.0); R[R .> 0.5] .= 0.5
G = replace!(g, NaN => 0.0); G[G .> 0.5] .= 0.5
B = replace!(b, NaN => 0.0); B[B .> 0.5] .= 0.5

fac = 2
RGB_img = RGB.(R*fac, G*fac, B*fac)

# Create the heatmap
fig = Figure(size=(1500, 900), fontsize=24, font = :bold)
ax1 = Axis(fig[1, 1], title="RGB image")
image!(ax1,RGB_img)
# Draw the outline of the rectangle using linesegments
poly!(ax1, rect_corners, color = (:red,0.2), linewidth = 2)
ax11 = Axis(fig[2, 1])
image!(ax11,RGB_img[zoomX,zoomY])

ax2 = Axis(fig[1, 2], title="O₂-A Band at 780nm")
heatmap!(ax2, o2, colorrange=(0.0,0.45),colormap=:grays)
poly!(ax2, rect_corners, color = (:red,0.2), linewidth = 2)
ax22 = Axis(fig[2, 2])
heatmap!(ax22, o2[zoomX,zoomY], colorrange=(0.0,0.45),colormap=:grays)

ax3 = Axis(fig[1, 3], title="Carbon-I band at 2.2µm")
c = to_color(:red)
ci = heatmap!(ax3, carbonI, colorrange=(0.0,0.45),colormap=:grays)#, highclip=cc)
poly!(ax3, rect_corners, color = (:red,0.2), linewidth = 2)
#Colorbar(fig[1, 4], ci)#, label="CH₄ error (ppb)")
ax33 = Axis(fig[2, 3])
c = to_color(:red)
ci2 = heatmap!(ax33, carbonI[zoomX, zoomY], colorrange=(0.0,0.45),colormap=:grays)#, highclip=cc)
Colorbar(fig[1:2, 4], ci)#, label="CH₄ error (ppb)")

for ax in (ax1,ax2, ax3, ax11, ax22, ax33)
    hidespines!(ax)
    hidexdecorations!(ax)
    hideydecorations!(ax)
end
fig
save("plots/GreekFire_CarbonI.png", fig)



f = Dataset("data/Sentinel2_Delhi_300m.nc")
f2 = Dataset("data/Sentinel2_Delhi_30m.nc")


greek300 = f2["data"][:]
carbonI = greek300[:,end:-1:1,1]
r   = greek300[:,end:-1:1,2]
g  = greek300[:,end:-1:1,3] 
b   = greek300[:,end:-1:1,4]
o2 = greek300[:,end:-1:1,5]
wco2 = greek300[:,end:-1:1,6]

zoomX = 2000:2400
zoomY = 2400:2800
# Define the coordinates of the rectangle
x_rect = zoomX[1]  # x-coordinate of the lower-left corner
y_rect = zoomY[1] # y-coordinate of the lower-left corner
width_rect = zoomX[end]-zoomX[1]  # Width of the rectangle
height_rect = zoomY[end]-zoomY[1]  # Height of the rectangle
# Calculate the corners of the rectangle
rect_corners = [
    Point(x_rect, y_rect),  # Bottom-left
    Point(x_rect + width_rect, y_rect),  # Bottom-right
    Point(x_rect + width_rect, y_rect + height_rect),  # Top-right
    Point(x_rect, y_rect + height_rect),  # Top-left
    Point(x_rect, y_rect)  # Closing the loop back to bottom-left
]

R = replace!(r, NaN => 0.0); R[R .> 0.5] .= 0.5
G = replace!(g, NaN => 0.0); G[G .> 0.5] .= 0.5
B = replace!(b, NaN => 0.0); B[B .> 0.5] .= 0.5

fac = 2
RGB_img = RGB.(R*fac, G*fac, B*fac)

# Create the heatmap
fig = Figure(size=(1900, 900), fontsize=24, font = :bold)
ax1 = Axis(fig[1, 1], title="RGB image")
image!(ax1,RGB_img)
# Draw the outline of the rectangle using linesegments
poly!(ax1, rect_corners, color = (:red,0.2), linewidth = 2)
ax11 = Axis(fig[2, 1])
image!(ax11,RGB_img[zoomX,zoomY])

ax2 = Axis(fig[1, 2], title="O₂-A Band at 780nm")
heatmap!(ax2, o2, colorrange=(0.0,0.45),colormap=:grays)
poly!(ax2, rect_corners, color = (:red,0.2), linewidth = 2)
ax22 = Axis(fig[2, 2])
heatmap!(ax22, o2[zoomX,zoomY], colorrange=(0.0,0.45),colormap=:grays)

ax3 = Axis(fig[1, 3], title="Carbon-I band at 2.2µm")
c = to_color(:red)
ci = heatmap!(ax3, carbonI, colorrange=(0.0,0.45),colormap=:grays)#, highclip=:red)
poly!(ax3, rect_corners, color = (:red,0.2), linewidth = 2)
#Colorbar(fig[1, 4], ci)#, label="CH₄ error (ppb)")
ax33 = Axis(fig[2, 3])
c = to_color(:red)
ci2 = heatmap!(ax33, carbonI[zoomX, zoomY], colorrange=(0.0,0.45),colormap=:grays)#, highclip=:red)
Colorbar(fig[1:2, 4], ci)#, label="CH₄ error (ppb)")

for ax in (ax1,ax2, ax3, ax11, ax22, ax33)
    hidespines!(ax)
    hidexdecorations!(ax)
    hideydecorations!(ax)
end
fig
save("plots/Delhi_CarbonI.png", fig)

using NCDatasets, CairoMakie, Colors

f = Dataset("data/Sentinel2_Amazon_30m.nc")

cc = to_color((:red,0.4))
greek300 = f["data"][:]
carbonI = greek300[1:3500,end:-1:1,1]
r   = greek300[:,end:-1:1,2]
g  = greek300[:,end:-1:1,3] 
b   = greek300[:,end:-1:1,4]
o2 = greek300[1:3500,end:-1:1,5]
wco2 = greek300[1:3500,end:-1:1,6]

zoomX = 1900:2300
zoomY = 1450:1850
# Define the coordinates of the rectangle
x_rect = zoomX[1]  # x-coordinate of the lower-left corner
y_rect = zoomY[1] # y-coordinate of the lower-left corner
width_rect = zoomX[end]-zoomX[1]  # Width of the rectangle
height_rect = zoomY[end]-zoomY[1]  # Height of the rectangle
# Calculate the corners of the rectangle
rect_corners = [
    Point(x_rect, y_rect),  # Bottom-left
    Point(x_rect + width_rect, y_rect),  # Bottom-right
    Point(x_rect + width_rect, y_rect + height_rect),  # Top-right
    Point(x_rect, y_rect + height_rect),  # Top-left
    Point(x_rect, y_rect)  # Closing the loop back to bottom-left
]

R = replace!(r, NaN => 0.0); R[R .> 0.5] .= 0.5
G = replace!(g, NaN => 0.0); G[G .> 0.5] .= 0.5
B = replace!(b, NaN => 0.0); B[B .> 0.5] .= 0.5

fac = 2
RGB_img = RGB.(R*fac, G*fac, B*fac)

# Create the heatmap
fig = Figure(size=(1500, 900), fontsize=24, font = :bold)
ax1 = Axis(fig[1, 1], title="RGB image")
image!(ax1,RGB_img)
# Draw the outline of the rectangle using linesegments
poly!(ax1, rect_corners, color = (:red,0.2), linewidth = 2)
ax11 = Axis(fig[2, 1])
image!(ax11,RGB_img[zoomX,zoomY])

ax2 = Axis(fig[1, 2], title="O₂-A Band at 780nm")
heatmap!(ax2, o2, colorrange=(0.0,0.45),colormap=:grays)
poly!(ax2, rect_corners, color = (:red,0.2), linewidth = 2)
ax22 = Axis(fig[2, 2])
heatmap!(ax22, o2[zoomX,zoomY], colorrange=(0.0,0.45),colormap=:grays)

ax3 = Axis(fig[1, 3], title="Carbon-I band at 2.2µm")
c = to_color(:red)
ci = heatmap!(ax3, carbonI, colorrange=(0.0,0.45),colormap=:grays)#, highclip=cc)
poly!(ax3, rect_corners, color = (:red,0.2), linewidth = 2)
#Colorbar(fig[1, 4], ci)#, label="CH₄ error (ppb)")
ax33 = Axis(fig[2, 3])
c = to_color(:red)
ci2 = heatmap!(ax33, carbonI[zoomX, zoomY], colorrange=(0.0,0.45),colormap=:grays)#, highclip=cc)
Colorbar(fig[1:2, 4], ci)#, label="CH₄ error (ppb)")

for ax in (ax1,ax2, ax3, ax11, ax22, ax33)
    hidespines!(ax)
    hidexdecorations!(ax)
    hideydecorations!(ax)
end
fig
save("plots/AmazonFire_CarbonI_RGB.png", fig)


function export_image(nc_file, name)
    f = Dataset(nc_file)
    greek300 = f["data"][:]
    nir = greek300[:,end:-1:1,1]
    r   = greek300[:,end:-1:1,2]
    g  = greek300[:,end:-1:1,3] 
    b   = greek300[:,end:-1:1,4]
    R = replace!(r, NaN => 0.0); #R[R .> 0.5] .= 0.5
    G = replace!(g, NaN => 0.0); #G[G .> 0.5] .= 0.5
    B = replace!(b, NaN => 0.0); #B[B .> 0.5] .= 0.5
    nir = replace!(nir, NaN => 0.0); nir[nir .> 1.0] .= 1.0
    RGB_img = RGB.(R, G, B)
    NIR_img = RGB.(nir, nir, nir)
    # Create a figure with the resolution matching the image size
    name_rgb = name * "_RGB.png"
    name_nir = name * "_NIR.png"

    fig = Figure(resolution = (size(RGB_img, 1), size(RGB_img, 2)))  # (width, height)

    # Add an axis without any borders, ticks, or labels
    ax = Axis(fig[1, 1])
    image!(ax, RGB_img)
    hidespines!(ax)
    hidexdecorations!(ax)
    hideydecorations!(ax)

    # Save the figure as PNG, with pixel-per-pixel resolution
    save(name_rgb, fig)

    fig = Figure(resolution = (size(NIR_img, 1), size(NIR_img, 2)))  # (width, height)

    # Add an axis without any borders, ticks, or labels
    ax = Axis(fig[1, 1])
    image!(ax, NIR_img)
    hidespines!(ax)
    hidexdecorations!(ax)
    hideydecorations!(ax)

    # Save the figure as PNG, with pixel-per-pixel resolution
    save(name_nir, fig)
end

