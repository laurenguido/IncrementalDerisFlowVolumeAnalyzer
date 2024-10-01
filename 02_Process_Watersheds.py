# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 12:48:02 2024

@author: lguido

The code in this file (Process_Watersheds.py) reads a DEM raster, processes the DEM to 
condition it for hydrological analysis, and then calculates flow direction and 
flow accumulation using the pysheds library. The results are saved as GeoTIFF files 
and visualized for user QC.
"""

from pysheds.grid import Grid   # For handling raster and flow direction/accumulation calculations
import rasterio                # For reading/writing GeoTIFF files
from rasterio.crs import CRS    # For handling coordinate reference systems
import numpy as np             # For numerical operations and array handling
import matplotlib.pyplot as plt # For visualization
import matplotlib.colors as colors # For color normalization in plots

"""
Step 1: Read elevation raster and initialize the grid
The DEM raster is read into the 'grid' object, which will be used for further processing
"""

grid = Grid.from_raster("C:/Users/path/to/DEM.tif")  # Path to the DEM file
dem = grid.read_raster("C:/Users/path/to/DEM.tif")  # Read the raster data

"""
Step 2: Condition the DEM (Digital Elevation Model)
Fill pits, depressions, and resolve flats to ensure proper flow analysis
Pits, depressions, and flats can cause errors in hydrological models if not addressed
"""
# Fill pits in DEM
pit_filled_dem = grid.fill_pits(dem)
print("DEM Filled")
# Fill depressions in DEM
flooded_dem = grid.fill_depressions(pit_filled_dem)
print("DEM Flooded")
# Resolve flats in DEM
inflated_dem = grid.resolve_flats(flooded_dem)
print("DEM Inflated")

"""
Step 3: Define D8 flow direction mapping and compute
The D8 method considers the flow direction to one of the 8 surrounding cells (cardinal and diagonal)
"""
# Specify directional mapping for d8 flow routing
dirmap = (64, 128, 1, 2, 4, 8, 16, 32)
# Compute flow directions
fdir_d8 = grid.flowdir(inflated_dem, dirmap=dirmap)
print("fdir_d8 done")

"""
Step 4: Compute flow accumulation using flow direction
"""
# Calculate flow accumulation
acc_d8 = grid.accumulation(fdir_d8, dirmap=dirmap)
print("acc_d8 done")

"""
Step 5: Save flow direction and accumulation data as a GeoTIFF
    in the same coordinate system as the original data (update CRS if needed)
"""
# Save flow direction (fdir_d8) as a TIFF file
fdir_output_file = "C:/Users/path/to/flow/direction.tif"
crs = CRS.from_proj4("+proj=utm +zone=13 +ellps=GRS80 +units=m +no_defs")
fdir_d8 = fdir_d8.astype(np.float32)

with rasterio.open(
    fdir_output_file,
    'w',
    driver='GTiff',
    height=fdir_d8.shape[0],
    width=fdir_d8.shape[1],
    count=1,
    dtype=fdir_d8.dtype,
    crs=crs,
    transform=grid.affine
) as dst:
    dst.write(fdir_d8, 1)

print("fdir_d8 saved as TIFF")

# Save flow accumulation (acc_d8) as a TIFF file
acc_output_file = "C:/Users/path/to/flow/accumulation.tif"
crs = CRS.from_proj4("+proj=utm +zone=13 +ellps=GRS80 +units=m +no_defs")

with rasterio.open(
    acc_output_file,
    'w',
    driver='GTiff',
    height=acc_d8.shape[0],
    width=acc_d8.shape[1],
    count=1,
    dtype=acc_d8.dtype,
    crs=crs,
    transform=grid.affine
) as dst:
    dst.write(acc_d8, 1)

print("acc_d8 saved as TIFF")


# Plot for QC
# Plotting Flow Accumulation
fig, ax = plt.subplots(figsize=(8,6))
fig.patch.set_alpha(0)
plt.grid('on', zorder=0)
im = ax.imshow(acc_d8, extent=grid.extent, zorder=2,
               cmap='cubehelix',
               norm=colors.LogNorm(1, acc_d8.max()),
               interpolation='bilinear')
plt.colorbar(im, ax=ax, label='Upstream Cells')
plt.title('Flow Accumulation', size=14)
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.tight_layout()
plt.show()

# Plotting Flow Direction (optional)
fig, ax = plt.subplots(figsize=(8,6))
fig.patch.set_alpha(0)
plt.grid('on', zorder=0)
im = ax.imshow(fdir_d8, extent=grid.extent, zorder=2,
               cmap='viridis',
               interpolation='bilinear')
plt.colorbar(im, ax=ax, label='Flow Direction')
plt.title('Flow Direction', size=14)
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.tight_layout()
plt.show()
