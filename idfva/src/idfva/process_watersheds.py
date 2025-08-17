# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 12:48:02 2024

@author: lguido

This module processes a DEM for hydrological analysis using pysheds. It fills pits, depressions, and resolves flats,
then calculates flow direction and flow accumulation, and saves the results as GeoTIFF files.
"""

from pysheds.grid import Grid
import rasterio
from rasterio.crs import CRS
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors


def read_dem(dem_path):
    """
    Read a DEM raster and initialize a pysheds Grid.

    Parameters
    ----------
    dem_path : str
        Path to the DEM raster file.

    Returns
    -------
    grid : pysheds.grid.Grid
        The initialized grid object.
    dem : np.ndarray
        The DEM data as an array.
    """
    grid = Grid.from_raster(dem_path)
    dem = grid.read_raster(dem_path)
    return grid, dem

def condition_dem(grid, dem):
    """
    Condition the DEM by filling pits, depressions, and resolving flats.

    Parameters
    ----------
    grid : pysheds.grid.Grid
        The grid object for raster analysis.
    dem : np.ndarray
        The DEM data.

    Returns
    -------
    inflated_dem : np.ndarray
        The conditioned DEM ready for hydrological processing.
    """
    pit_filled_dem = grid.fill_pits(dem)
    print("DEM Filled")
    flooded_dem = grid.fill_depressions(pit_filled_dem)
    print("DEM Flooded")
    inflated_dem = grid.resolve_flats(flooded_dem)
    print("DEM Inflated")
    return inflated_dem

def compute_flow_direction(grid, dem, dirmap=(64, 128, 1, 2, 4, 8, 16, 32)):
    """
    Compute D8 flow direction for the DEM.

    Parameters
    ----------
    grid : pysheds.grid.Grid
        The grid object.
    dem : np.ndarray
        The conditioned DEM.
    dirmap : tuple of int, optional
        Directional mapping for D8 flow routing.

    Returns
    -------
    fdir_d8 : np.ndarray
        Flow direction array.
    """
    fdir_d8 = grid.flowdir(dem, dirmap=dirmap)
    print("fdir_d8 done")
    return fdir_d8

def compute_flow_accumulation(grid, flowdir, dirmap=(64, 128, 1, 2, 4, 8, 16, 32)):
    """
    Compute flow accumulation from a flow direction array.

    Parameters
    ----------
    grid : pysheds.grid.Grid
        The grid object.
    flowdir : np.ndarray
        Flow direction array.
    dirmap : tuple of int, optional
        Directional mapping for D8 flow routing.

    Returns
    -------
    acc_d8 : np.ndarray
        Flow accumulation array.
    """
    acc_d8 = grid.accumulation(flowdir, dirmap=dirmap)
    print("acc_d8 done")
    return acc_d8

def save_raster(array, output_file, crs, transform, dtype=None):
    """
    Save a NumPy array as a GeoTIFF file.

    Parameters
    ----------
    array : np.ndarray
        The array to save.
    output_file : str
        Path to the output GeoTIFF file.
    crs : rasterio.crs.CRS
        Coordinate reference system for the output.
    transform : affine.Affine
        Affine transform for the output.
    dtype : np.dtype, optional
        Data type for output file.

    Returns
    -------
    None
    """
    if dtype is None:
        dtype = array.dtype
    with rasterio.open(
        output_file,
        'w',
        driver='GTiff',
        height=array.shape[0],
        width=array.shape[1],
        count=1,
        dtype=dtype,
        crs=crs,
        transform=transform
    ) as dst:
        dst.write(array, 1)
    print(f"{output_file} saved as TIFF")

def plot_array(array, grid, title, cmap, colorbar_label, log_norm=False):
    """
    Plot a raster array with matplotlib for quality control.

    Parameters
    ----------
    array : np.ndarray
        The array to plot.
    grid : pysheds.grid.Grid
        The grid object (for extent).
    title : str
        Title for the plot.
    cmap : str
        Colormap for the image.
    colorbar_label : str
        Label for the colorbar.
    log_norm : bool, optional
        Whether to use logarithmic normalization.

    Returns
    -------
    None
    """
    fig, ax = plt.subplots(figsize=(8, 6))
    fig.patch.set_alpha(0)
    plt.grid('on', zorder=0)
    norm = colors.LogNorm(1, array.max()) if log_norm else None
    im = ax.imshow(
        array,
        extent=grid.extent,
        zorder=2,
        cmap=cmap,
        norm=norm,
        interpolation='bilinear'
    )
    plt.colorbar(im, ax=ax, label=colorbar_label)
    plt.title(title, size=14)
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.tight_layout()
    plt.show()