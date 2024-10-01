# -*- coding: utf-8 -*-
"""
Created on Tue Feb 20 09:11:07 2024

@author: laurenguido

The code in this file (Las2Ras.py) processes LAS files to read scalar field values, prints a subset 
of these values for user QC, and converts the LAS data to a raster format. The raster is generated based on a 
specified scalar field from the LAS file.
"""

import laspy  # Library for reading and writing LAS files (LIDAR data).
import numpy as np  # For numerical operations and handling arrays.
import rasterio  # For reading and writing raster files (e.g., GeoTIFF).
from scipy.interpolate import griddata, Rbf, RectBivariateSpline  # Interpolation methods.
import matplotlib.tri as tri  # Used for triangulation and interpolation.


def print_scalar_field(las_file, scalar_field, num_values=10):
    """
    Reads a LAS file and prints the first `num_values` of the specified scalar field.

    Parameters:
        las_file (str): Path to the LAS file.
        scalar_field (str): The scalar field to print (e.g., intensity, elevation, etc.).
        num_values (int): Number of values to print (default is 10).
    """
    las_data = laspy.read(las_file)
    scalar_values = las_data[scalar_field]
    print(f"Scalar field values for '{scalar_field}':")
    for i in range(min(num_values, len(scalar_values))):
        print(f"Value {i+1}: {scalar_values[i]}")

def fill_nan_values(interpolated_values, x, y, scalar_values, grid_x, grid_y, method='nearest'):
    """
    Fills NaN values in the interpolated data using the specified interpolation method.

    Parameters:
        interpolated_values (ndarray): The array of interpolated values that may contain NaNs.
        x (ndarray): X-coordinates of the original LAS data points.
        y (ndarray): Y-coordinates of the original LAS data points.
        scalar_values (ndarray): Scalar values (e.g., intensity, elevation) from the LAS data.
        grid_x (ndarray): Grid of x-coordinates for the output raster.
        grid_y (ndarray): Grid of y-coordinates for the output raster.
        method (str): Interpolation method to use (default is 'nearest').

    Returns:
        filled_values (ndarray): The array of interpolated values with NaNs filled.
    """
    nan_indices = np.isnan(interpolated_values)
    valid_values = scalar_values[~np.isnan(scalar_values)]
    valid_x = x[~np.isnan(scalar_values)]
    valid_y = y[~np.isnan(scalar_values)]
    filled_values = interpolated_values.copy()
    filled_values[nan_indices] = griddata((valid_x, valid_y), valid_values, 
                                          (grid_x[nan_indices], grid_y[nan_indices]), 
                                          method=method)
    return filled_values

def las_to_raster(las_file, scalar_field, output_raster, cell_size, method):
    """
    Converts LAS data to a raster format using the specified interpolation method.

    Parameters:
        las_file (str): Path to the input LAS file.
        scalar_field (str): The scalar field to interpolate (e.g., intensity, elevation).
        output_raster (str): Path to the output raster file.
        cell_size (float): Size of each raster cell.
        method (str): Interpolation method ('nearest', 'linear', 'cubic', 'idw', 'spline', 'rbf', 'nn').

    This function reads the LAS data, interpolates it to create a grid, and writes the result as a raster.
    """
    las_data = laspy.read(las_file)
    x = las_data.x
    y = las_data.y
    scalar_values = las_data[scalar_field]
    
    print(f"Number of values: {len(scalar_values)}")


    min_value = np.nanmin(scalar_values)
    max_value = np.nanmax(scalar_values)
    print(f"Minimum value of {scalar_field}: {min_value}")
    print(f"Maximum value of {scalar_field}: {max_value}")

    min_x, min_y = np.min(x), np.min(y)
    max_x, max_y = np.max(x), np.max(y)
    width = int(np.ceil((max_x - min_x) / cell_size))
    height = int(np.ceil((max_y - min_y) / cell_size))
    print("Raster dimensions:", width, "x", height)

    transform = rasterio.transform.from_origin(min_x, max_y, cell_size, cell_size)

    profile = {
        'driver': 'GTiff',
        'count': 1,
        'dtype': rasterio.float32,
        'width': width,
        'height': height,
        'crs': '+init=epsg:26913',
        'transform': transform,
    }

    grid_x, grid_y = np.meshgrid(np.linspace(min_x, max_x, width), np.linspace(max_y, min_y, height))

    if method == 'nearest':
        interpolated_values = griddata((x, y), scalar_values, (grid_x, grid_y), method='nearest')
    elif method == 'linear':
        interpolated_values = griddata((x, y), scalar_values, (grid_x, grid_y), method='linear')
    elif method == 'cubic':
        interpolated_values = griddata((x, y), scalar_values, (grid_x, grid_y), method='cubic')
    elif method == 'idw':
        interpolated_values = idw_interpolation(x, y, scalar_values, grid_x, grid_y, power=2)
    elif method == 'spline':
        interpolated_values = rect_bivariate_spline_interpolation(x, y, scalar_values, grid_x, grid_y)
    elif method == 'rbf':
        interpolated_values = rbf_interpolation(x, y, scalar_values, grid_x, grid_y, function='multiquadric')
    elif method == 'nn':
        interpolated_values = natural_neighbor(x, y, scalar_values, grid_x, grid_y)
    else:
        raise ValueError(f"Unknown interpolation method: {method}")

    nan_count = np.isnan(interpolated_values).sum()
    print(f"Number of NaN values in the interpolated data: {nan_count}")

    # Fill NaN values
    interpolated_values = fill_nan_values(interpolated_values, x, y, scalar_values, grid_x, grid_y, method='nearest')

    # Check for NaN values in the interpolated data
    nan_count = np.isnan(interpolated_values).sum()
    print(f"Number of NaN values in the interpolated and filled data: {nan_count}")

    with rasterio.open(output_raster, 'w', **profile) as dst:
        print("Writing to raster...")
        dst.write(interpolated_values.astype(np.float32), 1)
        print(f"Raster created using {method} interpolation")

def idw_interpolation(x, y, values, grid_x, grid_y, power=2):
    """
    Performs Inverse Distance Weighting (IDW) interpolation.

    Parameters:
        x, y (ndarray): X and Y coordinates of the data points.
        values (ndarray): Scalar values at the data points.
        grid_x, grid_y (ndarray): Grids for interpolation.
        power (int): The power parameter for IDW (default is 2).

    Returns:
        interpolated_values (ndarray): Interpolated values on the grid.
    """
    interpolated_values = np.zeros_like(grid_x)
    for i in range(grid_x.shape[0]):
        for j in range(grid_x.shape[1]):
            distances = np.sqrt((x - grid_x[i, j])**2 + (y - grid_y[i, j])**2)
            if np.any(distances == 0):
                interpolated_values[i, j] = values[np.argmin(distances)]
            else:
                weights = 1 / distances**power
                interpolated_values[i, j] = np.sum(weights * values) / np.sum(weights)
    return interpolated_values

def rect_bivariate_spline_interpolation(x, y, values, grid_x, grid_y):
    """
    Performs interpolation using RectBivariateSpline.

    Parameters:
        x, y (ndarray): X and Y coordinates of the data points.
        values (ndarray): Scalar values at the data points.
        grid_x, grid_y (ndarray): Grids for interpolation.

    Returns:
        interpolated_values (ndarray): Interpolated values on the grid.
    """
    try:
        spline = RectBivariateSpline(x, y, values)
        result = spline.ev(grid_x, grid_y)
        print("Spline interpolation completed successfully.")
        return result
    except MemoryError as e:
        print("MemoryError:", e)
    except Exception as e:
        print("Exception during spline interpolation:", e)
    return np.full(grid_x.shape, np.nan)

def rbf_interpolation(x, y, values, grid_x, grid_y, function='multiquadric'):
    """
    Performs interpolation using Radial Basis Function (RBF).
    
    Parameters:
        x, y (ndarray): X and Y coordinates of the data points.
        values (ndarray): Scalar values at the data points.
        grid_x, grid_y (ndarray): Grids for interpolation.
        function (str): The type of RBF to use (default is 'multiquadric').
    
    Returns:
        interpolated_values (ndarray): Interpolated values on the grid.
    """
    rbf = Rbf(x, y, values, function=function)
    return rbf(grid_x, grid_y)

def natural_neighbor(x, y, values, grid_x, grid_y):
    """
    Performs Natural Neighbor interpolation using Delaunay triangulation.

    Parameters:
        x, y (ndarray): X and Y coordinates of the data points.
        values (ndarray): Scalar values at the data points.
        grid_x, grid_y (ndarray): Grids for interpolation.

    Returns:
        interpolated_values (ndarray): Interpolated values on the grid.
    """
    triangulation = tri.Triangulation(x, y)
    interpolator = tri.LinearTriInterpolator(triangulation, values)
    interpolated_values = interpolator(grid_x, grid_y)
    return interpolated_values

# Example usage
las_file = "C:/Users/path/to/input/las.las"
scalar_field = "M3C2 distance"
output_raster = "C:/Users/path/to/output/raster.tif"
cell_size = 1
method = 'linear'
print_scalar_field(las_file, scalar_field, num_values=10)
las_to_raster(las_file, scalar_field, output_raster, cell_size, method)
