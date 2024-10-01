# -*- coding: utf-8 -*-
"""
Created on Fri Aug 16 12:25:29 2024

@author: lguido

The code in this file (Generate_Query_points.py) locates the point of greatest accumulation 
(the outlet) of a user-defined debris flow path, and sets this as a pour point with
the distance to the outlet at this point equal to zero. Points are then generated moving 
upstream at a user-defined spacing (default 10 m) for future use centering voronoi 
polygons. A distance upstream and catchment area are assigned to each point.

This script has been developed to handle simple (non-branching) as well as complex,
multi-strahler-order branching debris flow paths. 
"""

import rasterio
# rasterio: For reading raster data files (e.g., flow direction, accumulation rasters).
import matplotlib.pyplot as plt
# matplotlib: For generating plots (e.g., scatter plot of points color-coded by distance).
from pysheds.grid import Grid
# pysheds: For creating grids and analyzing hydrological data (e.g., flow direction and distance to outlet).
from shapely.geometry import LineString, Point, mapping
from shapely.ops import unary_union, linemerge
# shapely: For handling geometric operations like creating/interpolating points along lines.
import numpy as np
# numpy: For numerical computations, including distance calculations and array manipulations.
from pyproj import Transformer
# pyproj: For converting between coordinate reference systems (CRS).
from shapely.ops import transform
# shapely.ops.transform: For transforming geometries between coordinate systems.
import fiona
from fiona.crs import from_epsg
# fiona: For reading/writing shapefiles and handling spatial data attributes (e.g., point geometries, attributes).
from osgeo import ogr
# osgeo: For accessing and manipulating shapefile features, layers, and attributes (e.g., polyline data).


# Function to find the point along a polyline with the maximum raster value
def find_max_value_point(polyline_path, raster_path):
    """
    Finds the point along a polyline that has the maximum raster value.
    
    Args:
        polyline_path (str): Path to the polyline shapefile.
        raster_path (str): Path to the raster file.
    
    Returns:
        tuple: Coordinates of the point with the maximum raster value (x, y), or None if not found.
    """
    # Open the polyline shapefile using the OGR driver
    driver = ogr.GetDriverByName('ESRI Shapefile')
    polyline_ds = driver.Open(polyline_path, 0)
    polyline_layer = polyline_ds.GetLayer()
    
    # Open the raster file using rasterio
    with rasterio.open(raster_path) as src:
        max_value = -float('inf')  # Initialize the max value as negative infinity
        max_point = None  # Initialize the max point as None
      
        # Loop through each feature (line) in the polyline shapefile
        for feature in polyline_layer:
            # Extract the geometry of the line as a LineString
            geometry = feature.GetGeometryRef()
            segment = LineString([geometry.GetPoint(i)[:2] for i in range(geometry.GetPointCount())])

            # Loop through each point (vertex) along the LineString
            for point in segment.coords:
                point_geom = Point(point)  # Create a Point geometry
                row, col = src.index(point_geom.x, point_geom.y)  # Get the raster row and column index for the point
                value = src.read(1, masked=True)[row, col]  # Read the raster value at this point
                
                # Update max_value and max_point if the current value is greater than the max_value
                if value > max_value:
                    max_value = value
                    max_point = point_geom.xy  # Store the coordinates of the point
        
        polyline_ds = None  # Close the shapefile
        return max_point  # Return the coordinates of the point with the maximum raster value

# Function to initialize the grid and read the accumulation and flow direction rasters
def initialize_grid(direction_path, accumulation_path):
    """
    Initializes the grid and loads the flow direction and accumulation rasters.
    
    Args:
        direction_path (str): Path to the flow direction raster.
        accumulation_path (str): Path to the accumulation raster.
    
    Returns:
        tuple: A tuple containing the grid, the accumulation array, and the flow direction array.
    """
    grid = Grid.from_raster(direction_path)  # Create a grid from the flow direction raster
    acc = grid.read_raster(accumulation_path)  # Read the accumulation raster into the grid
    fdir = grid.read_raster(direction_path)  # Read the flow direction raster into the grid
    return grid, acc, fdir  # Return the grid, accumulation, and flow direction rasters

# Function to calculate the distance to the outlet from a given point
def calculate_distance_to_outlet(grid, acc, fdir, x, y):
    """
    Calculates the distance to the outlet from a given point on the grid.
    
    Args:
        grid (pysheds.Grid): Grid object containing the flow direction and accumulation data.
        acc (numpy.ndarray): Flow accumulation array.
        fdir (numpy.ndarray): Flow direction array.
        x (float): X-coordinate of the point.
        y (float): Y-coordinate of the point.
    
    Returns:
        tuple: Distance to the outlet, snapped x-coordinate, and snapped y-coordinate.
    """
    dirmap = (64, 128, 1, 2, 4, 8, 16, 32)  # Direction mapping for flow directions
    x_snap, y_snap = grid.snap_to_mask(acc > 1000, (x, y))  # Snap the point to the nearest cell with accumulation > 1000
    dist = grid.distance_to_outlet(x=x_snap, y=y_snap, fdir=fdir, dirmap=dirmap, xytype='coordinate')  # Calculate the distance to the outlet
    return dist, x_snap, y_snap  # Return the distance to the outlet, and the snapped coordinates

# Function to convert coordinates from raster grid to geographic coordinates
def convert_coordinates(direction_path, x_snap, y_snap):
    """
    Converts raster grid coordinates to geographic coordinates.
    
    Args:
        direction_path (str): Path to the flow direction raster.
        x_snap (float): Snapped x-coordinate.
        y_snap (float): Snapped y-coordinate.
    
    Returns:
        tuple: Longitude and latitude of the point.
    """
    with rasterio.open(direction_path) as src:
        transformer = src.transform  # Get the affine transform from the raster
        lon, lat = transformer * (x_snap, y_snap)  # Transform the snapped coordinates to geographic coordinates
    return lon, lat  # Return the longitude and latitude

def merge_segments(layer):
    """
    Merges all line segments in a shapefile layer into a single continuous LineString.
    
    Args:
        layer (ogr.Layer): Layer containing line features.
    
    Returns:
        shapely.geometry.LineString: A merged LineString geometry.
    """
    lines = []
    for feature in layer:
        geom = feature.GetGeometryRef()
        line = LineString(geom.GetPoints())
        lines.append(line)
    
    merged_line = linemerge(unary_union(lines))  # Merge all lines into one continuous LineString
    return merged_line

# Function to interpolate points along a polyline and extract distance and catchment values
def interpolate_points_along_polyline(polyline_path, direction_path, accumulation_path, dist):
    """
    Interpolates points along a polyline and extracts flow distance and catchment values for each point.
    
    Args:
        polyline_path (str): Path to the polyline shapefile.
        direction_path (str): Path to the flow direction raster.
        accumulation_path (str): Path to the accumulation raster.
        dist (numpy.ndarray): Distance-to-outlet array.
    
    Returns:
        tuple: Lists of interpolated points, distances to outlet, and catchment values.
    """
    driver = ogr.GetDriverByName("ESRI Shapefile")
    dataSource = driver.Open(polyline_path, 0)
    layer = dataSource.GetLayer()

    # Merge all segments into one continuous LineString
    polyline = merge_segments(layer)

    all_points = []  # List to store all interpolated points
    distances = []  # List to store distances to the outlet for each point
    catchments = []  # List to store catchment values for each point

    spacing = 10  # Define the spacing between interpolated points
    distance = np.arange(0, polyline.length, spacing)  # Generate an array of distances along the polyline
    points = [polyline.interpolate(d) for d in distance]  # Interpolate points along the polyline

    with rasterio.open(direction_path) as src:
        crs = src.crs  # Get the CRS of the raster
    transformer = Transformer.from_crs(crs.to_string(), "EPSG:26913", always_xy=True)  # Define a coordinate transformer
    points = [transform(lambda x, y: transformer.transform(x, y), point) for point in points]  # Transform the points to the target CRS

    with rasterio.open(accumulation_path) as acc_src:
        for point in points:
            px, py = point.x, point.y
            row, col = src.index(px, py)  # Get the raster index for the point
            distance_value = float(dist[row, col])  # Get the distance to the outlet from the distance raster
            catchment_value = float(acc_src.read(1, masked=True)[row, col])  # Get the catchment value from the accumulation raster
            distances.append(distance_value)  # Append the distance value to the list
            catchments.append(catchment_value)  # Append the catchment value to the list
    
    all_points.extend(points)  # Append the points to the list of all points
    
    return all_points, distances, catchments  # Return the points, distances, and catchment values

# Function to save interpolated points and attributes to a shapefile
def save_points_to_shapefile(output_shapefile_path, all_points, distances, catchments):
    """
    Saves interpolated points and their attributes (distance and catchment) to a shapefile.
    
    Args:
        output_shapefile_path (str): Path to the output shapefile.
        all_points (list): List of interpolated points.
        distances (list): List of distances to outlet for each point.
        catchments (list): List of catchment values for each point.
    """
    schema = {
        'geometry': 'Point',  # Define the geometry type as Point
        'properties': {
            'id': 'int',  # ID field as integer
            'Distance': 'float',  # Distance field as float
            'Catch_m2': 'float'  # Catchment field as float
        },
    }

    # Open a new shapefile for writing with the defined schema
    with fiona.open(output_shapefile_path, 'w', driver='ESRI Shapefile', schema=schema, crs=from_epsg(26913)) as shp:
        for i, point in enumerate(all_points):
            shp.write({
                'geometry': mapping(point),  # Convert the point geometry to GeoJSON format
                'properties': {
                    'id': i,  # Assign an ID to each point
                    'Distance': distances[i],  # Assign the distance value to each point
                    'Catch_m2': catchments[i]  # Assign the catchment value to each point
                },
            })

# Function to plot flow distance and interpolated points
def plot_flow_distance(grid, dist, all_points, polyline_path):
    """
    Plots generated points with hydrologic data for QC

    """
    # Create a bounding box around the polyline
    driver = ogr.GetDriverByName("ESRI Shapefile")
    dataSource = driver.Open(polyline_path, 0)
    layer = dataSource.GetLayer()
    polyline = merge_segments(layer)
    minx, miny, maxx, maxy = polyline.bounds  # Get the bounding box of the polyline
    
    # Calculate the expansion amount
    x_range = maxx - minx
    y_range = maxy - miny
    x_expand = x_range * 0.1
    y_expand = y_range * 0.1
    
    # Expand the bounding box
    minx -= x_expand
    maxx += x_expand
    miny -= y_expand
    maxy += y_expand
    
    fig, ax = plt.subplots(figsize=(10, 10))  # Create a figure and axes for the plot
    fig.patch.set_alpha(0)  # Set the figure background to be transparent
    plt.grid('on', zorder=0)  # Enable grid lines
    
    # Display the distance raster with a colormap
    im = ax.imshow(dist, extent=grid.extent, zorder=2, cmap='cubehelix_r')
    plt.colorbar(im, ax=ax, label='Distance to outlet (cells)')  # Add a colorbar with a label
    plt.xlabel('Longitude')  # Label the x-axis
    plt.ylabel('Latitude')  # Label the y-axis
    plt.title('Flow Distance', size=14)  # Set the plot title
    
    # Set the plot limits to the bounding box of the polyline
    ax.set_xlim(minx, maxx)
    ax.set_ylim(miny, maxy)
    
    for point in all_points:
        ax.plot(point.x, point.y, 'ro')  # Plot each point on the map in red
    
    plt.show()  # Display the plot


# Main function to run all steps
def main():
    # File paths for input and output data
    accumulation_path = "C:/Users/path/to/flow/accumulation.tif"
    direction_path = "C:/Users/path/to/flow/direction.tif"
    debris_flow_path = "C:/Users/path/to/flowpath.shp"
    output_shapefile_path = "C:/Users/path/to/querypoints.shp"

    # Find the point on the stream network with the maximum accumulation value
    max_point = find_max_value_point(debris_flow_path, accumulation_path)
    
    if max_point:
        X, Y = max_point
        x = X[0]  # Extract the x-coordinate of the maximum point
        y = Y[0]  # Extract the y-coordinate of the maximum point

        # Initialize the grid and read the accumulation and flow direction rasters
        grid, acc, fdir = initialize_grid(direction_path, accumulation_path)

        # Calculate the distance to the outlet from the maximum accumulation point
        dist, x_snap, y_snap = calculate_distance_to_outlet(grid, acc, fdir, x, y)

        # Convert the snapped coordinates to geographic coordinates
        lon, lat = convert_coordinates(direction_path, x_snap, y_snap)

        # Interpolate points along the polyline and extract distances and catchment values
        all_points, distances, catchments = interpolate_points_along_polyline(debris_flow_path, direction_path, accumulation_path, dist)

        # Save the interpolated points and attributes to a shapefile
        save_points_to_shapefile(output_shapefile_path, all_points, distances, catchments)

        # Plot the flow distance and interpolated points
        plot_flow_distance(grid, dist, all_points, debris_flow_path)

if __name__ == "__main__":
    main()  # Run the main function
