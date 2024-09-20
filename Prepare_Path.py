# -*- coding: utf-8 -*-
"""
Created on Mon Aug 12 07:50:46 2024

@author: lguido

This script (Prepare_PAth.py) reads the user-defined debris flow path shapefile, dissolves the 
polylines into a single multiline object, creates a buffer around it at a user-specified distance, and 
saves the result as a new shapefile. The script also allows for plotting of both the original polylines 
and the resulting buffer for QC. Projection transformations are handled as needed, and the output can be
visualized using matplotlib.
"""

import fiona  # For reading and writing shapefiles.
from shapely.geometry import shape, mapping, Polygon, MultiPolygon  # Handles geometric operations.
from shapely.ops import unary_union, transform  # Unary_union combines multiple geometries into one.
import matplotlib.pyplot as plt  # For plotting shapefiles.
import pyproj  # Library for projection transformations.
from functools import partial  # Allows partial function application.
from matplotlib.patches import Polygon as mpl_Polygon  # For creating matplotlib polygons.


def create_dissolved_buffer(debris_flow_path, output_shapefile, buffer_width):
    """
    Reads a polyline shapefile, dissolves all polylines into a single object, and creates a buffer around the result.
    The buffer is then saved to a new shapefile.

    Parameters:
        polyline_shapefile (str): Path to the input polyline shapefile.
        output_shapefile (str): Path to the output shapefile where the buffer will be saved.
        buffer_width (float): The width of the buffer to create around the dissolved polyline.

    Returns:
        geometries (list): List of original polyline geometries from the input shapefile.
        buffer (Polygon): The resulting buffer around the dissolved polylines.
        crs (dict): Coordinate reference system (CRS) of the input shapefile.
    """
    # Read the polyline shapefile
    with fiona.open(debris_flow_path, 'r') as src:
        # Extract the geometries from the shapefile
        geometries = [shape(feature['geometry']) for feature in src]
        crs = src.crs
        
        # Dissolve (unify) all polylines into a single multiline object
        multiline = unary_union(geometries)
        
        # Create a buffer around the dissolved multiline
        buffer = multiline.buffer(buffer_width)
        
        # Define a schema for the output shapefile
        schema = {
            'geometry': 'Polygon',
            'properties': {key: val for key, val in src.schema['properties'].items()}
        }
        
        # Write the buffer to a new shapefile
        with fiona.open(output_shapefile, 'w', driver='ESRI Shapefile', crs=src.crs, schema=schema) as dst:
            # Create a new feature from the buffer and write it to the output shapefile
            dst.write({
                'geometry': mapping(buffer),
                'properties': {key: None for key in schema['properties']}  # No properties are carried over
            })
    
    return geometries, buffer, crs

def plot_shapefiles(original_geometries, buffer, crs):
    """
    Plots the original polyline geometries and the resulting buffer polygon.

    Parameters:
        original_geometries (list): List of geometries representing the original polylines.
        buffer (Polygon or MultiPolygon): The buffer polygon created around the dissolved polylines.
        crs (dict): Coordinate reference system (CRS) of the shapefile.
    
    This function creates a plot with the original geometries and the buffer overlaid for visualization purposes.
    """
    # Project the buffer based on CRS (if needed)
    if crs['init'].startswith('epsg:26913'):  
        project = partial(
            pyproj.transform,
            pyproj.Proj(init=crs['init']),  # Source CRS
            pyproj.Proj(proj='merc', datum='WGS84')  # Projected CRS (Mercator)
        )
    else:
        project = None
    
    # Project the original geometries
    projected_geometries = [transform(project, geom) for geom in original_geometries] if project else original_geometries
    projected_buffer = transform(project, buffer) if project else buffer
    
    # Debug: Check the buffer type and properties
    print(f"Buffer Type: {type(projected_buffer)}")
    print(f"Buffer is valid: {projected_buffer.is_valid if hasattr(projected_buffer, 'is_valid') else 'N/A'}")
    print(f"Buffer is empty: {projected_buffer.is_empty}")
    
    # Check if the buffer is valid and non-empty
    if projected_buffer.is_empty or not isinstance(projected_buffer, (Polygon, MultiPolygon)):
        print("The buffer geometry is empty or invalid. Please check the input geometries and buffer width.")
        return
    
    # Create a plot
    fig, ax = plt.subplots()
    
    # Plot the original polylines
    for geom in projected_geometries:
        if geom.is_empty:
            continue
        x, y = geom.xy
        ax.plot(x, y, color='blue', label='Original Polyline')
    
    # Plot the buffer
    if isinstance(projected_buffer, Polygon):
        # Handle single Polygon case
        exterior_coords = list(projected_buffer.exterior.coords)
        patch = mpl_Polygon(exterior_coords, facecolor='green', edgecolor='black', alpha=0.5, label='Dissolved Buffer')
        ax.add_patch(patch)
    elif isinstance(projected_buffer, MultiPolygon):
        # Handle MultiPolygon case
        for poly in projected_buffer:
            exterior_coords = list(poly.exterior.coords)
            patch = mpl_Polygon(exterior_coords, facecolor='green', edgecolor='black', alpha=0.5, label='Dissolved Buffer')
            ax.add_patch(patch)
    
    # Set aspect ratio to equal to avoid distortion
    ax.set_aspect('equal')
    
    # Add title and legend
    ax.set_title('Original Polyline and Dissolved Buffer (Projected)')
    ax.legend()
    
    # Display the plot
    plt.show()

# User-defined inputs
debris_flow_path = "C:/Users/path/to/flowpath.shp"
output_shapefile = "C:/Users/path/to/buffered/polygon.shp"
buffer_width = 30.0  # Replace with your desired buffer width

# Create the dissolved buffer and get the geometries
original_geometries, buffer, crs = create_dissolved_buffer(debris_flow_path, output_shapefile, buffer_width)

# Plot the original polyline and the new buffer
plot_shapefiles(original_geometries, buffer, crs)


