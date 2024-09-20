# -*- coding: utf-8 -*-
"""
Created on Thu Aug  8 11:35:05 2024

@author: lguido

This script (Apply_Modified_Voronoi.py) processes the query point and debris flow path shapefiles used previously to 
generate modified Voronoi polygons within the specified buffer (from Prepare_Path.py) of expected change from the debris 
flow, and plots the result on a map for visualization/QC.
"""

import matplotlib.pyplot as plt  # For plotting maps and visualizations.
from shapely.geometry import Point, Polygon, mapping, shape  # For handling geometric objects.
from shapely.ops import unary_union  # For geometric operations (e.g., union of polygons).
import fiona  # For reading/writing shapefiles.
import numpy as np  # For numerical operations.
from scipy.stats import gaussian_kde  # For kernel density estimation (not used directly here).
from scipy.spatial import Voronoi  # For generating Voronoi diagrams.

def read_points_from_shapefile(point_file):
    """
    Reads points from a shapefile.

    Parameters:
        point_file (str): Path to the input point shapefile.

    Returns:
        points (list): A list of Shapely Point objects.
    """
    with fiona.open(point_file, 'r') as src:
        return [Point(feature['geometry']['coordinates']) for feature in src]

def compute_voronoi_polygons(points, original_polygon, buffer_distance):
    """
    Computes Voronoi polygons from a set of points and clips them to a buffered region of an original polygon.

    Parameters:
        points (list): A list of Shapely Point objects.
        original_polygon (Polygon): The original polygon to clip Voronoi polygons to.
        buffer_distance (float): Distance for buffering the original polygon.

    Returns:
        voronoi_polygons (list): A list of clipped Voronoi polygons.
    """
    points_array = np.array([[p.x, p.y] for p in points])
    vor = Voronoi(points_array)

    buffered_polygon = original_polygon.buffer(buffer_distance)

    voronoi_polygons = []

    for region in vor.regions:
        if not -1 in region and region:
            polygon = Polygon([vor.vertices[i] for i in region])
            if polygon.is_valid:
                # Clip the Voronoi polygon to the buffered polygon
                clipped_polygon = polygon.intersection(buffered_polygon)
                # Further clip to the original polygon (unbuffered)
                final_polygon = clipped_polygon.intersection(original_polygon)
                if not final_polygon.is_empty:
                    voronoi_polygons.append(final_polygon)

    return voronoi_polygons

def plot_voronoi_on_map(points, extra_points, voronoi_polygons, debris_flow_path):
    """
    Plots Voronoi polygons, points, and a polyline on a map for visualization.

    Parameters:
        points (list): Points along the polyline (Shapely Point objects).
        extra_points (list): Extra points (Shapely Point objects) to plot.
        voronoi_polygons (list): List of Voronoi polygons.
        polyline_file (str): Path to the polyline shapefile.
    """
    with fiona.open(debris_flow_path, 'r') as src:
        polyline_coords = src[0]['geometry']['coordinates']
    
    x, y = zip(*polyline_coords)
    aspect_ratio = (max(x) - min(x)) / (max(y) - min(y))

    plt.figure(figsize=(10 * aspect_ratio, 10))
    plt.plot(x, y, color='blue', linewidth=2, label='Polyline')
    plt.scatter([p.x for p in points], [p.y for p in points], color='red', label='Points along Polyline')
    plt.scatter([p.x for p in extra_points], [p.y for p in extra_points], color='purple', label='Extra Points')

    for poly in voronoi_polygons:
        plt.plot(*poly.exterior.xy, color='green', linewidth=1, alpha=0.7)

    plt.title('Voronoi Polygons Generated Along Polyline')
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.legend()
    plt.grid(True)
    plt.show()

def save_voronoi_polygons_to_shapefile(voronoi_polygons, shapefile_path):
    """
    Save a list of Voronoi polygons to a shapefile.

    Parameters:
    - voronoi_polygons (list): A list of shapely Polygon objects representing the Voronoi polygons.
    - shapefile_path (str): The file path where the shapefile will be saved.

    This function creates a shapefile with each polygon's geometry and an associated ID and area in square meters.
    The schema for the shapefile includes 'id' as an integer and 'area_m2' as a float.
    """
    schema = {
        'geometry': 'Polygon',
        'properties': {
            'id': 'int',
            'area_m2': 'float'
        },
    }

    with fiona.open(shapefile_path, 'w', 'ESRI Shapefile', schema) as c:
        for idx, poly in enumerate(voronoi_polygons):
            c.write({
                'geometry': mapping(poly),
                'properties': {'id': idx, 'area_m2': poly.area},
            })

def plot_voronoi_polygon_areas(voronoi_polygons):
    """
    Plot a histogram and kernel density estimate (KDE) of the areas of Voronoi polygons.

    Parameters:
    - voronoi_polygons (list): A list of shapely Polygon objects representing the Voronoi polygons.

    This function generates a histogram of the areas of the provided Voronoi polygons,
    along with a KDE curve to illustrate the distribution of polygon areas.
    Vertical lines are drawn at the mean and standard deviations of the area.
    The x-axis represents the area, while the y-axis shows the density or frequency.
    """
    areas = [poly.area for poly in voronoi_polygons]
    
    plt.figure(figsize=(8, 6))
    hist_values, bins, _ = plt.hist(areas, bins=30, edgecolor='black', alpha=0.7, density=True, label='Histogram')
    
    kde = gaussian_kde(areas)
    x_values = np.linspace(min(areas), max(areas), 100)
    plt.plot(x_values, kde(x_values), color='red', linewidth=2, label='Density Function')
    
    mean_area, std_dev_area = np.mean(areas), np.std(areas)
    for i in range(1, 4):
        plt.axvline(x=mean_area + i * std_dev_area, color='blue', linestyle='--', linewidth=1, label=f'{i}$\sigma$')
        plt.axvline(x=mean_area - i * std_dev_area, color='blue', linestyle='--', linewidth=1)
    
    plt.xlabel('Area')
    plt.ylabel('Density / Frequency')
    plt.title('Histogram and Density Function of Voronoi Polygon Areas')
    plt.grid(True)
    plt.legend()
    plt.show()

def voronoi_polygon_qc(investigation_polygons_shapefile):
    """
    Perform a quality check on Voronoi polygons to ensure they touch exactly two other polygons.

    Parameters:
    - investigation_polygons_shapefile (str): The file path to the shapefile containing the investigation polygons.

    This function reads the specified shapefile and checks each polygon to count how many
    other polygons it touches. If a polygon touches a number of other polygons that is
    not equal to two, a message is printed indicating the polygon ID and the count.
    This quality check is useful for validating the Voronoi polygon geometry.
    """
    adjacent_counts = {}
    
    with fiona.open(investigation_polygons_shapefile, 'r') as src:
        features = [shape(feat['geometry']) for feat in src]
        ids = [feat['properties']['id'] for feat in src]

    for idx, polygon in enumerate(features):
        adjacent_counts[ids[idx]] = sum(1 for poly in features if polygon.touches(poly) and poly != polygon)

    for poly_id, count in adjacent_counts.items():
        if count != 2:
            print(f"Polygon {poly_id} touches {count} other polygons instead of 2.")
    
    print("Voronoi polygon quality check completed.")

# Parameters
debris_flow_path = "C:/Users/path/to/flowpath.shp"
point_file = "C:/Users/path/to/querypoints.shp"
buffered_path = "C:/Users/path/to/buffered/polygon.shp"
investigation_polygons_shapefile = "C:/Users/path/to/investigation/polygons.shp"
buffer_distance = 0 # Additional buffer that is ultimately masked out, sometimes needed to generate 
                    # appropriate polygons in complex flow networks

# Execution
points = read_points_from_shapefile(point_file)
points_array = np.array([[p.x, p.y] for p in points])
x_min, y_min = points_array.min(axis=0)
x_max, y_max = points_array.max(axis=0)
extra_points = [
    Point(x_min - buffer_distance, y_min - buffer_distance),
    Point(x_min - buffer_distance, y_max + buffer_distance),
    Point(x_max + buffer_distance, y_min - buffer_distance),
    Point(x_max + buffer_distance, y_max + buffer_distance)
]
points.extend(extra_points)

# Read the original (unbuffered) polygon
with fiona.open(buffered_path, 'r') as src:
    original_polygon = unary_union([Polygon(feat['geometry']['coordinates'][0]) for feat in src])

clipped_voronoi_polygons = compute_voronoi_polygons(points, original_polygon, buffer_distance)
save_voronoi_polygons_to_shapefile(clipped_voronoi_polygons, investigation_polygons_shapefile)
plot_voronoi_on_map(points, extra_points, clipped_voronoi_polygons, debris_flow_path)
plot_voronoi_polygon_areas(clipped_voronoi_polygons)
voronoi_polygon_qc(investigation_polygons_shapefile)
