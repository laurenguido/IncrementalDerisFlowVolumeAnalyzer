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
import matplotlib.pyplot as plt
from pysheds.grid import Grid
from shapely.geometry import LineString, Point, mapping
from shapely.ops import unary_union, linemerge
import numpy as np
from pyproj import Transformer
from shapely.ops import transform
import fiona
from fiona.crs import from_epsg
from osgeo import ogr


def find_max_value_point(polyline_path, raster_path):
    """
    Finds the point along a polyline that has the maximum raster value.

    Parameters
    ----------
    polyline_path : str
        Path to the polyline shapefile.
    raster_path : str
        Path to the raster file.

    Returns
    -------
    tuple
        Coordinates of the point with the maximum raster value (x, y), or None if not found.
    """
    driver = ogr.GetDriverByName('ESRI Shapefile')
    polyline_ds = driver.Open(polyline_path, 0)
    polyline_layer = polyline_ds.GetLayer()

    with rasterio.open(raster_path) as src:
        max_value = -float('inf')
        max_point = None
        for feature in polyline_layer:
            geometry = feature.GetGeometryRef()
            segment = LineString([geometry.GetPoint(i)[:2] for i in range(geometry.GetPointCount())])
            for point in segment.coords:
                point_geom = Point(point)
                row, col = src.index(point_geom.x, point_geom.y)
                value = src.read(1, masked=True)[row, col]
                if value > max_value:
                    max_value = value
                    max_point = point_geom.xy
        polyline_ds = None
        return max_point


def initialize_grid(direction_path, accumulation_path):
    """
    Initializes the grid and loads the flow direction and accumulation rasters.

    Parameters
    ----------
    direction_path : str
        Path to the flow direction raster.
    accumulation_path : str
        Path to the accumulation raster.

    Returns
    -------
    tuple
        (grid, acc, fdir): pysheds Grid, accumulation array, flow direction array.
    """
    grid = Grid.from_raster(direction_path)
    acc = grid.read_raster(accumulation_path)
    fdir = grid.read_raster(direction_path)
    return grid, acc, fdir


def calculate_distance_to_outlet(grid, acc, fdir, x, y):
    """
    Calculates the distance to the outlet from a given point on the grid.

    Parameters
    ----------
    grid : pysheds.grid.Grid
        Grid object containing the flow direction and accumulation data.
    acc : numpy.ndarray
        Flow accumulation array.
    fdir : numpy.ndarray
        Flow direction array.
    x : float
        X-coordinate of the point.
    y : float
        Y-coordinate of the point.

    Returns
    -------
    tuple
        (dist, x_snap, y_snap): distance to the outlet, snapped x/y coordinate.
    """
    dirmap = (64, 128, 1, 2, 4, 8, 16, 32)
    x_snap, y_snap = grid.snap_to_mask(acc > 1000, (x, y))
    dist = grid.distance_to_outlet(x=x_snap, y=y_snap, fdir=fdir, dirmap=dirmap, xytype='coordinate')
    return dist, x_snap, y_snap


def convert_coordinates(direction_path, x_snap, y_snap):
    """
    Converts raster grid coordinates to geographic coordinates.

    Parameters
    ----------
    direction_path : str
        Path to the flow direction raster.
    x_snap : float
        Snapped x-coordinate.
    y_snap : float
        Snapped y-coordinate.

    Returns
    -------
    tuple
        (lon, lat): Longitude and latitude of the point.
    """
    with rasterio.open(direction_path) as src:
        transformer = src.transform
        lon, lat = transformer * (x_snap, y_snap)
    return lon, lat


def merge_segments(layer):
    """
    Merges all line segments in a shapefile layer into a single continuous LineString.

    Parameters
    ----------
    layer : ogr.Layer
        Layer containing line features.

    Returns
    -------
    shapely.geometry.LineString
        A merged LineString geometry.
    """
    lines = []
    for feature in layer:
        geom = feature.GetGeometryRef()
        line = LineString(geom.GetPoints())
        lines.append(line)
    merged_line = linemerge(unary_union(lines))
    return merged_line


def interpolate_points_along_polyline(polyline_path, direction_path, accumulation_path, dist):
    """
    Interpolates points along a polyline and extracts flow distance and catchment values for each point.

    Parameters
    ----------
    polyline_path : str
        Path to the polyline shapefile.
    direction_path : str
        Path to the flow direction raster.
    accumulation_path : str
        Path to the accumulation raster.
    dist : numpy.ndarray
        Distance-to-outlet array.

    Returns
    -------
    tuple
        (all_points, distances, catchments): Lists of interpolated points, distances to outlet, and catchment values.
    """
    driver = ogr.GetDriverByName("ESRI Shapefile")
    dataSource = driver.Open(polyline_path, 0)
    layer = dataSource.GetLayer()

    polyline = merge_segments(layer)

    all_points = []
    distances = []
    catchments = []

    spacing = 10
    distance = np.arange(0, polyline.length, spacing)
    points = [polyline.interpolate(d) for d in distance]

    with rasterio.open(direction_path) as src:
        crs = src.crs
    transformer = Transformer.from_crs(crs.to_string(), "EPSG:26913", always_xy=True)
    points = [transform(lambda x, y: transformer.transform(x, y), point) for point in points]

    with rasterio.open(accumulation_path) as acc_src, rasterio.open(direction_path) as src:
        for point in points:
            px, py = point.x, point.y
            row, col = src.index(px, py)
            distance_value = float(dist[row, col])
            catchment_value = float(acc_src.read(1, masked=True)[row, col])
            distances.append(distance_value)
            catchments.append(catchment_value)
    all_points.extend(points)

    return all_points, distances, catchments


def save_points_to_shapefile(output_shapefile_path, all_points, distances, catchments):
    """
    Saves interpolated points and their attributes (distance and catchment) to a shapefile.

    Parameters
    ----------
    output_shapefile_path : str
        Path to the output shapefile.
    all_points : list of shapely.geometry.Point
        List of interpolated points.
    distances : list of float
        List of distances to outlet for each point.
    catchments : list of float
        List of catchment values for each point.

    Returns
    -------
    None
    """
    schema = {
        'geometry': 'Point',
        'properties': {
            'id': 'int',
            'Distance': 'float',
            'Catch_m2': 'float'
        },
    }

    with fiona.open(output_shapefile_path, 'w', driver='ESRI Shapefile', schema=schema, crs=from_epsg(26913)) as shp:
        for i, point in enumerate(all_points):
            shp.write({
                'geometry': mapping(point),
                'properties': {
                    'id': i,
                    'Distance': distances[i],
                    'Catch_m2': catchments[i]
                },
            })


def plot_flow_distance(grid, dist, all_points, polyline_path):
    """
    Plots generated points with hydrologic data for quality control.

    Parameters
    ----------
    grid : pysheds.grid.Grid
        Grid object for plotting extent.
    dist : numpy.ndarray
        Array of distance-to-outlet values.
    all_points : list of shapely.geometry.Point
        Points to plot.
    polyline_path : str
        Path to the polyline shapefile.

    Returns
    -------
    None
    """
    driver = ogr.GetDriverByName("ESRI Shapefile")
    dataSource = driver.Open(polyline_path, 0)
    layer = dataSource.GetLayer()
    polyline = merge_segments(layer)
    minx, miny, maxx, maxy = polyline.bounds

    x_range = maxx - minx
    y_range = maxy - miny
    x_expand = x_range * 0.1
    y_expand = y_range * 0.1

    minx -= x_expand
    maxx += x_expand
    miny -= y_expand
    maxy += y_expand

    fig, ax = plt.subplots(figsize=(10, 10))
    fig.patch.set_alpha(0)
    plt.grid('on', zorder=0)

    im = ax.imshow(dist, extent=grid.extent, zorder=2, cmap='cubehelix_r')
    plt.colorbar(im, ax=ax, label='Distance to outlet (cells)')
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.title('Flow Distance', size=14)

    ax.set_xlim(minx, maxx)
    ax.set_ylim(miny, maxy)

    for point in all_points:
        ax.plot(point.x, point.y, 'ro')

    plt.show()
