# -*- coding: utf-8 -*-
"""
Created on Tue Jul  9 10:09:19 2024

@author: lguido

This script (Estimate_Volume_and_Plot.py) combines data and outputs from the preceeding scripts to 
process change detection and query point data to estimate volumes of erosion and deposition within
each investigation polygon. This script also assigns the distance upstream and catchment area of the 
query point each investigation polygon is centered on to the attribute table of the polygon for ease 
of analysis and plotting. This script also applies a limit of detection (user-defined) before querying 
volume data. Finally, a series of plots are generated for preliminary data analysis and QC.The user 
may define the ID number of known tributary polygons to have them plotted as vertical lines. 
"""
import rasterio  # For reading and writing raster files (e.g., GeoTIFF).
from rasterio.mask import mask  # For masking raster data with polygon geometries.
import fiona  # For reading and writing shapefiles.
from shapely.geometry import shape, box, LineString, MultiLineString  # For geometric operations.
import numpy as np  # For numerical operations and handling arrays.
import matplotlib.pyplot as plt  # For plotting geometries and attributes.
from matplotlib.colors import Normalize  # For normalizing color maps.
from matplotlib.cm import ScalarMappable  # For creating a colormap based on values.
from matplotlib.patches import Polygon  # For adding polygons to plots.



# Function to calculate the length of the polyline within each polygon
def calculate_flow_length(debris_flow_path, polygon_geom):
    """
    Calculates the total length of the polyline within a specified polygon geometry.

    Parameters:
        debris_flow_path (str): Path to the polyline shapefile.
        polygon_geom (shapely.geometry): The polygon geometry to intersect with.

    Returns:
        float: Total length of the polyline within the polygon geometry.
    """
    total_length = 0.0
    with fiona.open(debris_flow_path, 'r') as polyline_shapefile:
        for line_feature in polyline_shapefile:
            line_geom = shape(line_feature['geometry'])
            intersection = line_geom.intersection(polygon_geom)
            if isinstance(intersection, LineString):
                total_length += intersection.length
            elif isinstance(intersection, MultiLineString):
                total_length += sum(line.length for line in intersection)
    return total_length


# Function to extract raster values and cell areas within polygons, with LoD
def extract_raster_values(raster_path, investigation_polygons_path, debris_flow_path, LoD):
    """
    Extracts raster values and calculates associated metrics within polygon geometries 
    from a specified shapefile. Updates polygon attributes with calculated volumes and lengths.

    Parameters:
        raster_path (str): Path to the input raster file.
        investigation_polygons_path (str): Path to the input shapefile containing polygons.
        debris_flow_path (str): Path to the polyline shapefile for flow length calculations.
        LoD (float or None): Level of Detection for setting cell values to zero.

    This function reads the raster data, masks it with the polygon geometries, 
    calculates volumes based on raster values, and updates the shapefile with new attributes.
    """
    with rasterio.open(raster_path) as src:
        nodata = src.nodata  # Get the NoData value from the raster
        transform = src.transform  # Get the raster's transform

        # Create an empty list to store feature modifications
        new_features = []
        

        with fiona.open(investigation_polygons_path, 'r') as shapefile:
            # Define schema for output shapefile (including new fields)
            schema = shapefile.schema.copy()
            schema['properties']['dep_vol'] = 'float'
            schema['properties']['ero_vol'] = 'float'
            schema['properties']['net_mob'] = 'float'
            schema['properties']['tot_mob'] = 'float'
            schema['properties']['flow_len'] = 'float'
            schema['properties']['dep_YR'] = 'float'
            schema['properties']['ero_YR'] = 'float'
            schema['properties']['net_YR'] = 'float'

            for feature in shapefile:
                geom = shape(feature['geometry'])  # Convert shapefile feature to shapely geometry
                
                # Calculate the length of the polyline within the polygon
                flow_len = calculate_flow_length(debris_flow_path, geom)

                # Mask the data with the polygon geometry
                mask_data, mask_transform = mask(src, [geom], crop=True)
                mask_data = mask_data[0]  # Extract the single band from the masked data

                # Apply the LoD to set cell values below LoD to zero
                if LoD is not None:
                    mask_data[np.abs(mask_data) < LoD] = 0

                # Create a mask based on the NoData value
                data_mask = mask_data != nodata

                # Calculate cell areas, volumes, and sum of positive/negative volumes
                cell_areas = np.zeros(mask_data.shape)
                cell_volumes = np.zeros(mask_data.shape)
                dep_vol = 0.0  # Initialize sum of positive volumes
                ero_vol = 0.0  # Initialize sum of negative volumes
                for row in range(mask_data.shape[0]):
                    for col in range(mask_data.shape[1]):
                        if data_mask[row, col]:  # Check if the cell is not NoData
                            # Define the raster cell as a shapely box
                            cell_polygon = box(
                                mask_transform[2] + col * mask_transform[0],
                                mask_transform[5] + row * mask_transform[4],
                                mask_transform[2] + (col + 1) * mask_transform[0],
                                mask_transform[5] + (row + 1) * mask_transform[4]
                            )
                            # Calculate the intersection area between the cell and the polygon
                            intersection = cell_polygon.intersection(geom)
                            cell_areas[row, col] = intersection.area  # Store the intersection area

                            # Calculate the volume (value * area)
                            cell_volumes[row, col] = mask_data[row, col] * cell_areas[row, col]

                            # Accumulate positive and negative volumes
                            if cell_volumes[row, col] > 0:
                                dep_vol += cell_volumes[row, col]
                            elif cell_volumes[row, col] < 0:
                                ero_vol += cell_volumes[row, col]

                # Calculate net_mob and tot_mob
                net_mob = dep_vol + ero_vol
                tot_mob = np.abs(dep_vol) + np.abs(ero_vol)
                dep_YR = dep_vol / flow_len
                ero_YR = ero_vol / flow_len
                net_YR = net_mob / flow_len


                # Update feature attributes with dep_vol, ero_vol, net_mob, and tot_mob
                feature['properties']['dep_vol'] = dep_vol
                feature['properties']['ero_vol'] = ero_vol
                feature['properties']['net_mob'] = net_mob
                feature['properties']['tot_mob'] = tot_mob
                feature['properties']['flow_len'] = flow_len
                feature['properties']['dep_YR'] = dep_YR
                feature['properties']['ero_YR'] = ero_YR
                feature['properties']['net_YR'] = net_YR
                new_features.append(feature)

                # Print information for QC purposes
                print(f"Polygon {feature['id']}: Total Positive Volume (dep_vol): {dep_vol}")
                print(f"Polygon {feature['id']}: Total Negative Volume (ero_vol): {ero_vol}")
                print(f"Polygon {feature['id']}: Total Cell Area: {np.sum(cell_areas)}")
                print(f"Polygon {feature['id']}: Net Mobility (net_mob): {net_mob}")
                print(f"Polygon {feature['id']}: Total Mobility (tot_mob): {tot_mob}")
                for row in range(mask_data.shape[0]):
                    for col in range(mask_data.shape[1]):
                        if cell_areas[row, col] > 0:  # Check if the cell intersects the polygon
                            print(f"  Cell value: {mask_data[row, col]}, Cell area: {cell_areas[row, col]}, Cell volume: {cell_volumes[row, col]}")

        # Write updated features back to shapefile
        with fiona.open(investigation_polygons_path, 'w', crs=shapefile.crs, driver='ESRI Shapefile', schema=schema) as output:
            for feature in new_features:
                output.write(feature)



def plot_polygons(investigation_polygons_path, attribute, cmap='Blues'):
    """
    Plots the polygons from the shapefile with colors based on a specified attribute.

    Parameters:
        investigation_polygons_path (str): Path to the shapefile containing polygons.
        attribute (str): Attribute to visualize with color mapping.
        cmap (str): Colormap to use for visualization (default is 'Blues').

    This function creates a plot of the polygons, coloring them according to the specified attribute,
    providing a visual representation of the data.
    """
    # Open the shapefile to get the extent
    with fiona.open(investigation_polygons_path, 'r') as shapefile:
        minx, miny, maxx, maxy = shapefile.bounds

    # Calculate the width and height of the plot
    width = maxx - minx
    height = maxy - miny

    # Calculate the aspect ratio
    aspect_ratio = width / height
    
    # Create the plot with adjusted aspect ratio
    fig, ax = plt.subplots(figsize=(12, 12 / aspect_ratio))  # Adjust figsize based on aspect ratio

    # Open the shapefile
    with fiona.open(investigation_polygons_path, 'r') as shapefile:
        # Extract attribute values and geometries
        values = []
        geometries = []
        for feature in shapefile:
            geom = shape(feature['geometry'])
            geometries.append(geom)
            values.append(feature['properties'][attribute])
        
        # Normalize the values for color mapping
        vmin = min(values)
        vmax = max(values)
        norm = Normalize(vmin=vmin, vmax=vmax)
        cmap = plt.get_cmap(cmap)
        sm = ScalarMappable(norm=norm, cmap=cmap)
        
        for geom, value in zip(geometries, values):
            color = sm.to_rgba(value)
            poly = Polygon(list(geom.exterior.coords), facecolor=color, alpha=1)
            ax.add_patch(poly)
        
        # Set plot title and colorbar
        ax.set_title(f'{attribute} Map')
        sm.set_array(values)
        cbar = plt.colorbar(sm, ax=ax, orientation='vertical', label=f'{attribute} (m³)')  # Add the (m³) label here
        
        # Set limits and aspect ratio
        ax.set_xlim([minx, maxx])
        ax.set_ylim([miny, maxy])
        ax.set_aspect('equal')  # Set aspect ratio to be equal
        
        # Increase plot margins to accommodate larger polygons
        ax.margins(0.1)
        
        # Show the plot
        plt.show()

"""
The following functions all generate fairly simple plots to visualize incremental erosion and 
deposition volumes along the path of a debris flow using standard matplotlib plotting.
"""

def plot_attributes_vs_distance(investigation_polygons_path):
    dep_vol = []
    ero_vol = []
    distances = []
    
    # Open the shapefile
    with fiona.open(investigation_polygons_path, 'r') as shapefile:
        for feature in shapefile:
            # Extract feature properties
            properties = feature['properties']
            
            # Append dep_vol, ero_vol, and Distance to lists
            dep_vol.append(properties['dep_vol'])
            ero_vol.append(properties['ero_vol'])
            #ero_vol.append(abs(properties['ero_vol']))# Take absolute value of ero_vol
            distances.append(properties['Distance'])  # Assuming 'Distance' is the attribute name
            
    # Create plot for dep_vol and ero_vol vs Distance
    fig, ax = plt.subplots(figsize=(10, 6))

    ax.set_xlabel('Distance Upstream (m)')
    ax.set_ylabel('Volume (m³)')
    ax.scatter(distances, dep_vol, color='tab:blue', label='Depositional Volume', alpha=0.4)
    ax.scatter(distances, ero_vol, color='tab:red', label='Erosional Volume', alpha=0.4)
    ax.invert_xaxis()
    
    # Plot vertical lines at specified distances
    if polygon_ids is not None:
        vertical_lines, _ = get_dist_catch_from_ids(investigation_polygons_path, polygon_ids)
        for distance in vertical_lines:
            ax.axvline(x=distance, color='black', linestyle='--', linewidth=0.5)
    
    ax.legend()
    plt.title('Erosion and Deposition along Path by Distance')
    plt.show()

def plot_attributes_vs_catchment(investigation_polygons_path):
    dep_vol = []
    ero_vol = []
    catchments = []
    
    # Open the shapefile
    with fiona.open(investigation_polygons_path, 'r') as shapefile:
        for feature in shapefile:
            # Extract feature properties
            properties = feature['properties']
            
            # Append dep_vol, ero_vol, and Catch_m2 to lists
            dep_vol.append(properties['dep_vol'])
            ero_vol.append(properties['ero_vol'])
            #ero_vol.append(abs(properties['ero_vol']))# Take absolute value of ero_vol
            catchments.append(properties['Catch_m2'])  # Assuming 'Catch_m2' is the attribute name
    
    # Create plot for dep_vol and ero_vol vs Catchment Area
    fig, ax = plt.subplots(figsize=(10, 6))

    ax.set_xlabel('Catchment Area (m²)')
    ax.set_ylabel('Volume (m³)')
    ax.scatter(catchments, dep_vol, color='tab:blue', label='Depositional Volume', alpha=0.4)
    ax.scatter(catchments, ero_vol, color='tab:red', label='Erosional Volume', alpha=0.4)
    
    if polygon_ids is not None:
        _ , vertical_lines = get_dist_catch_from_ids(investigation_polygons_path, polygon_ids)
        for catchment in vertical_lines:
            ax.axvline(x=catchment, color='black', linestyle='--', linewidth=0.5)
    
    ax.legend()
    plt.title('Erosion and Deposition along Path by Catchment Area')
    plt.show()
    
def plot_attributes_vs_distance_bar(investigation_polygons_path):
    dep_vol = []
    ero_vol = []
    distances = []
    
    # Open the shapefile
    with fiona.open(investigation_polygons_path, 'r') as shapefile:
        for feature in shapefile:
            # Extract feature properties
            properties = feature['properties']
            
            # Append dep_vol, ero_vol, and Distance to lists
            dep_vol.append(properties['dep_vol'])
            ero_vol.append(properties['ero_vol'])
            distances.append(properties['Distance'])  
            
    # Create bar plot for dep_vol and ero_vol vs Distance
    fig, ax = plt.subplots(figsize=(10, 6))

    ax.set_xlabel('Distance Upstream (m)')
    ax.set_ylabel('Volume (m³)')
    bar_width = (max(distances) - min(distances)) / len(distances) * 4  # Set bar width
    ax.bar(distances, dep_vol, width=bar_width, color='tab:blue', label='Depositional Volume', alpha=0.4)
    ax.bar(distances, ero_vol, width=bar_width, color='tab:red', label='Erosional Volume', alpha=0.4)
    ax.invert_xaxis()
    
    if polygon_ids is not None:
        vertical_lines, _ = get_dist_catch_from_ids(investigation_polygons_path, polygon_ids)
        for distance in vertical_lines:
            ax.axvline(x=distance, color='black', linestyle='--', linewidth=0.5)
        
    ax.legend()
    plt.title('Depositional and Erosional Volume vs Distance')
    plt.show()

def plot_attributes_vs_catchment_bar(investigation_polygons_path):
    dep_vol = []
    ero_vol = []
    catchments = []
    
    # Open the shapefile
    with fiona.open(investigation_polygons_path, 'r') as shapefile:
        for feature in shapefile:
            # Extract feature properties
            properties = feature['properties']
            
            # Append dep_vol, ero_vol, and Catch_m2 to lists
            dep_vol.append(properties['dep_vol'])
            ero_vol.append(properties['ero_vol'])
            catchments.append(properties['Catch_m2'])  # Assuming 'Catch_m2' is the attribute name
    
    # Create bar plot for dep_vol and ero_vol vs Catchment Area
    fig, ax = plt.subplots(figsize=(10, 6))

    ax.set_xlabel('Catchment Area (m²)')
    ax.set_ylabel('Volume (m³)')
    bar_width = (max(catchments) - min(catchments)) / len(catchments) * 10  # Set bar width
    ax.bar(catchments, dep_vol, width=bar_width, color='tab:blue', label='Depositional Volume', alpha=0.2)
    ax.bar(catchments, ero_vol, width=bar_width, color='tab:red', label='Erosional Volume', alpha=0.2)
    
    if polygon_ids is not None:
        _ , vertical_lines = get_dist_catch_from_ids(investigation_polygons_path, polygon_ids)
        for catchment in vertical_lines:
            ax.axvline(x=catchment, color='black', linestyle='--', linewidth=0.5)
    
    ax.legend()
    plt.title('Depositional and Erosional Volume vs Catchment Area')
    plt.show()
    
def plot_attributes_vs_catchment_log_bar(investigation_polygons_path):
    dep_vol = []
    ero_vol = []
    catchments = []
    
    # Open the shapefile
    with fiona.open(investigation_polygons_path, 'r') as shapefile:
        for feature in shapefile:
            # Extract feature properties
            properties = feature['properties']
            
            # Append dep_vol, ero_vol, and Catch_m2 to lists
            dep_vol.append(properties['dep_vol'])
            ero_vol.append(properties['ero_vol'])
            catchments.append(properties['Catch_m2'])  # Assuming 'Catch_m2' is the attribute name
    
    # Create bar plot for dep_vol and ero_vol vs Catchment Area
    fig, ax = plt.subplots(figsize=(10, 6))
    
    ax.set_xlabel('Catchment Area (m²)')
    ax.set_ylabel('Volume (m³)')
    ax.set_xscale('log')  # Set x-axis to logarithmic scale
    
    # Define bar width on logarithmic scale
    log_width = 0.02  # Adjust this value to change bar width
    
    # Plot the bars using plt.fill
    for i in range(len(catchments)):
        ax.fill(
            [10**(np.log10(catchments[i])-log_width), 10**(np.log10(catchments[i])-log_width),
             10**(np.log10(catchments[i])+log_width), 10**(np.log10(catchments[i])+log_width)],
            [0, dep_vol[i], dep_vol[i], 0],
            'tab:blue', alpha=0.2, label='Depositional Volume' if i == 0 else ""
        )
        ax.fill(
            [10**(np.log10(catchments[i])-log_width), 10**(np.log10(catchments[i])-log_width),
             10**(np.log10(catchments[i])+log_width), 10**(np.log10(catchments[i])+log_width)],
            [0, ero_vol[i], ero_vol[i], 0],
            'tab:red', alpha=0.2, label='Erosional Volume' if i == 0 else ""
        )
        
    if polygon_ids is not None:    
        _ , vertical_lines = get_dist_catch_from_ids(investigation_polygons_path, polygon_ids)
        for catchment in vertical_lines:
            ax.axvline(x=catchment, color='black', linestyle='--', linewidth=0.5)
    
    ax.legend()
    plt.title('Depositional and Erosional Volume vs Catchment Area')
    plt.show()
    
def plot_cumulative_vs_distance(investigation_polygons_path):
    dep_vol_cumulative = []
    ero_vol_cumulative = []
    distances = []
    
    # Open the shapefile and read features
    with fiona.open(investigation_polygons_path, 'r') as shapefile:
        features = list(shapefile)
        # Sort features by Distance
        features_sorted = sorted(features, key=lambda x: x['properties']['Distance'], reverse=True)
        
        cumulative_dep_vol = 0
        cumulative_ero_vol = 0
        
        for feature in features_sorted:
            # Extract feature properties
            properties = feature['properties']
            
            # Update cumulative dep_vol and ero_vol
            cumulative_dep_vol += properties['dep_vol']
            cumulative_ero_vol += abs(properties['ero_vol'])  # Take absolute value of ero_vol
            
            # Append cumulative values and Distance to lists
            dep_vol_cumulative.append(cumulative_dep_vol)
            ero_vol_cumulative.append(cumulative_ero_vol)
            distances.append(properties['Distance'])  # Assuming 'Distance' is the attribute name
            
    # Create plot for cumulative dep_vol and ero_vol vs Distance
    fig, ax = plt.subplots(figsize=(10, 6))

    ax.set_xlabel('Distance Upstream (m)')
    ax.set_ylabel('Cumulative Volume (m³)')
    ax.plot(distances, dep_vol_cumulative, color='tab:blue', label='Cumulative Depositional Volume')
    ax.plot(distances, ero_vol_cumulative, color='tab:red', label='Cumulative Erosional Volume (absolute)')
    ax.invert_xaxis()
    
    if polygon_ids is not None:
        vertical_lines, _ = get_dist_catch_from_ids(investigation_polygons_path, polygon_ids)
        for distance in vertical_lines:
            ax.axvline(x=distance, color='black', linestyle='--', linewidth=0.5)
    
    ax.legend()
    plt.title('Cumulative Deposition and Erosion Volumes vs Distance')
    plt.show()

def plot_cumulative_vs_catchment(investigation_polygons_path):
    dep_vol_cumulative = []
    ero_vol_cumulative = []
    catchments = []
    
    # Open the shapefile and read features
    with fiona.open(investigation_polygons_path, 'r') as shapefile:
        features = list(shapefile)
        # Sort features by Catch_m2
        features_sorted = sorted(features, key=lambda x: x['properties']['Catch_m2'])
        
        cumulative_dep_vol = 0
        cumulative_ero_vol = 0
        
        for feature in features_sorted:
            # Extract feature properties
            properties = feature['properties']
            
            # Update cumulative dep_vol and ero_vol
            cumulative_dep_vol += properties['dep_vol']
            cumulative_ero_vol += abs(properties['ero_vol'])  # Take absolute value of ero_vol
            
            # Append cumulative values and Catch_m2 to lists
            dep_vol_cumulative.append(cumulative_dep_vol)
            ero_vol_cumulative.append(cumulative_ero_vol)
            catchments.append(properties['Catch_m2'])  # Assuming 'Catch_m2' is the attribute name
    
    # Create plot for cumulative dep_vol and ero_vol vs Catchment Area
    fig, ax = plt.subplots(figsize=(10, 6))

    ax.set_xlabel('Catchment Area (m²)')
    ax.set_ylabel('Cumulative Volume (m³)')
    ax.plot(catchments, dep_vol_cumulative, color='tab:blue', label='Cumulative Depositional Volume')
    ax.plot(catchments, ero_vol_cumulative, color='tab:red', label='Cumulative Erosional Volume (absolute)')
    
    
    if polygon_ids is not None:
        _ , vertical_lines = get_dist_catch_from_ids(investigation_polygons_path, polygon_ids)
        for catchment in vertical_lines:
            ax.axvline(x=catchment, color='black', linestyle='--', linewidth=0.5)
    
    ax.legend()
    plt.title('Cumulative Deposition and Erosion Volumes vs Catchment Area')
    plt.show()

def plot_cumulative_vs_distance_scatter(investigation_polygons_path):
    dep_vol_cumulative = []
    ero_vol_cumulative = []
    distances = []
    
    # Open the shapefile and read features
    with fiona.open(investigation_polygons_path, 'r') as shapefile:
        features = list(shapefile)
        # Sort features by Distance
        features_sorted = sorted(features, key=lambda x: x['properties']['Distance'], reverse=True)
        
        cumulative_dep_vol = 0
        cumulative_ero_vol = 0
        
        for feature in features_sorted:
            # Extract feature properties
            properties = feature['properties']
            
            # Update cumulative dep_vol and ero_vol
            cumulative_dep_vol += properties['dep_vol']
            cumulative_ero_vol += abs(properties['ero_vol'])  # Take absolute value of ero_vol
            
            # Append cumulative values and Distance to lists
            dep_vol_cumulative.append(cumulative_dep_vol)
            ero_vol_cumulative.append(cumulative_ero_vol)
            distances.append(properties['Distance'])  # Assuming 'Distance' is the attribute name
            
            
    # Create scatter plot for cumulative dep_vol and ero_vol vs Distance
    fig, ax = plt.subplots(figsize=(10, 6))

    ax.set_xlabel('Distance Upstream (m)')
    ax.set_ylabel('Cumulative Volume (m³)')
    ax.scatter(distances, dep_vol_cumulative, color='tab:blue', label='Cumulative Depositional Volume')
    ax.scatter(distances, ero_vol_cumulative, color='tab:red', label='Cumulative Erosional Volume (absolute)')
    ax.invert_xaxis()
    
    if polygon_ids is not None:
        vertical_lines, _ = get_dist_catch_from_ids(investigation_polygons_path, polygon_ids)
        for distance in vertical_lines:
            ax.axvline(x=distance, color='black', linestyle='--', linewidth=0.5)
        
    
    ax.legend()
    plt.title('Cumulative dep_vol and abs(ero_vol) vs Distance (Scatter)')
    plt.show()

def plot_cumulative_vs_catchment_scatter(investigation_polygons_path):
    dep_vol_cumulative = []
    ero_vol_cumulative = []
    catchments = []
    
    # Open the shapefile and read features
    with fiona.open(investigation_polygons_path, 'r') as shapefile:
        features = list(shapefile)
        # Sort features by Catch_m2
        features_sorted = sorted(features, key=lambda x: x['properties']['Catch_m2'])
        
        cumulative_dep_vol = 0
        cumulative_ero_vol = 0
        
        for feature in features_sorted:
            # Extract feature properties
            properties = feature['properties']
            
            # Update cumulative dep_vol and ero_vol
            cumulative_dep_vol += properties['dep_vol']
            cumulative_ero_vol += abs(properties['ero_vol'])  # Take absolute value of ero_vol
            
            # Append cumulative values and Catch_m2 to lists
            dep_vol_cumulative.append(cumulative_dep_vol)
            ero_vol_cumulative.append(cumulative_ero_vol)
            catchments.append(properties['Catch_m2'])  # Assuming 'Catch_m2' is the attribute name
    
    # Create scatter plot for cumulative dep_vol and ero_vol vs Catchment Area
    fig, ax = plt.subplots(figsize=(10, 6))

    ax.set_xlabel('Catchment Area (m²)')
    ax.set_ylabel('Cumulative Volume (m³)')
    ax.scatter(catchments, dep_vol_cumulative, color='tab:blue', label='Cumulative Depositional Volume')
    ax.scatter(catchments, ero_vol_cumulative, color='tab:red', label='Cumulative Erosional Volume (absolute)')
    
    if polygon_ids is not None:
        _ , vertical_lines = get_dist_catch_from_ids(investigation_polygons_path, polygon_ids)
        for catchment in vertical_lines:
            ax.axvline(x=catchment, color='black', linestyle='--', linewidth=0.5)
    
    ax.legend()
    plt.title('Cumulative dep_vol and abs(ero_vol) vs Catchment Area (Scatter)')
    plt.show()

# Function to plot polygons color-coded by attribute using fiona and matplotlib
def plot_attributes_vs_catchment_log(investigation_polygons_path):
    dep_vol = []
    ero_vol = []
    catchments = []
    
    # Open the shapefile
    with fiona.open(investigation_polygons_path, 'r') as shapefile:
        for feature in shapefile:
            # Extract feature properties
            properties = feature['properties']
            
            # Append dep_vol, ero_vol, and Catch_m2 to lists
            dep_vol.append(properties['dep_vol'])
            ero_vol.append(properties['ero_vol'])
            #ero_vol.append(abs(properties['ero_vol']))# Take absolute value of ero_vol
            catchments.append(properties['Catch_m2'])  # Assuming 'Catch_m2' is the attribute name
    
    # Create plot for dep_vol and ero_vol vs Catchment Area
    fig, ax = plt.subplots(figsize=(10, 6))

    ax.set_xlabel('Catchment Area (m²)')
    ax.set_ylabel('Volume (m³)')
    ax.scatter(catchments, dep_vol, color='tab:blue', label='Cumulative Depositional Volume', alpha=0.4)
    ax.scatter(catchments, ero_vol, color='tab:red', label='Cumulative Erosional Volume (absolute)', alpha=0.4)
    
    ax.set_xscale('log')  # Set x-axis to logarithmic scale
    
    if polygon_ids is not None:
        _ , vertical_lines = get_dist_catch_from_ids(investigation_polygons_path, polygon_ids)
        for catchment in vertical_lines:
            ax.axvline(x=catchment, color='black', linestyle='--', linewidth=0.5)
    
    ax.legend()
    plt.title('dep_vol and abs(ero_vol) vs Catchment Area (Log Scale)')
    plt.show()

def plot_cumulative_vs_catchment_log(investigation_polygons_path):
    dep_vol_cumulative = []
    ero_vol_cumulative = []
    catchments = []
    
    # Open the shapefile and read features
    with fiona.open(investigation_polygons_path, 'r') as shapefile:
        features = list(shapefile)
        # Sort features by Catch_m2
        features_sorted = sorted(features, key=lambda x: x['properties']['Catch_m2'])
        
        cumulative_dep_vol = 0
        cumulative_ero_vol = 0
        
        for feature in features_sorted:
            # Extract feature properties
            properties = feature['properties']
            
            # Update cumulative dep_vol and ero_vol
            cumulative_dep_vol += properties['dep_vol']
            cumulative_ero_vol += abs(properties['ero_vol'])  # Take absolute value of ero_vol
            
            # Append cumulative values and Catch_m2 to lists
            dep_vol_cumulative.append(cumulative_dep_vol)
            ero_vol_cumulative.append(cumulative_ero_vol)
            catchments.append(properties['Catch_m2'])  # Assuming 'Catch_m2' is the attribute name
    
    # Create plot for cumulative dep_vol and ero_vol vs Catchment Area
    fig, ax = plt.subplots(figsize=(10, 6))

    ax.set_xlabel('Catchment Area (m²)')
    ax.set_ylabel('Cumulative Volume (m³)')
    ax.plot(catchments, dep_vol_cumulative, color='tab:blue', label='Cumulative Depositional Volume')
    ax.plot(catchments, ero_vol_cumulative, color='tab:red', label='Cumulative Erosional Volume (absolute)')
    
    ax.set_xscale('log')  # Set x-axis to logarithmic scale
    
    if polygon_ids is not None:
        _ , vertical_lines = get_dist_catch_from_ids(investigation_polygons_path, polygon_ids)
        for catchment in vertical_lines:
            ax.axvline(x=catchment, color='black', linestyle='--', linewidth=0.5)
    
    ax.legend()
    plt.title('Cumulative dep_vol and abs(ero_vol) vs Catchment Area (Log Scale)')
    plt.show()

def plot_cumulative_vs_catchment_scatter_log(investigation_polygons_path):
    dep_vol_cumulative = []
    ero_vol_cumulative = []
    catchments = []
    
    # Open the shapefile and read features
    with fiona.open(investigation_polygons_path, 'r') as shapefile:
        features = list(shapefile)
        # Sort features by Catch_m2
        features_sorted = sorted(features, key=lambda x: x['properties']['Catch_m2'])
        
        cumulative_dep_vol = 0
        cumulative_ero_vol = 0
        
        for feature in features_sorted:
            # Extract feature properties
            properties = feature['properties']
            
            # Update cumulative dep_vol and ero_vol
            cumulative_dep_vol += properties['dep_vol']
            cumulative_ero_vol += abs(properties['ero_vol'])  # Take absolute value of ero_vol
            
            # Append cumulative values and Catch_m2 to lists
            dep_vol_cumulative.append(cumulative_dep_vol)
            ero_vol_cumulative.append(cumulative_ero_vol)
            catchments.append(properties['Catch_m2'])  # Assuming 'Catch_m2' is the attribute name
    
    # Create scatter plot for cumulative dep_vol and ero_vol vs Catchment Area
    fig, ax = plt.subplots(figsize=(10, 6))

    ax.set_xlabel('Catchment Area (m²)')
    ax.set_ylabel('Cumulative Volume (m³)')
    ax.scatter(catchments, dep_vol_cumulative, color='tab:blue', label='Cumulative Depositional Volume')
    ax.scatter(catchments, ero_vol_cumulative, color='tab:red', label='Cumulative Erosional Volume (absolute)')
    
    if polygon_ids is not None:   
        _ , vertical_lines = get_dist_catch_from_ids(investigation_polygons_path, polygon_ids)
        for catchment in vertical_lines:
            ax.axvline(x=catchment, color='black', linestyle='--', linewidth=0.5)
    
    ax.set_xscale('log')  # Set x-axis to logarithmic scale
    
    ax.legend()
    plt.title('Cumulative dep_vol and abs(ero_vol) vs Catchment Area (Scatter)')
    plt.show()


def plot_YR_vs_distance(investigation_polygons_path):
    dep_YR = []
    ero_YR = []
    distances = []
    
    # Open the shapefile
    with fiona.open(investigation_polygons_path, 'r') as shapefile:
        features = list(shapefile)
        features_sorted = sorted(features, key=lambda x: x['properties']['Distance'])
        for feature in features_sorted:
            # Extract feature properties
            properties = feature['properties']
            
            # Append dep_vol, ero_vol, and Distance to lists
            dep_YR.append(properties['dep_YR'])
            ero_YR.append(properties['ero_YR'])
            #ero_vol.append(abs(properties['ero_vol']))# Take absolute value of ero_vol
            distances.append(properties['Distance'])  # Assuming 'Distance' is the attribute name
            
    # Create plot for dep_vol and ero_vol vs Distance
    fig, ax = plt.subplots(figsize=(10, 6))

    ax.set_xlabel('Distance Upstream (m)')
    ax.set_ylabel('YR (m³/m)')
    ax.plot(distances, dep_YR, color='tab:blue', label='Depositional YR', alpha=0.4)
    ax.plot(distances, ero_YR, color='tab:red', label='Erosional YR', alpha=0.4)
    ax.invert_xaxis()
    
    # Plot vertical lines at specified distances
    if polygon_ids is not None:
        vertical_lines, _ = get_dist_catch_from_ids(investigation_polygons_path, polygon_ids)
        for distance in vertical_lines:
            ax.axvline(x=distance, color='black', linestyle='--', linewidth=0.5)
    
    ax.legend()
    plt.title('YR vs Distance')
    plt.show()

def plot_YR_vs_catchment(investigation_polygons_path):
    dep_YR = []
    ero_YR = []
    catchments = []
    
    # Open the shapefile
    with fiona.open(investigation_polygons_path, 'r') as shapefile:
        features = list(shapefile)
        features_sorted = sorted(features, key=lambda x: x['properties']['Catch_m2'])
        for feature in features_sorted:
            # Extract feature properties
            properties = feature['properties']
            
            
            # Append dep_vol, ero_vol, and Catch_m2 to lists
            dep_YR.append(properties['dep_YR'])
            ero_YR.append(properties['ero_YR'])
            #ero_vol.append(abs(properties['ero_vol']))# Take absolute value of ero_vol
            catchments.append(properties['Catch_m2'])  # Assuming 'Catch_m2' is the attribute name
    
    # Create plot for dep_vol and ero_vol vs Catchment Area
    fig, ax = plt.subplots(figsize=(10, 6))

    ax.set_xlabel('Catchment Area (m²)')
    ax.set_ylabel('Volume (m³/M)')
    ax.plot(catchments, dep_YR, color='tab:blue', label='Depositional YR', alpha=0.4)
    ax.plot(catchments, ero_YR, color='tab:red', label='Erosional YR', alpha=0.4)
    
    if polygon_ids is not None:
        _ , vertical_lines = get_dist_catch_from_ids(investigation_polygons_path, polygon_ids)
        for catchment in vertical_lines:
            ax.axvline(x=catchment, color='black', linestyle='--', linewidth=0.5)
    
    ax.legend()
    plt.title('YR vs Catchment Area')
    plt.show()

def get_dist_catch_from_ids(investigation_polygons_path, polygon_ids):
    """
    Extracts distances and catchment areas from specified investigation polygons 
    based on their IDs.

    Parameters:
        investigation_polygons_path (str): Path to the shapefile containing investigation polygons.
        polygon_ids (list): List of IDs corresponding to the polygons of interest.

    Returns:
        tuple: A tuple containing two lists:
            - distances (list): A list of distances associated with the specified polygon IDs.
            - catchments (list): A list of catchment areas (in square meters) 
              associated with the specified polygon IDs.
    """
    distances = []
    catchments = []
    with fiona.open(investigation_polygons_path, 'r') as shapefile:
        for feature in shapefile:
            if feature['id'] in polygon_ids:
                distances.append(feature['properties']['Distance'])
                catchments.append(feature['properties']['Catch_m2'])
    return distances, catchments
    

raster_path = "C:/Users/path/to/change/det/raster.tif"
investigation_polygons_path = "C:/Users/path/to/investigation/polygons.shp"
debris_flow_path = "C:/Users/path/to/flowpath.shp"
# Define locations of tributaries based on ID
polygon_ids = ['803','830','766'] or none # Defaults to none
# Define the Limit of Detection (LoD)
LoD = 0.25  # Example value, set as needed


extract_raster_values(raster_path, investigation_polygons_path, debris_flow_path, LoD)
plot_polygons(investigation_polygons_path, 'dep_vol', cmap='Blues')
plot_polygons(investigation_polygons_path, 'ero_vol', cmap='Reds_r')
plot_polygons(investigation_polygons_path, 'net_mob', cmap='Purples_r')
plot_attributes_vs_distance_bar(investigation_polygons_path)
plot_YR_vs_distance(investigation_polygons_path)
plot_YR_vs_catchment(investigation_polygons_path)
plot_attributes_vs_catchment_bar(investigation_polygons_path)
plot_attributes_vs_catchment_log_bar(investigation_polygons_path)
plot_cumulative_vs_distance(investigation_polygons_path)
plot_cumulative_vs_catchment(investigation_polygons_path)
plot_cumulative_vs_distance_scatter(investigation_polygons_path)
plot_cumulative_vs_catchment_scatter(investigation_polygons_path)
