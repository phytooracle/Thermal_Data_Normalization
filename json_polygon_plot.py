#!/usr/bin/env python3
"""
Author : Sebastian Calleja
Date   : 2021-03-17
Purpose: Rock the Casbah
"""

import argparse
import os
import sys
from osgeo import gdal, osr
from osgeo import osr
import geopandas as gpd
from osgeo import ogr
from shapely.geometry import Polygon, Point
from geopandas.geoseries import *
import re
import tarfile
import pandas as pd
import datetime as dt
import os
import glob
import subprocess
import os.path
from os import path
import json
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import utm
# from terrautils.spatial import scanalyzer_to_latlon


# --------------------------------------------------
def get_args():
    """Get command-line arguments"""

    parser = argparse.ArgumentParser(
        description='Associate image metadata with agriculutural plots',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # parser.add_argument('positional',
    #                     metavar='str',
    #                     help='A positional argument')

    parser.add_argument('-d',
                        '--date',
                        help='Date to process',
                        metavar='str',
                        type=str,
                        default=None,
                        required=True)

    parser.add_argument('-s',
                        '--season',
                        help='Season to which the scan date belongs to',
                        metavar='str',
                        type=str,
                        default=None,
                        required=True)


    parser.add_argument('-g',
                        '--geojson',
                        help='Plot polygon GeoJSON',
                        metavar='str',
                        type=str,
                        required=True)


    parser.add_argument('-o',
                        '--output',
                        help='Output file path',
                        metavar='str',
                        type=str,
                        default='image_plot_association')

    return parser.parse_args()


# Grabs the thermal raw data tar file for the Gantry's metadata for the images taken at the specified date
class JSON:
    # def get_tar(season, date):
    #     command = f'iget -rKTPf -N 0 /iplant/home/shared/phytooracle/{season}/level_0/flirIrCamera/flirIrCamera-{date}.tar'
    #     subprocess.call(command, shell = True)
    #     command = f'tar -xvf flirIrCamera-{date}.tar'
    #     subprocess.call(command, shell = True)

# Finds the individual json files and adds them to the filepath (end up with a list of paths for each individual json file)
    def pathlist(season, date):
        # json_data = JSON.get_tar(season, date)
        pathlist = Path(f"./flirIrCamera/{date}/").glob('**/*.json')
        JSON_path_list = []
        for path in pathlist:
            path_str = str(path)
            JSON_path_list.append(path_str)
        return JSON_path_list

# Uses the pathlists and searches each one
# Gathers the time, image name, gantry_x and gantry_y information from the metadata of each individual image and adds to dictionary
    def time_dict(season, date):
        file_path_list = JSON.pathlist(season, date)
        JSON_dict = dict()
        for file in file_path_list:
            path_metadata = glob.glob(f'{file}')
            metadata = str(path_metadata)[2:-2]
            with open(metadata) as f:
                meta = json.load(f)['lemnatec_measurement_metadata']
                time = (meta['gantry_system_variable_metadata']['time'])
                gantry_x = float(meta['gantry_system_variable_metadata']['position x [m]'])
                gantry_y = float(meta['gantry_system_variable_metadata']['position y [m]'])
                gantry_z = float(meta['gantry_system_variable_metadata']['position z [m]'])
                filename = os.path.basename(metadata)
            if JSON is not JSON_dict:
                JSON_dict[time, filename, gantry_x, gantry_y, gantry_z] = "Date, Time, Gantry_x, Gantry_y, Gantry_z"
            else:
                print("JSON already in Dictionary")
        return sorted(JSON_dict)

# Searches through the dictionary created and creates a dataframe of the information in the dictionary
    def time_df(season, date):
        JSON_time_d = JSON.time_dict(season, date)
        JSON_time_df = pd.DataFrame.from_dict(JSON_time_d)
        JSON_time_df.columns = ['Date and Time', 'Image Name', 'Gantry_x', 'Gantry_y', "Gantry_z"]
        return JSON_time_df

# Converts the gantry coordinates into GPS coordinates (using 'terrautils')
# Used when we defined an image as a Point
    def GPS_coord (season, date):
        data = JSON.time_df(season, date)
        gantry_x_pos = data['Gantry_x']
        gantry_y_pos = data['Gantry_y']        
        GPS_latlon = scanalyzer_to_latlon(gantry_x_pos, gantry_y_pos)
        GPS_df = pd.DataFrame(GPS_latlon)
        GPS_latlon_df = GPS_df.transpose()
        GPS_latlon_df.columns = ['GPS_lat', 'GPS_lon']
        data['GPS_lat'] = GPS_latlon_df['GPS_lat']
        data['GPS_lon'] = GPS_latlon_df['GPS_lon']
        return data
    
# Takes Metadata for the image and creates a bounding box
# Used when we defined image as a polygon rather than using the center point
    def b_box(season, date):
        file_path_list = JSON.pathlist(season, date)
        JSON_dict = dict()
        for file in file_path_list:
            path_metadata = glob.glob(f'{file}')
            metadata = str(path_metadata)[2:-2]
            with open(metadata) as f:
                meta = json.load(f)['lemnatec_measurement_metadata']
                time = (meta['gantry_system_variable_metadata']['time'])
                loc_gantry_x = float(meta['sensor_fixed_metadata']['location in camera box x [m]'])
                loc_gantry_y = float(meta['sensor_fixed_metadata']['location in camera box y [m]'])
                loc_gantry_z = float(meta['sensor_fixed_metadata']['location in camera box z [m]'])
                z_offset = 0.76
                gantry_x = float(meta['gantry_system_variable_metadata']['position x [m]']) + loc_gantry_x
                gantry_y = float(meta['gantry_system_variable_metadata']['position y [m]']) + loc_gantry_y
                gantry_z = float(meta['gantry_system_variable_metadata']['position z [m]']) + z_offset + loc_gantry_z #offset in m
                fov_x, fov_y = float(meta['sensor_fixed_metadata']['field of view x [m]']), float(meta['sensor_fixed_metadata']['field of view y [m]'])
                B = gantry_z
                A_x = np.arctan((0.5*float(fov_x))/2)
                A_y = np.arctan((0.5*float(fov_y))/2)
                L_x = 2*B*np.tan(A_x)
                L_y = 2*B*np.tan(A_y)
                x_n = gantry_x + (L_x/2)
                x_s = gantry_x - (L_x/2)
                y_w = gantry_y + (L_y/2)
                y_e = gantry_y - (L_y/2)
                bbox_nw_latlon = scanalyzer_to_latlon(x_n, y_w)
                bbox_se_latlon = scanalyzer_to_latlon(x_s, y_e)

                # TERRA-REF
                lon_shift = 0.000020308287

                # Drone
                lat_shift = 0.000018292 #0.000015258894
                b_box =  ( bbox_se_latlon[0] - lat_shift,
                            bbox_nw_latlon[0] - lat_shift,
                            bbox_nw_latlon[1] + lon_shift,
                            bbox_se_latlon[1] + lon_shift)
                filename = os.path.basename(metadata)
            if JSON is not JSON_dict:
                JSON_dict[time, filename, gantry_x, gantry_y, gantry_z, b_box] = "Date, Time, Gantry_x, Gantry_y, Gantry_z, b_box"
            else:
                print("JSON already in Dictionary")
        return sorted(JSON_dict)
    
    def bbox_df(season, date):
        image_bbox = JSON.b_box(season, date)
        bbox_df = pd.DataFrame.from_dict(image_bbox)
        bbox_df.columns = ['Date and Time', 'Image Name', 'Gantry_x', 'Gantry_y', "Gantry_z", "b_box"]
        return bbox_df
    

# Create polygons using the two image points
# creates df from list returned above and adds it to original df with all the information
polygon_list = []
class Get_poly:
    def to_poly (image_bbox):
        for i, row in image_bbox.iterrows():
            bbox = image_bbox['b_box'].loc[i]
            polygon = Polygon([[bbox[2], bbox[1]], [bbox[3], bbox[1]], [bbox[3], bbox[0]], [bbox[2], bbox[0]]])
            polygon_list.append(polygon)
            polygon_df = pd.DataFrame(polygon_list)
            polygon_df.rename(columns={0: 'bbox_geometry'}, inplace=True)
            image_bbox['bbox_geometry'] = polygon_df
        return image_bbox
    
# Takes bbox polygon and iterates it through plot polygons to see where they intersect
def intersection(bbox_polygon, shp):
    intersects = bbox_polygon.intersects
    plot = None
    intersection_list = []
    for i, row in shp.iterrows():
        plot_polygon = row['geometry']
        intersection = intersects(plot_polygon)
        if intersection == True:
            plot = [row['ID']]
            intersection_list.append(plot)
    return intersection_list


def scanalyzer_to_latlon(gantry_x, gantry_y):
    """Convert coordinates from gantry to lat/lon"""
    utm_x, utm_y = scanalyzer_to_utm(gantry_x, gantry_y)
    return utm_to_latlon(utm_x, utm_y)


def scanalyzer_to_utm(gantry_x, gantry_y):
    """Convert coordinates from gantry to UTM 12N"""

    # TODO: Hard-coded
    # Linear transformation coefficients
    ay = 3659974.971; by = 1.0002; cy = 0.0078;
    ax = 409012.2032; bx = 0.009; cx = - 0.9986;

    utm_x = ax + (bx * gantry_x) + (cx * gantry_y)
    utm_y = ay + (by * gantry_x) + (cy * gantry_y)

    return utm_x, utm_y

def utm_to_latlon(utm_x, utm_y):
    """Convert coordinates from UTM 12N to lat/lon"""

    # Get UTM information from southeast corner of field
    SE_utm = utm.from_latlon(33.07451869, -111.97477775)
    utm_zone = SE_utm[2]
    utm_num  = SE_utm[3]

    return utm.to_latlon(utm_x, utm_y, utm_zone, utm_num)

# --------------------------------------------------
def main():
    """Make a jazz noise here"""

    args = get_args()
    image_bbox = JSON.bbox_df(args.season, args.date)

    # GeoJson file with lat and lon information
    shp = gpd.read_file(args.geojson)

    Get_poly.to_poly(image_bbox)

    # iterates through all of the bbox polygons and feeds them into the intersection function (Takes ~ 35 min)
    image_bbox["plot"] = None
    for i, row in image_bbox.iterrows():
        bbox_polygon = row['bbox_geometry']
        plot = intersection(bbox_polygon, shp)
        image_bbox.at[i,'plot'] = plot
        
    image_bbox.to_csv(f'{args.output}.csv')


# --------------------------------------------------
if __name__ == '__main__':
    main()
