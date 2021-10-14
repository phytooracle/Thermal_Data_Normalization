#!/usr/bin/env python3

import geopandas as gpd
from shapely.geometry import Polygon, Point
from geopandas.geoseries import *
import tarfile
import pandas as pd
import datetime as dt
import os
import glob
import subprocess
import os.path
import json
from pathlib import Path
import numpy as np
import urllib
import argparse
import utm
#----------------------------------

def get_args():
    """Get command-line arguments"""

    parser = argparse.ArgumentParser(
        description='Individual plant temperature extraction',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('dir',
                        metavar='dir',
                        help='Directory containing geoTIFFs')

    parser.add_argument('-g',
                        '--geojson',
                        help='GeoJSON containing plot boundaries',
                        metavar='geojson',
                        type=str,
                        default=None,
                        required=True)

    parser.add_argument('-s',
                        '--season',
                        help=' "season_10_lettuce_yr_2020" OR "season_11_sorghum_yr_2020" ',
                        metavar='season',
                        type=str,
                        default=None,
                        required=True)

    parser.add_argument('-d',
                        '--date',
                        help='date to process, usually the same as dir',
                        metavar='date',
                        type=str,
                        default=None,
                        required=True)

    parser.add_argument('-o',
                        '--outdir',
                        help='Output directory where resulting csv will be saved',
                        metavar='str',
                        type=str,
                        default='Thermal_Output')

    parser.add_argument('-c',
                        '--csv',
                        help=' "Individual Plant Temp CSV" ',
                        metavar='csv',
                        type=str,
                        default=None,
                        required=True)

    return parser.parse_args()

#----------------------------------
# UTM Functions
def utm_to_latlon(utm_x, utm_y):
    """Convert coordinates from UTM 12N to lat/lon"""

    # Get UTM information from southeast corner of field
    SE_utm = utm.from_latlon(33.07451869, -111.97477775)
    utm_zone = SE_utm[2]
    utm_num  = SE_utm[3]
    return utm.to_latlon(utm_x, utm_y, utm_zone, utm_num)

def scanalyzer_to_utm(gantry_x, gantry_y):
    """Convert coordinates from gantry to UTM 12N"""

    # TODO: Hard-coded
    # Linear transformation coefficients
    ay = 3659974.971; by = 1.0002; cy = 0.0078;
    ax = 409012.2032; bx = 0.009; cx = - 0.9986;

    utm_x = ax + (bx * gantry_x) + (cx * gantry_y)
    utm_y = ay + (by * gantry_x) + (cy * gantry_y)

    return utm_x, utm_y

def scanalyzer_to_latlon(gantry_x, gantry_y):
    """Convert coordinates from gantry to lat/lon"""
    utm_x, utm_y = scanalyzer_to_utm(gantry_x, gantry_y)
    return utm_to_latlon(utm_x, utm_y)


#----------------------------------
# Gathers gantry x, y, z, coordinates as well as time
# Converts coordinates to lat lon
def md_shp():
    
    args = get_args()

    # Load required files as well as files to process
    pathlist = glob.glob(f'{args.dir}/*/*.json')
    
    shp = gpd.read_file(f'{args.geojson}')

    JSON_path_list = []
    for path in pathlist:
        path_str = str(path)
        JSON_path_list.append(path_str)

    # Create dictionary and populates it
    JSON_dict = {}
    cnt = 0
    # JSON_dict[time, filename, gantry_x, gantry_y, gantry_z] = "Date, Time, Gantry_x, Gantry_y, Gantry_z"
    for i in JSON_path_list:
        with open(i) as f:
            cnt += 1
            meta = json.load(f)['lemnatec_measurement_metadata']
            time = (meta['gantry_system_variable_metadata']['time'])
            filename = i.split('/')[-1]
            # Gantry loc metadata
            gantry_x = float(meta['gantry_system_variable_metadata']['position x [m]'])
            gantry_y = float(meta['gantry_system_variable_metadata']['position y [m]'])
            gantry_z = float(meta['gantry_system_variable_metadata']['position z [m]'])
            
            # Sensor loc metadata
            sens_x = float(meta['sensor_fixed_metadata']['location in camera box x [m]'])
            sens_y = float(meta['sensor_fixed_metadata']['location in camera box y [m]'])
            sens_z = float(meta['sensor_fixed_metadata']['location in camera box z [m]'])
            #  gantry_x_pos = 
            z_offset = 0.76
            sens_loc_x = gantry_x + sens_x
            sens_loc_y = gantry_y + sens_y
            sens_loc_z = gantry_z + z_offset + sens_z #offset in m
            fov_x, fov_y = float(meta['sensor_fixed_metadata']['field of view x [m]']), float(meta['sensor_fixed_metadata']['field of view y [m]'])
            B = sens_loc_z
            A_x = np.arctan((0.5*float(fov_x))/2)
            A_y = np.arctan((0.5*float(fov_y))/2)
            L_x = 2*B*np.tan(A_x)
            L_y = 2*B*np.tan(A_y)
            x_n = sens_loc_x + (L_x/2)
            x_s = sens_loc_x - (L_x/2)
            y_w = sens_loc_y + (L_y/2)
            y_e = sens_loc_y - (L_y/2)
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
            print(b_box)
            
            JSON_dict[cnt] = {
                "time": time,
                "filename": filename,
                "gantry_x": sens_loc_x,
                "gantry_y": sens_loc_y,
                "gantry_z": sens_loc_z,
                "b_box": b_box}

    JSON_df = pd.DataFrame.from_dict(JSON_dict, orient='index', columns=['time','filename','gantry_x','gantry_y','gantry_z', 'b_box'])
    
    # Converts gantry/scanners location to lat lon
    GPS_latlon = scanalyzer_to_latlon(JSON_df['gantry_x'], JSON_df['gantry_y'])
    GPS_latlon_df = pd.DataFrame(GPS_latlon).transpose()
    GPS_latlon_df.columns = ['GPS_lon', 'GPS_lat']

    # Creates polygons for plots    
    polygon_list = []

    for i, row in JSON_df.iterrows():
        bbox = JSON_df['b_box'].loc[i]
        polygon = Polygon([[bbox[2], bbox[1]], [bbox[3], bbox[1]], [bbox[3], bbox[0]], [bbox[2], bbox[0]]])
        polygon_list.append(polygon)

    JSON_df['bbox_geometry'] = polygon_list

    JSON_df['time'] = pd.to_datetime(JSON_df.time)
    JSON_df = JSON_df.sort_values(by ='time')

    def intersection(bbox_polygon):
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

    JSON_df["plot"] = None
    for i, row in JSON_df.iterrows():
        bbox_polygon = row['bbox_geometry']
        plot = intersection(bbox_polygon)
        JSON_df.at[i,'plot'] = plot  
    
    return JSON_df

#----------------------------------

def Env_data():
    # Env logger data
    args = get_args()    
    #command = f'iget -rKTPf -N 0 /iplant/home/shared/phytooracle/{args.season}/level_1/EnvironmentLogger/{args.date}_clean.tar.gz'
    command = f'wget https://data.cyverse.org/dav-anon/iplant/projects/phytooracle/{args.season}/level_1/EnvironmentLogger/{args.date}_clean.tar.gz'
    subprocess.call(command, shell = True)
    command = f'tar -xvf {args.date}_clean.tar.gz'
    subprocess.call(command, shell = True)

    # Retrieve csv data and organize/clean up
    EnvL_data = pd.read_csv(f'./{args.date}_clean.csv')
    EnvL_data['Time'] = pd.to_datetime(EnvL_data['Time'])
    Envlog_clean = EnvL_data[['Time', 'Sun Direction', 'Temperature', 'Photosynthetically active radiation', 'Wind velocity']]
    #print("Env Data Retrieved")
    return Envlog_clean

#----------------------------------
def AZMget(JSON_df):
    
    args = get_args()
    season = args.season
    date = args.date

    year = date[2:4]

    AZmet_data = pd.read_csv(f"https://cals.arizona.edu/azmet/data/06{year}rh.txt", names = ["Year", "Day", "Hour", 
                                            "Air Temperature", "Relative Humidity", 
                                            "VPD", "Solar Radiation", "Precipitation", 
                                            "4 inch Soil T", "12 inch Soil T", 
                                            "Avg Wind Speed", "Wind Vector Magnitude", 
                                            "Wind Vector Direction", "Wind Direction STDEV", 
                                            "Max Wind Speed", "Reference Evapotranspiration", 
                                            "Actual Vapor Pressure", "Dewpoint"])
    print("Document downloaded, loaded")

    AZmet_df = pd.DataFrame(AZmet_data)

    AZmet_df['combined'] = AZmet_df["Year"]*1000 + AZmet_df["Day"]
    AZmet_df['date'] = pd.to_datetime(AZmet_df["combined"], format = "%Y%j")
    AZmet_df = AZmet_df.set_index('date')
    
    del AZmet_df['combined']
    
    date_of_interest = AZmet_df[AZmet_df.index == date]
    
    date_of_interest = date_of_interest.reset_index()
    date_of_interest['Hour'] = pd.to_timedelta(date_of_interest['Hour'], unit = 'h')
    date_of_interest['date'] = date_of_interest['date'] + date_of_interest['Hour']
    date_of_interest = date_of_interest.set_index('date')
    date_of_interest = date_of_interest.reset_index()
    print("* * AZmet Date of Interest Gathered * *")
    
    #image_file = JSON_df
    JSON_df['time'] = pd.to_datetime(JSON_df['time'])
    JSON_df = pd.DataFrame(JSON_df)
    
    print("* * JSON_df Formatted * *")
    
    EnvLog = Env_data()
    EnvLog = EnvLog.set_index('Time')
    EnvLog = EnvLog.reset_index()
    
    print("* * Environment Logger Date of Interest Gathered * *")
    
    return JSON_df, EnvLog, date_of_interest
#----------------------------------
def azmet_dict(JSON_df):
    AZmet_dict = {}
    JSON_df, EnvLog, date_of_interest = AZMget(JSON_df)
    cnt = 0
    for i, row in JSON_df.iterrows():
        cnt += 1
        time = row['time'].round('H')
        result_index = time.hour
        AZmet_temp = date_of_interest['Air Temperature'].iloc[result_index]
        AZmet_wind = date_of_interest['Avg Wind Speed'].iloc[result_index]
        AZmet_vpd = date_of_interest['VPD'].iloc[result_index]
        AZmet_solar = date_of_interest['Solar Radiation'].iloc[result_index]
        AZmet_rh = date_of_interest['Relative Humidity'].iloc[result_index]
        Env_temp = EnvLog['Temperature'].iloc[result_index]
        Env_wind = EnvLog['Wind velocity'].iloc[result_index]
        AZmet_dict[cnt] = {'azmet_atm_temp': AZmet_temp, 'azmet_wind_velocity': AZmet_wind, 'azmet_VPD': AZmet_vpd, 'azmet_solar_radiation':
                          AZmet_solar, 'relative_humidity': AZmet_rh, 'env_temp': Env_temp, 'env_wind': Env_wind}
    return pd.DataFrame.from_dict(AZmet_dict)

#----------------------------------

def compile_info(JSON_df, environmental_df):
    JSON_df['azmet_atm_temp'] = environmental_df['azmet_atm_temp']
    JSON_df['azmet_wind_velocity'] = environmental_df['azmet_wind_velocity']
    JSON_df['azmet_VPD'] = environmental_df['azmet_VPD']
    JSON_df['azmet_solar_radiation'] = environmental_df['azmet_solar_radiation']
    JSON_df['relative_humidity'] = environmental_df['relative_humidity']
    JSON_df['env_temp'] = environmental_df['env_temp']
    JSON_df['env_wind'] = environmental_df['env_wind']
    return JSON_df


#----------------------------------

def listToStringWithoutBrackets(list1):
    return str(list1).replace('[','').replace(']','')


#----------------------------------

# Takes each image (with multiple plots per row) and breaks it up so you have all of the plots listed out with associated data
def expand_plots(clean_file):
    file = clean_file
    file['plot'] = (file['plot'].apply(listToStringWithoutBrackets)).apply(eval)
    plot_expand = file['plot'].apply(pd.Series)
    plot_expand['time'] = file['time']
    plot_expand['Image Name'] = file['filename']
    plot_expand['env_temp'] = file['env_temp']
    plot_expand['env_wind'] = file['env_wind']
    plot_expand['azmet_atm_temp'] = file['azmet_atm_temp']
    plot_expand['azmet_wind_velocity'] = file['azmet_wind_velocity']
    plot_expand['azmet_VPD'] = file['azmet_VPD']
    plot_expand['azmet_solar_radiation'] = file['azmet_solar_radiation']
    plot_expand['relative_humidity'] = file['relative_humidity']
    stacked = plot_expand.set_index(['time', 'Image Name', 'env_temp', 'env_wind', 'azmet_atm_temp', 'azmet_wind_velocity',
                                     'azmet_VPD', 'azmet_solar_radiation', 'relative_humidity']).stack()
    stack_df = pd.DataFrame(stacked).reset_index()
    del stack_df['level_9']
    final_df = stack_df.rename(columns = {0:'Plot'})
    return final_df

# Takes a plot as paramter and finds how long it took to image that plot
def plot_scan_time(interest_plot):
    final_df = expand_plots(file)
    plot_duration = final_df.loc[final_df['Plot'] == interest_plot]
    start = pd.to_datetime(plot_duration['time'].iloc[(0)])
    end = pd.to_datetime(plot_duration['time'].iloc[(-1)])
    return end - start

# Takes a plot as parameter and finds how much the temperature changed during the time it took to image the plot
def plot_temp_change(interest_plot):
    final_df = expand_plots(file)
    plot_duration = final_df.loc[final_df['Plot'] == interest_plot]
    start = plot_duration['env_temp'].iloc[(0)]
    end = plot_duration['env_temp'].iloc[(-1)]
    return end - start

# Expands plots so each plot is represented
def images_making_plot(interest_plot):
    final_df = expand_plots(file)
    plot_list = final_df.loc[final_df['Plot'] == interest_plot]
    return plot_list


#----------------------------------

def main():
    args = get_args()
    season = args.season

    JSON_df = md_shp()
    
    environmental_df = azmet_dict(JSON_df).transpose()
    file = compile_info(JSON_df, environmental_df)
    finale_df = expand_plots(file)
    env_logger = JSON_df[['time', 'filename', 'azmet_atm_temp', 'azmet_wind_velocity', 'azmet_VPD', 
                                        'azmet_solar_radiation', 'relative_humidity', 'env_temp','env_wind']].set_index('filename')


    print("Im through!")

    img_plot = finale_df[['Image Name', 'Plot']]
    img_plot.columns = ['image', 'plot']
    ######

    img_plot['Date and Time'] = img_plot['azmet_atm_temp'] = img_plot['azmet_wind_velocity'] = img_plot['azmet_VPD'] = img_plot['azmet_solar_radiation'] = img_plot['relative_humidity'] = img_plot['env_temp'] = img_plot['env_wind'] = None

    for i, row in img_plot.iterrows():
        meta = row['image']
        try:
            datetime = env_logger.loc[meta, 'time']
            azmet_temp = env_logger.loc[meta, 'azmet_atm_temp']
            azmet_wind_vel = env_logger.loc[meta, 'azmet_wind_velocity']
            azmet_vpd = env_logger.loc[meta, 'azmet_VPD']
            sol_rad = env_logger.loc[meta, 'azmet_solar_radiation']
            temp = env_logger.loc[meta, 'env_temp']
            win_vel = env_logger.loc[meta, 'env_wind']
            rel_hum = env_logger.loc[meta, 'relative_humidity']
            
            img_plot.at[i, 'Date and Time'] = datetime
            img_plot.at[i, 'azmet_atm_temp'] = azmet_temp
            img_plot.at[i, 'azmet_wind_velocity'] = azmet_wind_vel
            img_plot.at[i, 'azmet_VPD'] = azmet_vpd
            img_plot.at[i, 'azmet_solar_radiation'] = sol_rad
            img_plot.at[i, 'env_temp'] = temp
            img_plot.at[i, 'env_wind'] = win_vel
            img_plot.at[i, 'relative_humidity'] = rel_hum
            
        except:
            pass

    df_agg = img_plot.drop(['image', 'Date and Time'], axis=1)#.set_index('plot')#.groupby('plot')

    temp_dict = {}
    cnt = 0

    for plot in df_agg['plot'].unique().tolist():
        try:
            cnt += 1 

            select_df = df_agg.set_index('plot').loc[plot]
            temp_median = select_df['azmet_atm_temp'].median()
            temp_mean = select_df['azmet_atm_temp'].mean()
            temp_std = select_df['azmet_atm_temp'].std()
            
            azmet_wind_vel = select_df['azmet_wind_velocity'].median()
            azmet_vpd = select_df['azmet_VPD'].median()
            sol_rad = select_df['azmet_solar_radiation'].median()
            temp = select_df['env_temp'].median()
            wind_vel = select_df['env_wind'].median()
            rel_hum = select_df['relative_humidity'].median()
            

            temp_dict[cnt] = {'plot': plot,
                              'median': temp_median,
                              'mean': temp_mean, 
                              'std_dev': temp_std, 
                              'azmet_wind_velocity': azmet_wind_vel, 
                              'azmet_VPD': azmet_vpd, 
                              'azmet_solar_radiation': sol_rad, 
                              'env_temp': temp, 
                              'env_wind': wind_vel,
                             'relative_humidity': rel_hum}
        except:
            pass

    result = pd.DataFrame.from_dict(temp_dict, orient='index').set_index('plot')

    # Reads in csv file of individual plant temperatures
    plant_detections = pd.read_csv(args.csv)

    # Adds the field information and PCT values to the already existing csv that is indexed by plot
    plant_detections['norm_temp'] = plant_detections['atm_temp'] = None

    for i, row in plant_detections.iterrows():

        try:
            plot = row['plot'].replace('_', ' ')
            plant_temp = row['median']
            
            temp_df = result.loc[plot]
            atm_temp = temp_df['median']
            norm_temp =  atm_temp - plant_temp
            
            azmet_wind_vel = temp_df['azmet_wind_velocity']
            azmet_vpd = temp_df['azmet_VPD']
            sol_rad = temp_df['azmet_solar_radiation']
            temp = temp_df['env_temp']
            wind_vel = temp_df['env_wind']
            rel_hum = temp_df['relative_humidity']
            
            plant_detections.at[i, 'norm_temp'] = norm_temp
            plant_detections.at[i, 'atm_temp'] = atm_temp
            
            plant_detections.at[i, 'azmet_wind_velocity'] = azmet_wind_vel
            plant_detections.at[i, 'azmet_VPD'] = azmet_vpd
            plant_detections.at[i, 'azmet_solar_radiation'] = sol_rad
            plant_detections.at[i, 'env_temp'] = temp
            plant_detections.at[i, 'env_wind'] = wind_vel
            plant_detections.at[i, 'relative_humidity'] = rel_hum
        except:
            pass

    plant_detections.to_csv(os.path.join(args.outdir, f'{args.date}_indiv_atm_temp.csv'))
    print(f'Done, see outputs in ./{args.outdir}.')

#----------------------------------
if __name__ == '__main__':
    main()
