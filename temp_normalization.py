#!/usr/bin/env python3

import pandas as pd
import datetime
from datetime import timedelta
import numpy as np
import scipy
from scipy.optimize import curve_fit
from itertools import cycle
from datetime import datetime
from datetime import timedelta
from scipy.interpolate import UnivariateSpline
import subprocess
from shapely.geometry import Polygon, Point, mapping
import shapely.wkt
import geopandas as gpd
from pyproj import Proj
import glob
import utm
import json
import argparse
import os
import glob


# ----------------------------------

def get_args():
    """Get command-line arguments"""

    parser = argparse.ArgumentParser(
        description="Individual plant temperature extraction",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument("dir", metavar="dir", help="Directory containing geoTIFFs")

    parser.add_argument(
        "-g",
        "--geojson",
        help="GeoJSON containing plot boundaries",
        metavar="geojson",
        type=str,
        default=None,
        required=True,
    )

    parser.add_argument(
        "-s",
        "--season",
        help=' "season_10_lettuce_yr_2020" OR "season_11_sorghum_yr_2020" ',
        metavar="season",
        type=str,
        default=None,
        required=True,
    )

    parser.add_argument(
        "-d",
        "--date",
        help="date to process, usually the same as dir",
        metavar="date",
        type=str,
        default=None,
        required=True,
    )

    parser.add_argument(
        "-o",
        "--outdir",
        help="Output directory where resulting csv will be saved",
        metavar="str",
        type=str,
        default="Thermal_Output",
    )

    parser.add_argument(
        "-c",
        "--csv",
        help=' "Individual Plant Temp detections CSV" ',
        metavar="csv",
        type=str,
        default=None,
        required=True,
    )

    return parser.parse_args()


# ----------------------------------

def utm_to_latlon(utm_x, utm_y):
    p = Proj(proj='utm',zone=12,ellps='WGS84')
    lon, lat = p(utm_x, utm_y, inverse=True)
    return lat, lon

# ----------------------------------

def scanalyzer_to_utm(gantry_x, gantry_y):
    """Convert coordinates from gantry to UTM 12N"""

    # TODO: Hard-coded
    # Linear transformation coefficients
    ay = 3659974.971
    by = 1.0002
    cy = 0.0078
    ax = 409012.2032
    bx = 0.009
    cx = -0.9986

    utm_x = ax + (bx * gantry_x) + (cx * gantry_y)
    utm_y = ay + (by * gantry_x) + (cy * gantry_y)

    return utm_x, utm_y

# ----------------------------------

def scanalyzer_to_latlon(gantry_x, gantry_y):
    """Convert coordinates from gantry to lat/lon"""
    utm_x, utm_y = scanalyzer_to_utm(gantry_x, gantry_y)
    return utm_to_latlon(utm_x, utm_y)

# ----------------------------------

def md_shp():

    args = get_args()

    # Load required files as well as files to process
    pathlist = glob.glob(f"{args.dir}/*/*.json")

    shp = gpd.read_file(f"{args.geojson}")

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
            meta = json.load(f)["lemnatec_measurement_metadata"]
            time = meta["gantry_system_variable_metadata"]["time"]
            filename = i.split("/")[-1]
            # Gantry loc metadata
            gantry_x = float(meta["gantry_system_variable_metadata"]["position x [m]"])
            gantry_y = float(meta["gantry_system_variable_metadata"]["position y [m]"])
            gantry_z = float(meta["gantry_system_variable_metadata"]["position z [m]"])

            # Sensor loc metadata
            sens_x = float(
                meta["sensor_fixed_metadata"]["location in camera box x [m]"]
            )
            sens_y = float(
                meta["sensor_fixed_metadata"]["location in camera box y [m]"]
            )
            sens_z = float(
                meta["sensor_fixed_metadata"]["location in camera box z [m]"]
            )
            #  gantry_x_pos =
            z_offset = 0.76
            sens_loc_x = gantry_x + sens_x
            sens_loc_y = gantry_y + sens_y
            sens_loc_z = gantry_z + z_offset + sens_z  # offset in m
            fov_x, fov_y = (
                float(meta["sensor_fixed_metadata"]["field of view x [m]"]),
                float(meta["sensor_fixed_metadata"]["field of view y [m]"]),
            )
            B = sens_loc_z
            A_x = np.arctan((0.5 * float(fov_x)) / 2)
            A_y = np.arctan((0.5 * float(fov_y)) / 2)
            L_x = 2 * B * np.tan(A_x)
            L_y = 2 * B * np.tan(A_y)
            x_n = sens_loc_x + (L_x / 2)
            x_s = sens_loc_x - (L_x / 2)
            y_w = sens_loc_y + (L_y / 2)
            y_e = sens_loc_y - (L_y / 2)
            bbox_nw_latlon = scanalyzer_to_latlon(x_n, y_w)
            bbox_se_latlon = scanalyzer_to_latlon(x_s, y_e)

            # TERRA-REF
            lon_shift = 0.000020308287

            # Drone
            lat_shift = 0.000018292  # 0.000015258894
            b_box = (
                bbox_se_latlon[0] - lat_shift,
                bbox_nw_latlon[0] - lat_shift,
                bbox_nw_latlon[1] + lon_shift,
                bbox_se_latlon[1] + lon_shift,
            )
            print(b_box)

            JSON_dict[cnt] = {
                "time": time,
                "filename": filename,
                "gantry_x": sens_loc_x,
                "gantry_y": sens_loc_y,
                "gantry_z": sens_loc_z,
                "b_box": b_box,
            }

    JSON_df = pd.DataFrame.from_dict(
        JSON_dict,
        orient="index",
        columns=["time", "filename", "gantry_x", "gantry_y", "gantry_z", "b_box"],
    )

    GPS_latlon = scanalyzer_to_latlon(JSON_df["gantry_x"], JSON_df["gantry_y"])
    GPS_latlon_df = pd.DataFrame(GPS_latlon).transpose()
    GPS_latlon_df.columns = ["GPS_lon", "GPS_lat"]

    # Creates polygons for plots
    polygon_list = []

    for i, row in JSON_df.iterrows():
        bbox = JSON_df["b_box"].loc[i]
        polygon = Polygon(
            [
                [bbox[2], bbox[1]],
                [bbox[3], bbox[1]],
                [bbox[3], bbox[0]],
                [bbox[2], bbox[0]],
            ]
        )
        polygon_list.append(polygon)

    JSON_df["bbox_geometry"] = polygon_list

    
    JSON_df["time"] = pd.to_datetime(JSON_df.time)
    JSON_df = JSON_df.sort_values(by="time")

    def intersection(bbox_polygon):
        intersects = bbox_polygon.intersects
        plot = None
        intersection_list = []
        for i, row in shp.iterrows():
            plot_polygon = row["geometry"]
            intersection = intersects(plot_polygon)
            if intersection == True:
                plot = [row["ID"]]
                intersection_list.append(plot)
        return intersection_list

    JSON_df["plot"] = None
    for i, row in JSON_df.iterrows():
        bbox_polygon = row["bbox_geometry"]
        plot = intersection(bbox_polygon)
        JSON_df.at[i, "plot"] = plot

    return JSON_df


# ----------------------------------

def AZmet(date):
    year = date[2:4]

    AZmet_data = pd.read_csv(
        f"https://cals.arizona.edu/azmet/data/06{year}rh.txt",
        names=[
            "Year",
            "Day",
            "Hour",
            "Air Temperature",
            "Relative Humidity",
            "VPD",
            "Solar Radiation",
            "Precipitation",
            "4 inch Soil T",
            "12 inch Soil T",
            "Avg Wind Speed",
            "Wind Vector Magnitude",
            "Wind Vector Direction",
            "Wind Direction STDEV",
            "Max Wind Speed",
            "Reference Evapotranspiration",
            "Actual Vapor Pressure",
            "Dewpoint",
        ],
    )
    print("Document downloaded, loaded")
    AZmet_df = pd.DataFrame(AZmet_data)

    AZmet_df["combined"] = AZmet_df["Year"] * 1000 + AZmet_df["Day"]
    AZmet_df["date"] = pd.to_datetime(AZmet_df["combined"], format="%Y%j")
    AZmet_df = AZmet_df.set_index("date")
    
    del AZmet_df["combined"]
    
    return AZmet_df

# ----------------------------------

def find_date(AZmet_df, date):

    previous_day = str((pd.to_datetime(date) - timedelta(days=1)).date())
    next_day = str((pd.to_datetime(date) + timedelta(days=1)).date())

    yesterday = AZmet_df[AZmet_df.index == previous_day]
    today = AZmet_df[AZmet_df.index == date]
    tomorrow = AZmet_df[AZmet_df.index == next_day]

    concat = pd.concat([yesterday, today]).reset_index()
    concat_all = pd.concat([concat, tomorrow.reset_index()]).reset_index()

    Hour_0_index = concat[concat["date"] == date].index[0] - 1
    Hour_25_index = concat[concat["date"] == date].index[-1] + 1

    Hour_0 = pd.DataFrame(concat.iloc[Hour_0_index]).transpose()
    Hour_25 = pd.DataFrame(concat_all.iloc[Hour_25_index]).transpose()

    date_of_interest_pre = pd.concat([Hour_0, today.reset_index()])
    date_of_interest = pd.concat([date_of_interest_pre, Hour_25])

    date_of_interest = date_of_interest.reset_index()
    del date_of_interest["index"]
    del date_of_interest["level_0"]

    date_of_interest["date"][0] = date
    date_of_interest["Hour"][0] = 0

    date_of_interest["date"][25] = date
    date_of_interest["Hour"][25] = 25

    return date_of_interest

# ----------------------------------

def Env_data(date):
    # Env logger data
    args = get_args()
    print("Downloading Environement Logger tarfile")
    # command = f'iget -rKTPf -N 0 /iplant/home/shared/phytooracle/{args.season}/level_1/EnvironmentLogger/{args.date}_clean.tar.gz'
    command = f"wget https://data.cyverse.org/dav-anon/iplant/projects/phytooracle/{args.season}/level_1/EnvironmentLogger/{date}_clean.tar.gz"
    subprocess.call(command, shell=True)
    command = f"tar -xvf {date}_clean.tar.gz"
    subprocess.call(command, shell=True)
    print("Environment Logger data has been downloaded and uncompressed")
    # Retrieve csv data and organize/clean up
    EnvL_data = pd.read_csv(f"./{date}_clean.csv")
    EnvL_data["Time"] = pd.to_datetime(EnvL_data["Time"])
    Envlog_clean = EnvL_data[
        [
            "Time",
            "Sun Direction",
            "Temperature",
            "Photosynthetically active radiation",
            "Wind velocity",
        ]
    ]
    # print("Env Data Retrieved")
    return Envlog_clean

# ----------------------------------

def splines(
    df, xvar, yvar
):  # xdata would be the information you want to use ex: df['Hour']
    xdata = df[xvar]
    ydata = df[yvar]
    x, y = xdata.values, ydata.values
    spl = UnivariateSpline(x, y)
    # spl.set_smoothing_factor(50)

    xrange = np.arange(0, 24, 0.01667)

    d = {"Minute": np.arange(len(xrange)), yvar: spl(xrange)}
    finer_df = pd.DataFrame(data=d)

    hour_list = np.arange(0, 24, 1)

    K = 60
    res = [ele for ele in hour_list for i in range(K)]

    finer_df["Hour"] = res

    minute_cycle = cycle(np.arange(0, 60, 1))
    finer_df["Minute"] = [next(minute_cycle) for cycle in range(len(finer_df))]

    year = df["Year"].unique()
    date = df["date"].unique()
    finer_df["Year"] = year[0]
    finer_df["date"] = date[0]

    finer_df["Hour"] = pd.to_timedelta(finer_df["Hour"], unit="h")
    finer_df["Minute"] = pd.to_timedelta(finer_df["Minute"], unit="m")

    finer_df["date"] = finer_df["date"] + finer_df["Hour"] + finer_df["Minute"]
    return finer_df

# ----------------------------------

def retrieve_splines(df):
    temp_df = splines(df, "Hour", "Air Temperature")
    temp_df["VPD"] = splines(df, "Hour", "VPD")["VPD"]
    temp_df["Relative Humidity"] = splines(df, "Hour", "Relative Humidity")[
        "Relative Humidity"
    ]
    temp_df["Avg Wind Speed"] = splines(df, "Hour", "Avg Wind Speed")["Avg Wind Speed"]
    temp_df["Solar Radiation"] = splines(df, "Hour", "Solar Radiation")[
        "Solar Radiation"
    ]
    return temp_df

# ----------------------------------

def AZMget(JSON_df):
    # ----------------
    args = get_args()

    season = args.season

#     date = args.date
    # ----------------

    # Finds unique dates
    date_list = []
    JSON_df['time'] = pd.to_datetime(JSON_df['time'])
    dates = JSON_df["time"].dt.date.unique()
    for date in dates:
        date_list.append(date)

    if len(date_list) > 1:
        date1 = str(date_list[0])
        date2 = str(date_list[1])

        AZmet_date1 = AZmet(date1)
        AZmet_date2 = AZmet(date2)

        if AZmet_date1.equals(AZmet_date2) == False:
            AZmet_df = AZmet_date1
            AZmet_df = pd.concat([AZmet_date1, AZmet_date2])
        else:
            AZmet_df = AZmet_date1

        date_of_interest1 = find_date(AZmet_df, date1)
        date_of_interest2 = find_date(AZmet_df, date2)

        env_date1 = Env_data(date1)
        env_date2 = Env_data(date2)

        temp_df1 = retrieve_splines(date_of_interest1)
        temp_df2 = retrieve_splines(date_of_interest2)

        temp_df = pd.concat([temp_df1, temp_df2])

        EnvLog = pd.concat([env_date1, env_date2])

    else:
        date1 = str(date_list[0])
        AZmet_df = AZmet(date1)

        date_of_interest = find_date(AZmet_df, date1)

        EnvLog = Env_data(date1)

        temp_df = retrieve_splines(date_of_interest)

    temp_df["Hour"] = temp_df["date"].dt.hour
    temp_df["Minute"] = temp_df["date"].dt.minute

    EnvLog = EnvLog.set_index("Time")
    EnvLog = EnvLog.reset_index() 

    date_of_interest = temp_df
    print("* * AZmet Date of Interest Gathered * *")

    print("* * Environment Logger Date of Interest Gathered * *")

    # image_file = JSON_df
    JSON_df["time"] = pd.to_datetime(JSON_df["time"])
    JSON_df = pd.DataFrame(JSON_df)

    print("* * JSON_df Formatted * *")

    return JSON_df, date_of_interest, EnvLog

# ----------------------------------

def azmet_dict(JSON_df):
    AZmet_dict = {}
    JSON_df, date_of_interest, EnvLog = AZMget(JSON_df)
    cnt = 0
    for i, row in JSON_df.iterrows():
        cnt += 1
        time = row["time"]
        timestamp = time.round("min")
        result_index = date_of_interest[
            date_of_interest["date"] == timestamp
        ].index.values[0]
        result_index_env = EnvLog[EnvLog["Time"] == timestamp].index.values[0]
        # result_index = date_of_interest["date"].sub(time).abs().idxmin()
        # result_index_env = EnvLog["Time"].sub(time).abs().idxmin()
        #         time = row['time'].round('H')
        #         result_index = time.hour
        AZmet_temp = date_of_interest["Air Temperature"].iloc[result_index]
        AZmet_wind = date_of_interest["Avg Wind Speed"].iloc[result_index]
        AZmet_vpd = date_of_interest["VPD"].iloc[result_index]
        AZmet_solar = date_of_interest["Solar Radiation"].iloc[result_index]
        AZmet_rh = date_of_interest["Relative Humidity"].iloc[result_index]
        Env_temp = EnvLog["Temperature"].iloc[result_index_env]
        Env_wind = EnvLog["Wind velocity"].iloc[result_index_env]
        AZmet_dict[cnt] = {
            "azmet_atm_temp": AZmet_temp,
            "azmet_wind_velocity": AZmet_wind,
            "azmet_VPD": AZmet_vpd,
            "azmet_solar_radiation": AZmet_solar,
            "relative_humidity": AZmet_rh,
            "env_temp": Env_temp,
            "env_wind": Env_wind,
        }
    return pd.DataFrame.from_dict(AZmet_dict)

# ----------------------------------

def compile_info(JSON_df, environmental_df):
    JSON_df["azmet_atm_temp"] = environmental_df["azmet_atm_temp"]
    JSON_df["azmet_wind_velocity"] = environmental_df["azmet_wind_velocity"]
    JSON_df["azmet_VPD"] = environmental_df["azmet_VPD"]
    JSON_df["azmet_solar_radiation"] = environmental_df["azmet_solar_radiation"]
    JSON_df["relative_humidity"] = environmental_df["relative_humidity"]
    JSON_df["env_temp"] = environmental_df["env_temp"]
    JSON_df["env_wind"] = environmental_df["env_wind"]
    return JSON_df

# ----------------------------------

def listToStringWithoutBrackets(list1):
    return str(list1).replace("[", "").replace("]", "")

# ----------------------------------

def expand_plots(clean_file):
    file = clean_file
    file["plot"] = (file["plot"].apply(listToStringWithoutBrackets)).apply(eval)
    plot_expand = file["plot"].apply(pd.Series)
    plot_expand["time"] = file["time"]
    plot_expand["Image Name"] = file["filename"]
    plot_expand["bbox_geometry"] = file["bbox_geometry"]
    plot_expand["env_temp"] = file["env_temp"]
    plot_expand["env_wind"] = file["env_wind"]
    plot_expand["azmet_atm_temp"] = file["azmet_atm_temp"]
    plot_expand["azmet_wind_velocity"] = file["azmet_wind_velocity"]
    plot_expand["azmet_VPD"] = file["azmet_VPD"]
    plot_expand["azmet_solar_radiation"] = file["azmet_solar_radiation"]
    plot_expand["relative_humidity"] = file["relative_humidity"]
    stacked = plot_expand.set_index(
        [
            "time",
            "Image Name",
            "bbox_geometry",
            "env_temp",
            "env_wind",
            "azmet_atm_temp",
            "azmet_wind_velocity",
            "azmet_VPD",
            "azmet_solar_radiation",
            "relative_humidity",
        ]
    ).stack()
    stack_df = pd.DataFrame(stacked).reset_index()
    del stack_df["level_10"]
    final_df = stack_df.rename(columns={0: "Plot"})
    return final_df

# ----------------------------------

def get_point(df):
    df['Point'] = None
    for i, row in df.iterrows():
        lon = row['lon']
        lat = row['lat']
        point = Point(lon, lat)

        df.at[i, 'Point'] = point
    return df

# ----------------------------------

## Results in a dictionary that finds which image a particular point has been imaged in
def find_images(polygon_df, point):
    intersection_dict = {}
    cnt = 0
    for i, row in polygon_df.iterrows():
        cnt += 1
        # polygon = shapely.wkt.loads(row['bbox_geometry'])
        polygon = row['bbox_geometry']
        intersection = point.intersects(polygon)
        if intersection == True:
            time = row['time']
            plot = row['plot']
            image = row["filename"]
            azmet_temp = row['azmet_atm_temp']
            azmet_wind = row['azmet_wind_velocity']
            azmet_VPD = row['azmet_VPD']
            azmet_solar_radiation = row['azmet_solar_radiation']
            relative_humidity = row['relative_humidity']
            env_temp = row['env_temp']
            env_wind = row['env_wind']
            intersection_dict[cnt] = {
                "Date and Time": time,
                "image": image,
                "plot": plot,
                "polygon": polygon,
                "azmet_atm_temp": azmet_temp,
                "azmet_wind_velocity": azmet_wind,
                "azmet_VPD": azmet_VPD,
                "azmet_solar_radiation": azmet_solar_radiation,
                "azmet_relative_humidity": relative_humidity,
                "env_temp": env_temp,
                "env_wind": env_wind
            }
        dict_df = pd.DataFrame.from_dict(intersection_dict).T

    return dict_df

# ----------------------------------

def main():
    # --------------
    args = get_args()
    season = args.season
    # --------------

    JSON_df = md_shp()

    environmental_df = azmet_dict(JSON_df).transpose()

    file = compile_info(JSON_df, environmental_df)

    # --------------
    # Reads in csv file of individual plant temperatures
    plant_detections = pd.read_csv(args.csv)

    plant_detections = get_point(plant_detections)
    # --------------


    plant_detections['azmet_atm_temp'] = plant_detections['azmet_wind_velocity'] = plant_detections['azmet_VPD'] = plant_detections['env_temp'] = plant_detections['env_wind'] = None
    for i, row in plant_detections.iterrows():
        
        try:
            point = row['Point']
            dict_df = find_images(file, point)

            temp_mean = dict_df['azmet_atm_temp'].mean()
            wind_mean = dict_df['azmet_wind_velocity'].mean()
            VPD_mean = dict_df['azmet_VPD'].mean()
            solar_mean = dict_df['azmet_solar_radiation'].mean()
            rh_mean = dict_df['azmet_relative_humidity'].mean()
            env_temp_mean = dict_df['env_temp'].mean()
            env_wind_mean = dict_df['env_wind'].mean()

            plant_detections.at[i, "azmet_atm_temp"] = temp_mean
            plant_detections.at[i, "azmet_wind_velocity"] = wind_mean
            plant_detections.at[i, "azmet_VPD"] = VPD_mean
            plant_detections.at[i, "azmet_solar_radiation"] = solar_mean
            plant_detections.at[i, "azmet_relative_humidity"] = rh_mean
            plant_detections.at[i, "env_temp"] = env_temp_mean
            plant_detections.at[i, "env_wind"] = env_wind_mean

        except:
            pass

    
    out_dir = args.outdir

    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)
    
    plant_detections.to_csv(
        os.path.join(out_dir, f"{args.date}_indiv_atm_temp.csv")
    )
    print(f"Done, see outputs in ./{out_dir}.")

# ----------------------------------

if __name__ == "__main__":
    main()
