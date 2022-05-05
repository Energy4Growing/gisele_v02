"""
GIS For Electrification (GISEle)
Developed by the Energy Department of Politecnico di Milano
Supporting Code

Group of supporting functions used inside all the process of GISEle algorithm
"""

import os
import requests
import pandas as pd
import geopandas as gpd
import numpy as np
import json
import shapely.ops
import iso8601
from scipy.spatial import distance_matrix
from scipy.spatial.distance import cdist
from shapely.geometry import Point, box, LineString, MultiPoint
from shapely.ops import split
from gisele.michele.michele import start
from gisele.data_import import import_pv_data, import_wind_data
from datetime import datetime


def l():
    """Print long separating lines."""
    print('-' * 100)


def s():
    """Print short separating lines."""
    print("-" * 40)


def nearest(row, df, src_column=None):
    """
    Find the nearest point and return the value from specified column.
    :param row: Iterative row of the first dataframe
    :param df: Second dataframe to be found the nearest value
    :param src_column: Column of the second dataframe that will be returned
    :return value: Value of the desired src_column of the second dataframe
    """

    # Find the geometry that is closest
    nearest_p = df['geometry'] == shapely.ops.nearest_points(row['geometry'],
                                                             df.unary_union)[1]
    # Get the corresponding value from df2 (matching is based on the geometry)
    value = df.loc[nearest_p, src_column].values[0]
    return value


def distance_2d(df1, df2, x, y):
    """
    Find the 2D distance matrix between two datasets of points.
    :param df1: first point dataframe
    :param df2: second point dataframe
    :param x: column representing the x coordinates (longitude)
    :param y: column representing the y coordinates (latitude)
    :return value: 2D Distance matrix between df1 and df2
    """

    d1_coordinates = {'x': df1[x], 'y': df1[y]}
    df1_loc = pd.DataFrame(data=d1_coordinates)
    df1_loc.index = df1['ID']


    d2_coordinates = {'x': df2[x], 'y': df2[y]}
    df2_loc = pd.DataFrame(data=d2_coordinates)
    df2_loc.index = df2['ID']

    value = distance_matrix(df1_loc, df2_loc)
    return value


def distance_3d(df1, df2, x, y, z):
    """
    Find the 3D euclidean distance matrix between two datasets of points.
    :param df1: first point dataframe
    :param df2: second point dataframe
    :param x: column representing the x coordinates (longitude)
    :param y: column representing the y coordinates (latitude)
    :param z: column representing the z coordinates (Elevation)
    :return value: 3D Distance matrix between df1 and df2
    """

    d1_coordinates = {'x': df1[x], 'y': df1[y], 'z': df1[z]}
    df1_loc = pd.DataFrame(data=d1_coordinates)
    df1_loc.index = df1['ID']

    d2_coordinates = {'x': df2[x], 'y': df2[y], 'z': df2[z]}
    df2_loc = pd.DataFrame(data=d2_coordinates)
    df2_loc.index = df2['ID']

    value = pd.DataFrame(cdist(df1_loc.values, df2_loc.values, 'euclidean'),
                         index=df1_loc.index, columns=df2_loc.index)
    return value


def river_intersection(gdf,resolution):
    '''
    Check which lines in the adjacency matrix cross the rivers and assign
     a very high weight to those lines
    :param gdf:
    :return:
    '''
    print('begin river intersection')
    n = gdf['X'].size
    weight_columns_x1 = np.repeat(gdf['X'].values[:,np.newaxis], n,1)
    weight_columns_y1 = np.repeat(gdf['Y'].values[:,np.newaxis], n,1)
    weight_columns_x2 = np.repeat(gdf['X'].values[np.newaxis,:], n,0)
    weight_columns_y2 = np.repeat(gdf['Y'].values[np.newaxis,:], n,0)

    weight_columns_x1_res = np.reshape(weight_columns_x1, (n*n, 1), order='F')
    weight_columns_x2_res = np.reshape(weight_columns_x2, (n*n, 1), order='F')
    weight_columns_y1_res = np.reshape(weight_columns_y1, (n*n, 1), order='F')
    weight_columns_y2_res = np.reshape(weight_columns_y2, (n*n, 1), order='F')

    df=pd.DataFrame()
    df['x1'] = weight_columns_x1_res[:,0]
    df['x2'] = weight_columns_x2_res[:,0]
    df['y1'] = weight_columns_y1_res[:,0]
    df['y2'] = weight_columns_y2_res[:,0]
    # todo-> very slow passage, need to speed it up, no sense to compute it each time,
    #todo -> take only short lines in a predefined buffer around rivers
    #import rivers, intersect them according to the area considered
    #create a buffer around them
    #filter geodf to take only point that are closer than 1.5*resolution to the river
    # it would be better to associate weights in advance,
    df['2d_dist'] = ((df['x1']-df['x2'])**2+(df['y1']-df['y2'])**2)**0.5
    # select only short lines and create linestrings
    df_short_lines = df.loc[(df['2d_dist']<resolution * 1.5) &(df['2d_dist']>0)]
    df_short_lines['geometry'] = df.apply(
        lambda x: LineString([(x['x1'], x['y1']), (x['x2'], x['y2'])]), axis=1)

    geodf = gpd.GeoDataFrame(df_short_lines, geometry='geometry')

    # todo -> automatize this step!!!

    geodf.crs = 'epsg:22287'
    case_study='test_3'
    dir='Case studies/'+case_study
    # intersect short lines
    test_inters = gpd.read_file('C:/Users/silvi/Progetti Github/Gisele_development/Case studies/test_3/Input/rivers.shp')
    a = np.empty(shape=(geodf['geometry'].size, test_inters['geometry'].size))
    for i, row in test_inters.iterrows():
        a[:, i] = geodf['geometry'].intersects(row['geometry'])
    a_tot = np.amax(a, 1)

    geodf['intersection'] = a_tot
    df['a_tot'] = 0
    df.loc[geodf.index,'a_tot'] = geodf['intersection']
    matrix = df['a_tot'].values.reshape((n, n), order='F') * 100

    # df['geometry']=df.apply(lambda x: LineString([(x['x1'], x['y1']),(x['x2'], x['y2']) ]),axis=1)
    # geodf= gpd.GeoDataFrame(df,geometry='geometry')
    # geodf.crs = 'epsg:32737'
    #
    # test_inters=gpd.read_file('test_inters.shp')
    # a=np.empty(shape =(geodf['geometry'].size,test_inters['geometry'].size))
    # for i, row in test_inters.iterrows():
    #       a[:,i] = geodf['geometry'].intersects(row['geometry'])
    # a_tot = np.amax(a, 1)
    # geodf['intersection']=a_tot
    # matrix=a_tot.reshape((n, n), order='F')*100
    return matrix

def river_intersection(graph_gdf,box,graph,rivers):
    #todo ->create the rivers geodf
    rivers_clipped = gpd.clip(rivers, box) # box needs to be gdf with same crs as rivers

    graph_gdf.loc[graph_gdf[rivers_clipped['geometry'].intersects(graph_gdf['geometry'])],'Cost'] = \
        graph_gdf.loc[graph_gdf[rivers_clipped['geometry'].intersects(graph_gdf['geometry'])],'Cost']*5

    graph_gdf['inters'] =''
    for i, row in graph_gdf.iterrows():
        graph_gdf.loc[i,'inters'] = rivers_clipped['geometry'].intersects(row['geometry'])
    graph_gdf.loc[graph_gdf['inters']==True,'Cost'] = graph_gdf.loc[graph_gdf['inters']==True,'Cost']*10
    graph_intersect =graph_gdf[graph_gdf['inters']==True]

    for i,row in graph_intersect.iterrows():
        graph[row['ID1']][row['ID2']]['weight'] = row['Cost']

    return (graph,graph_gdf)

def cost_matrix(gdf, dist_3d_matrix, line_bc,resolution,Rivers_option):
    """
    Creates the cost matrix in €/km by finding the average weight between
    two points and then multiplying by the distance and the line base cost.
    :param gdf: Geodataframe being analyzed
    :param dist_3d_matrix: 3D distance matrix of all points [meters]
    :param line_bc: line base cost for line deployment [€/km]
    :return value: Cost matrix of all the points present in the gdf [€]
    """
    # Altitude distance in meters
    weight = gdf['Weight'].values
    n = gdf['X'].size

    weight_columns = np.repeat(weight[:, np.newaxis], n, 1)
    weight_rows = np.repeat(weight[np.newaxis, :], n, 0)
    if Rivers_option:
        river_inters =river_intersection(gdf,resolution)
        total_weight = (weight_columns + weight_rows) / 2 + river_inters
    else:
        total_weight = (weight_columns + weight_rows) / 2

    # 3D distance
    value = (dist_3d_matrix * total_weight) * line_bc / 1000

    return value


def line_to_points(line, df):
    """
    Finds all the points of a linestring geodataframe correspondent to a
    point geodataframe.
    :param line: Linestring geodataframe being analyzed
    :param df: Point geodataframe where the linestring nodes will be referenced
    :return nodes_in_df: Point geodataframe containing all nodes of linestring
    """
    nodes = list(line.ID1.astype(int)) + list(line.ID2.astype(int))
    nodes = list(dict.fromkeys(nodes))
    nodes_in_df = gpd.GeoDataFrame(crs=df.crs, geometry=[])
    for i in nodes:
        nodes_in_df = nodes_in_df.append(df[df['ID'] == i], sort=False)
    nodes_in_df.reset_index(drop=True, inplace=True)

    return nodes_in_df


def create_roads(gdf_roads, geo_df):
    '''
    Creates two geodataframes
    :param gdf_roads: geodataframe with all roads in the area, as imported from OSM
    :param geo_df: geodataframe with the grid of points
    :return:line_gdf: point geodataframe containing verteces of the roads (in all the area)
            segments: geodataframe containing all the roads segments (in all the area)
    the GeoDataframes are also saves as shapefiles

    '''
    #w = geo_df.shape[0]
    if not geo_df.empty: #in case we are just processing the roads without pre-existing grid of points
        w = int(geo_df['ID'].max())+1 # this is because not all road points are included in the weighted grid of points. Basically,
    #it could be 10500,10501 and then 10504. df.shape[0] will give us a bad starting point in this case-> we want to start from 10505
    else:
        w=0
    line_vertices = pd.DataFrame(
        index=pd.Series(range(w, w + len(gdf_roads.index))),
        columns=['ID', 'X', 'Y', 'ID_line', 'Weight', 'Elevation'], dtype=int)
    # create geodataframe with all the segments that compose the road
    segments = gpd.GeoDataFrame(columns=['geometry', 'ID1', 'ID2'])
    k = 0
    gdf_roads.reset_index(inplace=True)
    x = 0
    for i, row in gdf_roads.iterrows():
        for j in list(row['geometry'].coords):
            line_vertices.loc[w, 'X'] = j[0]
            line_vertices.loc[w, 'Y'] = j[1]
            line_vertices.loc[w, 'ID_line'] = k
            line_vertices.loc[w, 'ID'] = w
            line_vertices.loc[w, 'Weight'] = 1
            w = w + 1
        k = k + 1
        points_to_split = MultiPoint(
            [Point(x, y) for x, y in row['geometry'].coords[1:]])
        splitted = split(row['geometry'], points_to_split)
        for j in splitted:
            segments.loc[x, 'geometry'] = j
            segments.loc[x, 'length'] = segments.loc[
                                            x, 'geometry'].length / 1000
            segments.loc[x, 'ID1'] = line_vertices[
                (line_vertices['X'] == j.coords[0][0]) & (
                        line_vertices['Y'] == j.coords[0][1])][
                'ID'].values[0]
            segments.loc[x, 'ID2'] = line_vertices[
                (line_vertices['X'] == j.coords[1][0]) & (
                        line_vertices['Y'] == j.coords[1][1])][
                'ID'].values[0]
            x = x + 1
        print('\r' + str(i) + '/' + str(gdf_roads.index.__len__()),
              sep=' ', end='', flush=True)
    if not geo_df.empty:
        line_vertices.loc[:, 'Elevation'] = geo_df[geo_df['Elevation']>0].Elevation.mean()
    else:
        line_vertices.loc[:, 'Elevation']=1000
        # line_vertices.loc[:, 'Elevation'] = 300
    geometry = [Point(xy) for xy in
                zip(line_vertices['X'], line_vertices['Y'])]
    line_gdf = gpd.GeoDataFrame(line_vertices, crs=geo_df.crs,
                                geometry=geometry)
    #line_gdf.to_file('Output/Datasets/Roads/gdf_roads.shp')
    #segments.to_file('Output/Datasets/Roads/roads_segments.shp')
    #segments.crs=22287
    # line_gdf.to_file('Testing_strategy/Points.shp')
    # segments.to_file('Testing_strategy/lines.shp')
    return line_gdf, segments

def create_roads2(gdf_roads, geo_df,crs):
    '''
    Creates two geodataframes
    :param gdf_roads: geodataframe with all roads in the area, as imported from OSM
    :param geo_df: geodataframe with the grid of points
    :return:line_gdf: point geodataframe containing verteces of the roads (in all the area)
            segments: geodataframe containing all the roads segments (in all the area)
    the GeoDataframes are also saves as shapefiles

    '''
    #w = geo_df.shape[0]
    if not geo_df.empty: #in case we are just processing the roads without pre-existing grid of points
        w = int(geo_df['ID'].max())+1 # this is because not all road points are included in the weighted grid of points. Basically,
    #it could be 10500,10501 and then 10504. df.shape[0] will give us a bad starting point in this case-> we want to start from 10505
    else:
        w=0
    line_vertices = pd.DataFrame(
        index=pd.Series(range(w, w + len(gdf_roads.index))),
        columns=['ID', 'X', 'Y', 'ID_line', 'Weight', 'Elevation'], dtype=int)
    # create geodataframe with all the segments that compose the road
    segments = gpd.GeoDataFrame(columns=['geometry', 'ID1', 'ID2'])
    k = 0
    gdf_roads.reset_index(inplace=True)
    x = 0
    for i, row in gdf_roads.iterrows():
        for j in list(row['geometry'].coords):
            if not (j[0] in line_vertices['X'].to_list() and j[1] in line_vertices['Y'].to_list()):
                line_vertices.loc[w, 'X'] = j[0]
                line_vertices.loc[w, 'Y'] = j[1]
                line_vertices.loc[w, 'ID_line'] = k
                line_vertices.loc[w, 'ID'] = w
                line_vertices.loc[w, 'Weight'] = 1
                w = w + 1
            else:
                pass
                #print('Double road point!')

        k = k + 1
        points_to_split = MultiPoint(
            [Point(x, y) for x, y in row['geometry'].coords[1:]])
        splitted = split(row['geometry'], points_to_split)
        for j in splitted:
            segments.loc[x, 'geometry'] = j
            segments.loc[x, 'length'] = segments.loc[
                                            x, 'geometry'].length / 1000
            segments.loc[x, 'ID1'] = line_vertices[
                (line_vertices['X'] == j.coords[0][0]) & (
                        line_vertices['Y'] == j.coords[0][1])][
                'ID'].values[0]
            segments.loc[x, 'ID2'] = line_vertices[
                (line_vertices['X'] == j.coords[1][0]) & (
                        line_vertices['Y'] == j.coords[1][1])][
                'ID'].values[0]
            x = x + 1
        print('\r' + str(i) + '/' + str(gdf_roads.index.__len__()),
              sep=' ', end='', flush=True)
    if not geo_df.empty:
        line_vertices.loc[:, 'Elevation'] = geo_df[geo_df['Elevation']>0].Elevation.mean()
    else:
        line_vertices.loc[:, 'Elevation']=1000
        # line_vertices.loc[:, 'Elevation'] = 300
    geometry = [Point(xy) for xy in
                zip(line_vertices['X'], line_vertices['Y'])]
    line_gdf = gpd.GeoDataFrame(line_vertices, crs=crs,
                                geometry=geometry)
    #line_gdf.to_file('Output/Datasets/Roads/gdf_roads.shp')
    #segments.to_file('Output/Datasets/Roads/roads_segments.shp')
    #segments.crs=22287
    # line_gdf.to_file('Testing_strategy/Points.shp')
    # segments.to_file('Testing_strategy/lines.shp')
    return line_gdf, segments

def create_box(limits, df):
    """
    Creates a delimiting box around a geodataframe.
    :param limits: Linestring geodataframe being analyzed
    :param df: Point geodataframe to be delimited
    :return df_box: All points of df that are inside the delimited box
    """
    x_min = min(limits.X)
    x_max = max(limits.X)
    y_min = min(limits.Y)
    y_max = max(limits.Y)

    dist = Point(x_min, y_min).distance(Point(x_max, y_max))
    if dist < 5000:
        extension = dist
    elif dist < 15000:
        extension = dist * 0.6
    else:
        extension = dist / 4

    bubble = box(minx=x_min - extension, maxx=x_max + extension,
                 miny=y_min - extension, maxy=y_max + extension)
    df_box = df[df.within(bubble)]
    df_box.index = pd.Series(range(0, len(df_box.index)))

    return df_box

def create_box(limits, df,resolution):
    """
    Creates a delimiting box around a geodataframe.
    :param limits: Linestring geodataframe being analyzed
    :param df: Point geodataframe to be delimited
    :return df_box: All points of df that are inside the delimited box
    """
    x_min = min(limits.X)
    x_max = max(limits.X)
    y_min = min(limits.Y)
    y_max = max(limits.Y)

    dist = Point(x_min, y_min).distance(Point(x_max, y_max))
    if dist < 5*resolution:
        extension = dist
    elif dist < 15*resolution:
        extension = dist * 0.6
    else:
        extension = dist / 4

    bubble = box(minx=x_min - extension, maxx=x_max + extension,
                 miny=y_min - extension, maxy=y_max + extension)
    df_box = df[df.within(bubble)]
    df_box.index = pd.Series(range(0, len(df_box.index)))

    return df_box

def edges_to_line(path, df, edges_matrix):
    """
    Transforms a list of NetworkX graph edges into a linestring geodataframe
    based on a input point geodataframe
    :param path: NetworkX graph edges sequence
    :param df: Point geodataframe to be used as reference
    :param edges_matrix: Matrix containing the cost to connect a pair of points
    :return line: Linestring geodataframe containing point IDs and its cost
    :return line_points: All points of df that are part of the linestring
    """
    steps = len(path)
    line = gpd.GeoDataFrame(index=range(0, steps),
                            columns=['ID1', 'ID2', 'Cost', 'geometry'],
                            crs=df.crs)
    line_points = []
    for h in range(0, steps):
        line.at[h, 'geometry'] = LineString(
            [(df.loc[path[h][0], 'X'],
              df.loc[path[h][0], 'Y']),
             (df.loc[path[h][1], 'X'],
              df.loc[path[h][1], 'Y'])])

        # int here is necessary to use the command .to_file
        line.at[h, 'ID1'] = int(df.loc[path[h][0], 'ID'])
        line.at[h, 'ID2'] = int(df.loc[path[h][1], 'ID'])
        line.at[h, 'Cost'] = int(edges_matrix.loc[df.loc[path[h][0], 'ID'],
                                                  df.loc[path[h][1], 'ID']])
        line_points.append(list(df.loc[path[h], 'ID']))
    line.drop(line[line['Cost'] == 0].index, inplace=True)
    line.Cost = line.Cost.astype(int)
    return line, line_points


def load(clusters_list, grid_lifetime, input_profile,gisele_folder, case_study):
    """
    Reads the input daily load profile from the input csv. Reads the number of
    years of the project and the demand growth from the data.dat file of
    Micehele. Then it multiplies the load profile by the Clusters' peak load
    and append values to create yearly profile composed of 12 representative
    days.
    :param grid_lifetime: Number of years the grid will operate
    :param clusters_list: List of clusters ID numbers
    :return load_profile: Cluster load profile for the whole period
    :return years: Number of years the microgrid will operate
    :return total_energy: Energy provided by the grid in its lifetime [kWh]
    """
    l()
    print("5. Microgrid Sizing")
    l()
    case_folder = gisele_folder + '/Case studies/' + case_study

    data_michele = pd.read_table(gisele_folder+"/gisele/michele/Inputs/data.dat", sep="=",
                                 header=None)
    print("Creating load profile for each cluster..")
    daily_profile = pd.DataFrame(index=range(1, 25),
                                 columns=clusters_list.Cluster)
    for column in daily_profile:
        daily_profile.loc[:, column] = \
            (input_profile.loc[:, 'Hourly Factor']
             * float(clusters_list.loc[clusters_list['Cluster']==column, 'Load [kW]'])).values
    rep_days = int(data_michele.loc[0, 1].split(';')[0])
    grid_energy = daily_profile.append([daily_profile] * 364,
                                       ignore_index=True)
    #  append 11 times since we are using 12 representative days in a year
    load_profile = daily_profile.append([daily_profile] * (rep_days - 1),
                                        ignore_index=True)

    years = int(data_michele.loc[1, 1].split(';')[0])
    demand_growth = float(data_michele.loc[87, 1].split(';')[0])
    daily_profile_new = daily_profile
    #  appending for all the years considering demand growth
    for i in range(grid_lifetime - 1):
        daily_profile_new = daily_profile_new.multiply(1 + demand_growth)
        if i < (years - 1):
            load_profile = load_profile.append([daily_profile_new] * rep_days,
                                               ignore_index=True)
        grid_energy = grid_energy.append([daily_profile_new] * 365,
                                         ignore_index=True)
    total_energy = pd.DataFrame(index=clusters_list.Cluster,
                                columns=['Energy'])
    for cluster in clusters_list.Cluster:
        total_energy.loc[cluster, 'Energy'] = \
            grid_energy.loc[:, cluster].sum().round(2)
    print("Load profile created")
    total_energy.to_csv(case_folder +'/Intermediate/Microgrid/Grid_energy.csv')
    return load_profile, years, total_energy


def shift_timezone(df, shift):
    """
    Move the values of a dataframe with DateTimeIndex to another UTC zone,
    adding or removing hours.
    :param df: Dataframe to be analyzed
    :param shift: Amount of hours to be shifted
    :return df: Input dataframe with values shifted in time
    """
    if shift > 0:
        add_hours = df.tail(shift)
        df = pd.concat([add_hours, df], ignore_index=True)
        df.drop(df.tail(shift).index, inplace=True)
    elif shift < 0:
        remove_hours = df.head(abs(shift))
        df = pd.concat([df, remove_hours], ignore_index=True)
        df.drop(df.head(abs(shift)).index, inplace=True)
    return df


def sizing(load_profile, clusters_list, geo_df_clustered, wt, mg_types, gisele_folder,case_study):
    """
    Imports the solar and wind production from the RenewablesNinja api and then
    Runs the optimization algorithm MicHEle to find the best microgrid
    configuration for each Cluster.
    :param load_profile: Load profile of all clusters during all years
    :param clusters_list: List of clusters ID numbers
    :param geo_df_clustered: Point geodataframe with Cluster identification
    :param wt: Wind turbine model used for computing the wind velocity
    :param mg_types: number of times to evaluate microgrids in each cluster.
                renewables fraction in michele changes accordingly
    :return mg: Dataframe containing the information of the Clusters' microgrid
    """
    speed_up = True
    case_folder = gisele_folder + '/Case studies/' + case_study
    geo_df_clustered = geo_df_clustered.to_crs(4326)
    mg = {}
    # mg = pd.DataFrame(index=clusters_list.index,
    #                   columns=['Cluster','PV [kW]', 'Wind [kW]', 'Hydro [kW]'
    #                            'Diesel [kW]',
    #                            'BESS [kWh]', 'Inverter [kW]',
    #                            'Investment Cost [kEUR]', 'OM Cost [kEUR]',
    #                            'Replace Cost [kEUR]', 'Total Cost [kEUR]',
    #                            'Energy Demand [MWh]', 'Energy Produced [MWh]',
    #                            'LCOE [EUR/kWh]','CO2 [kg]', 'Unavailability [MWh/y]'],
    #                   dtype=float)

    for i in range(mg_types):
        mg[i] = pd.DataFrame(index=clusters_list.index,
                          columns=['Cluster','Renewable fraction index', 'PV [kW]', 'Wind [kW]', 'Diesel [kW]',
                                   'BESS [kWh]', 'Inverter [kW]',
                                   'Investment Cost [kEUR]', 'OM Cost [kEUR]',
                                   'Replace Cost [kEUR]', 'Total Cost [kEUR]',
                                   'Energy Demand [MWh]', 'Energy Produced [MWh]',
                                   'LCOE [EUR/kWh]','CO2 [kg]', 'Unavailability [MWh/y]'],
                          dtype=float)

    #save useful values from michele input data
    with open(gisele_folder+'/gisele/michele/Inputs/data.json') as f:
        input_michele = json.load(f)
    proj_lifetime = input_michele['num_years']
    num_typ_days = input_michele['num_days']
    clusters = clusters_list.Cluster
    for index in range(len(clusters)):
        # try:
        cluster_n = clusters[index]
        l()
        print('Creating the optimal Microgrid for Cluster ' + str(cluster_n))
        l()
        load_profile_cluster = load_profile.loc[:, cluster_n]
        lat = geo_df_clustered[geo_df_clustered['Cluster']
                             == cluster_n].geometry.y.values[0]
        lon = geo_df_clustered[geo_df_clustered['Cluster']
                             == cluster_n].geometry.x.values[0]
        all_angles = pd.read_csv(gisele_folder+'/general_input/TiltAngles.csv')
        tilt_angle = abs(all_angles.loc[abs(all_angles['lat'] - lat).idxmin(),
                                      'opt_tilt'])
        if (speed_up==True and index==0) or speed_up==False:
            pv_prod = import_pv_data(lat, lon, tilt_angle)
            wt_prod = import_wind_data(lat, lon, wt)
        utc = pv_prod.local_time[0]
        if type(utc) is pd.Timestamp:
          time_shift = utc.hour
        else:
          utc = iso8601.parse_date(utc)
          time_shift = int(utc.tzinfo.tzname(utc).split(':')[0])
        div_round = 8760 // (num_typ_days * 24)
        length = num_typ_days * 24
        new_length= length *div_round
        # pv_avg = pv_prod.groupby([pv_prod.index.month,
        #                           pv_prod.index.hour]).mean()

        pv_avg_new=np.zeros(24*num_typ_days)
        pv_avg = pv_prod.values[0:new_length,1].reshape(24,div_round*num_typ_days,order='F')
        wt_avg_new = np.zeros(24 * num_typ_days)
        wt_avg = wt_prod.values[0:new_length].reshape(24,
                                                       div_round * num_typ_days,
                                                       order='F')
        for i in range(num_typ_days):
          pv_avg_new[i*24:(i+1)*24] = pv_avg[:,div_round*i:div_round*(i+1)].mean(axis=1)

          wt_avg_new[i * 24:(i + 1) * 24] = wt_avg[:, div_round * i:div_round * (
                  i + 1)].mean(axis=1)



        pv_avg = pd.DataFrame(pv_avg_new)
        pv_avg = pv_avg.append([pv_avg] * (proj_lifetime - 1), ignore_index=True)
        pv_avg.reset_index(drop=True, inplace=True)
        pv_avg = shift_timezone(pv_avg, time_shift)


        # wt_prod = import_wind_data(lat, lon, wt)
        # wt_avg = wt_prod.groupby([wt_prod.index.month,
        #                           wt_prod.index.hour]).mean()
        wt_avg = pd.DataFrame(wt_avg_new)
        wt_avg = wt_avg.append([wt_avg] * (proj_lifetime - 1), ignore_index=True)
        wt_avg.reset_index(drop=True, inplace=True)
        wt_avg = shift_timezone(wt_avg, time_shift)

        #todo ->implement hydro resource, for the moment creation of a fake input
        ht_avg = wt_avg

        results = start(load_profile_cluster, pv_avg, wt_avg,input_michele, ht_avg, mg_types)

        for i in range(mg_types):
          mg[i].loc[index, 'Cluster'] = 'C' + str(cluster_n)
          mg[i].loc[index, 'Renewable fraction index'] = str(i)
          mg[i].loc[index, 'PV [kW]'] = results[str(i)]['inst_pv']
          mg[i].loc[index, 'Wind [kW]'] = results[str(i)]['inst_wind']
          mg[i].loc[index, 'Diesel [kW]'] = results[str(i)]['inst_dg']
          mg[i].loc[index, 'BESS [kWh]'] = results[str(i)]['inst_bess']
          mg[i].loc[index, 'Inverter [kW]'] = results[str(i)]['inst_inv']
          mg[i].loc[index, 'Investment Cost [kEUR]'] = results[str(i)]['init_cost']
          mg[i].loc[index, 'OM Cost [kEUR]'] = results[str(i)]['om_cost']
          mg[i].loc[index, 'Replace Cost [kEUR]'] = results[str(i)]['rep_cost']
          mg[i].loc[index, 'Total Cost [kEUR]'] = results[str(i)]['npc']
          mg[i].loc[index, 'Energy Produced [MWh]'] = results[str(i)]['gen_energy']
          mg[i].loc[index, 'Energy Demand [MWh]'] = results[str(i)]['load_energy']
          mg[i].loc[index, 'LCOE [EUR/kWh]'] = results[str(i)]['npc'] / \
                                                   results[str(i)]['gen_energy']
          mg[i].loc[index, 'CO2 [kg]'] = results[str(i)]['emissions']
          mg[i].loc[index, 'Unavailability [MWh/y]'] = results[str(i)]['tot_unav']
          print(mg)
        # except:
        #   print('Region too large to compute the optimal microgrid.')

    microgrid = pd.DataFrame()
    for i in range(mg_types):
        microgrid = microgrid.append(mg[i].round(decimals=4))
    microgrid.to_csv(case_folder+'/Intermediate/Microgrid/microgrids.csv', index=False)

    return microgrid


def download_url(url, out_path, chunk_size=128):
    """
    Download zip file from specified url and save it into a defined folder
    Note: check chunk_size parameter
    """
    r = requests.get(url, stream=True)
    with open(out_path, 'wb') as fd:
        for chunk in r.iter_content(chunk_size=chunk_size):
            fd.write(chunk)


def download_tif(area, crs, scale, image, out_path):
    """
    Download data from Earth Engine
    :param area: GeoDataFrame with the polygon of interest area
    :param crs: str with crs of the project
    :param scale: int with pixel size in meters
    :param image: image from the wanted database in Earth Image
    :param out_path: str with output path
    :return:
    """
    min_x, min_y, max_x, max_y = area.geometry.total_bounds
    path = image.getDownloadUrl({
        'scale': scale,
        'crs': 'EPSG:' + str(crs),
        'region': [[min_x, min_y], [min_x, max_y], [max_x, min_y],
                   [max_x, max_y]]
    })
    print(path)
    download_url(path, out_path)
    return
def MultiLine_to_Line(lines_shapefile):
    lines=[]
    for index, row in lines_shapefile.iterrows():
        try:
            row['geometry'][0]
            multi_line = row['geometry']
            lines_shapefile.drop(index,axis=0,inplace=True)
            for i in multi_line:
                lines.append(i)
        except:
            a=1
    lines_shapefile=lines_shapefile.append(gpd.GeoDataFrame({'geometry': lines}))
    lines_shapefile = lines_shapefile.reset_index(drop=True)
    return lines_shapefile
# def lcoe_analysis(clusters_list, total_energy, grid_resume, mg, coe,
#                   grid_ir, grid_om, grid_lifetime):
#     """
#      Computes the LCOE of the on-grid and off-grid and compares both of them
#      to find the best solution.
#      :param clusters_list: List of clusters ID numbers
#      :param total_energy: Energy provided by the grid in its lifetime [kWh]
#      :param grid_resume: Table summarizing grid and connection costs
#      :param mg: Table summarizing all the microgrids and its costs.
#      :param coe: Cost of electricity in the market [€/kWh]
#      :param grid_ir: Inflation rate for the grid investments [%]
#      :param grid_om: Operation and maintenance cost of grid [% of invest cost]
#      :param grid_lifetime: Number of years the grid will operate
#      :return final_lcoe: Table summary of all LCOEs and the proposed solution
#      """
#
#     final_lcoe = pd.DataFrame(index=clusters_list.Cluster,
#                               columns=['Grid NPC [k€]', 'MG NPC [k€]',
#                                        'Grid Energy Consumption [MWh]',
#                                        'MG LCOE [€/kWh]',
#                                        'Grid LCOE [€/kWh]'],
#                               dtype=float)
#     for i in clusters_list.Cluster:
#         final_lcoe.at[i, 'Grid Energy Consumption [MWh]'] = \
#             total_energy.loc[i, 'Energy'] / 1000
#
#         # finding the npv of the cost of O&M for the whole grid lifetime
#
#         total_grid_om = grid_om * (grid_resume.loc[i, 'Grid Cost [k€]'] +
#                                    grid_resume.loc
#                                    [i, 'Connection Cost [k€]']) * 1000
#         total_grid_om = [total_grid_om] * grid_lifetime
#         total_grid_om = np.npv(grid_ir, total_grid_om)
#
#         cluster_grid_om = grid_om * grid_resume.loc[i, 'Grid Cost [k€]'] * 1000
#         cluster_grid_om = [cluster_grid_om] * grid_lifetime
#         cluster_grid_om = np.npv(grid_ir, cluster_grid_om)
#         cluster_grid_npc = cluster_grid_om / 1000 + \
#                            grid_resume.loc[i, 'Grid Cost [k€]']
#         cluster_grid_lcoe = \
#             cluster_grid_npc / final_lcoe.loc[
#                 i, 'Grid Energy Consumption [MWh]']
#         final_lcoe.at[i, 'Grid NPC [k€]'] = \
#             (grid_resume.loc[i, 'Grid Cost [k€]'] +
#              grid_resume.loc[i, 'Connection Cost [k€]']) \
#             + total_grid_om / 1000
#
#         final_lcoe.at[i, 'MG NPC [k€]'] = mg.loc[i, 'Total Cost [k€]']
#         final_lcoe.at[i, 'MG LCOE [€/kWh]'] = mg.loc[i, 'LCOE [€/kWh]'] + \
#                                               cluster_grid_lcoe
#
#         grid_lcoe = final_lcoe.loc[i, 'Grid NPC [k€]'] / final_lcoe.loc[
#             i, 'Grid Energy Consumption [MWh]'] + coe
#
#         final_lcoe.at[i, 'Grid LCOE [€/kWh]'] = grid_lcoe
#
#         if mg.loc[i, 'LCOE [€/kWh]'] + cluster_grid_lcoe > grid_lcoe:
#             ratio = grid_lcoe / (mg.loc[i, 'LCOE [€/kWh]'] + cluster_grid_lcoe)
#             if ratio < 0.95:
#                 proposal = 'ON-GRID'
#             else:
#                 proposal = 'BOTH'
#
#         else:
#             ratio = (mg.loc[i, 'LCOE [€/kWh]'] + cluster_grid_lcoe) / grid_lcoe
#             if ratio < 0.95:
#                 proposal = 'OFF-GRID'
#             else:
#                 proposal = 'BOTH'
#
#         final_lcoe.at[i, 'Best Solution'] = proposal
#     final_lcoe = final_lcoe.round(decimals=4)
#     l()
#     print(final_lcoe)
#
#     final_lcoe.to_csv(r'Output/Microgrids/LCOE_Analysis.csv')
#
#     return final_lcoe
