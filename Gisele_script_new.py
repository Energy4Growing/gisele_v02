import os
import base64
import io
import pandas as pd
import geopandas as gpd
from shapely.geometry import Polygon
import numpy as np
import plotly.graph_objs as go
from shapely.geometry import Point
from gisele.functions import *
from gisele.functions2 import *
from gisele import initialization, clustering, processing, collecting, \
    optimization, results, grid, branches
import pyutilib.subprocess.GlobalData
from gisele import QGIS_processing_polygon as qgis_process
from gisele import Local_area_optimization as LAO
from gisele import MILP_Input_creation,MILP_models,process_output,grid_routing,Secondary_substations
import time
import rasterio
############# INPUT ELECTRICAL PARAMETERS #############

# parameters for the MV cable that is used
# generic cable
# resistance = 0.42  # [ohm/km]
# reactance = 0.39  # [ohm/km]
# Pmax = 11  # [MVA]

#cable of 48mm2
resistance = 0.321  # [ohm/km]
reactance = 0.372  # [ohm/km]
Pmax = 4  # [MVA]
line_cost = 10000  # Cost of MV voltage feeders [€/km]
#second cable - this cable is used only in case the MILP option with 2 different lines is chosen
resistance2 = 0.35
reactance2 = 0.38
Pmax2 = 2
line_cost2 = 8000 # Cost of MV voltage feeders [€/km]

# #cable of 70 mm2
# resistance = 0.413  # [ohm/km]
# reactance = 0.36  # [ohm/km]
# Pmax = 5.5  # [MVA]

voltage = 11  # [kV]

LV_base_cost=10000 # for the LV we are not checking electrical parameters

###########STARTING THE SCRIPT###########

gisele_folder=os.getcwd()
villages_file = 'Villages_areas.geojson'
villages_file = 'test_large_village/test_large_village.shp'
local_database =False
country  = 'Lesotho'
case_study='lesotho_poster'

crs = 22287
resolution = 240
load_capita=0.7 #kW
pop_per_household=4
resolution_population = 30 # this should be automatic from the raster
#data used for the local area optimization
max_length_segment = resolution_population*1.5
#simplify_coef = 5
#crit_dist = simplify_coef/2


if local_database ==False:
    #database= r'C:\Users\silvi\OneDrive - Politecnico di Milano\Documents\2020-2021\Gisele shared\8.Case_Study'
    database= r'C:\Users\alekd\Politecnico di Milano\Silvia Corigliano - Gisele shared\8.Case_Study'
    cluster_folder = database+'\Lesotho\Villages_areas_46.geojson'
    #cluster_folder = database + '\Lesotho/test_large_village/test_large_village.shp'
    substations_folder = database+'\Lesotho\con_point.shp'
    study_area_folder = database + '/' + country + '/Study_area/Study_area_big.shp'
else:
    database =gisele_folder + '/Database'
    cluster_folder = database + '/' + country + '/'+villages_file
    substations_folder = database + '/' + country +'/con_point.shp'
    study_area_folder = database + '/' + country + '/Study_area/Study_area_test1.shp'





#for ind,rows in Clusters.iterrows():
#    Clusters.loc[ind,'geometry']  = rows['geometry'].buffer(resolution_population*0.9)




if not os.path.exists(r'Case studies/'+case_study): # if this is a new project, create the starting point for the analysis
    # Create new folders for the study case
    os.makedirs(r'Case studies/'+case_study)
    os.makedirs(r'Case studies/'+case_study+'/Input')
    os.makedirs(r'Case studies/'+case_study+'/Output')
    os.makedirs(r'Case studies/'+case_study+'/Intermediate')
    os.makedirs(r'Case studies/' + case_study + '/Intermediate/Communities  ')
    os.makedirs(r'Case studies/' + case_study + '/Intermediate/Microgrid')
    os.makedirs(r'Case studies/'+case_study+'/Intermediate/Optimization')
    os.makedirs(r'Case studies/' + case_study + '/Intermediate/Geospatial_Data')
    os.makedirs(r'Case studies/' + case_study + '/Intermediate/Optimization/MILP_output')
    os.makedirs(r'Case studies/' + case_study + '/Output/MILP_processed')
    os.makedirs(r'Case studies/' + case_study + '/Intermediate/Optimization/all_data')
    os.makedirs(r'Case studies/' + case_study + '/Intermediate/Optimization/MILP_input')
    os.makedirs(r'Case studies/' + case_study + '/Intermediate/Optimization/all_data/Lines_connections')
    os.makedirs(r'Case studies/' + case_study +'/Intermediate/Optimization/all_data/Lines_marked')
    # Copy the Configuration file from the general input
    pd.read_csv(r'general_input/Configuration.csv').to_csv(r'Case studies/'+case_study+'/Input/Configuration.csv')
    # Read the possible connection points and write them in the case study's folder
    Substations = gpd.read_file(substations_folder)
    Substations_crs = Substations.to_crs(crs)
    Substations['X'] = [Substations_crs['geometry'].values[i].xy[0][0] for i in range(Substations.shape[0])]
    Substations['Y'] = [Substations_crs['geometry'].values[i].xy[1][0] for i in range(Substations.shape[0])]
    Substations.to_file(r'Case studies/' + case_study + '/Input/substations')
    # Read the polygon of the study area and write it in the local database.
    study_area = gpd.read_file(study_area_folder)
    study_area.to_file(r'Case studies/'+case_study+'/Input/Study_area')
    # Read the communities and write them in the local database
    Clusters = gpd.read_file(cluster_folder)
    Clusters = Clusters.to_crs(crs)
    Clusters['cluster_ID'] = [*range(1, Clusters.shape[0] + 1)]
    for i, row in Clusters.iterrows(): # this is just in case one of the polygons is saved as a MP with just 1 polygon
        if row['geometry'].geom_type == 'MultiPolygon':
            Clusters.loc[i, 'geometry'] = row['geometry'][0]
    Clusters.to_file(r'Case studies/'+case_study+'/Input/Communities_boundaries')
else: # not a new project, just read the files from the local folder
    Clusters = gpd.read_file(r'Case studies/'+case_study+'/Input/Communities_boundaries/Communities_boundaries.shp')
    study_area = gpd.read_file(r'Case studies/'+case_study+'/Input/Study_area/Study_area.shp')
    Substations = gpd.read_file(r'Case studies/' + case_study + '/Input/substations/substations.shp')
# specific data on the case study
LV_distance=500 # Maximum length of the LV network.
ss_data = 'ss_data_evn.csv' # folder in which the costs for substations can be found
simplify_road_coef_inside = 5 # in meters, used for the routing inside the clusters.
simplify_road_coef_outside = 30 # in meters, used for creating connections among clusters/substations.
road_coef = 2
roads_weight=0.3
### USER OPTIONS
mg_option = True
mg_types =1 #if more than one mg for each cluster needs to be computed, with different reliability levels (needed for multiobjective: mg_types=3)
multi_objective_option = False
reliability_option = False
Roads_option=True # Boolean stating whether we want to use the roads for the steiner tree and dijkstra.
Rivers_option=False
n_line_type=1
run_genetic=False

# parameters for the economical factors
coe = 100# euro/MWh of electrical energy supplied
grid_lifetime = 40 #years
landcover_option='ESACCI'

# NEW - suggest connection points
'''Create the grid of points'''
print('1. CREATE A WEIGHTED GRID OF POINTS')
# df_weighted = qgis_process.create_input_csv(crs,resolution,resolution_population,landcover_option,country,case_study,database,study_area)
# accepted_road_types = ['living_street', 'pedestrian', 'primary', 'primary_link', 'secondary', 'secondary_link',
#                           'tertiary', 'tertiary_link', 'unclassified','residential']
# Road_nodes,Road_lines = create_roads_new(gisele_folder, case_study, Clusters,crs, accepted_road_types,resolution,resolution_population)
# Merge_Roads_GridOfPoints(gisele_folder,case_study)
''' CLUSTERING PROCEDURE'''
'''For each cluster, perform further aglomerative clustering, locate secondary substations and perform MV grid routing'''
print('2. LOCATE SECONDARY SUBSTATIONS INSIDE THE CLUSTERS.')
#Secondary_substations.locate_secondary_ss(crs, resolution_population, load_capita, pop_per_household, road_coef,
#                       Clusters, case_study, LV_distance, ss_data,landcover_option, gisele_folder)
#start = time.time()
#LAO.optimize(crs, resolution_population, load_capita, pop_per_household, road_coef, Clusters, case_study, LV_distance,
#            ss_data,landcover_option,gisele_folder,roads_weight,run_genetic, max_length_segment,simplify_coef, crit_dist,LV_base_cost)
#end = time.time()
#print(end-start)
print('3. ROUTE THE LV GRID FOR EACH CLUSTER')
# grid_routing.LV_routing(gisele_folder,case_study,crs,resolution_population,Roads_option,simplify_road_coef_inside,LV_base_cost,roads_weight)
print('4. ROUTE THE MV GRID INSIDE THE CLUSTERS.')
# grid_routing.routing(Clusters,gisele_folder,case_study,crs,resolution,Roads_option,roads_weight,simplify_road_coef_inside,Rivers_option,line_cost)
'''Create the input for the MILP'''
print('5. Create the input for the MILP.')
MILP_Input_creation.create_input(gisele_folder,case_study,crs,line_cost,resolution,reliability_option,Roads_option,
                         simplify_road_coef_outside, Rivers_option,mg_option,mg_types)

if multi_objective_option:
    calculate_mg_multiobjective(gisele_folder, case_study, crs)
elif mg_option:
     calculate_mg(gisele_folder, case_study, crs,mg_types)
if reliability_option:
    MILP_Input_creation.create_MILP_input(gisele_folder,case_study,crs,mg_option) # all the lines are from i-j and j-i, only positive powers to consider reliability
else:
    MILP_Input_creation.create_MILP_input_1way(gisele_folder,case_study,crs) # only i-j, without reliability


# '''Execute the desired MILP model'''
# # print('6. Execute the MILP according to the selected options.')
# n_clusters = Clusters.shape[0]
# start = time.time()
# if mg_option == False and reliability_option==True and n_line_type==1:
#     MILP_models.MILP_without_MG(gisele_folder,case_study,n_clusters,coe,voltage,resistance,reactance,Pmax,line_cost)
# elif mg_option == True and reliability_option==False and n_line_type==1:
#     MILP_models.MILP_MG_noRel(gisele_folder, case_study, n_clusters, coe, voltage, resistance, reactance, Pmax, line_cost)
# elif mg_option == True and reliability_option==False and n_line_type==2:
#     MILP_models.MILP_MG_2cables(gisele_folder, case_study, n_clusters, coe, voltage, resistance, reactance, Pmax, line_cost
#                                  ,resistance2,reactance2,Pmax2,line_cost2)
#
# # '''Process the output from the MILP'''
# # print('7. Process MILP output')
# process_output.process(gisele_folder,case_study,crs,mg_option,reliability_option)
# process_output.create_final_output(gisele_folder, case_study)
#
# end = time.time()
# print(end-start)