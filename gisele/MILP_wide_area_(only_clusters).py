from math import ceil
from gisele import dijkstra
from shapely.ops import nearest_points
from gisele.functions2 import *
from scipy.spatial import Delaunay
from shapely.geometry import Point
from gisele.Local_area_optimization import *
from gisele.MILP_Input_creation import *
from gisele.Secondary_substations import *
from math import *
from __future__ import division
from shapely import wkt

from pyomo.opt import SolverFactory
from pyomo.core import AbstractModel
from pyomo.dataportal.DataPortal import DataPortal
from pyomo.environ import *
import pandas as pd
from datetime import datetime
import os
def input(clusters,gisele_folder,case_study,crs,resolution_LV,load_capita,pop_per_household,triangulation_logic,
          line_bc,resolution,mg_option):
    case_folder = gisele_folder + '/Case studies/' + case_study
    data_folder = case_folder + '/Intermediate/Optimization/all_data'
    Roads_option=True
    Rivers_option = False
    tolerance_outside=30

    crs_str = "EPSG:" + str(crs)
    dir_input = gisele_folder+'/Case studies/' + case_study + '/Intermediate/Geospatial_Data'
    case_folder = gisele_folder + '/Case studies/' + case_study
    all_points = pd.read_csv(dir_input + '/weighted_grid_of_points_with_roads.csv')

    geometry = [Point(xy) for xy in zip(all_points.X, all_points.Y)]

    # This here is to make sure that the points of the internal MV grid are included in all_points
    all_points = gpd.GeoDataFrame(all_points, crs=crs_str, geometry=geometry)
    next_ID = int(all_points['ID'].max() + 1)

    # Create the file All_Nodes.csv
    Population = rasterio.open(dir_input + '/Population_' + str(crs) + '.tif')
    All_nodes = pd.DataFrame()
    for i,row in clusters.iterrows():
        geometry = row['geometry'].centroid
        grid_of_points = create_grid(crs, resolution_LV, row['geometry'])
        coords = [(x, y) for x, y in zip(grid_of_points.X, grid_of_points.Y)]
        grid_of_points['Population'] = [x[0] for x in Population.sample(coords)]
        population = grid_of_points['Population'].sum()
        cluster = row['cluster_ID']
        mv_load = population * load_capita * coincidence_factor(population, pop_per_household)
        All_nodes=All_nodes.append(pd.DataFrame({'ID':[next_ID],'Weight':[1],'Cluster':[cluster],'MV_Power':[mv_load],
                                                 'Population':[population],'Substation':[1],'geometry':[geometry],'Type':['Cluster']}))
        next_ID+=1

    All_nodes = gpd.GeoDataFrame(All_nodes, crs=crs_str, geometry=All_nodes.geometry)
    All_nodes['X'] = All_nodes.geometry.x
    All_nodes['Y'] = All_nodes.geometry.y
    All_nodes_with_elev=All_nodes.copy()
    All_nodes_with_elev['Elevation']=1000
    all_points_with_clusters = all_points.append(All_nodes_with_elev)
    all_points_with_clusters.reset_index(inplace=True)
    all_points_with_clusters.to_csv(
        case_folder + '/Intermediate/Geospatial_Data/weighted_grid_of_points_with_ss_and_roads.csv',index=False)


    Primary_substations = gpd.read_file(
        gisele_folder + '/Case studies/' + case_study + '/Input/substations/substations.shp')
    Primary_substations = Primary_substations.to_crs(crs)
    Primary_substations['ID'] = [*range(next_ID, next_ID + Primary_substations.shape[0])]
    Primary_substations.to_file(
        gisele_folder + '/Case studies/' + case_study + '/Input/substations/substations.shp')

    all_points_withPS = add_PS_to_grid_of_points(all_points_with_clusters, Primary_substations, next_ID, add_elevation=True)

    All_nodes_withPS = add_PS_to_grid_of_points(All_nodes, Primary_substations, next_ID, add_elevation=False)
    All_nodes_withPS.to_csv(case_folder + '/Intermediate/Optimization/all_data/All_Nodes.csv', index=False)

    if triangulation_logic:
        allowed_cluster_connections = delaunay_for_MILP_input_creation(gisele_folder,case_study,crs,All_nodes_withPS)
    else:
        allowed_cluster_connections=[]

    #Create connections
    Roads_points_all = gpd.read_file(case_folder + '/Intermediate/GeoSpatial_Data/Roads_points/Roads_points.shp')
    Roads_lines_all = gpd.read_file(case_folder + '/Intermediate/GeoSpatial_Data/Roads_lines/Roads_lines.shp')
    Lines_connections=pd.DataFrame()
    n_clusters = int(All_nodes['Cluster'].max())
    distance_limit_clusters = 30 * resolution
    c_grid_points=[] # no grid inside clusters
    for i in range(1,n_clusters+1):
        for j in range(i + 1, n_clusters+1, 1):
            if not triangulation_logic or (i,j) in allowed_cluster_connections:
                print('Creating the connection C'+str(i)+' and C' +str(j))
                p1 = All_nodes[All_nodes['Cluster'] == i]
                p2 = All_nodes[All_nodes['Cluster'] == j]
                #include roads here
                dist = p1['geometry'].values[0].distance(p2['geometry'].values[0])
                if Roads_option and dist< 30*resolution:

                    if dist < 1000:
                        extension = dist
                    elif dist < 5000:
                        extension = dist * 0.6
                    else:
                        extension = dist / 4
                    x_min = min(p1['X'].values[0],p2['X'].values[0])
                    x_max = max(p1['X'].values[0],p2['X'].values[0])
                    y_min = min(p1['Y'].values[0],p2['Y'].values[0])
                    y_max = max(p1['Y'].values[0],p2['Y'].values[0])
                    bubble = box(minx=x_min - extension, maxx=x_max + extension,
                                 miny=y_min - extension, maxy=y_max + extension)

                    gdf_roads = gpd.clip(Roads_points_all, bubble)
                    roads_segments = Roads_lines_all [(Roads_lines_all ['ID1'].isin(gdf_roads.ID.to_list()) &
                                             Roads_lines_all ['ID2'].isin(gdf_roads.ID.to_list()))]
                    connection, connection_cost, connection_length, pts = dijkstra.dijkstra_connection_roads_new(all_points_withPS, p1,
                                p2, c_grid_points,    line_bc,  resolution,gdf_roads, roads_segments,Rivers_option,1,1.5,distance_limit_clusters )
                if not connection_cost == 999999 :
                    Line = array_to_LineString(all_points_withPS, pts)
                    Data = {'ID1': pts[0][0], 'ID2': pts[len(pts)-1][1], 'C1': i, 'C2': j, 'Cost': connection_cost,
                            'Length': connection_length, 'geometry': [Line], 'Type': 'C-C'}
                    Lines_connections = Lines_connections.append(pd.DataFrame(Data))

    # Cluster to substation
    for i in range(1,n_clusters+1):
        cluster_nodes = All_nodes[All_nodes['Cluster'] == i]
        points_list = cluster_nodes.loc[cluster_nodes['ID'] > 0, 'geometry'].to_list()
        # points_list =  all_points['geometry'].to_list()
        grid_multipoint = MultiPoint(points_list)
        print('S-C '+str(i))
        for index, row in Primary_substations.iterrows():
            if not triangulation_logic or (i,str(row['ID'])) in allowed_cluster_connections:
                sub_node = pd.DataFrame({'X': [row['X']], 'Y': [row['Y']],
                                         'ID': [row['ID']]})  # this is to fix an issue with this being data series
                dist_matrix = pd.DataFrame(distance_2d(cluster_nodes, sub_node, 'X', 'Y'),
                                           index=cluster_nodes.ID.values,
                                           columns=[Primary_substations.ID.values[index]])
                p1 = All_nodes_withPS[All_nodes_withPS['ID'] == dist_matrix.min().idxmin()]
                p2 = All_nodes_withPS[All_nodes_withPS['ID'] ==
                               dist_matrix.min(axis=1).idxmin()]
                print(p1['ID'])
                print(p2['ID'])
                dist = p1['geometry'].values[0].distance(p2['geometry'].values[0])
                if Roads_option and dist<50*resolution:

                    if dist < 1000:
                        extension = dist
                    elif dist < 5000:
                        extension = dist * 0.6
                    else:
                        extension = dist / 4
                    x_min = min(p1['X'].values[0],p2['X'].values[0])
                    x_max = max(p1['X'].values[0],p2['X'].values[0])
                    y_min = min(p1['Y'].values[0],p2['Y'].values[0])
                    y_max = max(p1['Y'].values[0],p2['Y'].values[0])
                    bubble = box(minx=x_min - extension, maxx=x_max + extension,
                                 miny=y_min - extension, maxy=y_max + extension)
                    gdf_roads = gpd.clip(Roads_points_all, bubble)
                    roads_segments = Roads_lines_all[(Roads_lines_all['ID1'].isin(gdf_roads.ID.to_list()) &
                                                      Roads_lines_all['ID2'].isin(gdf_roads.ID.to_list()))]

                    connection, connection_cost, connection_length, pts = dijkstra.dijkstra_connection_roads_new(all_points_withPS, p1,
                                p2, c_grid_points,    line_bc,  resolution,gdf_roads, roads_segments,Rivers_option,1,1.5,resolution * 25 )

                    if not connection_cost == 999999:

                        Line = array_to_LineString(all_points_withPS, pts)
                        Data = {'ID1': pts[0][0], 'ID2': pts[len(pts) - 1][1], 'C1': i, 'C2':  'S' + str(index),
                                'Cost': connection_cost,
                                'Length': connection_length, 'geometry': [Line], 'Type': 'C-C'}
                        Lines_connections = Lines_connections.append(pd.DataFrame(Data))

    Lines_connections = gpd.GeoDataFrame(Lines_connections, geometry=Lines_connections['geometry'], crs=crs)
    Lines_connections.to_file(data_folder + '/Lines_connections')
    Lines_connections.to_csv(data_folder + '/Lines_connections.csv')

    MILP_input_folder = gisele_folder + '/Case studies/' + case_study + '/Intermediate/Optimization/MILP_Input'
    Nodes = All_nodes_withPS.copy()
    Nodes = Nodes.drop('geometry', axis=1)
    if mg_option:
        Microgrids = pd.read_csv(case_folder+'/Output/Microgrid.csv')
        n_clusters = int(Nodes['Cluster'].max())
        microgrid_nodes = [*range(1000000,1000000*(n_clusters+1),1000000)]
    ps_list = Primary_substations['ID'].to_list()
    Nodes.reset_index(inplace=True)
    for i in ps_list:
        Nodes = Nodes.drop(Nodes[Nodes['ID'] == i].index)
    if mg_option:
        all_nodes = Nodes['ID'].to_list() + ps_list + microgrid_nodes
    else:
        all_nodes = Nodes['ID'].to_list() + ps_list
    Primary_substations['ID'].to_csv(MILP_input_folder + '/nodes_PS.csv', index=False)
    pd.DataFrame({'ID': all_nodes}).to_csv(MILP_input_folder + '/nodes.csv', index=False)
    Nodes['ID'].to_csv(MILP_input_folder + '/nodes_clusters.csv', index=False)
    power_list = Nodes['MV_Power'].to_list()
    power_list = [i / 1000 for i in power_list]  # transfer to MW
    power_list = [ceil(i * 1000) / 1000 for i in power_list]  # make sure to round to the third decimal, basically 1 kW
    Nodes['MV_Power'] = power_list
    Nodes[['ID', 'MV_Power']].to_csv(MILP_input_folder + '/power_nodes.csv', index=False)
    # Write files about the PS
    Primary_substations[['ID', 'Power']].to_csv(MILP_input_folder + '/PS_power_max.csv', index=False)
    Primary_substations[['ID', 'Cost']].to_csv(MILP_input_folder + '/PS_costs.csv', index=False)
    Primary_substations[['ID', 'Voltage']].to_csv(
        MILP_input_folder + '/PS_voltage.csv', index=False)
    # Write files about the MGs
    if mg_option:
        mg_powers = [Nodes.loc[Nodes['Cluster'] == i, 'MV_Power'].sum() for i in range(1, n_clusters + 1)]
        pd.DataFrame({'ID': microgrid_nodes}).to_csv(MILP_input_folder + '/microgrids_nodes.csv', index=False)
        pd.DataFrame({'ID': microgrid_nodes, 'power': mg_powers}).to_csv(MILP_input_folder + '/microgrids_powers.csv',
                                                                         index=False)
        pd.DataFrame({'ID': microgrid_nodes, 'energy': Microgrids.loc[:, 'Energy Demand [MWh]']}).to_csv(
            MILP_input_folder + '/energy.csv', index=False)
        pd.DataFrame({'ID': microgrid_nodes, 'Cost': Microgrids.loc[:, 'Total Cost [kEUR]']}).to_csv(
            MILP_input_folder + '/microgrids_costs.csv', index=False)

    Lines_connections[['ID1', 'ID2']].to_csv(MILP_input_folder + '/links_decision.csv', index=False)
    Lines_connections[['ID1', 'ID2', 'Cost']].to_csv(MILP_input_folder + '/weights_decision_lines.csv', index=False)
    Lines_connections[['ID1', 'ID2', 'Length']].to_csv(MILP_input_folder + '/distances.csv', index=False)



    # MILP
def MILP_1_cluster_1_point(gisele_folder,case_study,n_clusters,coe,voltage,resistance,reactance,Pmax,line_cost):


    model = AbstractModel()
    data = DataPortal()
    MILP_input_folder = gisele_folder + '/Case studies/' + case_study + '/Intermediate/Optimization/MILP_input'
    MILP_output_folder = gisele_folder + '/Case studies/' + case_study + '/Intermediate/Optimization/MILP_output'
    os.chdir(MILP_input_folder)

    # Define some basic parameter for the per unit conversion and voltage limitation
    Abase = 1
    Vmin = 0.9

    ####################Define sets#####################
    model.N = Set()
    data.load(filename='nodes.csv', set=model.N)  # first row is not read
    model.N_clusters = Set()
    data.load(filename='nodes_clusters.csv', set=model.N_clusters)
    model.N_PS = Set()
    data.load(filename='nodes_PS.csv', set=model.N_PS)
    # Node corresponding to primary substation

    model.links_decision = Set(dimen=2)
    data.load(filename='links_decision.csv', set=model.links_decision)
    # Connection distance of all the edges
    model.dist = Param(model.links_decision)
    data.load(filename='distances.csv', param=model.dist)

    # Nodes are divided into two sets, as suggested in https://pyomo.readthedocs.io/en/stable/pyomo_modeling_components/Sets.html:
    # NodesOut[nodes] gives for each node all nodes that are connected to it via outgoing links
    # NodesIn[nodes] gives for each node all nodes that are connected to it via ingoing links

    def NodesOut_init(model, node):
        retval = []
        for (i, j) in model.links_decision:
            if i == node:
                retval.append(j)
        return retval

    model.NodesOut = Set(model.N, initialize=NodesOut_init)

    def NodesIn_init(model, node):
        retval = []
        for (i, j) in model.links_decision:
            if j == node:
                retval.append(i)
        return retval

    model.NodesIn = Set(model.N, initialize=NodesIn_init)

    #####################Define parameters#####################

    # Electric power in the nodes (injected (-) or absorbed (+))
    model.Psub = Param(model.N_clusters)
    data.load(filename='power_nodes.csv', param=model.Psub)

    model.ps_cost = Param(model.N_PS)
    data.load(filename='PS_costs.csv', param=model.ps_cost)

    model.PSmax = Param(model.N_PS)
    data.load(filename='PS_power_max.csv', param=model.PSmax)

    model.PS_voltage = Param(model.N_PS)
    data.load(filename='PS_voltage.csv', param=model.PS_voltage)

    model.weights = Param(model.links_decision)
    data.load(filename='weights_decision_lines.csv', param=model.weights)
    # Electrical parameters of all the cables
    model.V_ref = Param(initialize=voltage)
    model.A_ref = Param(initialize=Abase)
    model.E_min = Param(initialize=Vmin)

    model.R_ref = Param(initialize=resistance)
    model.X_ref = Param(initialize=reactance)
    model.P_max = Param(initialize=Pmax)
    model.cf = Param(initialize=line_cost)

    model.Z = Param(initialize=model.R_ref + model.X_ref * 0.5)
    model.Z_ref = Param(initialize=model.V_ref ** 2 / Abase)

    model.n_clusters = Param(initialize=n_clusters)
    model.coe = Param(initialize=coe)

    #####################Define variables#####################

    # binary variable x[i,j]: 1 if the connection i,j is present, 0 otherwise
    model.x = Var(model.links_decision, within=Binary)
    # power[i,j] is the power flow of connection i-j
    model.P = Var(model.links_decision)
    # positive variables E(i) is p.u. voltage at each node
    model.E = Var(model.N, within=NonNegativeReals)

    # binary variable k[i]: 1 if node i is a primary substation, 0 otherwise
    model.k = Var(model.N_PS, within=Binary)
    # Power output of Primary substation
    model.PPS = Var(model.N_PS, within=NonNegativeReals)
    model.cable_type = Var(model.links_decision)

    #####################Define constraints###############################

    # Radiality constraint
    def Radiality_rule(model):
        return summation(model.x) == model.n_clusters

    model.Radiality = Constraint(rule=Radiality_rule)

    # Power flow constraints
    def Power_flow_conservation_rule(model, node):
        return (sum(model.P[j, node] for j in model.NodesIn[node]) - sum(
            model.P[node, j] for j in model.NodesOut[node])) == model.Psub[node]

    model.Power_flow_conservation = Constraint(model.N_clusters, rule=Power_flow_conservation_rule)

    def Power_flow_conservation_rule3(model, node):
        return (sum(model.P[j, node] for j in model.NodesIn[node]) - sum(
            model.P[node, j] for j in model.NodesOut[node])) == - model.PPS[node]

    model.Power_flow_conservation3 = Constraint(model.N_PS, rule=Power_flow_conservation_rule3)

    def Power_upper_decision(model, i, j):
        return model.P[i, j] <= model.P_max * model.x[i, j]

    model.Power_upper_decision = Constraint(model.links_decision, rule=Power_upper_decision)

    def Power_lower_decision(model, i, j):
        return model.P[i, j] >= -model.P_max * model.x[i, j]

    model.Power_lower_decision = Constraint(model.links_decision, rule=Power_lower_decision)

    # Voltage constraints
    def Voltage_balance_rule(model, i, j):
        return (model.E[i] - model.E[j]) + model.x[i, j] - 1 <= model.dist[i, j] / 1000 * model.P[
            i, j] * model.Z / model.Z_ref

    model.Voltage_balance_rule = Constraint(model.links_decision, rule=Voltage_balance_rule)

    def Voltage_balance_rule2(model, i, j):
        return (model.E[i] - model.E[j]) - model.x[i, j] + 1 >= model.dist[i, j] / 1000 * model.P[
            i, j] * model.Z / model.Z_ref

    model.Voltage_balance_rule2 = Constraint(model.links_decision, rule=Voltage_balance_rule2)

    def Voltage_limit(model, i):
        return model.E[i] >= model.k[i] * (model.PS_voltage[i] - model.E_min) + model.E_min

    model.Voltage_limit = Constraint(model.N_PS, rule=Voltage_limit)

    def Voltage_PS2(model, i):
        return model.E[i] <= model.PS_voltage[i]

    model.Voltage_PS2 = Constraint(model.N_PS, rule=Voltage_PS2)

    def Voltage_limit_clusters2(model, i):
        return model.E[i] >= model.E_min

    model.Voltage_limit_clusters2 = Constraint(model.N_clusters, rule=Voltage_limit_clusters2)

    def PS_power_rule_upper(model, i):
        return model.PPS[i] <= model.PSmax[i] * model.k[i]

    model.PS_power_upper = Constraint(model.N_PS, rule=PS_power_rule_upper)

    def Balance_rule(model):
        return (sum(model.PPS[i] for i in model.N_PS) - sum(model.Psub[i] for i in model.N_clusters)) == 0

    model.Balance = Constraint(rule=Balance_rule)

    ####################Define objective function##########################

    ####################Define objective function##########################
    reliability_index = 1000

    def ObjectiveFunction(model):
        return summation(model.weights, model.x) * model.cf / 1000 + summation(model.ps_cost, model.k)

        # return summation(model.weights, model.x) * model.cf / 1000  + summation(model.ps_cost,model.k) - sum(model.Psub[i]*model.Distance[i] for i in model.N_clusters)*reliability_index

    model.Obj = Objective(rule=ObjectiveFunction, sense=minimize)

    #############Solve model##################

    instance = model.create_instance(data)
    print('Instance is constructed:', instance.is_constructed())
    # opt = SolverFactory('cbc',executable=r'C:\Users\Asus\Desktop\POLIMI\Thesis\GISELE\Gisele_MILP\cbc')
    opt = SolverFactory('gurobi')
    opt.options['TimeLimit'] = 600
    # opt.options['numericfocus']=0
    # opt.options['mipgap'] = 0.0002
    # opt.options['presolve']=2
    # opt.options['mipfocus']=2
    # opt = SolverFactory('cbc',executable=r'C:\Users\Asus\Desktop\POLIMI\Thesis\GISELE\New folder\cbc')
    print('Starting optimization process')
    time_i = datetime.now()
    opt.solve(instance, tee=True, symbolic_solver_labels=True)
    time_f = datetime.now()
    print('Time required for optimization is', time_f - time_i)
    links = instance.x
    power = instance.P
    subs = instance.k
    voltage = instance.E
    PS = instance.PPS


    # voltage_drop=instance.z
    connections_output = pd.DataFrame(columns=[['id1', 'id2', 'power']])
    PrSubstation = pd.DataFrame(columns=[['index', 'power']])
    all_lines = pd.DataFrame(columns=[['id1', 'id2', 'power']])
    Voltages = pd.DataFrame(columns=[['index', 'voltage [p.u]']])

    k = 0
    for index in links:
        if int(round(value(links[index]))) == 1:
            connections_output.loc[k, 'id1'] = index[0]
            connections_output.loc[k, 'id2'] = index[1]
            connections_output.loc[k, 'power'] = value(power[index])
            k = k + 1
    k = 0
    for index in subs:
        if int(round(value(subs[index]))) == 1:
            PrSubstation.loc[k, 'index'] = index
            PrSubstation.loc[k, 'power'] = value(PS[index])
            print(((value(PS[index]))))
            k = k + 1
    k = 0
    for v in voltage:
        Voltages.loc[k, 'index'] = v
        Voltages.loc[k, 'voltage [p.u]'] = value(voltage[v])
        k = k + 1
    k = 0
    for index in power:
        all_lines.loc[k, 'id1'] = index[0]
        all_lines.loc[k, 'id2'] = index[1]
        all_lines.loc[k, 'power'] = value(power[index])
        k = k + 1

    connections_output.to_csv(MILP_output_folder + '/connections_output.csv', index=False)
    PrSubstation.to_csv(MILP_output_folder + '/PrimarySubstations.csv', index=False)
    Voltages.to_csv(MILP_output_folder + '/Voltages.csv', index=False)
    all_lines.to_csv(MILP_output_folder + '/all_lines.csv', index=False)

def process_output_1cluster_1point(gisele_folder, case_study):
    # This part is to calculate the cost for connecting to the primary substations / MV network
    read_input = gisele_folder + '\Case studies/' + case_study + '\Intermediate/Optimization\MILP_input/'
    read_output = gisele_folder + '\Case studies/' + case_study + '/Intermediate/Optimization\MILP_output/'
    read_general_data = gisele_folder + '\Case studies/' + case_study + '/Intermediate/Optimization/all_data/'
    write = gisele_folder + '\Case studies/' + case_study + '\Output\MILP_processed/'
    Connections = pd.read_csv(read_output + 'connections_output.csv')
    Voltages = pd.read_csv(read_output + 'Voltages.csv')
    Nodes = pd.read_csv(read_general_data + 'All_Nodes.csv')
    Lines = gpd.read_file(read_general_data+'/Lines_connections/Lines_connections.shp')
    Lines['geometry'] = Lines.geometry.apply(wkt.loads)
    for i,row in Connections.iterrows():
        id1 = row['id1']
        id2 = row['id2']
        geom = Lines.loc[(Lines['ID1']==id1) & (Lines['ID2']==id2),'geometry'].values
        Connections.loc[i,'geometry']=geom[0]