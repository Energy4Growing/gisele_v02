import pandas as pd
from shapely.geometry import Point, box, LineString, MultiPoint,MultiLineString
import geopandas as gpd
import os
import numpy as np
def process(gisele_folder,case_study,crs,mg_option,reliability_option):
    read_input = gisele_folder+'\Case studies/'+ case_study+ '\Intermediate/Optimization\MILP_input/'
    read_output=gisele_folder+'\Case studies/' + case_study+'/Intermediate/Optimization\MILP_output/'
    read_general_data = gisele_folder+'\Case studies/' + case_study+'/Intermediate/Optimization/all_data/'
    write = gisele_folder+'\Case studies/'+ case_study+'\Output\MILP_processed/'

    Connections=pd.read_csv(read_output+'connections_output.csv')
    All_connections_dist=pd.read_csv(read_input+'distances.csv')
    Connections_costs=pd.read_csv(read_input+'weights_decision_lines.csv')
    Voltages=pd.read_csv(read_output+'Voltages.csv')
    Nodes=pd.read_csv(read_general_data+'All_Nodes.csv')
    Lines_clusters = pd.read_csv(read_output+'links_clusters.csv')
    Lines_clusters_original = pd.read_csv(read_general_data+'Lines_clusters.csv')
    Lines_clusters_marked = gpd.read_file(read_general_data+'Lines_marked/Lines_marked.shp')
    Lines_clusters_marked['Conn_param'] = Lines_clusters_marked['Conn_param'].astype(int)

    #Lines=pd.DataFrame(columns=[['id1','id2','Length','Cost','geom']])
    if mg_option:
        Microgrid=pd.read_csv(read_output+'Microgrid.csv')
    else:
        print('Only on-grid')
    if reliability_option:
        dist_from_PS=pd.read_csv(read_output+'Distances.csv')
        case = 'reliability'
    else:
        case = 'no_reliability'

    if 'Type' in Lines_clusters.columns:
        case2 = 'Multiple_Types'
    else:
        case2 = 'One_Type'
    id1=[]
    id2=[]
    Length=[]
    Cost=[]
    geom=[]
    Loading=[]
    Type=[]
    Segments_connections = gpd.read_file(read_general_data +'Lines_connections/Lines_connections.shp')
    Segments_connections['ID2']=Segments_connections['ID2'].astype(int)
    Segments_connections['ID1']=Segments_connections['ID1'].astype(int)
    for index, row in Connections.iterrows():
        try:
            #geom.append(LineString(
            #    [(Nodes.loc[Nodes['ID'] == int(row['id1']), 'X'], Nodes.loc[Nodes['ID'] == int(row['id1']), 'Y']),
            #     (Nodes.loc[Nodes['ID'] == int(row['id2']), 'X'], Nodes.loc[Nodes['ID'] == int(row['id2']), 'Y'])]))
            linestring  = Segments_connections.loc[(Segments_connections['ID1']==row['id1']) & (Segments_connections['ID2']==row['id2']),'geometry']
            try:
                linestring=linestring.values[0]
            except:
                linestring = Segments_connections.loc[(Segments_connections['ID1'] == row['id2']) & (Segments_connections['ID2'] == row['id1']), 'geometry'].values[0]
            geom.append(linestring)
            id1.append(row['id1'])
            id2.append(row['id2'])
            if case2=='Multiple_Types':
                Type.append(row['Type'])
            Loading.append(round(row['power'],3))
            distance=All_connections_dist.loc[(All_connections_dist['ID1']==row['id1']) & (All_connections_dist['ID2']==row['id2']),'Length']
            cost= Connections_costs.loc[(Connections_costs['ID1']==row['id1']) & (Connections_costs['ID2']==row['id2']),'Cost']
            Length.append(round(abs(float(distance)),3))
            Cost.append(float(cost))


        except:
            print('MG connection - fake node')
    #Data = {'id1': id1,'id2': id2, 'Length': Length, 'Cost': Cost, 'Loading': Loading, 'Geom': geom}

    try:
        if case2 == 'Multiple_Types':
            Data = {'id1': id1, 'id2': id2, 'Length': Length, 'Loading': Loading,'Cost': Cost, 'Geom': geom ,'Type': Type}
            Lines = pd.DataFrame(Data)
            Lines.to_csv(write+'output_connections.csv',index=False)
        else:
            Data = {'id1': id1, 'id2': id2, 'Length': Length, 'Loading': Loading,'Cost': Cost, 'Geom': geom}
            Lines = pd.DataFrame(Data)
            Lines.to_csv(write + 'output_connections.csv', index=False)
    except:
        print('empty file - all Microgrids')
    if not Lines.empty:
        Lines=gpd.GeoDataFrame(Lines, geometry = 'Geom')
        Lines.crs=crs
        Lines.to_file(write+'output_connections')
    x=[]
    y=[]
    Power=[]
    voltage=[]
    ID=[]
    mg=[]
    cluster=[]
    dist=[]
    if not case == 'reliability':
        for index,row in Voltages.iterrows():
            try:
                id=row['index']

                x.append(float(Nodes.loc[(Nodes['ID']==id),'X']))
                y.append(float(Nodes.loc[(Nodes['ID']==id),'Y']))
                cluster.append(float(Nodes.loc[(Nodes['ID']==id),'Cluster']))
                ID.append(id)
                voltage.append(round(row['voltage [p.u]'],3))
                Power.append(round(float(Nodes.loc[(Nodes['ID']==id),'MV_Power']),3))
            except:
                print('microgrid fake node')
                print(row['index'])
        Data = {'ID': ID, 'X': x, 'Y': y, 'Voltage': voltage, 'Power': Power, 'Cluster': cluster}
    else:
        for index, row in Voltages.iterrows():
            try:
                id = row['index']

                x.append(float(Nodes.loc[(Nodes['ID'] == id), 'X']))
                y.append(float(Nodes.loc[(Nodes['ID'] == id), 'Y']))
                cluster.append(float(Nodes.loc[(Nodes['ID'] == id), 'Cluster']))
                ID.append(id)
                voltage.append(round(row['voltage [p.u]'],3))
                Power.append(round(float(Nodes.loc[(Nodes['ID']==id),'MV_Power']),3))
                dist.append(round(-float(dist_from_PS.loc[(dist_from_PS['index'] == id), 'length[m]']),3))
            except:
                print('microgrid fake node')
                print(row['index'])
        Data = {'ID': ID, 'X': x, 'Y': y, 'Voltage': voltage, 'Power': Power, 'Cluster': cluster,'Distance':dist}

    #Data={'ID': ID, 'X': x, 'Y': y, 'Voltage': voltage,'Microgrid': mg,'Power': Power }
    #Data={'ID': ID, 'X': x, 'Y': y, 'Voltage': voltage,'Microgrid': mg,'Power': Power, 'dist_from_PS': dist}
    Nodes=pd.DataFrame(Data)
    Nodes.to_csv(write+'output_nodes.csv',index=False)
    Nodes=gpd.GeoDataFrame(Nodes, geometry=gpd.points_from_xy(Nodes.X, Nodes.Y))
    Nodes.crs=crs
    Nodes.to_file(write+'output_nodes')
    id1=[]
    id2=[]
    Length=[]
    Cost=[]
    geom=[]
    Loading=[]
    Type=[]
    Lines_clusters=Lines_clusters[Lines_clusters['power']!=0]
    for index, row in Lines_clusters.iterrows():
        id1.append(row['id1'])
        id2.append(row['id2'])
        Loading.append(round(row['power'],3))
        distance=All_connections_dist.loc[(All_connections_dist['ID1']==row['id1']) & (All_connections_dist['ID2']==row['id2']),'Length']
        if case2 == 'Multiple_Types':
            Type.append(row['Type'])
        #cost= All_connections_cost.loc[(All_connections_cost['id1']==row['id1']) & (All_connections_cost['id2']==row['id2']),'distance']
        Length.append(round(abs(float(distance)),3))
        #Cost.append(float(cost))

        # locate the conn param that connects the output of the MILP with the actual line following the roads
        conn_param = Lines_clusters_original.loc[(Lines_clusters_original['ID1'] == row['id1']) &
                                            (Lines_clusters_original['ID2'] == row['id2']), 'Conn_param']
        cost = Lines_clusters_original.loc[(Lines_clusters_original['ID1'] == row['id1']) &
                                            (Lines_clusters_original['ID2'] == row['id2']), 'Cost']
        try:
            conn_param=int(conn_param.values[0])
            cost = float(cost.values[0])
        except:
            conn_param = int(Lines_clusters_original.loc[(Lines_clusters_original['ID1'] == row['id2']) &
                                                     (Lines_clusters_original['ID2'] == row['id1']), 'Conn_param'].values[0])
            cost = float(Lines_clusters_original.loc[(Lines_clusters_original['ID1'] == row['id2']) &
                                               (Lines_clusters_original['ID2'] == row['id1']), 'Cost'].values[0])
        Cost.append(float(cost))
        # now we have the conn_param, we should access the file Lines_clusters_marked
        LiNeS = Lines_clusters_marked.loc[Lines_clusters_marked['Conn_param'] == conn_param, 'geometry']
        geometry = MultiLineString([line for line in LiNeS])
        geom.append(geometry)
    if case2 == 'Multiple_Types':
        Data = {'id1': id1, 'id2': id2, 'Length': Length, 'Loading': Loading,'Cost':Cost, 'Geom': geom,'Type':Type}
    else:
        Data = {'id1': id1, 'id2': id2, 'Length': Length, 'Loading': Loading, 'Cost': Cost, 'Geom': geom}
    Lines2 = pd.DataFrame(Data)
    Lines2.to_csv(write + 'output_cluster_lines.csv', index=False)
    Lines2=gpd.GeoDataFrame(Lines2, geometry = 'Geom')
    Lines2.crs=crs
    Lines2.to_file(write+'output_cluster_lines')
def create_final_output(gisele_folder, case_study):
    # This part is to calculate the cost for connecting to the primary substations / MV network
    read_output = gisele_folder + '\Case studies/' + case_study + '/Intermediate/Optimization\MILP_output/'
    write = gisele_folder + '\Case studies/' + case_study + '\OUtput\MILP_processed/'

    Lines2 = pd.read_csv(write + 'output_cluster_lines.csv')
    Lines = pd.read_csv(write + 'output_connections.csv')
    Primary_substations = gpd.read_file(
        gisele_folder + '/Case studies/' + case_study + '/Input/substations/substations.shp')
    PS_chosen = gpd.read_file(read_output+'PrimarySubstations.csv')
    LV_resume = pd.read_csv(gisele_folder + '/Case studies/' + case_study + '/Output/LV_resume.csv')
    PS_list = PS_chosen['index'].to_list()
    costs_of_connection=[]
    for i in PS_list:
        cost = Primary_substations.loc[Primary_substations['ID']==int(i),'Cost']
        costs_of_connection.append(float(cost))
    connection_costs = sum(costs_of_connection)
    number_connection_points = PS_chosen.shape[0]
    cost_of_MV_lines = Lines['Cost'].sum() + Lines2['Cost'].sum()
    length_of_MV_lines= (Lines['Length'].sum() + Lines2['Length'].sum())/1000
    total_population = LV_resume['Population'].sum()
    cost_secondary_subs = LV_resume['Cost[euro]'].sum()
    length_LV_network = LV_resume['Grid_Length [km]'].sum()
    cost_LV_network = LV_resume['Grid Cost [euro]'].sum()
    number_clusters = LV_resume['Cluster'].max()
    number_secondary_substations = LV_resume.shape[0]
    total_costs = cost_of_MV_lines + cost_secondary_subs+cost_LV_network + connection_costs
    # number of villages, number of people,
    Data = {'Number_clusters': number_clusters,'Number_connection_points': number_connection_points,'Total population':total_population,
            'Number_substations':number_secondary_substations,'Length_LV_network': length_LV_network,'Length_MV_network':length_of_MV_lines,
            'Cost_substations': cost_secondary_subs,'Cost_LV_network[euro]': cost_LV_network,'Cost_MV_network[euro]': cost_of_MV_lines,
            'Cost_connection_points[euro]': connection_costs,'Total_costs[euro]': total_costs}

    grid_resume = pd.DataFrame(Data,index=[0])
    grid_resume.transpose().to_csv(gisele_folder + '/Case studies/' + case_study + '/Output/Final_resume.csv')