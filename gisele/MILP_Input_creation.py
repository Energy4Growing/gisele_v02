from math import ceil
from gisele import dijkstra
from shapely.ops import nearest_points
from gisele.functions2 import *


#import SS_Clustering
def create_node_lines(Nodes,gisele_folder,case_study,crs):
    ''' The goal of this function is to create a csv with the Nodes and Lines (inside clusters) that are supposed to be considered
    in the MILP. However, further processing is required to get the exact files that are needed'''
    case_folder = gisele_folder+'/Case studies/'+case_study+'/'
    folder_clusters=case_folder+ 'Intermediate/Communities/'
    clusters=os.listdir(folder_clusters)
    conn_param=0
    New_Lines = pd.DataFrame()
    All_lines_marked=pd.DataFrame()

    for n_cluster in range(len(clusters)):
        cluster=clusters[n_cluster]
        print(cluster)
        #Nodes = gpd.read_file(folder_clusters + cluster + '/points500m.shp')
        Nodes_clusters = Nodes[Nodes['Cluster'] == int(cluster)]
        try:
            Lines = gpd.read_file(folder_clusters+cluster+'/Grid_substations.shp')
            Lines['ID1'] = Lines['ID1'].astype(int)
            Lines['ID2'] = Lines['ID2'].astype(int)
        except:
            print('Just 1 substation, no internal grid')
            if n_cluster == 0:
                all_substations = Nodes_clusters[Nodes_clusters.Substation == 1]
            else:
                add_nodes = Nodes_clusters[Nodes_clusters.Substation == 1]
                all_substations = all_substations.append(add_nodes)

            continue

        Lines_marked=Lines.copy() # this file is to allow after the MILP, to go back to the original lines through the grid of points
        Lines_marked['Conn_param']=''
        #Now we also need to add intersection nodes
        nodes = Lines.ID1.to_list() + Lines.ID2.to_list()
        occurence = Counter(nodes)
        terminal_nodes=[]
        for i in occurence:
            if occurence[i]>2:
                Nodes.loc[Nodes['ID']==int(i),'Substation']=1
                Nodes.loc[Nodes['ID']==int(i),'Population']=0
                Nodes.loc[Nodes['ID'] == int(i), 'MV_Power'] = 0
                Nodes.loc[Nodes['ID'] == int(i), 'Cluster'] = int(cluster)
                print(i)
        Nodes_clusters = Nodes[Nodes['Cluster'] == int(cluster)]
        Substations = Nodes_clusters.loc[Nodes_clusters.Substation == 1, 'ID'].to_list()
        if n_cluster==0:
            all_substations=Nodes_clusters[Nodes_clusters.Substation==1]
        else:
            add_nodes=Nodes_clusters[Nodes_clusters.Substation == 1]
            all_substations=all_substations.append(add_nodes)
        while not Lines.empty:
            nodes = Lines.ID1.to_list() + Lines.ID2.to_list()
            print(nodes)
            nodes=[int(i) for i in nodes]
            Substations=list(set(Substations) & set(nodes))
            current_node=int(Substations[0])
            print(current_node)
            no_lines=False
            tot_length=0
            tot_cost=0
            id1=current_node
            while not no_lines:
                next_index=Lines.index[Lines['ID1']== current_node].to_list()
                # check if there actually is a next node
                if next_index:
                    next_index=next_index[0] # i only care about the first line if there are many
                    next_node=Lines.loc[next_index,'ID2']
                    print(next_node)
                    tot_length = tot_length + Lines.loc[next_index, 'Length']
                    tot_cost = tot_cost + Lines.loc[next_index, 'Cost']
                    Lines.drop(index=next_index,inplace=True)
                else:
                    next_index=Lines.index[Lines['ID2']== current_node].to_list()
                    if next_index:
                        next_index = next_index[0]  # i only care about the first line if there are many
                        next_node = Lines.loc[next_index, 'ID1']
                        print(next_node)
                        tot_length = tot_length + Lines.loc[next_index, 'Length']
                        tot_cost = tot_cost + Lines.loc[next_index, 'Cost']
                        Lines.drop(index=next_index, inplace=True)
                    else:
                        no_lines=True #there are no lines starting from this node
                if not no_lines:
                    is_substation = Nodes.loc[Nodes.ID==int(next_node),'Substation']==1
                    Lines_marked.loc[next_index,'Conn_param']=conn_param
                    if is_substation.values[0]:
                        cond=False
                        Point1 = Nodes.loc[Nodes['ID']==int(id1),'geometry'].values[0]
                        Point2 = Nodes.loc[Nodes['ID']==int(next_node),'geometry'].values[0]
                        geom=LineString([Point1,Point2])
                        Data={'ID1':id1, 'ID2': next_node, 'Cost': tot_cost,'Length':tot_length, 'geometry': geom,'Conn_param':conn_param}
                        New_Lines=New_Lines.append(Data,ignore_index=True)
                        current_node=next_node
                        tot_length = 0
                        tot_cost = 0
                        id1 = current_node
                        conn_param=conn_param+1
                    else:
                        current_node=next_node
        All_lines_marked=All_lines_marked.append(Lines_marked)

    New_Lines.to_csv(case_folder+'Intermediate/Optimization/all_data/Lines_clusters.csv',index=False)
    All_lines_marked.to_csv(case_folder + 'Intermediate/Optimization/all_data/Lines_marked.csv',index=False)
    All_lines_marked=gpd.GeoDataFrame(All_lines_marked,geometry=All_lines_marked.geometry,crs=crs)
    All_lines_marked.to_file(case_folder + 'Intermediate/Optimization/all_data/Lines_'
                                           'marked',index=False)
    all_substations['Population']=[ceil(i) for i in all_substations['Population'].to_list()]
    all_substations[['X','Y','ID','Weight','Cluster','MV_Power','Population','Substation','geometry']].to_csv(
        case_folder+'Intermediate/Optimization/all_data/All_Nodes.csv')

def create_connections(all_points,Primary_substations,gisele_folder,case_study,line_bc,resolution,crs,Roads_option,tolerance_outside, Rivers_option):
    print('START CONNECTIONS')
    crs_string = "EPSG:"+str(crs)
    case_folder = gisele_folder + '/Case studies/' + case_study
    data_folder = case_folder + '/Intermediate/Optimization/all_data'
    MV_grid_inside_clusters = gpd.read_file(case_folder + '/Output/MV_grid/MV_grid.shp')
    Nodes = pd.read_csv(data_folder + '/All_Nodes.csv')
    Lines_clusters = pd.read_csv(data_folder + '/Lines_clusters.csv')
    n_clusters = int(Nodes['Cluster'].max())
    # start creating links between clusters with Djikstra
    geometry = [Point(xy) for xy in zip(Nodes.X, Nodes.Y)]
    Nodes_gdf = gpd.GeoDataFrame(Nodes, crs=crs_string, geometry=geometry)
    c_grid_points = []
    for i,row in MV_grid_inside_clusters.iterrows():
        c_grid_points.append([row.ID1,row.ID2])
    Lines_connections = pd.DataFrame(columns=Lines_clusters.columns)
    Lines_connections = Lines_connections.drop(['Conn_param'], axis=1)
    Lines_connections['Type'] = ''

    if Roads_option:
        Roads_points_all = gpd.read_file(case_folder + '/Intermediate/GeoSpatial_Data/Roads_points/Roads_points.shp')
        Roads_lines_all = gpd.read_file(case_folder + '/Intermediate/GeoSpatial_Data/Roads_lines/Roads_lines.shp')
    for i in range(1,n_clusters+1):
        for j in range(i + 1, n_clusters+1, 1):
            print('C-C '+str(i)+str(j))
            grid_1 = Nodes[Nodes['Cluster'] == i]
            grid_2 = Nodes[Nodes['Cluster'] == j]
            dist_2d = pd.DataFrame(distance_2d(grid_1, grid_2, 'X', 'Y'),
                                   index=grid_1.ID.values,
                                   columns=grid_2.ID.values)

            p1 = Nodes_gdf[Nodes_gdf['ID'] == dist_2d.min().idxmin()]
            p2 = Nodes_gdf[Nodes_gdf['ID'] ==
                           dist_2d.min(axis=1).idxmin()]
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
                connection, connection_cost, connection_length, pts = dijkstra.dijkstra_connection_roads_new(all_points, p1,
                            p2, c_grid_points,    line_bc,  resolution,gdf_roads, roads_segments,Rivers_option,1,1.5,resolution * 20 )


            elif not Roads_option and dist<30*resolution :
                connection, connection_cost, connection_length, pts = dijkstra.dijkstra_connection(all_points, p1, p2,
                                                                                               c_grid_points, line_bc,
                                                                                               resolution,resolution*30, Rivers_option)

            else:
                connection_cost = 999999
            if not connection_cost == 999999 :
                # THIS SEGMENT HERE IS TO MAKE SURE THAT THE CONNECTIONS ARE NOT ILLOGICAL AND FAR INTO THE INTERNAL
                # CLUSTER GRID
                MV_grid = MV_grid_inside_clusters[MV_grid_inside_clusters['Cluster'].isin([i,j])]
                Line, connection_length, connection_cost,term_pts = improve_connection\
                    (MV_grid, all_points, pts, grid_1, grid_2, line_bc)
                Data = {'ID1': term_pts[0], 'ID2': term_pts[1], 'C1': i, 'C2': j, 'Cost': connection_cost,
                        'Length': connection_length, 'geometry': [Line], 'Type': 'C-C'}
                Lines_connections = Lines_connections.append(pd.DataFrame(Data))

    print('Something')

    for i in range(1,n_clusters+1):
        cluster_nodes = Nodes_gdf[Nodes_gdf['Cluster'] == i]
        points_list = cluster_nodes.loc[cluster_nodes['ID'] > 0, 'geometry'].to_list()
        # points_list =  all_points['geometry'].to_list()
        grid_multipoint = MultiPoint(points_list)
        print('S-C '+str(i))
        for index, row in Primary_substations.iterrows():
            sub_node = pd.DataFrame({'X': [row['X']], 'Y': [row['Y']],
                                     'ID': [row['ID']]})  # this is to fix an issue with this being data series
            dist_matrix = pd.DataFrame(distance_2d(cluster_nodes, sub_node, 'X', 'Y'),
                                       index=cluster_nodes.ID.values,
                                       columns=[Primary_substations.ID.values[index]])
            p1 = Nodes_gdf[Nodes_gdf['ID'] == dist_matrix.min().idxmin()]
            p2 = Nodes_gdf[Nodes_gdf['ID'] ==
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

                connection, connection_cost, connection_length, pts = dijkstra.dijkstra_connection_roads_new(all_points, p1,
                            p2, c_grid_points,    line_bc,  resolution,gdf_roads, roads_segments,Rivers_option,1,1.5,resolution * 25 )

            elif not Roads_option and dist<50*resolution:
            #try:  # in case the PS is too close to the cluster
                connection, connection_cost, connection_length, pts = dijkstra.dijkstra_connection(all_points, p1, p2,
                         c_grid_points,  line_bc,   resolution, Rivers_option,resolution*50)

            else:
                connection_cost = 999999
            if not connection_cost == 999999:
                MV_grid = MV_grid_inside_clusters[MV_grid_inside_clusters['Cluster'].isin([i])]
                Line, connection_length, connection_cost, term_pts = improve_connection \
                    (MV_grid, all_points, pts, cluster_nodes, sub_node, line_bc)
                Data = {'ID1': term_pts[0], 'ID2': term_pts[1], 'C1': i, 'C2':  'S' + str(index), 'Cost': connection_cost,
                        'Length': connection_length, 'geometry': [Line], 'Type': 'S-C'}
                Lines_connections = Lines_connections.append(pd.DataFrame(Data))
            #except:
            #    nearest_geoms = nearest_points(row['geometry'], grid_multipoint)
            #    connection_length = nearest_geoms[0].distance(nearest_geoms[1])
            #   print('Losho')
            #    connection_cost = all_points.loc[all_points['geometry'] == nearest_geoms[1], 'Weight'].values[
            #                          0] * connection_length

            #    pt1 = row['geometry']
            #    pt2 = all_points.loc[all_points['geometry'] == nearest_geoms[1], 'geometry'].values[0]
            #    Line = LineString([pt1, pt2])
            #    Data = {'ID1': p2['ID'].values[0], 'ID2': p1['ID'].values[0], 'C1': i, 'C2': 'S' + str(index),
            #            'Cost': connection_cost,
            #            'Length': connection_length, 'geometry': [Line], 'Type': 'S-C'}
            #    Lines_connections = Lines_connections.append(pd.DataFrame(Data))

    Lines_connections = gpd.GeoDataFrame(Lines_connections, geometry=Lines_connections['geometry'], crs=crs)
    Lines_connections.to_file(data_folder + '/Lines_connections')
    Lines_connections.to_csv(data_folder + '/Lines_connections.csv')

def create_MILP_input(all_points,gisele_folder,case_study,line_bc,resolution,crs,mg_option):
    crs_string = "EPSG:"+str(crs)
    case_folder = gisele_folder + '/Case studies/' + case_study
    data_folder = case_folder + '/Intermediate/Optimization/all_data'
    MILP_input_folder = gisele_folder + '/Case studies/' + case_study + '/Intermediate/Optimization/MILP_Input'


    # read everything else
    Lines_clusters = pd.read_csv(data_folder+'/Lines_clusters.csv')
    Lines_connections = pd.read_csv(data_folder+'/Lines_connections.csv')
    Nodes = pd.read_csv(data_folder+'/All_Nodes.csv')
    Nodes=Nodes.drop('geometry',axis=1)
    Primary_substations = gpd.read_file(gisele_folder + '/Case studies/' + case_study+'/Input/substations/substations.shp')

    n_substations=Primary_substations.shape[0]


    n_clusters = int(Nodes['Cluster'].max())
    if mg_option:
        Microgrids = pd.read_csv(case_folder + '/Output/Microgrid.csv')
        microgrid_nodes = [*range(1000000, 1000000 * (n_clusters + 1), 1000000)]
    ps_list=Primary_substations['ID'].to_list()
    for i in ps_list:
        Nodes=Nodes.drop(Nodes[Nodes['ID']==i].index)
    if mg_option:
        all_nodes = Nodes['ID'].to_list()+ps_list+microgrid_nodes
    else:
        all_nodes = Nodes['ID'].to_list()+ps_list

    #Write all different types of nodes to csv files


    Primary_substations['ID'].to_csv(MILP_input_folder+'/nodes_PS.csv',index=False)
    pd.DataFrame({'ID':all_nodes}).to_csv(MILP_input_folder+'/nodes.csv',index=False)
    Nodes['ID'].to_csv(MILP_input_folder + '/nodes_clusters.csv', index=False)
    power_list = Nodes['MV_Power'].to_list()
    power_list = [ceil(i*10)/10 for i in power_list]
    Nodes['MV_Power'] = power_list
    Nodes['MV_Power']=Nodes['MV_Power']/1000
    Nodes[['ID','MV_Power']].to_csv(MILP_input_folder + '/power_nodes.csv', index=False)
    # Write files about the PS
    Primary_substations[['ID','Power']].to_csv(MILP_input_folder +'/PS_power_max.csv',index=False)
    Primary_substations[['ID', 'Cost']].to_csv(MILP_input_folder + '/PS_costs.csv', index=False)
    Primary_substations[['ID', 'Voltage']].to_csv(MILP_input_folder + '/PS_voltage.csv', index=False)
    Primary_substations[['ID', 'Distance']].to_csv(MILP_input_folder + '/PS_distance.csv', index=False)
    # Write files about the MGs
    if mg_option:
        mg_powers = [Nodes.loc[Nodes['Cluster']== i,'MV_Power'].sum() for i in range(n_clusters)]
        pd.DataFrame({'ID': microgrid_nodes}).to_csv(MILP_input_folder + '/microgrids_nodes.csv', index=False)
        pd.DataFrame({'ID': microgrid_nodes,'power':mg_powers}).to_csv(MILP_input_folder +'/microgrids_powers.csv',index=False)
        pd.DataFrame({'ID': microgrid_nodes,'energy': Microgrids.loc[:,'Energy Consumed [MWh]']}).to_csv(MILP_input_folder +'/energy.csv',index=False)
        pd.DataFrame({'ID': microgrid_nodes, 'Cost': Microgrids.loc[:,'Total Cost [k€]']}).to_csv(MILP_input_folder +'/microgrids_costs.csv',index=False)


    #Start processing lines - make sure to put both paths 1-2 and 2-1
    Lines_clusters_duplicate=Lines_clusters.copy()
    Lines_clusters_duplicate['ID1']=Lines_clusters['ID2']
    Lines_clusters_duplicate['ID2'] = Lines_clusters['ID1']
    Lines_clusters=Lines_clusters.append(Lines_clusters_duplicate)
    Lines_clusters[['ID1','ID2']].to_csv(MILP_input_folder + '/links_clusters.csv', index=False)
    Lines_clusters[['ID1', 'ID2','Cost']].to_csv(MILP_input_folder + '/weights_clusters.csv', index=False)

    #Links decision - here we also add the "Fake" microgrid points
    if mg_option:
        for i in range(n_clusters):
            node_mg = microgrid_nodes[i]
            node_cluster = Nodes.loc[Nodes['Cluster']==i,'ID'].values[0]
            Data = {'ID1':[node_mg],'ID2': [node_cluster], 'Cost': [0],'Length': [1],'Type':'MG'}
            Lines_connections=Lines_connections.append(pd.DataFrame(Data))

    Lines_connections_duplicate = Lines_connections.copy()
    Lines_connections_duplicate['ID1'] = Lines_connections['ID2']
    Lines_connections_duplicate['ID2'] = Lines_connections['ID1']
    Lines_connections = Lines_connections.append(Lines_connections_duplicate)


    Lines_connections[['ID1','ID2']].to_csv(MILP_input_folder + '/links_decision.csv', index=False)
    Lines_connections[['ID1', 'ID2','Cost']].to_csv(MILP_input_folder + '/weights_decision_lines.csv', index=False)

    All_lines = Lines_clusters.copy()
    All_lines = All_lines.append(Lines_connections)
    All_lines[['ID1','ID2']].to_csv(MILP_input_folder + '/links_all.csv', index=False)
    All_lines[['ID1', 'ID2','Length']].to_csv(MILP_input_folder + '/distances.csv', index=False)
def create_MILP_input_1way(all_points,gisele_folder,case_study,line_bc,resolution,crs):
    crs_string = "EPSG:"+str(crs)
    case_folder = gisele_folder + '/Case studies/' + case_study
    data_folder = case_folder + '/Intermediate/Optimization/all_data'
    MILP_input_folder = gisele_folder + '/Case studies/' + case_study + '/Intermediate/Optimization/MILP_Input'


    # read everything else
    Lines_clusters = pd.read_csv(data_folder+'/Lines_clusters.csv')
    Lines_connections = pd.read_csv(data_folder+'/Lines_connections.csv')
    Nodes = pd.read_csv(data_folder+'/All_Nodes.csv')
    Nodes=Nodes.drop('geometry',axis=1)

    Primary_substations = gpd.read_file(gisele_folder + '/Case studies/' + case_study+'/Input/substations/substations.shp')

    n_substations=Primary_substations.shape[0]
    Microgrids = pd.read_csv(case_folder+'/Output/Microgrid.csv')

    n_clusters = int(Nodes['Cluster'].max())

    microgrid_nodes = [*range(1000000,1000000*(n_clusters+1),1000000)]
    ps_list=Primary_substations['ID'].to_list()
    for i in ps_list:
        Nodes=Nodes.drop(Nodes[Nodes['ID']==i].index)
    all_nodes = Nodes['ID'].to_list()+ps_list+microgrid_nodes

    #Write all different types of nodes to csv files


    Primary_substations['ID'].to_csv(MILP_input_folder+'/nodes_PS.csv',index=False)
    pd.DataFrame({'ID':all_nodes}).to_csv(MILP_input_folder+'/nodes.csv',index=False)
    Nodes['ID'].to_csv(MILP_input_folder + '/nodes_clusters.csv', index=False)
    power_list = Nodes['MV_Power'].to_list()
    power_list = [ceil(i*100)/100 for i in power_list]
    Nodes['MV_Power'] = power_list
    Nodes['MV_Power']=Nodes['MV_Power']/1000
    Nodes[['ID','MV_Power']].to_csv(MILP_input_folder + '/power_nodes.csv', index=False)
    # Write files about the PS
    Primary_substations[['ID','Power']].to_csv(MILP_input_folder +'/PS_power_max.csv',index=False)
    Primary_substations[['ID', 'Cost']].to_csv(MILP_input_folder + '/PS_costs.csv', index=False)
    Primary_substations[['ID', 'Voltage']].to_csv(
        MILP_input_folder + '/PS_voltage.csv', index=False)
    # Write files about the MGs
    mg_powers = [Nodes.loc[Nodes['Cluster']== i,'MV_Power'].sum() for i in range(n_clusters)]
    pd.DataFrame({'ID': microgrid_nodes}).to_csv(MILP_input_folder + '/microgrids_nodes.csv', index=False)
    pd.DataFrame({'ID': microgrid_nodes,'power':mg_powers}).to_csv(MILP_input_folder +'/microgrids_powers.csv',index=False)
    pd.DataFrame({'ID': microgrid_nodes,'energy': Microgrids.loc[:,'Energy Demand [MWh]']}).to_csv(MILP_input_folder +'/energy.csv',index=False)
    pd.DataFrame({'ID': microgrid_nodes, 'Cost': Microgrids.loc[:,'Total Cost [kEUR]']}).to_csv(MILP_input_folder +'/microgrids_costs.csv',index=False)


    #Start processing lines

    Lines_clusters[['ID1','ID2']].to_csv(MILP_input_folder + '/links_clusters.csv', index=False)
    Lines_clusters[['ID1', 'ID2','Cost']].to_csv(MILP_input_folder + '/weights_clusters.csv', index=False)

    #Links decision - here we also add the "Fake" microgrid points
    for i in range(n_clusters):
        node_mg = microgrid_nodes[i]
        node_cluster = Nodes.loc[Nodes['Cluster']==i+1,'ID'].values[0]
        Data = {'ID1':[node_mg],'ID2': [node_cluster], 'Cost': [0],'Length': [1],'Type':'MG'}
        Lines_connections=Lines_connections.append(pd.DataFrame(Data))

    Lines_connections[['ID1','ID2']].to_csv(MILP_input_folder + '/links_decision.csv', index=False)
    Lines_connections[['ID1', 'ID2','Cost']].to_csv(MILP_input_folder + '/weights_decision_lines.csv', index=False)

    All_lines = Lines_clusters.copy()
    All_lines = All_lines.append(Lines_connections)
    All_lines[['ID1','ID2']].to_csv(MILP_input_folder + '/links_all.csv', index=False)
    All_lines[['ID1', 'ID2','Length']].to_csv(MILP_input_folder + '/distances.csv', index=False)


def array_to_LineString(all_points,pts):
    pts_new=[]
    i = 0
    for pt in pts:
        if i==0:
            pts_new.append(all_points.loc[all_points['ID']==pt[0],'geometry'].values[0])
            pts_new.append(all_points.loc[all_points['ID'] == pt[1], 'geometry'].values[0])
        else:
            pts_new.append(all_points.loc[all_points['ID']==pt[1],'geometry'].values[0])
        print(i)
        i+=1
    Line = LineString(pts_new)
    return Line

def add_PS_to_grid_of_points(all_points,Primary_substations,next_ID,add_elevation):
    '''ADDING A PRIMARY SUBSTATION POINT TO THE ALREADY EXISTING FILE WITH POINTS. THE WEIGHT WILL BE EQUAL
    TO THE ONE OF THE CLOSEST POINT. THIS CAN BE CHANGED BY CALCULATING THE AVERAGE OF THE 4 CLOSEST. BOTH
    INPUT ELEMENTS NEED TO BE GEODATAFRAMES'''
    points_list = all_points['geometry'].to_list()

    grid_multipoint = MultiPoint(points_list)
    for index,row in Primary_substations.iterrows():
        nearest_geoms = nearest_points(row['geometry'], grid_multipoint)
        weight_of_closest = all_points.loc[all_points['geometry']==nearest_geoms[1],'Weight'].values[0]
        if add_elevation:
            elevation_of_closest = all_points.loc[all_points['geometry']==nearest_geoms[1],'Elevation'].values[0]
            Data={'X':[row['X']],'Y': [row['Y']], 'ID': [next_ID], 'Substation':[1],'geometry': [row['geometry']],
              'Weight': [weight_of_closest], 'Cluster': [-1],'MV_Power': [0],'Elevation': [elevation_of_closest],'Type':['Primary substation']}
        else:
            Data={'X':[row['X']],'Y': [row['Y']], 'ID': [next_ID], 'Substation':[1],'geometry': [row['geometry']],
              'Weight': [weight_of_closest], 'Cluster': [-1],'MV_Power': [0]}
        all_points = all_points.append(gpd.GeoDataFrame(Data))
        next_ID+=1
    return all_points
def create_input(gisele_folder,case_study,crs,line_bc,resolution,reliability_option,Roads_option,tolerance_outside,Rivers_option,mg_option,mg_types):

    crs_str = "EPSG:"+str(crs)
    case_folder = gisele_folder+'/Case studies/'+case_study

    if Roads_option:
        all_points = pd.read_csv(case_folder+'/Intermediate/Geospatial_Data/weighted_grid_of_points_with_ss_and_roads.csv')
    else:
        all_points = pd.read_csv(case_folder + '/Intermediate/Geospatial_Data/weighted_grid_of_points_with_ss.csv')

    Primary_substations = gpd.read_file(
        gisele_folder + '/Case studies/' + case_study + '/Input/substations/substations.shp')
    Primary_substations=Primary_substations.to_crs(crs)
    geometry = [Point(xy) for xy in zip(all_points.X, all_points.Y)]


    # This here is to make sure that the points of the internal MV grid are included in all_points
    all_points = gpd.GeoDataFrame(all_points, crs=crs_str, geometry=geometry)
    final_houses = gpd.read_file(case_folder+'/Output/final_users/final_users.shp')
    MV_grid = gpd.read_file(case_folder+'/Output/MV_grid/MV_grid.shp')
    MV_grid_nodes = MV_grid.ID1.to_list()+MV_grid.ID2.to_list()
    MV_grid_nodes=list(set(MV_grid_nodes))
    final_houses = final_houses[(final_houses['ID'].isin(MV_grid_nodes) &
                                 ~final_houses['ID'].isin(all_points['ID'].to_list()))]

    final_houses['Type'] = 'LV_Users'
    final_houses['Weight'] = 2
    all_points=all_points.append(final_houses)

    next_ID = int(all_points['ID'].max()+1)
    all_points_withPS = add_PS_to_grid_of_points(all_points,Primary_substations,next_ID,add_elevation=True)

    #all_points_withPS.to_csv(case_folder+'/Input/weighted_grid_of_points_with_ss.csv',index=False)
    Primary_substations['ID'] = [*range(next_ID,next_ID+Primary_substations.shape[0])]
    Primary_substations.to_file(
        gisele_folder + '/Case studies/' + case_study + '/Input/substations/substations.shp')
    create_node_lines(all_points,gisele_folder,case_study,crs)
    # add the PS also to the file All_Nodes.csv
    All_nodes = pd.read_csv(case_folder+'/Intermediate/Optimization/all_data/All_Nodes.csv')
    geometry = [Point(xy) for xy in zip(All_nodes.X, All_nodes.Y)]
    All_nodes = gpd.GeoDataFrame(All_nodes, crs=crs_str, geometry=geometry)
    All_nodes_withPS = add_PS_to_grid_of_points(All_nodes,Primary_substations,next_ID,add_elevation=False)
    All_nodes_withPS.to_csv(case_folder+'/Intermediate/Optimization/all_data/All_Nodes.csv',index=False)
    #different functions to better organize the code
    # make sure the copy the important files in the Input folder ( configuration, Load profile etc)
    #
    create_connections(all_points_withPS,Primary_substations,gisele_folder,case_study,line_bc,resolution,crs,Roads_option,tolerance_outside,Rivers_option)
    #calculate_NPC = blablabla
    #calculate_NPC.to_csv
    if mg_option:
        calculate_mg(gisele_folder, case_study, crs,mg_types)
    if reliability_option:
        create_MILP_input(all_points_withPS,gisele_folder,case_study,line_bc,resolution,crs,mg_option) # all the lines are from i-j and j-i, only positive powers to consider reliability
    else:
        create_MILP_input_1way(all_points_withPS,gisele_folder,case_study,line_bc,resolution,crs) # only i-j, without reliability