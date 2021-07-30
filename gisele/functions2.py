from shapely.geometry import MultiLineString,MultiPolygon
from collections import Counter
import networkx as nx
from gisele.Steiner_tree_code import *
from gisele.functions import *
def create_roads_new(gisele_folder, case_study, Clusters,crs, accepted_road_types,resolution_MV,resolution_LV):
    weighted_grid_of_points = pd.read_csv(gisele_folder+'/Case studies/'+case_study+'/Intermediate/Geospatial_Data/weighted_grid_of_points.csv')
    starting_ID = weighted_grid_of_points['ID'].max()+1
    #if not os.path.exists(gisele_folder+'/Case studies/'+case_study+'/Input/Roads_points'):
    ROADS_unfiltered = gpd.read_file(gisele_folder+'/Case studies/'+case_study+'/Intermediate/Geospatial_Data/Roads.shp')
    ROADS_unfiltered = ROADS_unfiltered.to_crs(crs)

    ROADS = ROADS_unfiltered[ROADS_unfiltered['highway'].isin(accepted_road_types)]
    #ROADS = ROADS_unfiltered
    ROADS = MultiLine_to_Line(ROADS)
    all_points = gpd.GeoDataFrame()
    gdf_ROADS, ROADS_segments = create_roads2(ROADS, all_points, crs)
    gdf_ROADS.crs = crs
    ROADS_segments.crs = crs
    #else:
        #ROADS_segments=gpd.read_file(gisele_folder + '/Case studies/' + case_study + '/Input/Roads_lines/Roads_lines.shp')
        #gdf_ROADS=gpd.read_file(gisele_folder + '/Case studies/' + case_study + '/Input/Roads_points/Roads_points.shp')
    ### PROBLEM WITH gdf_ROADS -> THERE ARE MULTI_POINTS

    print(len(gdf_ROADS))
    MP = MultiPolygon([p for p in Clusters['geometry']])
    nodes = ROADS_segments.ID1.to_list() + ROADS_segments.ID2.to_list()
    nodes = [int(i) for i in nodes]
    occurence = Counter(nodes)

    intersection_IDs = []
    terminal_IDs = []
    for i in occurence:
        if occurence[i] == 1:
            terminal_IDs.append(i)
        elif occurence[i] > 2:
            intersection_IDs.append(i)
    new_nodes = terminal_IDs + intersection_IDs

    Substations = new_nodes
    Nodes = gdf_ROADS.copy()
    Nodes.loc[Nodes['ID'].isin(new_nodes), 'Substation'] = 1

    Nodes['inside_clusters'] = [1 if MP.contains(row['geometry']) else 0 for i,row in Nodes.iterrows()]
    Lines = ROADS_segments.copy()
    Lines.ID1 = Lines.ID1.astype(int)
    Lines.ID2 = Lines.ID2.astype(int)

    Lines_marked = Lines.copy()
    conn_param = 0
    New_Lines = gpd.GeoDataFrame()
    while not Lines.empty:
        nodes = Lines.ID1.to_list() + Lines.ID2.to_list()
        # print(nodes)
        nodes = [int(i) for i in nodes]
        Substations = list(set(Substations) & set(nodes))
        current_node = int(Substations[0])
        #print(Substations)
        #print(current_node)
        no_lines = False
        tot_length = 0
        tot_cost = 0
        id1 = current_node
        while not no_lines:

            next_index = Lines.index[Lines['ID1'] == current_node].to_list()
            # check if there actually is a next node
            if next_index:
                next_index = next_index[0]  # i only care about the first line if there are many
                next_node = Lines.loc[next_index, 'ID2']

                tot_length = tot_length + Lines.loc[next_index, 'length']
                tot_cost = tot_cost + Lines.loc[next_index, 'length']
                Lines.drop(index=next_index, inplace=True)
            else:
                next_index = Lines.index[Lines['ID2'] == current_node].to_list()
                if next_index:
                    next_index = next_index[0]  # i only care about the first line if there are many
                    next_node = Lines.loc[next_index, 'ID1']
                    tot_length = tot_length + Lines.loc[next_index, 'length']
                    tot_cost = tot_cost + Lines.loc[next_index, 'length']
                    Lines.drop(index=next_index, inplace=True)
                else:
                    no_lines = True  # there are no lines starting from this node

            if not no_lines:
                is_substation = Nodes.loc[Nodes.ID == int(next_node), 'Substation'] == 1
                is_inside = int(Nodes.loc[Nodes.ID==int(next_node),'inside_clusters'])
                if is_inside == 1:
                    max_tot_length = resolution_LV/1000
                else:
                    max_tot_length = resolution_MV/1000
                Lines_marked.loc[next_index, 'Conn_param'] = conn_param
                if is_substation.values[0]:
                    cond = False
                    Point1 = Nodes.loc[Nodes['ID'] == int(id1), 'geometry'].values[0]
                    Point2 = Nodes.loc[Nodes['ID'] == int(next_node), 'geometry'].values[0]
                    geom = LineString([Point1, Point2])
                    Data = {'ID1': id1, 'ID2': next_node, 'Cost': tot_cost, 'length': tot_length, 'geometry': geom,
                            'Conn_param': conn_param}
                    New_Lines = New_Lines.append(Data, ignore_index=True)
                    current_node = next_node
                    tot_length = 0
                    tot_cost = 0
                    id1 = current_node
                    conn_param = conn_param + 1
                elif tot_length > max_tot_length:
                    Point1 = Nodes.loc[Nodes['ID'] == int(id1), 'geometry'].values[0]
                    Point2 = Nodes.loc[Nodes['ID'] == int(next_node), 'geometry'].values[0]
                    geom = LineString([Point1, Point2])
                    Data = {'ID1': id1, 'ID2': next_node, 'Cost': tot_cost, 'length': tot_length, 'geometry': geom,
                            'Conn_param': conn_param}
                    New_Lines = New_Lines.append(Data, ignore_index=True)
                    current_node = next_node
                    tot_length = 0
                    tot_cost = 0
                    id1 = current_node
                    conn_param = conn_param + 1
                else:
                    current_node = next_node

    New_Lines.crs = Lines.crs
    new_lines = []
    for i, row in New_Lines.iterrows():
        actual_Lines = Lines_marked.loc[Lines_marked['Conn_param'] == row['Conn_param'], 'geometry']
        new_line = MultiLineString([actual_Lines.values[i] for i in range(len(actual_Lines))])
        new_lines.append(new_line)

    New_Lines.geometry = new_lines
    # New_Lines.to_file(r'New_Lines')

    new_nodes = New_Lines.ID1.to_list() + New_Lines.ID2.to_list()
    New_Nodes = gdf_ROADS[gdf_ROADS['ID'].isin(new_nodes)]
    New_Nodes.reset_index(inplace=True)
    for i, row in New_Nodes.iterrows():
        id = int(i)
        New_Nodes.loc[i, 'ID'] = id
        New_Lines.loc[New_Lines['ID1'] == row['ID'], 'ID1'] = id
        New_Lines.loc[New_Lines['ID2'] == row['ID'], 'ID2'] = id

    New_Nodes.ID+=starting_ID
    New_Lines.ID1+=starting_ID
    New_Lines.ID2+=starting_ID

    drop = New_Lines.loc[New_Lines['ID1'] == New_Lines['ID2'], :]
    if not len(drop)==0:
        New_Lines.drop(index=drop.index, inplace=True)
    print(len(New_Nodes))
    New_Lines.to_file(gisele_folder + '/Case studies/' + case_study + '/Intermediate/Geospatial_Data/Roads_lines')
    New_Nodes.to_file(gisele_folder + '/Case studies/' + case_study + '/Intermediate/Geospatial_Data/Roads_points')

    return New_Nodes, New_Lines

def Merge_Roads_GridOfPoints(gisele_folder,case_study):
    road_points = gpd.read_file(gisele_folder+'/Case studies/'+case_study+'/Intermediate/Geospatial_Data/Roads_points/Roads_points.shp')
    weighted_grid_points = pd.read_csv(gisele_folder+'/Case studies/'+case_study+'/Intermediate/Geospatial_Data/weighted_grid_of_points.csv')
    weighted_grid_points['Type'] = 'Standard'
    road_points['Type'] = 'Road'
    road_points.drop(columns=['geometry'],inplace=True)
    weighted_grid_points_with_roads = weighted_grid_points.append(road_points)
    weighted_grid_points_with_roads[['X','Y','ID','Elevation','Type','Population','Weight','Elevation']].\
        to_csv(gisele_folder+'/Case studies/'+case_study+'/Intermediate/Geospatial_Data/weighted_grid_of_points_with_roads.csv')

def improve_connection(MV_grid,all_points,pts,grid_1,grid_2,line_bc):
    graph = nx.Graph()
    for i, row in MV_grid.iterrows():
        graph.add_edge(int(row.ID1), int(row.ID2),length = row.Length, weight=row.Cost)

    for i in pts:
        if not graph.has_edge(i[0], i[1]):
            p1 = all_points.loc[all_points.ID==i[0],'geometry'].values[0]
            p2 = all_points.loc[all_points.ID == i[1], 'geometry'].values[0]
            av_weight = (all_points.loc[all_points.ID==i[0],'Weight'].values[0]+all_points.loc[all_points.ID==i[1],'Weight']
                         .values[0])/2
            graph.add_edge(int(i[0]), int(i[1]), length=p1.distance(p2)/1000, weight = p1.distance(p2)/1000*av_weight*line_bc)
    T_metric1 = metric_closure(graph, weight='weight')
    T_metric = metric_closure(graph, weight='length')
    min_dist = 100000
    terminal_points = []
    for ind1, row1 in grid_1.iterrows():
        for ind2, row2 in grid_2.iterrows():
            id1 = int(row1['ID'])
            id2 = int(row2['ID'])
            dist = T_metric[id1][id2]['distance']
            if dist < min_dist:
                min_dist = dist
                terminal_points = [id1, id2]
    points = T_metric[terminal_points[0]][terminal_points[1]]['path']
    steps = len(points)
    new_path = []
    for i in range(0, steps - 1):
        new_path.append(points[i + 1])

    pts = list(zip(points, new_path))
    Line = array_to_LineString(all_points, pts)
    length = min_dist*1000
    cost = T_metric1[terminal_points[0]][terminal_points[1]]['distance']
    return Line,length,cost,terminal_points

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

def calculate_mg(gisele_folder,case_study,crs):
    case_folder = gisele_folder + '/Case studies/' + case_study
    data_folder = case_folder + '/Input/MILP/all_data'
    Nodes = pd.read_csv(data_folder + '/All_Nodes.csv')
    n_clusters = Nodes['Cluster'].max()
    clusters_list = [*range(1,n_clusters+1)]
    cluster_powers = [Nodes.loc[Nodes['Cluster'] == i, 'MV_Power'].sum() for i in range(1,n_clusters+1)]
    cluster_population = [Nodes.loc[Nodes['Cluster'] == i, 'Population'].sum() for i in range(1,n_clusters+1)]
    clusters_list=pd.DataFrame({'Cluster':clusters_list,'Population': cluster_population,'Load [kW]': cluster_powers})
    clusters_list.to_csv(case_folder+'/Output/clusters_list.csv')
    input_profile = pd.read_csv(gisele_folder+'/general_input/Load Profile.csv').round(4)
    config = pd.read_csv(case_folder+'/Input/Configuration.csv',index_col='Parameter')
    wt=config.loc['wt','Value']
    grid_lifetime = int(config.loc['grid_lifetime','Value'])
    Nodes_gdf = gpd.GeoDataFrame(Nodes, geometry=gpd.points_from_xy(Nodes.X, Nodes.Y),
                              crs=crs)

    yearly_profile, years, total_energy = load(clusters_list,
                                               grid_lifetime,
                                               input_profile)
    mg = sizing(yearly_profile, clusters_list, Nodes_gdf, wt,n_mg=n_clusters)

    mg.to_csv(case_folder+'/Output/Microgrid.csv')

def calculate_mg_multiobjective(gisele_folder, case_study, crs):
    pass