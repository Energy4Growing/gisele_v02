from collections import Counter
from itertools import combinations
from scipy.spatial import Delaunay
from gisele.geneticalgorithm_github import geneticalgorithm as ga
from gisele.Secondary_substations import *
from shapely.geometry import Point, MultiPoint,LineString,MultiLineString
from shapely.ops import split,nearest_points
import networkx as nx
from gisele.Steiner_tree_code import *
def genetic2(clustered_points,points_new_graph,distance_matrix,n_clusters,graph):
    clustered_points.reset_index(drop=True,inplace=True)
    lookup_edges = [i for i in graph.edges]
    dim = len(lookup_edges)-1
    dist_matrix_df = pd.DataFrame(distance_matrix,columns = [i for i in points_new_graph['ID']],index = [i for i in points_new_graph['ID']])
    #initial_solution = np.array(clustered_points['Cluster'].to_list())
    varbound=np.array([[0,dim]]*(n_clusters-1))
    #lookup_edges = [i for i in graph.edges([190,184,29,171,202,201,206,205,209,210,22,221,231,127,235,244,230,229,228,220,210,227,215,216,226,234,204,198,197,56,194,191,179])]
    #dim = len(lookup_edges) - 1
    def fitness(X):
        T = graph.copy()
        length_deleted_lines = 0
        count=Counter(X)
        for i in count:
            if count[i]>1:
                return 1000000 # this is in case it is trying to cut the same branch more than once
        for i in X:
            delete_edge = lookup_edges[int(i)]
            length_deleted_lines += graph[lookup_edges[int(i)][0]][lookup_edges[int(i)][1]]['weight']['distance']
            T.remove_edge(*delete_edge)
        islands = [c for c in nx.connected_components(T)]
        cost = 0
        penalty = 0
        for i in range(len(islands)):
            subgraph = T.subgraph(islands[i])
            subset_IDs = [i for i in subgraph.nodes]
            population =points_new_graph[points_new_graph['ID'].isin(subset_IDs)]['Population'].sum()
            power = population*0.7*0.3
            if power < 25:
                cost += 1500
            elif power < 50:
                cost += 2300
            elif power < 100:
                cost += 3500
            else:
                cost += 100000
            sub_dist_matrix = dist_matrix_df.loc[subset_IDs, subset_IDs]
            max_dist = sub_dist_matrix.max().max()
            if max_dist >1000:
                penalty = penalty+ 50000+ (max_dist-500)*25

        cost = cost - length_deleted_lines/1000*10000 # divided by 1000 for m->km and the multiplied by 10000euro/km

        if penalty>0:
            return penalty
        else:
            return cost

    algorithm_param = {'max_num_iteration': 1000, 'population_size': 40, 'mutation_probability': 0.1,
                       'elit_ratio': 0.025, 'crossover_probability': 0.6, 'parents_portion': 0.25,
                       'crossover_type': 'one_point', 'max_iteration_without_improv': 100}
    model = ga(function=fitness, dimension=n_clusters-1, variable_type='int', variable_boundaries=varbound,
             function_timeout=10000,algorithm_parameters=algorithm_param)
    model.run()
    cut_edges = model.best_variable
    T=graph.copy()
    for i in cut_edges:
        delete_edge = lookup_edges[int(i)]
        T.remove_edge(*delete_edge)

    islands = [c for c in nx.connected_components(T)]
    for i in range(len(islands)):
        subgraph = T.subgraph(islands[i])
        subset_IDs = [i for i in subgraph.nodes]
        clustered_points.loc[clustered_points['ID'].isin(subset_IDs),'Cluster']=i

    return clustered_points, cut_edges
def delaunay_test(graph,new_points,new_lines):
    tocki = new_points['geometry'].values
    number_points = new_points.shape[0]
    arr = np.zeros([number_points,2])
    counter=0
    for i in tocki:
        x = i.xy[0][0]
        y=i.xy[1][0]
        arr[counter,0] = x
        arr[counter,1] = y
        counter+=1
    tri = Delaunay(arr)
    triangle_sides = tri.simplices
    final_sides = []
    for i in triangle_sides:
        a=i[0]
        b=i[1]
        c=i[2]
        if a>b:
            final_sides.append((i[0],i[1]))
        else:
            final_sides.append((i[1], i[0]))
        if b>c:
            final_sides.append((i[1],i[2]))
        else:
            final_sides.append((i[2], i[1]))
        if a>c:
            final_sides.append((i[0],i[2]))
        else:
            final_sides.append((i[2], i[0]))
    final_sides2 = list(set(final_sides))
    new_lines_old=new_lines.copy() # dataframe without the new possible connections

    if not nx.is_empty(graph): # this is for the standard case with roads in the cluster
        for i,j in final_sides2:
            point1 = new_points.loc[new_points['order']==i,'geometry'].values[0]
            point2 = new_points.loc[new_points['order'] == j, 'geometry'].values[0]
            id1 = int(new_points.loc[new_points['order'] == i, 'ID'])
            id2 = int(new_points.loc[new_points['order'] == j, 'ID'])
            length = point1.distance(point2)
            line = LineString([point1, point2])
            if length<500 and not graph.has_edge(id1,id2) and ((sum([line.intersects(line1) for line1 in new_lines_old.geometry]) == 0) or
                    ((new_points.loc[new_points['ID'] == id1, 'pop_bool'] == 0).values[0]) or
                    (new_points.loc[new_points['ID'] == id2, 'pop_bool'] == 0).values[0]):
                graph.add_edge(id1,id2 , weight=length, length=length)

                data_segment = {'ID1': [id1], 'ID2': [id2], 'length': [point1.distance(point2) / 1000],
                                'geometry': [line], 'Type': ['Colateral']}
                new_lines = new_lines.append(gpd.GeoDataFrame(data_segment))
    else: # this is for the case without roads in the cluster, just create the lines in a straightforward way
        new_points = new_points.reset_index()
        for i, j in final_sides2:
            point1 = new_points.loc[new_points.index == i, 'geometry'].values[0]
            point2 = new_points.loc[new_points.index== j, 'geometry'].values[0]
            id1 = int(new_points.loc[new_points.index == i, 'ID'].values[0])
            id2 = int(new_points.loc[new_points.index== j, 'ID'].values[0])
            length = point1.distance(point2)
            line = LineString([point1, point2])
            graph.add_edge(id1, id2, weight=length, length=length)

            data_segment = {'ID1': [id1], 'ID2': [id2], 'length': [point1.distance(point2) / 1000],
                            'geometry': [line], 'Type': ['Colateral']}
            new_lines = new_lines.append(gpd.GeoDataFrame(data_segment))
    return graph,new_lines

def create_clean_graph(graph,points,terminal_points,T_metric,crs):
    '''This function returns a graph that is composed only of population nodes(translated on the roads) and the intersection points(points
    which are present in the existing graph more than 2 times. The idea is to start cutting the highest cost lines as the path
    to a much better clustering that includes the actual electrical distances.'''
    #WORKS
    #
    # STEP 1. Take all the terminal nodes + the intersection nodes

    terminal_IDs = terminal_points['ID'].to_list()
    edges_tuples = [i for i in graph.edges]
    nodes = [edges_tuples[i][0] for i in range(len(edges_tuples))]
    nodes+=[edges_tuples[i][1] for i in range(len(edges_tuples))]
    occurence = Counter(nodes)
    intersection_IDs=[]
    for i in occurence:
        if occurence[i]>2 and not i in terminal_IDs:
            intersection_IDs.append(i)
    new_nodes = terminal_IDs + intersection_IDs

    # STEP 2. Create the new graph
    start_node = new_nodes[0]
    #start_node = 154
    current_node = start_node
    graph_copy = graph.copy()
    new_graph=nx.Graph()
    terminal_IDs_2 = terminal_IDs.copy()
    unique_nodes=new_nodes.copy()
    while True:
        try:
            next_node = [i for i in graph_copy[current_node]][0]
            #print(next_node)
        except:
            print('A terminal node has been reached, back to the set of points')
            if current_node in terminal_IDs_2:
                terminal_IDs_2.remove(current_node)
                #print('Node ' + str(current_node) + ' was deleted.')
            #print('Next point is '+str(unique_nodes[0]))
            start_node = unique_nodes[0]
            current_node = start_node
            next_node = [i for i in graph_copy[start_node]][0]
        if next_node in new_nodes:
            new_graph.add_edge(start_node, next_node, weight=T_metric[start_node][next_node])
            #print('add ' + str(start_node) + ' and ' + str(next_node))
            graph_copy.remove_edge(current_node,next_node)
            #print('remove '+str(current_node)+' and ' + str(next_node))
            if start_node in terminal_IDs_2:
                terminal_IDs_2.remove(start_node)
                print('Node '+ str(start_node)+' was deleted.')
            start_node = next_node
            current_node = start_node
        else:
            graph_copy.remove_edge(current_node, next_node)
            #print('remove ' + str(current_node) + ' and ' + str(next_node))
            current_node = next_node

        if nx.is_empty(graph_copy):
            break
        new_edges = [i for i in graph_copy.edges]
        unique_nodes = list(set([new_edges[i][0] for i in range(len(new_edges))] + [new_edges[i][1] for i in range(len(new_edges))]))
        unique_nodes = list(set(unique_nodes) & set(new_nodes))
    new_edges = [i for i in new_graph.edges]
    new_lines=gpd.GeoDataFrame()
    for j in new_edges:
        point1 = points.loc[points['ID'] == j[0],'geometry'].values[0]
        point2 = points.loc[points['ID'] == j[1],'geometry'].values[0]
        length = new_graph[j[0]][j[1]]['weight']['distance']
        new_lines=new_lines.append({'geometry': LineString([point1,point2]),'length': length},ignore_index=True)
    new_lines.crs = crs


    new_points = points[points['ID'].isin(new_nodes)]
    return new_lines, new_points, new_graph
def connect_unconnected_graph(graph,lines,points,weight):
    if nx.is_connected(graph):
        return graph,lines
    else:
        islands = [c for c in nx.connected_components(graph)]
        for i in range(len(islands)):
            for j in range(i+1,len(islands)):
                subgraph_1 = [val for val in islands[i]]
                subgraph_2 = [val for val in islands[j]]
                points_s1 = points.loc[points['ID'].isin(subgraph_1),:]
                points_s2 = points.loc[points['ID'].isin(subgraph_2), :]
                multi_point1= MultiPoint([row['geometry'] for i, row in points_s1.iterrows()])
                multi_point2 = MultiPoint([row['geometry'] for i, row in points_s2.iterrows()])
                closest_points = nearest_points(multi_point1,multi_point2)
                distance = multi_point1.distance(multi_point2)#in km
                id_point1 = int(points.loc[points['geometry']==closest_points[0],'ID'])
                id_point2 = int(points.loc[points['geometry'] == closest_points[1], 'ID'])
                lines=lines.append(gpd.GeoDataFrame({'ID1':[id_point1],'ID2':[id_point2],'length':[distance]
                                            , 'geometry':[LineString([closest_points[0],closest_points[1]])]}))
                graph.add_edge(id_point1,id_point2,weight = distance*weight,length = distance)
    return graph,lines

def add_roads_points_to_gdf(gdf,gdf_roads,c_grid,cluster):
    nodes_in_lines = c_grid.ID1.to_list() + c_grid.ID2.to_list()
    nodes_in_lines = list(set(nodes_in_lines)) # have the unique values
    gdf_roads['Cluster'] = cluster


    for i in nodes_in_lines:
        if i in gdf_roads.ID.to_list():
            a = gdf_roads.loc[gdf_roads['ID']==i,:]

            a.loc[a['Weight'] == i, 'ID'] = 1 # change it in the future
            gdf = gdf.append(a)
    return gdf

def fix_lines_intersecting(lines,points):
   next_node = points.shape[0]
   multi_line = MultiLineString([row['geometry'] for i,row in lines.iterrows()])
   index_to_drop=[]
   for line1, line2 in combinations([line for line in multi_line], 2):
       if line1.intersects(line2):
           intersection_point = line1.intersection(line2)
           if Point(line1.coords[0][0],line1.coords[0][1]).distance(intersection_point)>0.0001 and Point(line1.coords[1][0],line1.coords[1][1]).distance(intersection_point)>0.0001:
                print(intersection_point)
                line1_1 = Point(line1.coords[0][0],line1.coords[0][1])
                line1_2 = Point(line1.coords[1][0],line1.coords[1][1])
                line2_1 = Point(line2.coords[0][0],line2.coords[0][1])
                line2_2 = Point(line2.coords[1][0],line2.coords[1][1])
                point1_ID = int(points.loc[points['geometry']==line1_1,'ID'])
                point2_ID = int(points.loc[points['geometry'] == line1_2, 'ID'])
                point3_ID = int(points.loc[points['geometry'] == line2_1, 'ID'])
                point4_ID = int(points.loc[points['geometry'] == line2_2, 'ID'])
                dist1 = intersection_point.distance(line1_1)/1000
                dist2 = intersection_point.distance(line1_2)/1000
                dist3 = intersection_point.distance(line2_1)/1000
                dist4 = intersection_point.distance(line2_2)/1000 # to km
                new_line1 = LineString([intersection_point,line1_1])
                new_line2 = LineString([intersection_point, line1_2])
                new_line3 = LineString([intersection_point, line2_1])
                new_line4 = LineString([intersection_point, line2_2])
                points=points.append(gpd.GeoDataFrame({'ID':[next_node],'X':[intersection_point.coords[0][0]],
                                    'Y':[intersection_point.coords[0][1]],'Weight':[1],'Elevation':[1000],'geometry':[intersection_point]}))

                # the indices of the 2 lines that need to be deleted
                index_line1 = lines.loc[lines['geometry']==line1,:].index.values[0]
                index_line2 = lines.loc[lines['geometry']==line2,:].index.values[0]
                # add the 4 new lines
                Data = {'ID1':[next_node]*4,'ID2':[point1_ID,point2_ID,point3_ID,point4_ID],'length':[dist1,dist2,dist3,dist4],
                        'geometry':[new_line1,new_line2,new_line3,new_line4]}
                lines = lines.append(gpd.GeoDataFrame(Data))
                index_to_drop.append([index_line1,index_line2])
                next_node += 1
   for i in index_to_drop:
       lines.drop(index=i, inplace=True)
                # the new point is added, now just delete the old 2 lines and create 4 new lines.
   return lines,points
def fix_roads(lines,points,critical_points,critdist):
    ''' The goal of this function is to take points that are critical, actually points that are very close to a road,
    but not really on a road. Then, those points are translated on to the given road, creating a much more realistic mash
    that represents the roads with lines and points.'''

    for index,row in critical_points.iterrows():
        id1 = int(row['ID_point']) # this is the critical point
        if not row['Repeating_line']:
            id3 = int(row['ID1_line']) # this is the first point of the line on which we project
            id4 = int(row['ID2_line']) # this is the second point of the line on which we project
        else:
            multi_point = MultiPoint([row['geometry']])
            multi_line = MultiLineString([row['geometry'] for i, row in lines.iterrows()])
            lt = [each for each in list(map(lambda x: plt(x, multi_line, critdist), multi_point)) if each != False]
            id3 = lines.loc[lines['geometry'] == lt[0][1], 'ID1'].values[0]
            id4 = lines.loc[lines['geometry'] == lt[0][1], 'ID2'].values[0]
        point = points.loc[points['ID'] == id1, 'geometry']
        line = lines.loc[(lines['ID1'] == id3) & (lines['ID2'] == id4), 'geometry']
        #line2 = lines.loc[(lines['ID1'] == id1) & (lines['ID2'] == id2), 'geometry']
        point_geom = point.values[0]
        line_geom = line.values[0]
        #line2_geom = line2.values[0]
        #intersect_point = line_geom.intersection(line2_geom)
        #if intersect_point.is_empty: # it's empty if there is no intersection, then i find the closest distance.
        intersect_point = line_geom.interpolate(line_geom.project(point_geom))
        points.loc[points['ID'] == id1, 'geometry'] = intersect_point
        try:
            lines.loc[lines['ID1']==id1,'geometry'] = LineString([intersect_point, points.loc[points['ID'] ==
                                                            int(lines.loc[lines['ID1']==id1,'ID2']), 'geometry'].values[0]])
        except:
            pass
        try:
            lines.loc[lines['ID2'] == id1, 'geometry'] = LineString([intersect_point, points.loc[points['ID'] ==
                                                    int(lines.loc[lines['ID2'] == id1, 'ID1']), 'geometry'].values[0]])
        except:
            pass

        point1=points.loc[points['ID'] == id4, 'geometry'].values[0]
        dist1=point1.distance(intersect_point)/1000
        point2 = points.loc[points['ID'] == id3, 'geometry'].values[0]
        dist2 = point2.distance(intersect_point)/1000
        lines=lines.append(pd.Series({'ID1':id4,'ID2':id1,'length':dist1,'geometry':LineString([point1,intersect_point])},name=lines.index.max()+1))
        lines = lines.append(pd.Series({'ID1': id1, 'ID2': id3, 'length': dist2, 'geometry': LineString([intersect_point, point2])},name = lines.index.max()+1))


        index_to_drop = lines[(lines['ID1'] == id3) & (lines['ID2'] == id4)].index[0]
        lines.drop(index=index_to_drop, inplace=True)

    lines = lines.reset_index(drop=True)
    #for index,row in critical_points.iterrows():
    #    line_id1 = int(row['ID1_line'])
    #    line_id2 = int(row['ID2_line'])
    #    try:
    #        index_to_drop = lines[(lines['ID1']==line_id1) & (lines['ID2']==line_id2)].index[0]
    #        lines.drop(index=index_to_drop, inplace=True)
    #    except:
    #        try:
    #            index_to_drop = lines[(lines['ID1'] == line_id2) & (lines['ID2'] == line_id1)].index[0]
    #            lines.drop(index=index_to_drop, inplace=True)
    #        except:
    #            print('Line is already dropped ( there were 2 critical points on this one')

    return lines,points
def plt(point, multiLine, threshold):
    '''Function that is used in distance_point_to_roads. It return the distance if a point is too close to a road and needs to be
    translated ( less than the threshold ).'''
    for line in multiLine:
        dist = point.distance(line)
        if dist<threshold and dist>0:
            return (point,line, dist)
    return False
def distance_point_to_roads(segments,line_gdf,critdist,crs):
    '''Returns "critical_points", which is a geodataframe containing points which are too close to the existing roads,
    which means they should be translated on to the roads. Also, it contains data on which road they are close to.'''
    multi_point=MultiPoint([row['geometry'] for i,row in line_gdf.iterrows()])
    multi_line = MultiLineString([row['geometry'] for i,row in segments.iterrows()])
    lt = [each for each in list(map(lambda x: plt(x, multi_line, critdist), multi_point)) if each != False]
    points = [tup[0] for tup in lt]
    lines = [tup[1] for tup in lt]
    points_ID=[]
    lines_ID1 = []
    lines_ID2 = []
    repeating_line = []
    for point in points:
        id = line_gdf.loc[line_gdf['geometry']==point,'ID'].values[0]
        points_ID.append(id)
    #make sure that it's not just a coincidence that a regular point is close to a road


    for line in lines:
        id1 = segments.loc[segments['geometry'] == line, 'ID1'].values[0]
        id2 = segments.loc[segments['geometry'] == line, 'ID2'].values[0]
        repeat = False
        if id1 in lines_ID1 and id2 in lines_ID2:
            if lines_ID1.index(id1) == lines_ID2.index(id2):
                repeat = True
        lines_ID1.append(id1)
        lines_ID2.append(id2)
        repeating_line.append(repeat)
    length = np.array([tup[2] for tup in lt])
    critical_points=pd.DataFrame({'ID_point':points_ID,'ID1_line':lines_ID1,'ID2_line':lines_ID2,'length':length,'Repeating_line': repeating_line})
    critical_points=gpd.GeoDataFrame(critical_points,geometry = points)
    critical_points.crs = crs
    return critical_points

def process_roads(Roads,roads_weight,crs,directory,max_length_segment,simplify_coef,crit_dist,starting_ID=0):
    crs_str = 'EPSG:' + str(crs)
    Roads = Roads[Roads['highway'] != 'path']
    Roads = MultiLine_to_Line(Roads)
    roads = Roads.simplify(simplify_coef)
    gdf_roads = gpd.GeoDataFrame(geometry=roads)
    os.chdir('..')
    gdf_roads.to_file(directory + '/roads_simplified')
    w = starting_ID # ID of the lines

    line_vertices = pd.DataFrame(
        index=pd.Series(range(w, w + len(gdf_roads.index))),
        columns=['ID', 'X', 'Y', 'ID_line', 'Weight', 'Elevation'], dtype=int)
    # create geodataframe with all the segments that compose the road
    segments = gpd.GeoDataFrame(columns=['geometry', 'ID1', 'ID2'])
    k = 0

    x = 0
    # Create a set of points and lines from the simplified roads.
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
    geometry = [Point(xy) for xy in
                zip(line_vertices['X'], line_vertices['Y'])]
    line_gdf = gpd.GeoDataFrame(line_vertices, crs=crs_str,
                                geometry=geometry)
    line_vertices.loc[:, 'Elevation'] = 1000
    segments.crs = crs
    #print('test')
    line_gdf.to_file(directory + '/points_first')
    segments.to_file(directory+'/lines_first')
    #line_gdf.to_file('Testing_strategy' + dir + '/Points_first')
    #segments.to_file('Testing_strategy' + dir + '/Lines_first')
    critical_points = distance_point_to_roads(segments, line_gdf, crit_dist,crs)
    #print('critical points found')
    segments, line_gdf = fix_roads(segments, line_gdf, critical_points, crit_dist)
    #print('roads fixed')
    #line_gdf.to_file('Testing_strategy' + dir + '/Points_fixed1')
    #segments.to_file('Testing_strategy' + dir + '/lines_fixed1')
    segments, line_gdf = fix_lines_intersecting(segments, line_gdf)
    #print('create points in intersections')
    #line_gdf.to_file('Testing_strategy' + dir + '/points_fixed2')
    #segments.to_file('Testing_strategy' + dir + '/lines_fixed2')
    segments = segments.reset_index(drop=True)
    new_lines = gpd.GeoDataFrame()

    next_node = line_gdf.ID.max()+1  # just logic stuff
    new_points = line_gdf.copy()
    # here make sure that when we find the points in the middle of a longer line, translate them on to the original Roads file ( non-simplified)
    # and use that point - to get a more realistic representation of the roads via the grid of points
    for index, segment in segments.iterrows():
        if segment['length'] > max_length_segment/1000:
            #print('da')
            number_segments = ceil(segment['length'] / (max_length_segment/1000))
            num_points = number_segments - 1
            new_lengths = segment['length'] / number_segments
            splitter = MultiPoint([segment['geometry'].interpolate((i / number_segments), normalized=True) for i in
                                   range(1, number_segments)])
            for i in range(num_points):
                if i == 0:
                    id1 = int(segment['ID1'])
                else:
                    id1 = next_node
                next_node += 1
                data_point = {'ID': [next_node], 'X': [splitter[i].xy[0][0]], 'Y': [splitter[i].xy[1][0]],
                              'Weight': [1], 'Elevation': [1000], 'geometry': [splitter[i]]}
                new_points = new_points.append(gpd.GeoDataFrame(data_point))
                line = LineString([new_points.loc[new_points['ID'] == id1, 'geometry'].values[0], splitter[i]])
                data_segment = {'ID1': [id1], 'ID2': [next_node], 'length': [new_lengths], 'geometry': [line]}
                new_lines = new_lines.append(gpd.GeoDataFrame(data_segment))
            line = LineString([new_points.loc[new_points['ID'] == next_node, 'geometry'].values[0],
                               new_points.loc[new_points['ID'] == int(segment['ID2']), 'geometry'].values[0]])
            data_segment = {'ID1': [next_node], 'ID2': [int(segment['ID2'])], 'length': [new_lengths],
                            'geometry': [line]}
            new_lines = new_lines.append(gpd.GeoDataFrame(data_segment))
        else:
            new_lines = new_lines.append(segments.loc[index, :])

    new_lines.crs = crs
    #new_lines.to_file(directory+ '/new_lines', index=False)
    return new_points,new_lines
def optimize(crs, resolution, load_capita, pop_per_household, road_coef,Clusters, case_study, LV_distance, ss_data,
             landcover_option,gisele_dir,roads_weight,run_genetic,max_length_segment,simplify_coef,crit_dist,LV_base_cost):
    dir_input = r'Case studies/' + case_study + '/Intermediate/Geospatial_Data'

    dir_output = '/Case studies/' + case_study + '/Output'
    grid_of_points = pd.read_csv(dir_input + '/weighted_grid_of_points_with_roads.csv')
    grid_of_points_GDF = gpd.GeoDataFrame(grid_of_points,
                                     geometry=gpd.points_from_xy(grid_of_points.X, grid_of_points.Y), crs=crs)

    Starting_node = int(grid_of_points_GDF['ID'].max()+1)
    # Create a clusters.exe file
    LV_resume = pd.DataFrame()
    LV_grid = gpd.GeoDataFrame()
    MV_grid = gpd.GeoDataFrame()
    secondary_substations=gpd.GeoDataFrame()
    all_houses = gpd.GeoDataFrame()
    #Clusters=Clusters[Clusters['cluster_ID']==18]
    for index, row in Clusters.iterrows():
        os.chdir(gisele_dir)
        print('WORKING ON CLUSTER '+str(row['cluster_ID']))
        dir = gisele_dir + '/Case studies/' + case_study +'/Intermediate/Communities/' + str(row['cluster_ID'])
        clus = row['cluster_ID']
        if not os.path.exists(dir):
            os.makedirs(dir)
            os.makedirs(dir + '/grids')
        area = row['geometry']
        area_buffered = area
        # area_buffered = row['geometry'].buffer((resolution_MV * 0.1 / 11250) / 2)
        area_list = [area_buffered]
        grid_of_points = create_grid(crs, resolution, area)
        grid_of_points.to_file(dir + '/points.shp')
        #min_x, min_y, max_x, max_y = area.bounds
        #area_for_roads = geometry.Polygon(
        #    [geometry.Point(min_x, min_y), geometry.Point(min_x, max_y), geometry.Point(max_x, max_y),
        #     geometry.Point(max_x, min_y)])
        #streets = gpd.read_file(dir_input + '/Roads.shp')
        area_for_roads = row['geometry'].buffer(resolution)
        road_points = gpd.read_file(dir_input+'/Roads_points/Roads_points.shp')
        road_points = road_points[['X', 'Y', 'ID', 'Weight', 'Elevation','geometry']]
        road_points = gpd.clip(road_points,area_for_roads)
        road_lines = gpd.read_file(dir_input+'/Roads_lines/Roads_lines.shp')
        road_lines = road_lines[(road_lines['ID1'].isin(road_points.ID.to_list()) &
                                                 road_lines['ID2'].isin(road_points.ID.to_list()))]
        # OPEN THE RASTERS FOR THE SPECIFIC REGION WHERE OUR CLUSTER IS
        Population = rasterio.open(dir_input + '/Population_' + str(crs) + '.tif')
        Elevation = rasterio.open(dir_input + '/Elevation_' + str(crs) + '.tif')
        Slope = rasterio.open(dir_input + '/Slope_' + str(crs) + '.tif')
        LandCover = rasterio.open(dir_input + '/LandCover_' + str(crs) + '.tif')

        # POPULATE THE GRID OF POINTS
        coords = [(x, y) for x, y in zip(grid_of_points.X, grid_of_points.Y)]
        grid_of_points = grid_of_points.reset_index(drop=True)
        grid_of_points['ID'] = grid_of_points.index

        grid_of_points['Population'] = [x[0] for x in Population.sample(coords)]
        grid_of_points['Elevation'] = [x[0] for x in Elevation.sample(coords)]
        grid_of_points['Slope'] = [x[0] for x in Slope.sample(coords)]
        grid_of_points['Land_cover'] = [x[0] for x in LandCover.sample(coords)]
        # THIS IS JUST A PROXY, NEEDS TO BE PROPERLY SET
        grid_of_points['Protected_area'] = ['FALSE' for x in LandCover.sample(coords)]
        print('Sampling rasters finished')
        grid_of_points.to_file(dir + '/points.shp')
        # perhaps change this in the future
        grid_of_points['Weight'] = 1
        #road_points,road_lines = process_roads(streets_clipped,roads_weight,crs,dir,max_length_segment,simplify_coef,crit_dist,Starting_node)
        #Starting_node += road_points.shape[0] +1
        # STEP 1 -> Translate the points on to the roads to find the "backbone".
        print('FINDING THE BACKBONE')
        Population = grid_of_points[grid_of_points['Population']>0]
        Population['ID'] = [*range(Starting_node, Starting_node + Population.shape[0])]
        Population['pop_bool'] = 1
        Starting_node += Population.shape[0]
        Population.to_file(dir + '/population', index=False)
        if not len(road_points)<2: # normal procedure in case there are roads in the cluster, actually more
            #than 1 road point is needed
            roads_multipoint = MultiPoint([point for point in road_points['geometry']])
            road_points['Population'] = 0
            road_points['pop_bool'] = 0
            road_points_populated=road_points.copy()
            for i, pop in Population.iterrows():
                point = pop['geometry']
                nearest_geoms = nearest_points(point, roads_multipoint)
                closest_road_point = nearest_geoms[1]
                road_points_populated.loc[road_points_populated['geometry'] == closest_road_point, 'Population'] = pop['Population']
                road_points_populated.loc[road_points_populated['geometry'] == closest_road_point, 'pop_bool'] = 1
            road_points_populated.crs = crs
            road_points_populated.to_file(dir + '/road_points_populated', index=False)
            road_points_populated = road_points_populated.set_index('ID', drop=False)
            road_lines = road_lines.set_index(pd.Index([*range(road_lines.shape[0])]))

            graph = nx.Graph()
            for index, ROW in road_lines.iterrows():
                id1 = ROW['ID1']
                id2 = ROW['ID2']
                #print(id1, id2)
                graph.add_edge(id1, id2, weight=ROW['length'] * 1000, length = ROW['length']*1000)
                #print(ROW['length'] * 1000 * roads_weight)
            # the next function is useful if the roads inside the village are not connected. Something smarter is required here.
            graph, new_lines = connect_unconnected_graph(graph, road_lines, road_points_populated,weight=5)
            road_lines.to_file(dir + '/road_lines', index=False)

            populated_points = road_points_populated[road_points_populated['pop_bool'] == 1]
            terminal_nodes = list(populated_points['ID'])
            road_points_populated[road_points_populated['ID'].isin(terminal_nodes)].to_file(dir + '/road_terminal_points',
                                                                      index=False)
            tree = steiner_tree(graph, terminal_nodes)
            path = list(tree.edges)
            grid_routing = gpd.GeoDataFrame()
            counter = 0
            for i in path:
                point1 = road_points_populated.loc[road_points_populated['ID'] == i[0], 'geometry'].values[0]
                point2 = road_points_populated.loc[road_points_populated['ID'] == i[1], 'geometry'].values[0]
                geom = LineString([point1, point2])
                grid_routing = grid_routing.append(gpd.GeoDataFrame({'ID': [counter], 'geometry': [geom]}))
                counter += 1
            grid_routing.crs = crs
            grid_routing.to_file(dir + '/LV_backbone', index=False)

            # STEP 2 -> Do the collaterals, the grid off roads ( new from 17th June)
            print('Connect the houses')
            road_points_backbone = road_points_populated[road_points_populated['ID'].isin(list(tree.nodes))]
            road_points_backbone['Population'] = 0
            road_points_backbone['pop_bool'] = 0
            index = [*range(road_points_backbone['ID'].astype(int).max() + 1,
                            road_points_backbone['ID'].astype(int).max() + 1 + Population.shape[0])]
            Population['ind'] = index
            Population.set_index('ind', inplace=True, drop=True)

            all_points = road_points_backbone.append(Population)
            new_graph = tree.copy()
            for n in new_graph.edges:
                new_graph[n[0]][n[1]]['weight'] = new_graph[n[0]][n[1]]['weight'] * 0.03

            road_lines_copy = road_lines.copy()
            road_lines['Type'] = 'Road'
            all_points['order'] = [*range(all_points.shape[0])]

            new_graph, all_lines = delaunay_test(new_graph, all_points, road_lines)
            new_graph, new_lines = connect_unconnected_graph(new_graph, new_lines, all_points,weight=3)
            terminal_nodes = all_points.loc[all_points['pop_bool'] == 1, 'ID'].to_list()

        else:  # in case there are no roads
            dist_2d_matrix = distance_2d(Population, Population, 'X', 'Y')
            dist_2d_matrix = pd.DataFrame(dist_2d_matrix, columns=Population.ID, index=Population.ID)
            terminal_nodes = Population['ID'].to_list()
            graph = nx.Graph()
            new_graph, all_lines = delaunay_test(graph, Population, road_lines)
            all_points = Population.copy() # in this case, there is no backbone - all_points is just the final houses
        tree_final = steiner_tree(new_graph, terminal_nodes)
        grid_final = gpd.GeoDataFrame()
        path = list(tree_final.edges)
        counter = 0
        for i in path:
            point1 = all_points.loc[all_points['ID'] == i[0], 'geometry'].values[0]
            point2 = all_points.loc[all_points['ID'] == i[1], 'geometry'].values[0]
            geom = LineString([point1, point2])
            grid_final = grid_final.append(gpd.GeoDataFrame({'ID': [counter], 'geometry': [geom]}))
            counter += 1
        grid_final.crs = crs
        grid_final.to_file(dir + '/grid_final', index=False)
        print('Clustering..')
        #end of step 2.
        T_metric = metric_closure(tree_final, weight='length')
        populated_points = all_points[all_points['pop_bool'] == 1]
        lines_new_graph, points_new_graph, new_graph = create_clean_graph(tree_final, all_points, populated_points,
                                                                          T_metric, crs)
        #T_metric = metric_closure(new_graph,weight='distance')
        points_set = all_points.loc[all_points['pop_bool'] == 1, 'ID'].values
        dist_matrix = np.zeros((len(points_set), len(points_set)))
        for i in range(len(points_set)):
            for j in range(len(points_set)):
                if not i == j:
                    dist_matrix[i, j] = T_metric[points_set[i]][points_set[j]]['distance']

        clustering = AgglomerativeClustering(n_clusters=None, affinity='precomputed', linkage='complete',
                                             distance_threshold=2*LV_distance).fit(dist_matrix)


        populated_points.loc[:,'Cluster'] = clustering.labels_
        clustered_points=populated_points.copy()
        clustered_points.to_file( dir + '/Clustered_points', index=False)
        #populated_points['Population'] = [ceil(i) for i in populated_points['Population']]
        populated_points['Population'] = 4
        number_clusters = populated_points['Cluster'].max() + 1
        if number_clusters>1:
            points_new_graph.loc[points_new_graph['Population'] > 0, 'Population'] = 4
            lookup_edges = [i for i in new_graph.edges]

            if run_genetic:
                points_set = points_new_graph.ID.to_list()
                dist_matrix2 = np.zeros((len(points_set), len(points_set)))
                for i in range(len(points_set)):
                    for j in range(len(points_set)):
                        if not i == j:
                            dist_matrix2[i, j] = T_metric[points_set[i]][points_set[j]]['distance']

                clustered_points, cut_edges= genetic2(populated_points, points_new_graph, dist_matrix2, number_clusters, new_graph)
                clustered_points.to_file(dir + '/Clustered_points_after_genetic', index=False)
            else:
                #determine which "not terminal" nodes belong to which cluster
                for number_clus in range(clustered_points['Cluster'].max()+1):
                    subset = clustered_points[clustered_points['Cluster']==number_clus]
                    if len(subset)==1: # if it's just a single house, it's easy.
                        points_new_graph.loc[points_new_graph['ID']==int(clustered_points.loc[clustered_points['Cluster']==number_clus,'ID']), 'Cluster'] =number_clus
                    #elif len(subset)==2:
                        #points_new_graph.loc[points_new_graph['ID'].isin(clustered_points.loc[clustered_points['Cluster']==number_clus,'ID'].to_list()),'Cluster'] = number_clus
                    else: #else, we need to do the procedure
                        edges = nx.edges(new_graph,subset.ID.to_list())
                        edges_nodes = [node for tuple in edges for node in tuple]
                        count=Counter(edges_nodes)
                        terminal_node = [i for i in count if count[i]==1]
                        terminal_node = [i for i in terminal_node if int(points_new_graph.loc[points_new_graph['ID']==i,'pop_bool'])==1] # filter for just the populated nodes
                        for i in range(len(terminal_node)-1):
                            for j in range(i+1,len(terminal_node)):
                                path = T_metric[terminal_node[i]][terminal_node[j]]['path']
                                #TODO fix the issue if the graph can't be simply cut. The following 6 lines of code are locating
                                # those nodes, we should find a way to create another node instead of deleting lines. We need
                                # to add a node in clustered_points and add lines in grid_final. Then, the temprary fix in 757
                                # will not be needed.
                                #a=points_new_graph.loc[points_new_graph['ID'].isin(path),'Cluster'].to_list()
                                #for i in range(len(a)):
                                #    if not isnan(a[i]) and a[i] != number_clus:
                                #        print('dont change')
                                #    else:
                                #        a[i] = number_clus
                                points_new_graph.loc[points_new_graph['ID'].isin(path),'Cluster']=number_clus
                for ind,row in clustered_points.iterrows(): # to fix a possible issue
                    points_new_graph.loc[points_new_graph['ID']==row['ID'],'Cluster']=row['Cluster']
                #cut the edges that are between clusters to form separate LV networks
                #points_new_graph.to_file(dir + '/Testing', index=False)
                cut_edges=[]
                for i in range(len(lookup_edges)):
                    try: # the only way for this to happen is if one of the points is an intersection between various clusters. In that case, we automatically delete those lines.
                        line = lookup_edges[i]
                        point1_cluster = int(points_new_graph.loc[points_new_graph['ID'] == line[0], 'Cluster'])
                        point2_cluster = int(points_new_graph.loc[points_new_graph['ID'] == line[1], 'Cluster'])
                        if not point1_cluster==point2_cluster:
                            cut_edges.append(i)
                    except:
                        cut_edges.append(i)
            # find the weighted centroid of the each cluster
            tree_final=nx.Graph(tree_final) # to unfreeze the graph
            tree_final_copy = tree_final.copy() # to save the full version of the tree, before cutting it.
            for i in cut_edges:
                edge = lookup_edges[int(i)]
                edge_path = nx.dijkstra_path(tree_final,edge[0],edge[1])
                for j in range(len(edge_path)-1):
                    tree_final.remove_edge(*(edge_path[j],edge_path[j+1]))
                    #print('path deleted'+str(edge_path[j])+'-'+str(edge_path[j+1]))
        islands = [c for c in nx.connected_components(tree_final)]
        islands = [i for i in islands if len(i)>1] #otherwise, many empty "islands" will be present
        # There is a specific problem if the aglomerative clustering gives an output that's not coherent with the graph.
        # Basically, it can happen that you can not cut just 1 line and properly divide the LV networks because 1 node overlaps.
        # This is a temporary solution
        if len(islands)>clustered_points['Cluster'].max()+1:
            for i in range(len(islands)):
                subgraph_IDs = list(islands[i])
                clustered_points.loc[clustered_points['ID'].isin(subgraph_IDs),'Cluster']=i
            number_clusters=len(islands)
        for i in range(len(islands)): # for each low voltage network
            print(i)
            subgraph = tree_final.subgraph(islands[i])
            LV_grid_length = 0
            for i in subgraph.edges:
                LV_grid_length += subgraph[i[0]][i[1]]['length']
            check_cluster = True
            #a = nx.algorithms.shortest_paths.dense.floyd_warshall_numpy(subgraph, nodelist=None, weight='weight') # next step, easy calculation of distances
            subset_IDs = [i for i in subgraph.nodes]
            all_points_subset = all_points.loc[all_points['ID'].isin(subset_IDs),:]
            all_points_subset['total_distance'] = 10000
            all_points_subset['feasible'] = True
            for index,row in all_points_subset.iterrows(): # cycle through the nodes and find the distance to all the others
                if check_cluster and row['pop_bool']==1:
                    check_id = row['ID']
                    cluster = int(clustered_points.loc[clustered_points['ID']==check_id,'Cluster'])
                total_weighted_distance = 0
                max_dist = 0
                if all_points_subset.loc[index,'feasible'] == True:
                    for index1,row1 in all_points_subset.loc[all_points_subset['pop_bool']==1,:].iterrows():

                        if index==index1:
                            total_distance = 0
                        else:
                            total_distance = T_metric[int(row['ID'])][int(row1['ID'])]['distance']
                            #print('distance between '+str(row['ID'] + ' and '+ str(row1['ID'] +' is '+str(total_distance))))
                        if total_distance >LV_grid_length*1.3:
                            all_points_subset.loc[index,'feasible'] = False
                            all_points_subset.loc[index1, 'feasible'] = False
                            continue
                        elif not total_distance==0:
                            total_weighted_distance += total_distance
                            if total_distance>max_dist:
                                max_dist = total_distance
                    all_points_subset.loc[index, 'av_distance'] = total_weighted_distance / len(all_points_subset)
                    all_points_subset.loc[index, 'max_distance'] =  max_dist
                    all_points_subset.loc[index,'final_distance'] = total_weighted_distance/len(all_points_subset)*0.9 + max_dist*0.1
            feasible_sites = all_points_subset.loc[all_points_subset['feasible']==True,:]
            best_site_ID = int(feasible_sites.loc[feasible_sites['final_distance']==feasible_sites['final_distance'].min(),'ID'].values[0])
            all_points.loc[all_points['ID']==best_site_ID,'substations']=True
            all_points.loc[all_points['ID']==best_site_ID,'Cluster']=cluster
            all_points.loc[all_points['ID']==best_site_ID,'LV_length'] = LV_grid_length
            all_points.loc[all_points['ID'] == best_site_ID, 'max_distance'] = float(feasible_sites.loc[feasible_sites['ID']==best_site_ID,'max_distance'])
        MV_LV_substations = all_points.loc[all_points['substations'] == True, :]

        grid_final = gpd.GeoDataFrame()
        path = list(tree_final.edges)
        counter = 0
        for i in path:
            point1 = all_points.loc[all_points['ID'] == i[0], 'geometry'].values[0]
            point2 = all_points.loc[all_points['ID'] == i[1], 'geometry'].values[0]
            length = T_metric[i[0]][i[1]]['distance']/1000
            cost = length*LV_base_cost
            geom = LineString([point1, point2])
            grid_final = grid_final.append(gpd.GeoDataFrame({'ID': [counter], 'geometry': [geom],'Length [km]':[length],'Cost [euro]':[cost]}))
            counter += 1
        grid_final.crs = crs
        grid_final.to_file(dir + '/grid_final_cut', index=False)
        LV_grid = LV_grid.append(grid_final)

        all_points[all_points['substations']==True].to_file(dir+'/secondary_substations',index=False)
        #start creating the final files and resumes




        # for the MV grid, create a copy of tree_final before cutting the branches for the LV network, and then find the steiner tree
        # considering only the secondary substations. When calculating the costs for the MV network, apply discounted rates for the lines in which
        # there is already a LV line.
        clusters_list = pd.DataFrame(columns=['Cluster', 'Sub_cluster', 'Population', 'Load [kW]'])
        for i in range(int(number_clusters)):
            subset = clustered_points[clustered_points['Cluster'] == i]
            if len(subset)==1: #specific case where its only 1 house in a cluster - we will add a substation here in that cluster, but in the future,
            # we might decide to simply assign a stand alone system
                LV_grid_length = 0
                LV_grid_cost = 0
                max_length = 0
                sum_pop = subset['Population'].sum()
                load = sum_pop*load_capita
                MV_LV_substations = MV_LV_substations.append(subset)
                MV_LV_substations.loc[MV_LV_substations['Cluster'] == i, 'LV_length'] = 0
                MV_LV_substations.loc[MV_LV_substations['Cluster'] == i, 'max_distance'] = 0
                MV_LV_substations.loc[MV_LV_substations['Cluster'] == i, 'MV_Power'] = load
            else:
                LV_grid_length = float(MV_LV_substations.loc[MV_LV_substations['Cluster'] == i ,'LV_length'])/1000
                LV_grid_cost = float(MV_LV_substations.loc[MV_LV_substations['Cluster'] == i, 'LV_length'])*LV_base_cost/1000 #fix
                max_length = float(MV_LV_substations.loc[MV_LV_substations['Cluster'] == i, 'max_distance'])/1000
                sum_pop = subset['Population'].sum()
                load = sum_pop * load_capita * coincidence_factor(sum_pop, pop_per_household)
            ID_substation = int(MV_LV_substations.loc[MV_LV_substations['Cluster'] == i ,'ID'])
            data = np.array([[int(clus), int(i), sum_pop, load,LV_grid_length, LV_grid_cost,max_length]])
            df2 = pd.DataFrame(data, columns=['Cluster', 'Sub_cluster', 'Population', 'Load [kW]','Grid_Length [km]', 'Grid Cost [euro]','Max length [km]'])
            clusters_list = clusters_list.append(df2)
            MV_LV_substations.loc[MV_LV_substations['Cluster']==i,'MV_Power' ] = load
            MV_LV_substations.loc[MV_LV_substations['Cluster'] == i, 'Population'] = sum_pop
            MV_LV_substations.to_file(dir+'/secondary_substations')
        MV_LV_substations['Cluster2']=MV_LV_substations['Cluster']
        MV_LV_substations['Cluster'] = clus
        secondary_substations = secondary_substations.append(MV_LV_substations)
        substation_data = pd.read_csv(gisele_dir + '/general_input/' + ss_data)
        clusters_list = categorize_substation(clusters_list, substation_data)
        clusters_list['Population'] = [ceil(i) for i in clusters_list['Population']]
        clusters_list.to_csv(dir + '/LV_networks_resume.csv', index=False)
        LV_resume = LV_resume.append(clusters_list)
        #total_costs = sum(clusters_list['Cost[euro]'].to_list())
        all_houses=all_houses.append(clustered_points)
        LV_grid = LV_grid.append(grid_final)
        terminal_MV_nodes = MV_LV_substations['ID'].to_list()
        ### ROUTE THE MV NETWORK BY ASSIGNING A LOWER COST TO THE EXISTING LV NETWORK
        if len(terminal_MV_nodes)>1:
            for i in tree_final_copy.edges:
                if tree_final.has_edge(*i):
                    # HERE PUT THIS AS AN ADDITIONAL DISCOUNT, BUT THE GRID HAS TO GO ACCROSS THE ROADS AND NOT PASS THROUGH SOMEONE's HOUSE
                    tree_final_copy[i[0]][i[1]]['weight'] *= 0.5

            tree_MV = steiner_tree(tree_final_copy, terminal_MV_nodes)
            grid_MV = gpd.GeoDataFrame()
            path = list(tree_MV.edges)
            counter = 0
            for i in path:
                point1 = all_points.loc[all_points['ID'] == i[0], 'geometry'].values[0]
                point2 = all_points.loc[all_points['ID'] == i[1], 'geometry'].values[0]
                id1 = int(all_points.loc[all_points['ID'] == i[0], 'ID'])
                id2 = int(all_points.loc[all_points['ID'] == i[1], 'ID'])
                length = T_metric[i[0]][i[1]]['distance']/1000
                cost = length*LV_base_cost # change this
                geom = LineString([point1, point2])
                grid_MV = grid_MV.append(gpd.GeoDataFrame({'ID1': [id1],'ID2':[id2],'geometry': [geom],'Length':[length],'Cost': [cost]}))
                counter += 1
            grid_MV.crs = crs
            grid_MV['Cluster'] = clus
            grid_MV.to_file(dir + '/grid_MV', index=False)
            #grid_of_points_GDF = add_roads_points_to_gdf(grid_of_points_GDF,all_points,grid_MV,clus)
            grid_MV.to_file(dir+'/Grid_substations.shp')
            MV_grid=MV_grid.append(grid_MV)

    LV_resume.to_csv(gisele_dir + '/' + dir_output + '/LV_resume.csv')
    LV_grid.to_file(gisele_dir + '/' + dir_output + '/LV_grid')
    secondary_substations.to_file(gisele_dir + '/' + dir_output + '/secondary_substations')
    all_houses.to_file(gisele_dir + '/' + dir_output + '/final_users')
    if not MV_grid.empty:
        MV_grid.to_file(gisele_dir + '/' + dir_output + '/MV_grid',index=False)

    #number_total_points = len(grid_of_points_GDF)
    #umber_substations = len(secondary_substations)
    #substations_new_IDs = [*range(number_total_points,number_total_points+number_substations)]
    #secondary_substations['ID']=substations_new_IDs
    secondary_substations['Substation'] = 1
    secondary_substations['Weight'] = 3 # for the logic with the routing, if it's 1 - the algorithm jumps the road and just
    #connects further point with the substation
    secondary_substations['Type'] = 'Secondary Substation'
    terminal_MV_nodes = secondary_substations.ID.to_list()
    grid_of_points_GDF.drop(grid_of_points_GDF[grid_of_points_GDF['ID'].isin(terminal_MV_nodes)].index, axis=0, inplace=True)
    grid_of_points_GDF=grid_of_points_GDF.append(secondary_substations)
    grid_of_points_GDF[['X','Y','ID','Population','Elevation','Weight','geometry','Land_cover','Cluster','MV_Power','Substation','Type']].\
        to_csv(gisele_dir + '/' + dir_input + '/weighted_grid_of_points_with_ss_and_roads.csv',index=False)
    #add the roads as well. Previously it was done in routing.py in row 212

    