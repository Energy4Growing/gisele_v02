from shapely.geometry import MultiLineString
from shapely.ops import nearest_points
from gisele.functions import *
from math import *
from scipy import sparse
import networkx as nx
from gisele.Steiner_tree_code import *
from itertools import combinations
import numpy as np
import time
from sklearn.cluster import AgglomerativeClustering
from collections import Counter
from gisele.geneticalgorithm_github import geneticalgorithm as ga
from scipy.spatial import Delaunay
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

    dir = '/grid_to_house'
    dir = '/houses_translated_on_roads'

    clustered_points.to_file('Testing_strategy' + dir + '/Clustered_points_after_genetic', index=False)
    return clustered_points
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
    for i,j in final_sides2:
        point1 = new_points.loc[new_points['order']==i,'geometry'].values[0]
        point2 = new_points.loc[new_points['order'] == j, 'geometry'].values[0]
        id1 = int(new_points.loc[new_points['order'] == i, 'ID'])
        id2 = int(new_points.loc[new_points['order'] == j, 'ID'])
        length = point1.distance(point2)
        line = LineString([point1, point2])
        if length<320 and not graph.has_edge(id1,id2) and ((sum([line.intersects(line1) for line1 in new_lines_old.geometry]) == 0) or
                ((new_points.loc[new_points['ID'] == id1, 'pop_bool'] == 0).values[0]) or
                (new_points.loc[new_points['ID'] == id2, 'pop_bool'] == 0).values[0]):
            graph.add_edge(id1,id2 , weight=length)

            data_segment = {'ID1': [id1], 'ID2': [id2], 'length': [point1.distance(point2) / 1000],
                            'geometry': [line], 'Type': ['Colateral']}
            new_lines = new_lines.append(gpd.GeoDataFrame(data_segment))

    return graph,new_lines

def create_clean_graph(graph,points,terminal_points,T_metric):
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
                print('Node ' + str(current_node) + ' was deleted.')
            print('Next point is '+str(unique_nodes[0]))
            start_node = unique_nodes[0]
            current_node = start_node
            next_node = [i for i in graph_copy[start_node]][0]
        if current_node == 155 and next_node == 156:
            print('tuka')
        if next_node in new_nodes:
            new_graph.add_edge(start_node, next_node, weight=T_metric[start_node][next_node])
            print('add ' + str(start_node) + ' and ' + str(next_node))
            graph_copy.remove_edge(current_node,next_node)
            print('remove '+str(current_node)+' and ' + str(next_node))
            if start_node in terminal_IDs_2:
                terminal_IDs_2.remove(start_node)
                print('Node '+ str(start_node)+' was deleted.')
            start_node = next_node
            current_node = start_node
        else:
            graph_copy.remove_edge(current_node, next_node)
            print('remove ' + str(current_node) + ' and ' + str(next_node))
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
    new_lines.crs = 22287
    dir = '/grid_to_house'
    new_lines.to_file('Testing_strategy' +dir + '/new_graph', index=False)
    new_points = points[points['ID'].isin(new_nodes)]
    return new_lines, new_points, new_graph
def connect_unconnected_graph(graph,lines,points):
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
                graph.add_edge(id_point1,id_point2,weight = distance*3)
    return graph,lines
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
    that represents the roads with lines and poinnts.'''
    id1=74
    id2=75
    id3=4
    id4=5
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
def distance_point_to_roads(segments,line_gdf,critdist):
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
    critical_points.crs = 22287
    return critical_points

    return critical_points

def houses_translated_on_roads():
    dir = '/houses_translated_on_roads'
    roads_weight = 1
    Roads = gpd.read_file(r'C:\Users\alekd\PycharmProjects\Gisele\Case studies\Thuso-10_clusters-180m_attempt2\Output\Clusters\10\Roads.shp')
    #perhaps add a function that separates the roads and puts a weight based on the type ( if it's a track - higher weight)
    Roads=Roads[Roads['highway']!='path']
    Roads = MultiLine_to_Line(Roads)
    roads=Roads.simplify(5)
    gdf_roads= gpd.GeoDataFrame(geometry=roads)
    os.chdir('..')
    gdf_roads.to_file('Testing_strategy'+dir+'/roads_simplified')
    w=0

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
                print('Double road point!')

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
    line_gdf = gpd.GeoDataFrame(line_vertices, crs='EPSG:22287',
                                    geometry=geometry)
    line_vertices.loc[:, 'Elevation']=1000
    segments.crs=22287
    line_gdf.to_file('Testing_strategy'+dir+'/Points_first')
    segments.to_file('Testing_strategy'+dir+'/Lines_first')
    critical_points = distance_point_to_roads(segments, line_gdf, critdist=7)
    segments, line_gdf = fix_roads(segments, line_gdf, critical_points,critdist = 7)
    line_gdf.to_file('Testing_strategy'+dir+'/Points_fixed1')
    segments.to_file('Testing_strategy'+dir+'/lines_fixed1')
    segments,line_gdf= fix_lines_intersecting(segments,line_gdf)
    line_gdf.to_file('Testing_strategy'+dir+'/points_fixed2')
    segments.to_file('Testing_strategy'+dir+'/lines_fixed2')
    segments = segments.reset_index(drop=True)
    new_lines = gpd.GeoDataFrame()

    next_node = line_gdf.shape[0]-1 # just logic stuff
    new_points = line_gdf.copy()
    #here make sure that when we find the points in the middle of a longer line, translate them on to the original Roads file ( non-simplified)
    # and use that point - to get a more realistic representation of the roads via the grid of points
    for index, segment in segments.iterrows():
        if segment['length']>0.045:
            print('da')
            number_segments = ceil(segment['length']/0.045)
            num_points = number_segments-1
            new_lengths = segment['length']/number_segments
            splitter = MultiPoint([segment['geometry'].interpolate((i / number_segments), normalized=True) for i in range(1, number_segments)])
            for i in range(num_points):
                if i==0:
                    id1=int(segment['ID1'])
                else:
                    id1=next_node
                next_node += 1
                data_point = {'ID':[next_node],'X':[splitter[i].xy[0][0]],'Y':[splitter[i].xy[1][0]],'Weight': [1], 'Elevation':[1000],'geometry':[splitter[i]]}
                new_points=new_points.append(gpd.GeoDataFrame(data_point))
                line = LineString([new_points.loc[new_points['ID'] == id1, 'geometry'].values[0], splitter[i]])
                data_segment={'ID1':[id1],'ID2':[next_node],'length':[new_lengths],'geometry':[line]}
                new_lines = new_lines.append(gpd.GeoDataFrame(data_segment))
            line = LineString([new_points.loc[new_points['ID'] == next_node, 'geometry'].values[0], new_points.loc[new_points['ID'] == int(segment['ID2']), 'geometry'].values[0]])
            data_segment = {'ID1': [next_node], 'ID2': [int(segment['ID2'])], 'length': [new_lengths], 'geometry': [line]}
            new_lines = new_lines.append(gpd.GeoDataFrame(data_segment))
        else:
            new_lines = new_lines.append(segments.loc[index,:])

    new_lines.crs=22287
    new_lines.to_file('Testing_strategy'+dir+'/new_lines.shp',index=False)
    print('Kraj')

    Population = gpd.read_file(r'C:\Users\alekd\PycharmProjects\Gisele\Case studies\Thuso-10_clusters-180m_attempt2\Output\Clusters\10\points.shp')
    Population=Population[Population['Population']>0]
    number_points = new_points.shape[0]
    roads_points = MultiPoint([point for point in new_points['geometry']])
    Population['ID'] = [*range(number_points,number_points+Population.shape[0])]
    Population.to_file('Testing_strategy'+dir+'/population.shp',index=False)
    new_points['Population']=0
    new_points['pop_bool']=0
    for i,pop in Population.iterrows():
        point = pop['geometry']
        nearest_geoms = nearest_points(point, roads_points)
        closest_road_point = nearest_geoms[1]
        new_points.loc[new_points['geometry']==closest_road_point,'Population']= pop['Population']
        new_points.loc[new_points['geometry'] == closest_road_point, 'pop_bool'] = 1
    new_points.crs=22287
    new_points.to_file('Testing_strategy'+dir+'/new_points',index=False)

    new_points = new_points.set_index('ID',drop=False)
    new_lines = new_lines.set_index(pd.Index([*range(new_lines.shape[0])]))
    start = time.time()
    strategy = 'Old without manipulating graph, just based on distance between points'
    strategy = 'NEW'
    if strategy == 'Old without manipulating graph, just based on distance between points':
        dist_2d_matrix = distance_2d(new_points, new_points, 'X', 'Y')
        dist_3d_matrix = distance_3d(new_points, new_points, 'X', 'Y', 'Elevation')
        line_cost = 1000
        edges_matrix=dist_3d_matrix.copy()
        edges_matrix[dist_2d_matrix > ceil(45*1.1)] = 0
        edges_matrix_sparse = sparse.csr_matrix(edges_matrix)
        graph = nx.from_scipy_sparse_matrix(edges_matrix_sparse)
        graph,new_lines = connect_unconnected_graph(graph,new_lines,new_points)
        new_lines.to_file('Testing_strategy'+dir+'/new_lines.shp',index=False)
        for index,row in new_lines.iterrows():
            id1= row['ID1']
            id2=row['ID2']
            print(id1,id2)
            graph[id1][id2]['weight'] = row['length']*1000*roads_weight
            print(row['length']*1000*roads_weight)
    elif strategy == 'NEW':
        graph=nx.Graph()
        for index,row in new_lines.iterrows():
            id1= row['ID1']
            id2=row['ID2']
            print(id1,id2)
            graph.add_edge(id1, id2, weight=row['length']*1000*roads_weight)
            print(row['length']*1000*roads_weight)
        graph, new_lines = connect_unconnected_graph(graph, new_lines, new_points)
        new_lines.to_file('Testing_strategy'+dir+'/new_lines', index=False)
    populated_points = new_points[new_points['pop_bool']==1]
    terminal_nodes = list(populated_points['ID'])
    new_points[new_points['ID'].isin(terminal_nodes)].to_file('Testing_strategy'+dir+'/terminal_points', index=False)
    tree = steiner_tree(graph, terminal_nodes)

    path = list(tree.edges)
    grid_routing = gpd.GeoDataFrame()
    counter = 0
    for i in path:
        point1 = new_points.loc[new_points['ID'] == i[0], 'geometry'].values[0]
        point2 = new_points.loc[new_points['ID'] == i[1], 'geometry'].values[0]
        geom = LineString([point1, point2])
        grid_routing = grid_routing.append(gpd.GeoDataFrame({'ID': [counter], 'geometry': [geom]}))
        counter += 1
    grid_routing.crs = 22287
    grid_routing.to_file('Testing_strategy' + dir + '/grid_routing', index=False)
    #new addition from 17 June
    new_points=new_points[new_points['ID'].isin(list(tree.nodes))]
    new_points['Population'] = 0
    new_points['pop_bool'] = 0
    index = [*range(new_points['ID'].astype(int).max()+1,new_points['ID'].astype(int).max()+1+Population.shape[0])]
    Population['ind'] = index
    Population.set_index('ind',inplace=True,drop=True)
    Population['Population'] = 4
    Population['pop_bool'] = 1
    new_points = new_points.append(Population)
    new_graph = tree.copy()
    for n in new_graph.edges:
         new_graph[n[0]][n[1]]['weight'] =new_graph[n[0]][n[1]]['weight']*0.05


    new_lines_old = new_lines.copy()
    new_lines['Type'] = 'Road'
    new_points['order'] = [*range(new_points.shape[0])]

    tactic='NEW'
    timestart = time.time()
    if tactic=='OLD':
        dist_2d_matrix = distance_2d(new_points, new_points, 'X', 'Y')
        dist_2d_matrix = pd.DataFrame(dist_2d_matrix)
        dist_2d_matrix.columns = new_points['ID'].to_list()
        dist_2d_matrix.index = new_points['ID'].to_list()
        for i, house in Population.iterrows():
            LIST = dist_2d_matrix[i].values
            idx = np.argpartition(LIST,
                                  5)  # find the 4 closest points and connect to them, 1 will be the connection to self.
            for j in range(5):
                length = LIST[idx[j]]
                Point1 = new_points.loc[new_points['ID'] == i, 'geometry'].values[0]
                Point2 = new_points.loc[new_points['order'] == idx[j], 'geometry'].values[0]
                line = LineString([Point1, Point2])
                pom = new_points.loc[new_points['order'] == idx[j], 'ID'].values[0]
                # now this part also makes sure that the houses are not connected jumping over a road
                try:
                    if i != pom and ((sum([line.intersects(line1) for line1 in new_lines_old.geometry]) == 0) or
                                        (new_points.loc[new_points['ID'] == pom, 'pop_bool'] == 0).values[0]):
                        new_graph.add_edge(i, pom, weight=length)  # *1.5
                        line = LineString([Point1, Point2])
                        data_segment = {'ID1': [i], 'ID2': [pom], 'length': [Point1.distance(Point2) / 1000],
                                        'geometry': [line], 'Type': ['Colateral']}
                        new_lines = new_lines.append(gpd.GeoDataFrame(data_segment))
                except:
                    print(pom)
    else:
        print('aowasda')
        new_graph, new_lines = delaunay_test(new_graph, new_points, new_lines)
    new_graph, new_lines = connect_unconnected_graph(new_graph, new_lines, new_points)
    new_lines.to_file('Testing_strategy' + dir + '/new_lines_final', index=False)
    terminal_nodes = new_points.loc[new_points['pop_bool']==1,'ID'].to_list()
    tree_final = steiner_tree(new_graph, terminal_nodes)
    grid_final = gpd.GeoDataFrame()
    path = list(tree_final.edges)
    counter = 0
    for i in path:
        point1 = new_points.loc[new_points['ID'] == i[0], 'geometry'].values[0]
        point2 = new_points.loc[new_points['ID'] == i[1], 'geometry'].values[0]
        geom = LineString([point1, point2])
        grid_final = grid_final.append(gpd.GeoDataFrame({'ID': [counter], 'geometry': [geom]}))
        counter += 1
    grid_final.crs = 22287
    grid_final.to_file('Testing_strategy' + dir + '/grid_final', index=False)
    for n in tree.edges:
        try:
            tree_final[n[0]][n[1]]['weight'] = tree[n[0]][n[1]]['weight']
        except:
            print('Edge is not in the final tree')
    # end of new addition

    T_metric = metric_closure(tree_final ,weight='weight')
    points_set = new_points.loc[new_points['pop_bool']==1,'ID'].values
    dist_matrix = np.zeros((len(points_set),len(points_set)))
    for i in range(len(points_set)):
        for j in range(len(points_set)):
            if not i==j:
                dist_matrix[i,j] = T_metric[points_set[i]][points_set[j]]['distance']


    clustering = AgglomerativeClustering(n_clusters=None, affinity = 'precomputed', linkage='complete',
                                                 distance_threshold=1000).fit(dist_matrix)
    populated_points = new_points[new_points['pop_bool']==1]
    populated_points['Cluster'] = clustering.labels_
    populated_points.to_file('Testing_strategy'+dir+'/Clustered_points',index=False)


    populated_points['Population'] = 4
    number_clusters = populated_points['Cluster'].max() +1
    lines_new_graph, points_new_graph, new_graph = create_clean_graph(tree_final, new_points, populated_points,
                                                                      T_metric)  # finish this function
    points_set = points_new_graph['ID'].values
    dist_matrix2 = np.zeros((len(points_set), len(points_set)))
    for i in range(len(points_set)):
        for j in range(len(points_set)):
            if not i == j:
                dist_matrix2[i, j] = T_metric[points_set[i]][points_set[j]]['distance']
    points_new_graph.loc[points_new_graph['Population'] > 0, 'Population'] = 4
    genetic2(populated_points, points_new_graph, dist_matrix2, number_clusters, new_graph)

    end=time.time()
    print(r'The time required for data processing + steiner is:'+str(end-start))



def houses_separate_points2():
    dir = '/grid_to_house'
    roads_weight = 1
    Roads = gpd.read_file(
        r'C:\Users\alekd\PycharmProjects\Gisele\Case studies\Thuso-10_clusters-180m_attempt2\Output\Clusters\8\Roads.shp')
    # perhaps add a function that separates the roads and puts a weight based on the type ( if it's a track - higher weight)
    Roads = Roads[Roads['highway'] != 'path']
    Roads = MultiLine_to_Line(Roads)
    roads = Roads.simplify(5)
    gdf_roads = gpd.GeoDataFrame(geometry=roads)
    os.chdir('..')
    gdf_roads.to_file('Testing_strategy' +dir + '/roads_simplified')
    w = 0

    line_vertices = pd.DataFrame(
        index=pd.Series(range(w, w + len(gdf_roads.index))),
        columns=['ID', 'X', 'Y', 'ID_line', 'Weight', 'Elevation'], dtype=int)
    # create geodataframe with all the segments that compose the road
    segments = gpd.GeoDataFrame(columns=['geometry', 'ID1', 'ID2'])
    k = 0

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
                print('Double road point!')

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
    line_gdf = gpd.GeoDataFrame(line_vertices, crs='EPSG:22287',
                                geometry=geometry)
    line_vertices.loc[:, 'Elevation'] = 1000
    segments.crs = 22287
    line_gdf.to_file('Testing_strategy' +dir + '/Points_first')
    segments.to_file('Testing_strategy' +dir + '/Lines_first')
    critical_points = distance_point_to_roads(segments, line_gdf, critdist=7)
    segments, line_gdf = fix_roads(segments, line_gdf, critical_points,critdist=7)
    line_gdf.to_file('Testing_strategy' +dir + '/Points_fixed1')
    segments.to_file('Testing_strategy' +dir + '/lines_fixexd1')
    segments, line_gdf = fix_lines_intersecting(segments, line_gdf)
    line_gdf.to_file('Testing_strategy' +dir + '/points_fixed2')
    segments.to_file('Testing_strategy' +dir + '/lines_fixed2')
    segments = segments.reset_index(drop=True)
    new_lines = gpd.GeoDataFrame()

    next_node = line_gdf.shape[0] - 1  # just logic stuff
    new_points = line_gdf.copy()
    # here make sure that when we find the points in the middle of a longer line, translate them on to the original Roads file ( non-simplified)
    # and use that point - to get a more realistic representation of the roads via the grid of points
    for index, segment in segments.iterrows():
        if segment['length'] > 0.045:
            print('da')
            number_segments = ceil(segment['length'] / 0.045)
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

    new_lines.crs = 22287
    new_lines.to_file('Testing_strategy' +dir + '/new_lines.shp', index=False)
    print('Kraj')

    Population = gpd.read_file(
        r'C:\Users\alekd\PycharmProjects\Gisele\Case studies\Thuso-10_clusters-180m_attempt2\Output\Clusters\8\points.shp')
    Population = Population[Population['Population'] > 0]
    number_points = new_points.shape[0]
    roads_points = MultiPoint([point for point in new_points['geometry']])
    Population['ID'] = [*range(number_points, number_points + Population.shape[0])]
    Population.to_file('Testing_strategy' +dir + '/population.shp', index=False)
    new_points['Population'] = 0
    new_points['pop_bool'] = 0
    Population = Population.set_index('ID', drop=False)
    Population['pop_bool'] = 1
    new_points = new_points.append(Population)
    new_points.crs = 22287
    new_points.to_file('Testing_strategy' +dir + '/new_points', index=False)
    new_points = new_points.set_index('ID', drop=False)
    new_lines = new_lines.set_index(pd.Index([*range(new_lines.shape[0])]))
    start = time.time()
    strategy = 'Old without manipulating graph, just based on distance between points'
    strategy = 'NEW'
    if strategy == 'Old without manipulating graph, just based on distance between points':
        dist_2d_matrix = distance_2d(new_points, new_points, 'X', 'Y')
        dist_3d_matrix = distance_3d(new_points, new_points, 'X', 'Y', 'Elevation')
        line_cost = 1000
        edges_matrix = dist_3d_matrix.copy()
        edges_matrix[dist_2d_matrix > ceil(45 * 1.1)] = 0
        edges_matrix_sparse = sparse.csr_matrix(edges_matrix)
        graph = nx.from_scipy_sparse_matrix(edges_matrix_sparse)
        graph, new_lines = connect_unconnected_graph(graph, new_lines, new_points)
        new_lines.to_file('Testing_strategy' +dir + 'new_lines.shp', index=False)
        for index, row in new_lines.iterrows():
            id1 = row['ID1']
            id2 = row['ID2']
            print(id1, id2)
            graph[id1][id2]['weight'] = row['length'] * 1000 * roads_weight
            print(row['length'] * 1000 * roads_weight)
    elif strategy == 'NEW':
        graph = nx.Graph()
        for index, row in new_lines.iterrows():
            id1 = row['ID1']
            id2 = row['ID2']
            print(id1, id2)
            graph.add_edge(id1, id2, weight=row['length'] * 1000 * roads_weight)
            print(row['length'] * 1000 * roads_weight)

        dist_2d_matrix = distance_2d(new_points, new_points, 'X', 'Y')
        dist_2d_matrix = pd.DataFrame(dist_2d_matrix)
        new_lines_old = new_lines.copy()
        new_lines['Type'] = 'Road'
        for i, house in Population.iterrows():
            LIST = dist_2d_matrix[i].values
            idx = np.argpartition(LIST,5)  # find the 4 closest points and connect to them, 1 will be the connection to self.
            for j in range(5):
                length = LIST[idx[j]]
                Point1 = new_points.loc[new_points['ID'] == i, 'geometry'].values[0]
                Point2 = new_points.loc[new_points['ID'] == idx[j], 'geometry'].values[0]
                line = LineString([Point1, Point2])
                # now this part also makes sure that the houses are not connected jumping over a road
                if i != idx[j] and ((sum([line.intersects(line1) for line1 in new_lines_old.geometry])==0) or
                                    (new_points.loc[new_points['ID']==idx[j],'pop_bool']==0).values[0]):
                    graph.add_edge(i, idx[j], weight=length) #*1.5
                    line = LineString([Point1,Point2])
                    data_segment = {'ID1': [i], 'ID2': [idx[j]], 'length': [Point1.distance(Point2)/1000], 'geometry': [line],'Type':['Colateral']}
                    new_lines = new_lines.append(gpd.GeoDataFrame(data_segment))
        graph, new_lines = connect_unconnected_graph(graph, new_lines, new_points)
        new_lines.to_file('Testing_strategy' +dir + '/new_lines_final', index=False)
    populated_points = new_points[new_points['pop_bool'] == 1]
    terminal_nodes = populated_points['ID'].to_list()
    tree = steiner_tree(graph, terminal_nodes)

    T_metric = metric_closure(tree, weight='weight')
    points_set = populated_points['ID'].values
    dist_matrix = np.zeros((len(points_set), len(points_set)))
    for i in range(len(points_set)):
        for j in range(len(points_set)):
            if not i == j:
                dist_matrix[i, j] = T_metric[points_set[i]][points_set[j]]['distance']
    clustering = AgglomerativeClustering(n_clusters=None, affinity='precomputed', linkage='complete',
                                         distance_threshold=1000).fit(dist_matrix)
    populated_points['Cluster'] = clustering.labels_
    populated_points['Population'] = 4
    number_clusters=populated_points['Cluster'].max()+1
    lines_new_graph, points_new_graph, new_graph = create_clean_graph(tree, new_points, populated_points,T_metric) # finish this function
    points_set = points_new_graph['ID'].values
    dist_matrix2 = np.zeros((len(points_set), len(points_set)))
    for i in range(len(points_set)):
        for j in range(len(points_set)):
            if not i == j:
                dist_matrix2[i, j] = T_metric[points_set[i]][points_set[j]]['distance']
    points_new_graph.loc[points_new_graph['Population']>0,'Population']  = 4
    genetic2(populated_points,points_new_graph,dist_matrix2, number_clusters, new_graph)

    Clusters=pd.DataFrame()
    for i in range(int(number_clusters)):
        subset = populated_points[populated_points['Cluster'] == i]
        Clusters.loc[i,'Population'] = subset['Population'].sum()
        Clusters.loc[i,'Load'] = subset['Population'].sum()*0.7*0.3
    populated_points.to_file('Testing_strategy' +dir + '/Clustered_points', index=False)

    end = time.time()
    print(r'The time required for data processing + steiner is:' + str(end - start))
    #path = list(tree.edges)
    grid_routing = gpd.GeoDataFrame()
    counter = 0
    for i in tree.edges:
        point1 = new_points.loc[new_points['ID'] == i[0], 'geometry'].values[0]
        point2 = new_points.loc[new_points['ID'] == i[1], 'geometry'].values[0]
        geom = LineString([point1, point2])
        grid_routing = grid_routing.append(gpd.GeoDataFrame({'ID': [counter], 'geometry': [geom]}))
        counter += 1
    grid_routing.crs = 22287
    grid_routing.to_file('Testing_strategy' +dir + '/grid_routing', index=False)

#houses_separate_points2()
houses_translated_on_roads()


