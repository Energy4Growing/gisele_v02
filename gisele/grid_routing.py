from gisele.functions import *
from shapely.geometry import Point,box
from gisele import Steinerman, Spiderman, dijkstra,initialization

def routing(Clusters,gisele_dir,case_study,crs,resolution_MV,Roads_option,roads_weight,simplify_road_coef_inside,Rivers_option,line_bc):
    dir_output = '/Case studies/' + case_study
    dir_input = r'Case studies/' + case_study + '/Input'
    grid_of_points_MV = pd.read_csv(dir_input + '/weighted_grid_of_points_with_ss.csv')
    grid_of_points_MV = gpd.GeoDataFrame(grid_of_points_MV,
                                     geometry=gpd.points_from_xy(grid_of_points_MV.X, grid_of_points_MV.Y), crs=crs)
    clus=1
    entire_MV_grid = gpd.GeoDataFrame()
    grid_of_points_MV2 = grid_of_points_MV.copy()
    for index, row in Clusters.iterrows():
        print('Routing cluster number ' + str(row['cluster_ID']))
        dir = gisele_dir + dir_output + '/Output/Clusters/' + str(row['cluster_ID'])
        os.chdir(gisele_dir)
        clusters_list = pd.read_csv(dir + '/clusters_list.csv')
        df = pd.read_csv(dir + '/Input.csv')
        df_weighted = initialization.weighting(df, 30, 'ESACCI')
        grid_of_points=gpd.read_file(dir + '/points.shp')
        grid_of_points['Weight'] = df_weighted['Weight']
        grid_of_points.to_file(dir + '/points.shp')
        area=row['geometry']
        # THIS IS SPECIFIC FOR THE CASE OF THUSO - ISSUE NEEDS TO BE FIXED WITH CLUSTERS THAT ARE TOO NARROW. If in one line
        # of the grid points there are no points -> the graph for the steiner tree will not be connected
        area_buffered=area
        # area_buffered = row['geometry'].buffer((resolution_MV * 0.1 / 11250) / 2)



        # STUPID USELESS PARAMETERS
        input_csv = 'imported_csv'
        input_sub = 'imported_subs'
        pop_thresh = 1

        sub_cost_hv = 100000
        sub_cost_mv = 10000


        study_area = area_buffered
        area = study_area

        area.crs = crs

        grid_500m_clip = gpd.clip(grid_of_points_MV2, area)  # this is our 500m resolution grid of points
        grid_500m_clip.to_file(dir + '/points500m.shp')

        try:
            c_grid, grid_resume, grid_of_points_MV, grid_500m_clip = routing_secondary_substations(gisele_dir, dir,
                 grid_500m_clip, grid_of_points_MV,clusters_list,resolution_MV,pop_thresh,input_sub, line_bc,
                        sub_cost_hv, sub_cost_mv,  Roads_option,roads_weight, clus,simplify_road_coef_inside, Rivers_option)


            c_grid.to_file(dir + '/grids/Grid_substations.shp')
            grid_resume.to_csv(dir + '/grids/grid_resume.csv')
            grid_500m_clip = create_box(grid_500m_clip, grid_of_points_MV, resolution_MV)
            grid_MV = grid_500m_clip.to_crs(crs)
            grid_MV.to_file(dir + '/points500m.shp')

            entire_MV_grid = entire_MV_grid.append(c_grid)
        except:
            print('Cant find the MV routing inside this cluster because its just 1 substation')
        clus+=1
    if Roads_option:
        grid_of_points_MV.to_csv(gisele_dir + '/' + dir_input + '/weighted_grid_of_points_with_ss_and_roads.csv')
    entire_MV_grid = entire_MV_grid.to_crs(crs)
    entire_MV_grid.to_file(gisele_dir + '/' + dir_input + '/clusters_MV_grid')

def LV_routing(gisele_dir,case_study,crs,resolution,Roads_option,simplify_road_coef_inside,line_LV_cost,roads_weight):
    dir_output = '/Case studies/' + case_study+'/Input'
    dir_input = r'/Case studies/' + case_study + '/Output/Clusters'
    cluster_folders = os.listdir(gisele_dir+dir_input)

    LV_resume = pd.read_csv(gisele_dir+dir_output+'/LV_resume.csv')
    Total_LV_grid = gpd.GeoDataFrame()

    for cluster in cluster_folders:
        directory = gisele_dir+dir_input+'/' + cluster +'/'
        points=gpd.read_file(directory + 'points.shp')
        # STUPID USELESS PARAMETERS
        pop_thresh = 1

        number_sub_clusters = int(points['Cluster'].max()+1)
        print('Cluster '+ cluster +' has a total of ' + str(number_sub_clusters)+ ' secondary substations.')
        for sub_cluster in range(number_sub_clusters):

            Total_LV_grid,LV_resume = routing_LV_network(directory, points,
                    pop_thresh, line_LV_cost,Roads_option,simplify_road_coef_inside,resolution,Total_LV_grid,LV_resume,cluster
                                                ,sub_cluster,crs,roads_weight)

    Total_LV_grid.crs = crs
    Total_LV_grid.to_file(gisele_dir+dir_output+'/LV_Network',index=False)
    LV_resume.to_csv(gisele_dir+dir_output+'/LV_resume.csv',index=False)
def routing_LV_network(dir,geo_df,pop_thresh, line_bc,Roads_option,simplify_road_coef_inside, resolution,Total_LV_grid, LV_resume,
                       cluster, sub_cluster,crs,roads_weight):
    ''' If we are considering the Roads, we use the houses + the roads points, otherwise, we use the grid of points'''
    gdf_cluster_pop = geo_df[geo_df['Cluster']==sub_cluster]
    n_terminal_nodes = int(len(gdf_cluster_pop))

    x_min = min(gdf_cluster_pop.X)
    x_max = max(gdf_cluster_pop.X)
    y_min = min(gdf_cluster_pop.Y)
    y_max = max(gdf_cluster_pop.Y)

    dist = Point(x_min, y_min).distance(Point(x_max, y_max))
    if dist < resolution*10:
        extension = resolution * 2
    elif dist < resolution*30:
        extension = resolution*3
    else:
        extension = resolution*5
    area = box(minx=x_min - extension, maxx=x_max + extension,
             miny=y_min - extension, maxy=y_max + extension)
    geo_df=gpd.clip(geo_df,area)
    if Roads_option:
        x_min = min(gdf_cluster_pop.X)
        x_max = max(gdf_cluster_pop.X)
        y_min = min(gdf_cluster_pop.Y)
        y_max = max(gdf_cluster_pop.Y)
        area_cluster = box(minx=x_min , maxx=x_max,
                   miny=y_min , maxy=y_max )
        roads = gpd.read_file(dir+'/Roads.shp')
        roads = gpd.clip(roads,area_cluster)
        #roads = roads.geometry.simplify(tolerance=0.0005)
        print('create points along roads..')
        roads = MultiLine_to_Line(roads)
        roads= roads.geometry.simplify(tolerance=simplify_road_coef_inside)
        roads = gpd.GeoDataFrame(geometry=roads)
        gdf_roads, roads_segments = create_roads(roads, geo_df)
        roads_segments.crs = crs
        gdf_roads.crs = crs
    else:
        gdf_roads=gpd.GeoDataFrame()
        roads_segments=gpd.GeoDataFrame()
    if not os.path.exists(dir+'/grids'):
        os.makedirs(dir+'/grids')


    if n_terminal_nodes > 1:
        l()
        print("Creating LV grid for cluster " + str(cluster)+' - sub_cluster' + str (sub_cluster))

        c_grid, c_grid_cost, c_grid_length, c_grid_points = \
            cluster_grid(gdf_cluster_pop, gdf_cluster_pop, resolution,
                         line_bc, n_terminal_nodes, gdf_cluster_pop, gdf_roads,
                         roads_segments, Roads_option, roads_weight,Rivers_option=False,length_max=10)
        print("Grid created")

        if Total_LV_grid.empty:
            Total_LV_grid = c_grid
        else:
            Total_LV_grid = Total_LV_grid.append(c_grid)
        c_grid.to_file(dir+'/grids/'+str(sub_cluster))

        LV_resume.loc[(LV_resume['Cluster']==int(cluster)) & (LV_resume['Sub_cluster']==sub_cluster),'Grid_Length [km]'] = c_grid_length/1000
        LV_resume.loc[(LV_resume['Cluster'] == int(cluster)) & (LV_resume['Sub_cluster'] == sub_cluster), 'Grid Cost [euro]'] = c_grid_cost

        return Total_LV_grid,LV_resume

def routing_secondary_substations(input_dir_dir,dir,geo_df_clustered, geo_df, clusters_list, resolution,
            pop_thresh, input_sub, line_bc, sub_cost_hv, sub_cost_mv,Roads_option,roads_weight,clus,simplify_road_coef_inside, Rivers_option):
    gdf_cluster = geo_df_clustered
    gdf_cluster_pop = gdf_cluster[gdf_cluster['Substation'] == 1]

    n_terminal_nodes = int(len(gdf_cluster_pop))
    if Roads_option:
        x_min = min(geo_df_clustered.X)
        x_max = max(geo_df_clustered.X)
        y_min = min(geo_df_clustered.Y)
        y_max = max(geo_df_clustered.Y)

        dist = Point(x_min, y_min).distance(Point(x_max, y_max))
        if dist < resolution * 10:
            extension = resolution * 2
        elif dist < resolution * 30:
            extension = resolution * 3
        else:
            extension = resolution * 5
        area = box(minx=x_min - extension, maxx=x_max + extension,
                 miny=y_min - extension, maxy=y_max + extension)
        roads = gpd.read_file(dir+'/Roads.shp')
        roads = gpd.clip(roads,area)

        print('create points along roads..')
        roads = MultiLine_to_Line(roads)
        roads= roads.geometry.simplify(tolerance=simplify_road_coef_inside)
        roads = gpd.GeoDataFrame(geometry=roads)
        gdf_roads, roads_segments = create_roads(roads, geo_df)
    else:
        gdf_roads=gpd.GeoDataFrame()
        roads_segments=gpd.GeoDataFrame()
    if not os.path.exists(dir+'/grids'):
        os.makedirs(dir+'/grids')
    os.chdir(dir+'/grids')

    if n_terminal_nodes > 1:
        l()
        print("Creating grid for the secondary substations")

        c_grid, c_grid_cost, c_grid_length, c_grid_points = \
            cluster_grid(geo_df, gdf_cluster_pop, resolution,
                         line_bc, n_terminal_nodes, gdf_cluster, gdf_roads,
                         roads_segments, Roads_option, roads_weight, Rivers_option,length_max=1.5)
        print("Grid created")

        if not gdf_roads.empty:
            gdf_roads['Cluster']=clus
            gdf_cluster = add_roads_points_to_gdf(gdf_cluster,gdf_roads,c_grid)
            geo_df = add_roads_points_to_gdf(geo_df,gdf_roads,c_grid)
        #os.chdir('..')
        c_grid.to_file('Grid_substations.shp')
        grid_resume = pd.DataFrame(index=clusters_list.index,
                                   columns=['Grid Length [km]', 'Grid Cost [k€]',
                                            'Connection Length [km]',
                                            'Connection Cost [k€]',
                                            'Connection Type', 'Substation ID'])

        grid_resume.loc[0, 'Grid Length [km]'] = \
            c_grid_length / 1000
        grid_resume.loc[0, 'Grid Cost [k€]'] = c_grid_cost / 1000
        grid_resume.to_csv('MV-grid_resume.csv')
        return c_grid,grid_resume,geo_df,gdf_cluster


def cluster_grid(geo_df, gdf_cluster_pop, resolution, line_bc,
                 n_terminal_nodes, gdf_cluster, gdf_roads, roads_segments, Roads_option,roads_weight, Rivers_option,length_max):
    c_grid = gdf_cluster_pop
    c_grid_cost = 0
    c_grid_length = 0
    c_grid_points = []

    c_grid2, c_grid_cost2, c_grid_length2, c_grid_points2 = Spiderman. \
        spider(geo_df, gdf_cluster_pop, line_bc, resolution, gdf_roads,
               roads_segments, Roads_option, Rivers_option,roads_weight)

    if n_terminal_nodes < resolution*2 and gdf_cluster.length.size < 1000:

        if gdf_roads.empty:
            c_grid1, c_grid_cost1, c_grid_length1, c_grid_points1 = Steinerman. \
                steiner(geo_df, gdf_cluster_pop, line_bc, resolution, Rivers_option)
        else:
            c_grid1, c_grid_cost1, c_grid_length1, c_grid_points1 = Steinerman. \
                steiner_roads(geo_df, gdf_cluster_pop, line_bc, resolution,
                        gdf_roads, roads_segments, Rivers_option,roads_weight,length_max)

        # c_grid2, c_grid_cost2, c_grid_length2, c_grid_points2 = Spiderman. \
        #     spider(geo_df, gdf_cluster_pop, line_bc, resolution, gdf_roads,
        #            roads_segments)
        #


        if c_grid_cost1 <= c_grid_cost2:
            print("Steiner algorithm has the better cost")
            c_grid = c_grid1
            c_grid_length = c_grid_length1
            c_grid_cost = c_grid_cost1
            c_grid_points = c_grid_points1

        elif c_grid_cost1 > c_grid_cost2:
            print("Spider algorithm has the better cost")
            c_grid = c_grid2
            c_grid_length = c_grid_length2
            c_grid_cost = c_grid_cost2
            c_grid_points = c_grid_points2

    else:
        print("Too many points to use Steiner, running Spider.")
        c_grid = c_grid2
        c_grid_length = c_grid_length2
        c_grid_cost = c_grid_cost2
        c_grid_points = c_grid_points2

    return c_grid, c_grid_cost, c_grid_length, c_grid_points



def links(geo_df_clustered, geo_df, all_collateral, resolution, line_bc,
          grid_resume, gdf_roads, roads_segments):
    all_link = pd.DataFrame()
    l()
    print('Connecting all the people outside the clustered area..')
    l()

    gdf_clusters_out = geo_df_clustered[geo_df_clustered['Cluster'] == -1]
    gdf_clusters_out = gdf_clusters_out[gdf_clusters_out['Population'] >= 1]

    grid_points = list(zip(all_collateral.ID1.astype(int),
                           all_collateral.ID2.astype(int)))
    all_points = line_to_points(all_collateral, geo_df)
    gdf_clusters_out = gdf_clusters_out.assign(
        nearest_id=gdf_clusters_out.apply(nearest, df=all_points,
                                          src_column='ID', axis=1))

    for row in gdf_clusters_out.iterrows():
        p1 = geo_df_clustered[geo_df_clustered['ID'] == row[1].ID]
        p2 = geo_df_clustered[geo_df_clustered['ID'] == row[1].nearest_id]

        link, link_cost, link_length, _ = \
            dijkstra.dijkstra_connection_roads(geo_df, p1, p2, grid_points,
                                               line_bc, resolution, gdf_roads,
                                               roads_segments)

        all_link = gpd.GeoDataFrame(pd.concat([all_link, link], sort=True))

    #  REMOVING DUPLICATIONS
    nodes = list(zip(all_link.ID1.astype(int), all_link.ID2.astype(int)))
    nodes = list(set(nodes))
    all_link.reset_index(drop=True, inplace=True)
    all_link_no_dup = pd.DataFrame()

    for pair in nodes:
        unique_line = all_link.loc[(all_link['ID1'] == pair[0]) & (
                all_link['ID2'] == pair[1])]
        unique_line = all_link[all_link.index ==
                               unique_line.first_valid_index()]
        all_link_no_dup = gpd.GeoDataFrame(pd.concat([all_link_no_dup,
                                                      unique_line], sort=True))

    grid_resume.loc[
        0, 'Link Length [km]'] = all_link_no_dup.Length.sum() / 1000
    grid_resume.loc[0, 'Link Cost [k€]'] = all_link_no_dup.Cost.sum() / 1000
    all_link_no_dup.crs = geo_df.crs
    all_link_no_dup.to_file('all_link')
    grid_resume.to_csv('grid_resume.csv', index=False)

    print('100% electrification achieved')

    return grid_resume

def add_roads_points_to_gdf(gdf,gdf_roads,c_grid):
    nodes_in_lines = c_grid.ID1.to_list() + c_grid.ID2.to_list()
    nodes_in_lines = list(set(nodes_in_lines)) # have the unique values
    for i in nodes_in_lines:
        if i in gdf_roads.ID.to_list():
            a = gdf_roads.loc[gdf_roads['ID']==i,:]
            gdf = gdf.append(a)
    return gdf