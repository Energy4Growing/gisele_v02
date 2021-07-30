import rasterio
from sklearn.cluster import AgglomerativeClustering
from shapely.geometry import Polygon
from shapely.ops import nearest_points
from math import sqrt
from gisele import initialization
from shapely import geometry
from gisele.functions import *
from math import *
#from gisele import LV_routing_new_strategy as new_strategy
def points_region_to_casestudy(points_region, polygon_casestudy):
    points_CaseStudy = points_region.clip(polygon_casestudy)
    return points_CaseStudy

def create_polygon_from_clusters(final_clus, clusters_file, substations_file, resolution, crs):
    Clusters = gpd.read_file(clusters_file)
    Substations = gpd.read_file(substations_file)
    Clusters = Clusters[Clusters['final_clus'] == final_clus]
    k = 0
    for index, row in Clusters.iterrows():
        area = row['geometry']
        minx = area.bounds[0]
        miny = area.bounds[1]
        maxx = area.bounds[2]
        maxy = area.bounds[3]
        if k == 0:
            min_x = minx
            min_y = miny
            max_x = maxx
            max_y = maxy
            k = 1
        else:
            if minx < min_x:
                min_x = minx
            if miny < min_y:
                min_y = miny
            if maxx > max_x:
                max_x = maxx
            if maxy > max_y:
                max_y = maxy
    for index, row in Substations.iterrows():
        substation = row['geometry']
        if substation.x < min_x:
            min_x = substation.x
        if substation.y < min_y:
            min_y = substation.y
        if substation.x > max_x:
            max_x = substation.x
        if substation.y > max_y:
            max_y = substation.y

    study_area = Polygon([Point(min_x, min_y), Point(min_x, max_y), Point(max_x, max_y), Point(max_x, min_y)])
    study_area_buffered = study_area.buffer(4 * resolution)
    polygon = gpd.GeoDataFrame({'ID': [0], 'geometry': study_area_buffered})
    polygon.crs = crs
    polygon.to_file('area_polygon', index=False)
    return polygon  # geodataframe with the polygon
    ''' The goal of this function is to create a polygon starting from a cluster of clusters'''


def create_grid(crs, resolution, study_area):
    # crs and resolution should be a numbers, while the study area is a polygon
    df = pd.DataFrame(columns=['X', 'Y'])
    min_x, min_y, max_x, max_y = study_area.bounds
    # create one-dimensional arrays for x and y
    lon = np.arange(min_x, max_x, resolution)
    lat = np.arange(min_y, max_y, resolution)
    lon, lat = np.meshgrid(lon, lat)
    df['X'] = lon.reshape((np.prod(lon.shape),))
    df['Y'] = lat.reshape((np.prod(lat.shape),))
    geo_df = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df.X, df.Y),
                              crs=crs)
    geo_df_clipped = gpd.clip(geo_df, study_area)
    # geo_df_clipped.to_file(r'Test\grid_of_points.shp')
    return geo_df_clipped


def street_to_points(streets):
    streets_points = []
    for line in streets['geometry']:
        # print(line.geometryType)
        if line.geometryType() == 'MultiLineString':
            for line1 in line:
                for x in list(zip(line1.xy[0], line1.xy[1])):
                    # print(line1)
                    streets_points.append(x)
        else:
            for x in list(zip(line.xy[0], line.xy[1])):
                # print(line)
                streets_points.append(x)
    streets_multipoint = MultiPoint(streets_points)
    return streets_multipoint


def coincidence_factor(population, pop_per_household):
    return 0.35 + 0.65 / sqrt(population / pop_per_household)


def categorize_substation(clusters_list, substations):
    values = []
    costs = []
    substations['Rated_power [kVA]'] = substations['Rated_power [kVA]']
    substations1 = substations['Rated_power [kVA]'].to_list()
    for index, row in clusters_list.iterrows():
        load_kVA = row.loc['Load [kW]'] / 0.9  # considering a power factor of 0.9
        substations2 = [i - load_kVA if i - load_kVA > 0 else 10000 for i in substations1]
        power = int(min(substations2) + load_kVA)
        values.append(power)
        locate_cost = substations[substations['Rated_power [kVA]'] == power]['Cost[euro]']
        costs.append(locate_cost.values[0])

    clusters_list['Transformer_rated_power [kVA]'] = values
    clusters_list['Cost[euro]'] = costs

    return clusters_list


def locate_secondary_ss(crs, resolution, load_capita, pop_per_household, road_coef,
                        Clusters, case_study, LV_distance, ss_data,landcover_option,gisele_dir):
    dir_input = r'Case studies/' + case_study + '/Input'


    dir_output = '/Case studies/' + case_study
    grid_500m_weighted = pd.read_csv(dir_input + '/weighted_grid_of_points.csv')
    grid_500m_with_ss = grid_500m_weighted.copy()
    grid_500m_gdf = gpd.GeoDataFrame(grid_500m_weighted,
                                     geometry=gpd.points_from_xy(grid_500m_weighted.X, grid_500m_weighted.Y), crs=crs)
    grid_500m_with_ss = gpd.GeoDataFrame(grid_500m_with_ss,
                                         geometry=gpd.points_from_xy(grid_500m_with_ss.X, grid_500m_with_ss.Y), crs=crs)

    # Create a clusters.exe file
    LV_resume = pd.DataFrame()
    for index, row in Clusters.iterrows():
        os.chdir(gisele_dir)

        dir = gisele_dir + dir_output + '/Output/Clusters/' + str(row['cluster_ID'])
        clus = row['cluster_ID']
        if not os.path.exists(dir):
            os.makedirs(dir)
            os.makedirs(dir + '/grids')
        area = row['geometry']
        # THIS IS SPECIFIC FOR THE CASE OF THUSO - ISSUE NEEDS TO BE FIXED WITH CLUSTERS THAT ARE TOO NARROW. If in one line
        # of the grid points there are no points -> the graph for the steiner tree will not be connected

        area_buffered = area
        # area_buffered = row['geometry'].buffer((resolution_MV * 0.1 / 11250) / 2)
        area_list = [area_buffered]
        # Create grid of points with a 30m resolution
        grid_of_points = create_grid(crs, resolution, area)
        grid_of_points.to_file(dir + '/points.shp')
        # Load and clip protected araes and streets
        # protected_areas = gpd.read_file(global_dir+db_dir+case_study+'/Protected_areas/Protected_area-Uganda.shp')
        # protected_areas = protected_areas.to_crs(crs)
        # protected_areas_clipped=gpd.clip(protected_areas,area)
        # protected_areas_clipped.to_file(dir + '/protected.shp)

        # To make sure the roads are not cut
        min_x, min_y, max_x, max_y = area.bounds
        area_for_roads = geometry.Polygon(
            [geometry.Point(min_x, min_y), geometry.Point(min_x, max_y), geometry.Point(max_x, max_y),
             geometry.Point(max_x, min_y)])
        streets = gpd.read_file(dir_input + '/Roads.shp')
        streets = streets.to_crs(crs)
        streets_clipped = gpd.clip(streets, area_for_roads)
        if not streets_clipped.empty:
            streets_clipped.to_file(dir + '/Roads.shp')
            Street_Multipoint = street_to_points(streets_clipped)

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

        # AGLOMERATIVE CLUSTERING
        scale_factor = 10
        populated_grid = grid_of_points[grid_of_points['Population'] > 0]
        populated_grid['Population'] = populated_grid['Population'].div(scale_factor).round(0) + 0.51
        populated_grid['Population'] = populated_grid['Population'].round(0)
        populated_grid = populated_grid.loc[populated_grid.index.repeat(populated_grid.Population.astype(int))]
        loc = {'x': populated_grid['X'], 'y': populated_grid['Y']}
        pop_points = pd.DataFrame(data=loc).values

        clustering = AgglomerativeClustering(distance_threshold=LV_distance * 1.8, linkage='complete',
                                             n_clusters=None).fit(pop_points)
        geo_df_clustered = gpd.GeoDataFrame(data=pop_points, columns=['X', 'Y'],
                                            geometry=gpd.points_from_xy(populated_grid['X'], populated_grid['Y']))
        geo_df_clustered['Cluster'] = clustering.labels_
        geo_df_clustered = geo_df_clustered.drop_duplicates()
        geo_df_clustered.set_crs(epsg=crs, inplace=True)
        geo_df_clustered.to_file(dir + '/pointsCluster.shp')
        # number_clusters=max(geo_df_clustered.Cluster)
        # calculate the total distances from one point to another inside the clusters
        for cluster in range(max(geo_df_clustered['Cluster']) + 1):
            total_distances = []
            geo_df_clustered_slice = geo_df_clustered[geo_df_clustered['Cluster'] == cluster]
            print(cluster)
            for index, row2 in geo_df_clustered_slice.iterrows():
                tot_dist = 0
                for index1, row1 in geo_df_clustered_slice.iterrows():
                    tot_dist += sqrt((row2.X - row1.X) ** 2 + (row2.Y - row1.Y) ** 2)
                total_distances.append(tot_dist)
            geo_df_clustered.loc[geo_df_clustered['Cluster'] == cluster, 'tot_distance'] = total_distances

        joinDF = gpd.sjoin(grid_of_points, geo_df_clustered, how='left', op="contains")
        grid_of_points['Cluster'] = joinDF['Cluster']
        grid_of_points['Cluster'] = grid_of_points['Cluster'].fillna(-1)

        grid_of_points['tot_distance'] = joinDF['tot_distance']

        # ASSIGN ROAD DISTANCE
        if not streets_clipped.empty:
            road_distances = []
            for index, point in grid_of_points.iterrows():
                x = point['geometry'].xy[0][0]
                y = point['geometry'].xy[1][0]
                nearest_geoms = nearest_points(Point(x, y), Street_Multipoint)
                road_distance = nearest_geoms[0].distance(nearest_geoms[1])
                road_distances.append(road_distance)

            grid_of_points['Road_dist'] = road_distances

        # CHOOSE A SUBSTATION

        number_clusters = max(grid_of_points.Cluster)
        substations = []
        # Create a clusters.exe file
        clusters_list = pd.DataFrame(columns=['Cluster','Sub_cluster', 'Population', 'Load [kW]'])
        # ASSIGN MV POWER FOR THE SECONDARY SUBSTATIONS
        grid_of_points['MV_Power'] = 0
        for i in range(int(number_clusters) + 1):
            subset = grid_of_points[grid_of_points['Cluster'] == i]
            sum_pop = subset['Population'].sum()
            load = sum_pop * load_capita * coincidence_factor(sum_pop, pop_per_household)
            data = np.array([[int(row['cluster_ID']),int(i), sum_pop, load]])
            df2 = pd.DataFrame(data, columns=['Cluster','Sub_cluster', 'Population', 'Load [kW]'])

            clusters_list = clusters_list.append(df2)
            average_distance = subset['tot_distance'] / len(subset)
            min_weight = 100000
            ID = 0
            for index, ROW in subset.iterrows():
                if not streets_clipped.empty:
                    weight = ROW['Road_dist'] * road_coef + average_distance[index]
                else:
                    weight = average_distance[index]
                if weight < min_weight:
                    min_weight = weight
                    ID = index
            substations.append(ID)
        grid_of_points['Substation'] = 0
        for substation in substations:
            grid_of_points.loc[grid_of_points['ID'] == substation, 'Substation'] = 1

        for i in range(int(number_clusters) + 1):
            population = clusters_list.loc[clusters_list['Sub_cluster'] == i, 'Population'][0]
            load = clusters_list.loc[clusters_list['Sub_cluster'] == i, 'Load [kW]'][0]
            grid_of_points.loc[
                (grid_of_points['Substation'] == 1) & (grid_of_points['Cluster'] == i), 'Population'] = population
            grid_of_points.loc[
                (grid_of_points['Substation'] == 1) & (grid_of_points['Cluster'] == i), 'MV_Power'] = load

        substation_data = pd.read_csv(gisele_dir + '/general_input/' + ss_data)
        clusters_list = categorize_substation(clusters_list, substation_data)

        clusters_list['Population']=[ceil(i) for i in clusters_list['Population']]
        clusters_list.to_csv(dir + '/clusters_list.csv',index=False)
        LV_resume=LV_resume.append(clusters_list)
        weights_grid = initialization.weighting(grid_of_points,resolution,landcover_option)
        grid_of_points['Weight'] = weights_grid['Weight']
        grid_of_points.crs = crs
        grid_of_points.to_csv(dir + '/Input.csv')
        grid_of_points.to_file(dir + '/points.shp')
        secondary_substations = grid_of_points[grid_of_points['Substation'] == 1]
        secondary_substations.to_file(dir + '/substations.shp')

        total_costs = sum(clusters_list['Cost[euro]'].to_list())
        print('The total costs for substations are ' + str(total_costs / 1000) + ' thousand euros')
        # cluster_polygons_gpd.to_file(dir + '/clusters_polygons.shp')
        print('The maximum loading of a cluster is ' + str(max(clusters_list['Load [kW]'])))
        print('The minimum loading of a cluster is ' + str(min(clusters_list['Load [kW]'])))



        # study_area = Clusters.loc[row['cluster_ID']]['geometry']
        study_area = area_buffered
        area = study_area
        area.crs = crs
        grid_500m_clip = gpd.clip(grid_500m_gdf, area)  # this is our 500m resolution grid of points
        substations = grid_of_points[grid_of_points['Substation'] == 1]
        number_substations = substations.shape[0]
        ss_starting_id = int(grid_500m_with_ss['ID'].max()) + 1
        substations['ID'] = [i for i in range(ss_starting_id, ss_starting_id + number_substations)]
        # ss_starting_id+=number_substations
        substations.to_file(dir + '/substations.shp')
        substations['Cluster'] = clus
        substations_weighted =  initialization.weighting(substations, resolution, landcover_option)
        substations['Weight'] = substations_weighted['Weight']
        grid_500m_clip = grid_500m_clip.append(substations)
        grid_500m_with_ss = grid_500m_with_ss.append(substations)
        #geo_df = initialization.weighting(grid_500m_clip, resolution_MV, landcover_option)
        grid_500m_clip.to_file(dir + '/points500m.shp')
       # geo_df['Substation'] = grid_500m_clip['Substation']
       # geo_df['geometry'] = grid_500m_clip['geometry']

        clus += 1

    grid_500m_with_ss.to_csv(gisele_dir + '/' + dir_input + '/weighted_grid_of_points_with_ss.csv')
    LV_resume.to_csv(gisele_dir + '/' + dir_input +'/LV_resume.csv')
