import geopandas as gpd
import rasterio
import numpy as np
from shapely.ops import unary_union
import pandas as pd
from sklearn.linear_model import LinearRegression
from matplotlib import pyplot as plt


river_discharge = gpd.read_file(r'C:\Users\silvi\Progetti Github\Gisele_development\Case studies/test\Intermediate files\Study area\Geodata/rchresults.shp')
elevation = rasterio.open(r'C:\Users\silvi\Progetti Github\Gisele_development\Case studies\test\Intermediate files\Study area\Geodata/DEM.tif')
osm = gpd.read_file(r'C:\Users\silvi\Progetti Github\Gisele_development\Case studies/test\Intermediate files\Study area\Geodata/osm_clipped.shp')

osm_line=osm.geometry.unary_union
sampled_points=river_discharge.geometry.copy(deep=True)
river_discharge['Power kW'] =0
river_discharge['aver_head m'] =0

#create points at a predefined distance along lines
for i, row in river_discharge.iterrows():
    line=row.geometry
    distances = np.arange(0, line.length, 200)
    points = [line.interpolate(distance) for distance in distances] + [line.boundary[1]]
    multipoint = unary_union(points)
    sampled_points[i] = (multipoint)
    # gpd.GeoSeries(multipoint).plot()
    # gpd.GeoSeries(line).plot()

#sample elevation layer
#reproject points on osm layer to avoid big errors in DEM estimation
#actually the error remains even with reprojection, it is difficult to estimate head between couple of points
#take an average of head along the river

    coords = [(osm_line.interpolate(osm_line.project(point)).x, osm_line.interpolate(osm_line.project(point)).y) for point in points]
    coords_old=[(line.interpolate(line.project(point)).x, line.interpolate(line.project(point)).y) for point in points]
    data = [x[0] for x in elevation.sample(coords_old)]
    data_old = [x[0] for x in elevation.sample(coords_old)]
    X = [x[0] for x in coords]
    Y = [x[1] for x in coords]
    X_old = [x[0] for x in coords_old]
    Y_old = [x[1] for x in coords_old]
    a = pd.DataFrame(data={'X':X,'Y':Y,'DEM':data})
    b=pd.DataFrame(data={'X':X_old,'Y':Y_old,'DEM':data_old})
    #check direction of river by looking at hom many delta h are in one direction and how many in the other
    diff = np.array(data[1:])-np.array(data[:-1])
    x=np.arange(len(data)).reshape((-1, 1))
    model = LinearRegression().fit(x, np.array(data))
    r_sq = model.score(x, np.array(data))
    print(r_sq)
    y_pred = model.predict(x)
    # fig,ax =plt.subplots()
    # ax.plot(x,np.array(data),'o')
    # ax.plot(x, y_pred, '-')
    # plt.show()
    river_discharge.loc[i, 'aver_head m'] = float(abs(model.coef_))
    river_discharge.loc[i,'Power kW'] =float(abs(model.coef_)*row.FLOW_INcms*9.8)
# a.to_csv(r'C:\Users\silvi\Progetti Github\Gisele_development\Case studies/test\Intermediate files\Study area\Geodata/osm_coords.csv')
# b.to_csv(r'C:\Users\silvi\Progetti Github\Gisele_development\Case studies/test\Intermediate files\Study area\Geodata/coords.csv')

river_discharge.to_file(r'C:\Users\silvi\Progetti Github\Gisele_development\Case studies/test\Intermediate files\Study area\Geodata/river_power_old.json')
a=0