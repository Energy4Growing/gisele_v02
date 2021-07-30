import pandas as pd
import geopandas as gpd
from math import *
import os
os.chdir('..')
points=gpd.read_file('SS_Clustering/Katakwi2956/points500m.shp')
for index,point in points.iterrows():
    if isnan(point['MV_Power']):
        points.loc[index,'MV_Power']=0


for index,point in points.iterrows():
    if isnan(point['Substation']):
        points.loc[index,'Substation']=0
    else:
        points.loc[index,'Substation']=1

points=points.drop(columns=['Elevation','Population','Slope','Land_cover','Protected_','tot_distan','Road_dist','Cluster','geometry'],axis=1)
points.to_csv('SS_Clustering/Adjumani2564/points_MILP.csv')
