import pandapower as pp
import geopandas as gpd
import pandas as pd
from pandapower.plotting.plotly.mapbox_plot import set_mapbox_token
from pandapower.plotting import simple_plotly
from shapely.geometry import Point
import shapely
import os
import plotly
import pandapower as pp
vn=13.8
linetypes = {"Alumoweld 4 AWG - 25 mm2": {"r_ohm_per_km": 0.01, "x_ohm_per_km": 0.02, "c_nf_per_km": 5, "max_i_ka": 0.145}}

set_mapbox_token('pk.eyJ1IjoibmExMjAzMTMyIiwiYSI6ImNrYW9qYThxMzFvb3cyc3A2cWJyaHdhdTMifQ.bfRDDy-DV4-VuVWjDNzodg')
gisele_folder = r'C:\Users\alekd\PycharmProjects\gisele_v02'
case_study = 'Isola_Giglio2'
output_folder = os.path.join(gisele_folder,'Case studies',case_study,'Output','MILP_processed')
Network_nodes = gpd.read_file(os.path.join(output_folder,'output_nodes','output_nodes.shp'))
Network_nodes=Network_nodes.to_crs(4326)
vn=15
#Cluster -1 means that it is a substation
network = pp.create_empty_network()
# Create nodes, main branches and collaterals
for col_name, node in Network_nodes.iterrows():
    try:
        x= node['geometry'].xy[0][0]
        y= node['geometry'].xy[1][0]
        pp.create_bus(net=network,vn_kv=vn,name=str(int(node.ID)),index=int(node.ID),geodata=(node.X,node.Y))
    except:
        print(node.ID)
        print('Node exists')

Branches_clusters = gpd.read_file(os.path.join(output_folder,'output_cluster_lines','output_cluster_lines.shp'))
for col_name, main_branch in Branches_clusters.iterrows():
    pp.create_line(net=network,from_bus=main_branch.id1,to_bus=main_branch.id2,length_km=main_branch.Length/1000
               ,std_type='94-AL1/15-ST1A 20.0') #Inom=140A
Branches_connections = gpd.read_file(os.path.join(output_folder,'output_connections','output_connections.shp'))
for col_name, main_branch in Branches_connections.iterrows():
    pp.create_line(net=network,from_bus=main_branch.id1,to_bus=main_branch.id2,length_km=main_branch.Length/1000
               ,std_type='94-AL1/15-ST1A 20.0') #Inom=140A



for col_name, node in Network_nodes.iterrows():
    pp.create_load(net=network, bus=int(node.ID),p_mw=node.Power/1000,
                   q_mvar=0.5*node.Power/1000)

#change in case of multiple nodes
Substations = Network_nodes[Network_nodes['Cluster']==-1]
for col_name, node in Substations.iterrows():
    try:
        pp.create_ext_grid(net=network, bus=int(node.ID), s_sc_max_mva=10000)
    except:
        print(node.ID)
        print('Node exists')




print(network)

#     pp.runpp(network, algorithm='nr')
# #transformer = Transformer.from_crs(4326, 32723)
#     plot=simple_plotly(network,on_map=True,projection='epsg:32723')
#     plotly.offline.plot(plot, filename='Cluster'+str(j)+'.html')
#     pp.to_excel(network, 'Resuts_cluster'+str(j)+'.xlsx')
#     #simple_plot(network)
#     pp.plotting.plot_voltage_profile(network)
#     pp.plotting.plotly.pf_res_plotly(network,
#                                         filename='Results_cluster'+str(j)+'.html',climits_volt=(0.85, 1.15))
    # plt.show()
    #pf_res_plotly(network

pp.runpp(network, algorithm='nr')
#plot=simple_plotly(network,on_map=True,projection='epsg:32632')
pp.plotting.plot_voltage_profile(network)
pp.to_excel(network, "results_powerflow_all_clusters.xlsx")
#plotly.offline.plot(plot, filename='All Clusters.html')
# map sttyle can be satelite, streets, light, dark, bright
pp.plotting.plotly.pf_res_plotly(network,on_map=True,map_style='streets',climits_volt=(0.92,1.08),climits_load=(0,60),bus_size=7,line_width=3)
    #Another way to plot results
    #bt = pp.plotting.plotly.create_bus_trace(network, cmap=True)
    #lt = pp.plotting.plotly.create_line_trace(network, cmap=True)
    #pp.plotting.plotly.draw_traces(lt + bt, showlegend=False)
