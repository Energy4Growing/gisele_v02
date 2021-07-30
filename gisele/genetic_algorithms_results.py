import numpy as np
import math
from gisele.geneticalgorithm_github import geneticalgorithm as ga
import pandas as pd
import pandapower as pp
import pandapower.topology
import pandapower.plotting
import geopandas as gpd
from datetime import datetime
import matplotlib.pyplot as plt

#input data
vn=33 # voltage level [kV]
vmax=1.1 #[pu]
vmin=0.9 #[pu]
linetype='184-AL1/30-ST1A 20.0'
# todo -> adapt results to the output of the MILP
def genetic_result(folder,iter):
    links=pd.read_csv('Output/LCOE/possible_links_complete.csv')
    c_links=pd.read_csv('Output/LCOE/cost_links_complete.csv')
    len_links=pd.read_csv('Output/LCOE/len_links_complete.csv')
    c_mg=pd.read_csv('Output/LCOE/c_npc.csv',index_col='ID')
    c_load=pd.read_csv('Output/LCOE/c_power.csv',index_col='ID')
    s_power=pd.read_csv('Output/LCOE/sub_power.csv',index_col='ID')
    if folder =='Aleks':
        linetype = 'new_line'
        c_centers = pd.read_csv('Silvia/' + str(folder) + '/' +'clusters.csv')
    else:
        c_centers = \
            gpd.read_file('Silvia/'+str(folder)+'/'+str(folder)+'.shp')
    substations = gpd.read_file('Silvia/' + str(folder) + '/subs.shp')
    substations.to_crs('epsg:21095', inplace=True)
    c_centers['ID'] = ['C' + str(i) for i in c_centers.index]
    substations['ID'] = ['S' + str(i) for i in substations.index]
    #solutions = pd.read_csv('Silvia/' + str(folder) + '/ga_solutions.csv',header=None)
    #solutions = solutions.values
    #solutions=solutions.reshape(10,15)

    #create list of all the nodes with unique id
    cluster_nodes=pd.Series(c_mg.index).sort_values()
    subs_nodes=pd.Series(s_power.index).sort_values()
    nodes=cluster_nodes.append(subs_nodes)
    nodes.index=range(1,len(nodes.index)+1)
    #remove lines where first column is a secondary substation
    links=links[~links['0'].str.contains("S")]
    c_links=c_links[~c_links['0'].str.contains("S")]
    len_links=len_links[~len_links['0'].str.contains("S")]

    links.drop_duplicates(['0', '1'], inplace=True)
    c_links.drop_duplicates(['0', '1'], inplace=True)
    len_links.drop_duplicates(['0', '1'], inplace=True)


    bounds=links.groupby('0').count()+1
    bounds.sort_index(inplace=True)

    #number of variables equal to number of clusters
    # (useful to reduce numb of variables and speed up the process)
    dim=np.shape(bounds)[0]
    #each variables is integer and a ID that varies from 1 to the

    #number of possible connections of the cluster
    varbound=np.ones((dim,2))
    varbound[:,1]=bounds['1'].values

    #add one line to the links representing the microgrid electrification
    for cluster in bounds.index:
        links = links.append(pd.DataFrame(data=[[cluster,'mg']],columns=['0','1']),ignore_index=True)
        c_links = c_links.append(pd.DataFrame(data=[[cluster, 'mg', c_mg.loc[cluster,'mg_npc']]], columns=['0', '1','Cost']),
                             ignore_index=True)
        len_links = len_links.append(
            pd.DataFrame(data=[[cluster, 'mg', 0]],
                         columns=['0', '1', 'Length']),
            ignore_index=True)

    #sort values of dataframe by cluster number
    links.sort_values(['0','1'],inplace=True)
    c_links.sort_values(['0','1'],inplace=True)
    len_links.sort_values(['0','1'],inplace=True)


    links.index=range(len(links.index))
    c_links.index=range(len(c_links.index))
    len_links.index=range(len(len_links.index))
    #need to associate value of variables to the specific link
    links['ID']=links['0']
    links['ID']=1
    links_id=links.groupby('0').cumsum()
    links['ID']=links_id['ID']
    c_links['ID']=links_id['ID']
    len_links['ID']=links_id['ID']
    #establish the conversion factor to pass from variable to the link
    #todo -> use double index notation to mak this part more readable
    conversion_factor = bounds['1'].cumsum()
    conversion_factor = conversion_factor[:-1].values
    conversion_factor =np.insert(conversion_factor,0,0)
    #find coordinates of cluster centers (useful to plot)
    # geo_df_clustered=geo_df_clustered[geo_df_clustered['Cluster']!=-1]
    # c_centers=geo_df_clustered.groupby('Cluster').mean()


    #create pandapower model

    net=pp.create_empty_network()
    id=1
    #create the loads
    for i in bounds.index:
        pp.create_bus(net,vn,index=id,name=i,
                      geodata=(c_centers.loc[int(i.split('C')[1]),'x'],
                               c_centers.loc[int(i.split('C')[1]),'y']))
        pp.create_load(net,id,c_load.loc[i,'Load [kW]']/1000,q_mvar=c_load.loc[i,'Load [kW]']/1000*0.5)
        id=id+1
    #create substations
    for i in subs_nodes:
        pp.create_bus(net,vn,index=id,name=i,
                      geodata=(substations.loc[int(i.split('S')[1]),'geometry'].x,
                               substations.loc[int(i.split('S')[1]),'geometry'].y))
        pp.create_ext_grid(net,id,max_p_mw=s_power.loc[i,'PowerAvailable']/1000)
        id=id+1

    #X= solutions[iter,:]
    #X=[3. ,1. ,3., 5., 3., 1., 2., 1., 4., 5., 7., 4., 1., 2., 4.]
    #X=[5.,5.,8.,1.,4.,4.,1.,4.,1.,1.,4.,1.,6.,5.,2.,3.,2.,8.,7.,2.,5.,2.,2.,2.
#,2.,3.,3.,5.,2.,2.,2.,2.,3.,2.,5.,3.,3.,3.,2.,4.,3.,1.,3.,4.,5.,3.,7.,1.
#,4.,6.]
    #X=[3., 1., 3., 6., 3., 4., 1., 1., 4., 5., 2., 4. ,4., 5., 1.]
    #X=[3., 1., 7., 1., 3., 1., 2. ,1. ,4., 7., 5., 4., 1., 2., 4.]
    #X=[9., 6., 4., 1., 2., 4., 5., 1., 4., 7., 5., 4., 4. ,5., 2.]
    #X=[9. ,3., 4., 1., 5., 4., 5., 1., 5., 5., 4., 4., 4., 5., 6.]
    #X=[9., 6., 4., 1. ,2., 4. ,5., 1. ,1., 5., 6., 7., 3. ,5., 4.]
    #X=[9., 3., 4., 7., 5., 1., 5., 2., 4., 7., 5., 4., 6., 3., 4.]
    #X= [9., 3., 3., 7., 5., 5., 2., 1., 4., 5., 7., 4., 4., 2., 6.]
    #X=[9., 4., 4., 2. ,2., 4., 5., 1., 4., 7., 5. ,4., 4., 5., 6.]
    #X=[9., 3., 4., 1., 5., 4., 5., 1., 4., 5., 7., 4., 4., 5., 6.]
    #X=[9., 4., 6., 2., 3., 5., 7., 4., 7., 4., 5., 5., 2., 2., 1., 1., 9., 2., 1., 2.]
#     X= [2.,7.,5.,6.,3.,1.,4.,4.,5.,2.,6.,1.,9.,3.,2.,6.,2.,2.,7.,7.,2.,2.,4.,3.
# ,1.,1.,1.,4.,3.,5.,5.,2.,4.,5.,7.,7.,4.,4.,6.,3.,3.,6.,3.,4.,8.,3.,8.,2.
# ,3.,9.]

    # X=[7., 7., 3., 1., 3., 1., 6., 3., 2., 1., 5., 4., 9., 3., 5., 1., 2., 4.,
    #  4., 7., 5., 2., 2., 2.
    #     , 6., 5., 5., 4., 1., 5., 6., 2., 4., 2., 4., 3., 3., 3., 4., 4., 2.,
    #  6., 7., 4., 7., 3., 2., 3.
    #     , 3., 6.]
    #X=[9., 2., 5., 2., 6., 3., 4., 4., 4., 4., 5., 5., 6., 2., 1., 2., 9., 2., 4., 7.]
    #[9. 6. 6. 7. 3. 3. 3. 4. 4. 4. 5. 5. 2. 3. 2. 6. 9. 2. 1. 7.]
    #X=[2., 3., 3., 7., 5., 6., 2., 7., 4.]
    X=[2., 3., 6., 7., 5., 2., 2., 7., 4.]
    C_max = c_links.groupby('0').max()['Cost'].sum()
    line_data = {"c_nf_per_km": 0, "r_ohm_per_km": 0.42,
                 "x_ohm_per_km": 0.39, "max_i_ka": 0.213}
    pandapower.create_std_type(net, line_data, "new_line", element='line')

    def f(X):
        #create integers (only for Nic function)
        #X = np.round(varbound[:, 0]+X*(varbound[:, 1]-varbound[:, 0]))

        sel_links = conversion_factor+X-1
        C_real=np.sum(c_links.loc[sel_links,'Cost'])
        line_id = 1
        gen_id = 0
        penalty= 0
    # create the links from the results of the iterations
        for i in range(len(X)):
            sel_link=sel_links[i]
            from_bus = links.loc[sel_link, '0']
            to_bus = links.loc[sel_link, '1']
            from_bus_id =nodes[(nodes==from_bus)].index.tolist()[0]
            if to_bus =='mg':
                pp.create_sgen(net,from_bus_id,c_load.loc[from_bus,'Load [kW]']/1000,
                               q_mvar=c_load.loc[from_bus,'Load [kW]']/1000*0.5)
                gen_id=gen_id+1
            else:
                to_bus_id = nodes[nodes == to_bus].index.tolist()[0]
                pp.create_line(net,from_bus_id,to_bus_id,len_links.loc[sel_link,'Length'],std_type=linetype,index=line_id)
                line_id=line_id+1
        uns_buses = pandapower.topology.unsupplied_buses(net)
        mg_buses = set(net.sgen.bus.values)
        #penalty if there are no unserved buses (need to increase it with the number of buses)
        diff= uns_buses -mg_buses
        # if uns buses are more than the microgrids, put each unserved bus as a microgrid
        if diff !=set():
            for from_bus_id in diff:
                pp.create_sgen(net, from_bus_id,
                           c_load.loc[nodes[from_bus_id], 'Load [kW]'] / 1000,
                           q_mvar=c_load.loc[
                                      nodes[from_bus_id], 'Load [kW]'] / 1000 * 0.5)
                gen_id = gen_id + 1
                C_real = C_real+c_mg.loc[nodes[from_bus_id], 'mg_npc']

             #penalty=C_real*10**(abs(len(uns_buses)-len(mg_buses))) #potrebbe essere alto (mettere 200*20 + 10*(numero di bus non supplied)
        pp.runpp(net, algorithm='nr')
        pp.plotting.plotly.pf_res_plotly(net)
        pp.plotting.simple_plotly(net)




    f(X)
    return
folder='Aleks'
for i in range(1):
    genetic_result(folder, i)