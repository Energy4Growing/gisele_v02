import numpy as np
import math
from gisele.geneticalgorithm_github import geneticalgorithm as ga
import pandas as pd
import pandapower as pp
import pandapower.topology
import pandapower.plotting
import geopandas as gpd
from datetime import datetime
from gisele import myFirstGA
import matplotlib.pyplot as plt

#input data
vn=33 # voltage level [kV]
vmax=1 #[pu]
vmin=0.9 #[pu]
int_nodes=True


def genetic_test(folder,iter,coe):
    links=pd.read_csv('Output/LCOE/possible_links_complete.csv')
    c_links=pd.read_csv('Output/LCOE/cost_links_complete.csv')
    len_links=pd.read_csv('Output/LCOE/len_links_complete.csv')
    c_mg=pd.read_csv('Output/LCOE/c_npc.csv',index_col='ID')
    c_load=pd.read_csv('Output/LCOE/c_power.csv',index_col='ID')
    c_energy = pd.read_csv('Output/LCOE/energy.csv', index_col='ID')
    s_power=pd.read_csv('Output/LCOE/sub_power.csv',index_col='ID')

    if folder =='Aleks':
        linetype = 'new_line'
        c_centers = pd.read_csv('Silvia/' + str(folder) + '/' +'clusters.csv')
        if int_nodes:
            c_nodes = pd.read_csv('Silvia/' + str(folder) + '/' +'output_nodes.csv',index_col='ID')
            c_internal_conn = pd.read_csv('Silvia/' + str(folder) + '/' +'output_cluster_lines.csv')
            info_links = pd.read_csv('Silvia/' + str(folder) + '/' +'links.csv')
        substations = gpd.read_file('Silvia/' + str(folder) + '/subs.shp',index_col='ID')
        substations.to_crs('epsg:21095', inplace=True)

    else:
        linetype = '184-AL1/30-ST1A 20.0'
        c_centers = \
            gpd.read_file('Silvia/'+str(folder)+'/'+str(folder)+'.shp')
        substations = gpd.read_file('Silvia/' + str(folder) + '/subs.shp')
        substations.to_crs('epsg:21095', inplace=True)

    c_centers['ID'] = ['C' + str(i) for i in c_centers.index]
    substations['ID'] = ['S' + str(i) for i in substations.index]


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

    #create minimum spanning tree
    import networkx as nx
    G = nx.Graph()
    for i,row in c_links.iterrows():
        G.add_edge(row['0'],row['1'],weight=row['Cost'])
    T=nx.minimum_spanning_tree(G)
    edges_min_span_tree = np.asarray(list(T.edges))
    degrees = np.asarray(list(T.degree))
    sorted_array = degrees[np.argsort(degrees[:, 1])]

    #find max number of connections for each cluster

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
    #add cost of energy aquired from the grid
    for i, row in c_links.iterrows():
        if row['1'] !='mg':
            c_links.loc[i,'Cost']+= coe* c_energy.loc[row['0'],'energy']


    conversion_factor = bounds['1'].cumsum()
    conversion_factor = conversion_factor[:-1].values
    conversion_factor =np.insert(conversion_factor,0,0)
    #find coordinates of cluster centers (useful to plot)
    # geo_df_clustered=geo_df_clustered[geo_df_clustered['Cluster']!=-1]
    # c_centers=geo_df_clustered.groupby('Cluster').mean()


    #create pandapower model
    #check which is the total load power
    Load_max= c_nodes ['Power'].sum()/1000 #['MW']
    # s_power['id_network'] =
    # s_power.loc[j, 'PowerAvailable'] / 1000

    net=pp.create_empty_network()
    id=1
    global line_id
    line_id=1
    line_data = {"c_nf_per_km": 4, "r_ohm_per_km": 0.42,
                 "x_ohm_per_km": 0.39, "max_i_ka": 0.213}
    pandapower.create_std_type(net, line_data, "new_line", element='line')
    #create the loads

    #create internal cluster nodes if the option is selected (to compute internal voltage drop)
    if int_nodes:
        for ii in c_nodes.index:
            if 'S' not in c_nodes.loc[ii,'Cluster']:
                pp.create_bus(net, vn, index=id, name=ii,
                              geodata=(c_nodes.loc[ii, 'X'],
                                       c_nodes.loc[ii, 'Y']))
                pp.create_load(net, id, c_nodes.loc[ii, 'Power'] / 1000,
                               q_mvar=c_nodes.loc[ii, 'Power'] / 1000 * 0.5)
                id = id + 1
        for iii, row in c_internal_conn.iterrows():
            from_bus_id = net.bus[net.bus['name']==row['id1']].index.values[0]
            to_bus_id = net.bus[net.bus['name']==row['id2']].index.values[0]
            pp.create_line(net, from_bus_id, to_bus_id,
                           row['Length']/1000,
                           std_type=linetype, index=line_id)
            line_id = line_id + 1
        #put an external grid that has a high power availability, not to create loadflow problems
        for j in subs_nodes:
            sub_name =c_nodes.loc[c_nodes['Cluster']==j].index.values[0]
            pp.create_bus(net,vn,index=id,name=sub_name,
                          geodata=(substations.loc[int(j.split('S')[1]),'geometry'].x,
                                   substations.loc[int(j.split('S')[1]),'geometry'].y))
            pp.create_ext_grid(net,id,max_p_mw=Load_max*2)
            id=id+1


    else:
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

    #create initial solution as the one from minimum spanning tree
    #start creating links from the clusters with only one connection
    '''
    presel_links=np.zeros(dim)
    for w in sorted_array[:,0]:
        if 'C' in w: #exclude substations
            start_node = w
            start_node_id = nodes[nodes==w].index.values[0]-1
            sel_edge = np.where(edges_min_span_tree==w)
            if sel_edge[1][0] ==0:
                end_node = edges_min_span_tree[sel_edge[0][0],1]
            else:
                end_node = edges_min_span_tree[sel_edge[0][0], 0]
            presel_link_id = links[(links['0']==start_node) & (links['1']==end_node)]['ID'].values[0]
            #presel_links[start_node_id-1] =presel_link_id

            presel_links[start_node_id]=(presel_link_id - (varbound[start_node_id, 0])) / (varbound[start_node_id, 1] - varbound[start_node_id, 0])
            #remove line from edges
            edges_min_span_tree =np.delete(edges_min_span_tree,sel_edge[0][0],axis=0)
        '''

    C_max = c_links.groupby('0').max()['Cost'].sum()

    #X=presel_links
    #X=[2., 3., 6., 7., 5., 2., 2., 7., 4.] #optimal solution, equal to the one of Aleks
    # X= [2., 1., 5., 3., 7., 2., 8., 3., 6.] #load flow not converging
  #   X= [14. , 5.,  1.,  1.,  2.,  5. , 5. , 4. , 7.,  9. , 4.,  8.,  3.,  8. , 3.,  2. , 7.,  1.,
  # 4. , 1.] cost 9560
  #   X = [14. , 1. , 5. , 1.,  5.,  5.,  5.,  1. , 4. ,14. , 8. , 8. , 2. , 7. , 6.,  2.,  1.,  1.,
  #    3. , 1.] cost 9077
  #   X=[14. , 5.,  3. , 1. , 1. , 4. , 5. , 4. , 4. , 9. , 2.,  8.,  1.,  7. , 6. , 5. , 7. , 1.,
  #    4. , 1.] cost 9010
  #   X=[14.  ,1. , 5. , 1. , 1. , 4.,  5.,  4. , 4. ,10. , 2. , 8. , 2. , 4. , 6. , 2. , 7. , 1.,
  #    4. , 1.] #cost 9070
  #   X= [ 9.,  1. , 3. , 1. , 1. , 4. , 5. , 4. , 3. ,13.,  2.,  8.,  2.,  4. , 3. , 5. , 7. , 1.,
  # 4. , 1.]
   # X=[14.,  5.,  5.,  1.,  1.,  4.,  5.,  4.,  4., 13. , 2.,  8.,  1.,  4. , 3.,  2. , 5. , 1.,
    # 4. , 1.]

    def f1(X):
        #create integers (only for Nic function)
        #X = np.round(varbound[:, 0]+X*(varbound[:, 1]-varbound[:, 0]))

        sel_links = conversion_factor+X-1
        C_real=np.sum(c_links.loc[sel_links,'Cost'])
        line_id_post=line_id+1
        gen_id = 0
        penalty= 0
        mg_buses=[]


    # create the links from the results of the iterations
        for i in range(len(X)):
            sel_link=sel_links[i]
            from_clust = links.loc[sel_link, '0']
            to_clust = links.loc[sel_link, '1']

            if to_clust =='mg':
                mg_buses.append(from_clust)
                from_bus = c_nodes[c_nodes['Cluster']==from_clust].index.values[0] #take first internal node in the cluster
                from_bus_id = net.bus[net.bus.name==from_bus].index.values[0]
                pp.create_sgen(net,from_bus_id,c_load.loc[from_clust,'Load [kW]']/1000,
                               q_mvar=c_load.loc[from_clust,'Load [kW]']/1000*0.5)
                gen_id=gen_id+1
            else:
                row_link = info_links.loc[(info_links.Cluster == from_clust)&(info_links.ID == to_clust)|(info_links.Cluster == to_clust)&(info_links.ID == from_clust)]
                from_bus = row_link['ID1'].values[0]
                to_bus = row_link['ID2'].values[0]
                from_bus_id =net.bus[net.bus.name==from_bus].index.values[0]
                to_bus_id = net.bus[net.bus.name == to_bus].index.values[0]
                pp.create_line(net,from_bus_id,to_bus_id,len_links.loc[sel_link,'Length'],std_type=linetype,index=line_id_post)
                line_id_post=line_id_post+1
        #long step, need to simplify in this case, no need to consider all internal nodes
        uns_buses = pandapower.topology.unsupplied_buses(net)
        uns_buses_names = net.bus.loc[uns_buses,'name']
        uns_buses=set(c_nodes.loc[uns_buses_names,'Cluster'].unique())
        #penalty if there are no unserved buses (need to increase it with the number of buses)
        diff= uns_buses -set(mg_buses)
        # if uns buses are more than the microgrids, put each unserved bus as a microgrid
        if diff !=set():
            for from_clust in diff:
                from_bus = \
                c_nodes[c_nodes['Cluster'] == from_clust].index.values[
                    0]  # take first internal node in the cluster
                from_bus_id = net.bus[net.bus.name == from_bus].index.values[0]
                pp.create_sgen(net, from_bus_id,
                               c_load.loc[from_clust, 'Load [kW]'] / 1000,
                               q_mvar=c_load.loc[
                                          from_clust, 'Load [kW]'] / 1000 * 0.5)
                gen_id = gen_id + 1
                C_real = C_real+c_mg.loc[from_clust, 'mg_npc']

             #penalty=C_real*10**(abs(len(uns_buses)-len(mg_buses))) #potrebbe essere alto (mettere 200*20 + 10*(numero di bus non supplied)

        # diag=pandapower.diagnostic(net,warnings_only=False,report_style='None')
        # try:
        #     diag['overload']
        #     if diag['overload']['load']:
        #         print('overload')
        #         con_comp= list(pandapower.topology.connected_components(
        #             pandapower.topology.create_nxgraph(net)))
        #     #check power not balanced
        #         for j in con_comp:
        #             if len(j)>1:
        #                 power=0
        #                 for sset in j:
        #                     if sset in net.sgen['bus'].values:
        #                         power = power - net.sgen[
        #                             net.sgen['bus'] == sset].p_mw.values[0]                       #if connected bus a microgrid, no overload problem
        #                     if sset in net.load['bus'].values: #check if connected bus is a load
        #                         power=power+net.load[net.load['bus']==sset].p_mw.values[0]
        #                     else: #bus is a substation
        #                         power = power - net.ext_grid[net.ext_grid['bus']==sset].max_p_mw.values[0]
        #                 if power>0:
        #                     penalty=penalty+power*10
        #         penalty = penalty +C_max*2


            #print('try loadflow')
        try:
            pp.runpp(net, algorithm='nr')
            #if substations are overloaded
            if (net.res_ext_grid.p_mw.values>s_power['PowerAvailable'].values / 1000).any():
                print('overload')
                penalty=C_max*2 + (net.res_ext_grid.p_mw.values-s_power['PowerAvailable'].values / 1000).sum()*1000


            #penalties if voltage is below limits or current is above limits
            elif (sum((net.res_bus.vm_pu>vmax)| (net.res_bus.vm_pu<vmin))
                            +sum((net.res_line.loading_percent>100))
                            +sum(net.res_ext_grid.p_mw>net.ext_grid.max_p_mw))>0:

                penalty=C_max+10*(sum((net.res_bus.vm_pu>vmax)| (net.res_bus.vm_pu<vmin))
                            +sum((net.res_line.loading_percent>100))
                            +sum(net.res_ext_grid.p_mw>net.ext_grid.max_p_mw))

        except :
            print('LoadflowNotConverged')
            print(X)
            penalty=C_max*3

        #remove al the lines after running load flow
        pp.drop_lines(net,range(line_id+1,line_id_post))
        # #how to remove generators?
        net.sgen.drop(range(gen_id),inplace=True)

        if penalty >0:
            #print('Penalty =' + str(penalty))
            return penalty
        else:
            #print('Cost ='+ str(C_real))
            return C_real

    def f(X):
        #create integers (only for Nic function)
        #X = np.round(varbound[:, 0]+X*(varbound[:, 1]-varbound[:, 0]))

        sel_links = conversion_factor+X-1
        C_real=np.sum(c_links.loc[sel_links,'Cost'])
        line_id=1
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

        diag=pandapower.diagnostic(net,warnings_only=False,report_style='None')
        try:
            diag['overload']
            if diag['overload']['load']:
                print('overload')
                con_comp= list(pandapower.topology.connected_components(
                    pandapower.topology.create_nxgraph(net)))
            #check power not balanced
                for j in con_comp:
                    if len(j)>1:
                        power=0
                        for sset in j:
                            if sset in net.sgen['bus'].values: #if connected bus a microgrid, no overload problem
                                continue
                            elif sset <dim+1: #check if connected bus is a load
                                power=power+net.load[net.load['bus']==sset].p_mw.values[0]
                            else: #bus is a substation
                                power = power - net.ext_grid[net.ext_grid['bus']==sset].max_p_mw.values[0]
                        if power>0:
                            penalty=penalty+power*10
                penalty = penalty +C_max*2


        except:
            #print('try loadflow')
            try:
                pp.runpp(net, algorithm='nr')
                #penalties if voltage is below limits or current is above limits
                if (sum((net.res_bus.vm_pu>vmax)| (net.res_bus.vm_pu<vmin))
                                +sum((net.res_line.loading_percent>100))
                                +sum(net.res_ext_grid.p_mw>net.ext_grid.max_p_mw))>0:

                    penalty=C_max+10*(sum((net.res_bus.vm_pu>vmax)| (net.res_bus.vm_pu<vmin))
                                +sum((net.res_line.loading_percent>100))
                                +sum(net.res_ext_grid.p_mw>net.ext_grid.max_p_mw))

            except :
                print('LoadflowNotConverged')
                penalty=C_max*3

        #remove al the lines after running load flow
        pp.drop_lines(net,range(1,line_id))
        # #how to remove generators?
        net.sgen.drop(range(gen_id),inplace=True)

        if penalty >0:
            #print('Penalty =' + str(penalty))
            return penalty
        else:
            #print('Cost ='+ str(C_real))
            return C_real

    #### pandapower backward forward model #####
    #net =pandapower.create_empty_network()
    #create buses (outside the loop)
    #create lines (inside the loop), making a for loop for each cluster,substations
    #are external grid with p.u. voltage

    #if the line is a microgrid, add external grid connection, if the line
    #if power flow is not respected because grid is not connected it will raise an error,
    #it should be inserted in the cost function

    #pandapower.topology.unsupplied_buses(net, mg=None, in_service_only=False, slacks=None, respect_switches=True)
    # pandapower.runpp(net,algorithm='bfsw',init='flat',max_iteration=1)
    # is it possible to run pandapower with two separate networks? it seems so
    #check if there are unsupplied buses: pandapower.topology.unsupplied_buses(net)

    #crossover probability =0.95, dovrebbe essere alto
    #mutation_probability default= 0.01

    algorithm_param = {'max_num_iteration': 300,\
                       'population_size':30,\
                       'mutation_probability':0.1,\
                       'elit_ratio': 0.02,\
                       'crossover_probability': 0.9,\
                       'parents_portion': 0.3,\
                       'crossover_type':'uniform',\
                       'max_iteration_without_improv':None}
    #f1(X)
    if int_nodes==False:
        model=ga(function=f,dimension=dim,variable_type='int',variable_boundaries=varbound,
                 function_timeout=20000,algorithm_parameters=algorithm_param)
    else:
        model = ga(function=f1, dimension=dim, variable_type='int',
                   variable_boundaries=varbound,
                   function_timeout=20000,
                   algorithm_parameters=algorithm_param)
    #
    model.run()


    #Ale Niccolai
    #f(X)
    '''
    minValue, avgValue, bestSolution =myFirstGA.GA(f,dim,presel_links,varbound)
    #
    plt.figure()
    plt.plot(minValue)
    plt.plot(avgValue)
    plt.grid(True)
    plt.xlabel('Iterations')
    plt.ylabel('Cost value')
    plt.legend(['Minimum value','Average value'])
    plt.savefig('Silvia/'+str(folder)+'/'+str(iter)+'.png')
    print(minValue)
    return minValue, bestSolution '''
