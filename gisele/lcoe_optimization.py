from __future__ import division

from pyomo.opt import SolverFactory
from pyomo.core import AbstractModel
from pyomo.dataportal.DataPortal import DataPortal
from pyomo.environ import *

import pandas as pd
import numpy as np
from datetime import datetime


def cost_optimization(p_max_lines, coe, nation_emis,nation_rel, line_rel, input_michele):
    #useful parameters from michele
    proj_lifetime =input_michele['num_years']

    # Initialize model
    model = AbstractModel()
    data = DataPortal()

    # Define sets

#    model.of = Set(initialize=[1,2])  # Set of objective functions
    model.of = Set(initialize=['cost', 'emis','rel'])  # Set of objective functions

    model.N = Set()  # Set of all nodes, clusters and substations
    data.load(filename=r'Output/LCOE/set.csv', set=model.N)

    model.clusters = Set()  # Set of clusters
    #data.load(filename='Output/LCOE/clusters.csv', set=model.clusters)

    model.substations = Set()  # Set of substations
    data.load(filename='Output/LCOE/subs.csv', set=model.substations)

    model.links = Set(dimen=2,
                      within=model.N * model.N)  # in the csv the values must be delimited by commas
    data.load(filename='Output/LCOE/possible_links_complete.csv', set=model.links)

    # Nodes are divided into two sets, as suggested in https://pyomo.readthedocs.io/en/stable/pyomo_modeling_components/Sets.html:
    # NodesOut[nodes] gives for each node all nodes that are connected to it via outgoing links
    # NodesIn[nodes] gives for each node all nodes that are connected to it via ingoing links

    def NodesOut_init(model, node):
        retval = []
        for (i, j) in model.links:
            if i == node:
                retval.append(j)
        return retval

    model.NodesOut = Set(model.N, initialize=NodesOut_init)

    def NodesIn_init(model, node):
        retval = []
        for (i, j) in model.links:
            if j == node:
                retval.append(i)
        return retval

    model.NodesIn = Set(model.N, initialize=NodesIn_init)

    def NodesOutSub_init(model, node):
        retval = []
        for (i, j) in model.links:
            if i == node:
                retval.append(j)
        return retval

    model.NodesOutSub = Set(model.substations, initialize=NodesOutSub_init)

    def NodesInSub_init(model, node):
        retval = []
        for (i, j) in model.links:
            if j == node:
                retval.append(i)
        return retval

    model.NodesInSub = Set(model.substations, initialize=NodesInSub_init)

    #####################Define parameters#####################

    # Direction of optimization for each objective function: -1 to minimize, +1 to maximize
    model.dir = Param(model.of)
    # Weight of objective functions in multi-objective optimization, range 0-1, sum over model.of has to be 1
    model.weight = Param(model.of)


    data.load(filename='Output/LCOE/data_MO.dat')

    # Parameters identifying the range of variation of each objective function (needed for normalization)
    model.min_obj = Param(model.of, initialize=0, mutable=True)
    model.max_obj = Param(model.of, initialize=1, mutable=True)

    # Electric power in the nodes (injected (-) or absorbed (+))
    model.p_clusters = Param(model.clusters)
    data.load(filename='Output/LCOE/c_power.csv', param=model.p_clusters)

    # Maximum power supplied by substations
    model.p_max_substations = Param(model.substations)
    data.load(filename='Output/LCOE/sub_power.csv',
              param=model.p_max_substations)

    # Total net present cost of microgrid to supply each cluster
    model.c_microgrids = Param(model.clusters)
    #data.load(filename='Output/LCOE/c_npc.csv', param=model.c_microgrids)

    # Total net present cost of substations
    model.c_substations = Param(model.substations)
    data.load(filename='Output/LCOE/sub_npc.csv',
              param=model.c_substations)

    # Connection cost of the possible links
    model.c_links = Param(model.links)
    data.load(filename='Output/LCOE/cost_links_complete.csv', param=model.c_links)

    # Energy consumed by each cluster in microgrid lifetime
    model.energy = Param(model.clusters)
    #data.load(filename='Output/LCOE/energy.csv', param=model.energy)

    # CO2 emission produced by each cluster in microgrid lifetime
    model.emission = Param(model.clusters)
    #data.load(filename='Output/LCOE/emissions.csv', param=model.emission)

    # CO2 emission related to construction of power infrastructure
    model.em_links = Param(model.links)
    data.load(filename='Output/LCOE/em_links.csv', param=model.em_links)

    #lol due to microgrid components_failure
    model.rel_mg = Param(model.clusters)
    # data.load(filename='Output/LCOE/mg_rel.csv', param=model.rel_mg)

    # Connection length associated to the nodes
    model.d_nodes = Param(model.N)
    data.load(filename='Output/LCOE/len_nodes.csv', param=model.d_nodes)

    # Connection length of the possible links
    model.d_links = Param(model.links)
    data.load(filename='Output/LCOE/len_links_complete.csv', param=model.d_links)

    # poximum power flowing on lines
    # model.p_max_lines = Param()  # max power flowing on MV lines
    # data.load(filename='Input/data_procedure2.dat')

    # M_max and M_min, values required to linearize the problem
    model.M_max = Param(initialize=10000)
    model.M_min = Param(initialize=-10000)

    data.load(filename='Output/Microgrids/microgrids.csv',
              select=('Cluster','Total Cost [kEUR]','Energy Demand [MWh]','CO2 [kg]','Unavailability [MWh/y]'),
              param=(model.c_microgrids,model.energy,model.emission,model.rel_mg), index=model.clusters)

    #####################Define variables#####################

    # objective function variables
    model.obj = Var(model.of, within=NonNegativeReals)
    # auxiliary variables for normalization step
    model.aux = Var(model.of)
    # normalized objective functions
    model.norm_obj = Var(model.of, within=NonNegativeReals)
    # binary variable x[i,j]: 1 if the connection i,j is present, 0 otherwise,initialize=x_rule
    model.x = Var(model.links, within=Binary)
    # binary variable y[i]: 1 if a substation is installed in node i, 0 otherwise,initialize=y_rule
    model.y = Var(model.substations, within=Binary)
    # binary variable z[i]: 1 if a microgrid is installed in node i, 0 otherwise,initialize=z_rule
    model.z = Var(model.clusters, within=Binary)
    # power[i,j] is the power flow of connection i-j
    model.P = Var(model.links, within=NonNegativeReals)
    # power[i] is the power provided by substation i
    model.p_substations = Var(model.substations, within=NonNegativeReals)
    # # variables k(i,j) is the variable necessary to linearize
    model.k = Var(model.links)
    #distance of cluster from substation
    model.dist = Var(model.N, within=NonNegativeReals)
    #lol due to MV lines
    model.lol_line = Var (model.clusters, within=NonNegativeReals)

    #####################Define constraints###############################

    def Radiality_rule(model):
        return summation(model.x) == len(model.clusters) - summation(model.z)

    model.Radiality = Constraint(
        rule=Radiality_rule)  # all the clusters are either connected to the MV grid or powered by microgrid

    def Power_flow_conservation_rule(model, node):
        return (sum(model.P[j, node] for j in model.NodesIn[node]) - sum(
            model.P[node, j]
            for j in model.NodesOut[node])) == +model.p_clusters[node] * (
                           1 - model.z[node])

    model.Power_flow_conservation = Constraint(model.clusters,
                                               rule=Power_flow_conservation_rule)  # when the node is powered by SHS all the power is transferred to the outgoing link

    def PS_Power_flow_conservation_rule(model, node):
        return (sum(model.P[j, node] for j in model.NodesIn[node]) - sum(
            model.P[node, j]
            for j in model.NodesOut[node])) == -model.p_substations[node]

    model.PS_Power_flow_conservation = Constraint(model.substations,
                                                  rule=PS_Power_flow_conservation_rule)  # outgoing power from PS to connected links

    def Power_upper_bounds_rule(model, i, j):
        return model.P[i, j] <= p_max_lines * model.x[i, j]

    model.upper_Power_limits = Constraint(model.links,
                                          rule=Power_upper_bounds_rule)  # limit to power flowing on MV lines

    def Power_lower_bounds_rule(model, i, j):
        return model.P[i, j] >= -p_max_lines * model.x[i, j]

    model.lower_Power_limits = Constraint(model.links,
                                          rule=Power_lower_bounds_rule)

    def Primary_substation_upper_bound_rule(model, i):
        return model.p_substations[i] <= model.p_max_substations[i] * model.y[
            i]

    model.Primary_substation_upper_bound = Constraint(model.substations,
                                                      rule=Primary_substation_upper_bound_rule)  # limit to power of PS

    def Number_substations_rule(model,node):
        return model.y[node] <= sum(model.x[j, node] for j in model.NodesInSub[node]) + \
            sum(model.x[node,j] for j in model.NodesOutSub[node])
    model.Number_substation = Constraint(model.substations,
                                                     rule=Number_substations_rule)

    #### Distance constraints #####
    def distance_balance_rule(model, i, j):
        return model.k[i, j] == (-model.d_links[i, j] -model.d_nodes[i])*model.x[i, j]

    model.distance_balance_rule = Constraint(model.links,
                                            rule=distance_balance_rule)

    def distance_linearization_rule_1(model, i, j):
        return model.k[i, j] <= model.M_max * model.x[i, j]

    model.distance_linearization_rule_1 = Constraint(model.links,
                                                    rule=distance_linearization_rule_1)
    #
    def distance_linearization_rule_2(model, i, j):
        return model.k[i, j] >= model.M_min * model.x[i, j]

    model.distance_linearization_rule_2 = Constraint(model.links,
                                                    rule=distance_linearization_rule_2)
    #
    def distance_linearization_rule_3(model, i, j):
        return model.dist[i] - model.dist[j] - (1 - model.x[i, j]) * model.M_max <= \
               model.k[i, j]

    model.distance_linearization_rule_3 = Constraint(model.links,
                                                    rule=distance_linearization_rule_3)

    def distance_linearization_rule_4(model, i, j):
        return model.dist[i] - model.dist[j] - (1 - model.x[i, j]) * model.M_min >= \
               model.k[i, j]

    model.distance_linearization_rule_4 = Constraint(model.links,
                                                    rule=distance_linearization_rule_4)

    def distance_linearization_rule_5(model, i, j):
        return model.dist[i] - model.dist[j] + (1 - model.x[i, j]) * model.M_max >= \
               model.k[i, j]

    model.distance_linearization_rule_5 = Constraint(model.links,
                                                    rule=distance_linearization_rule_5)

    # if a cluster is electrified with a microgrid, its distance is 0,
    # otherwise it must be less than a max treshold
    def distance_upper_bound_rule(model, i):
        return 50 *(1-model.z[i])>= model.dist[i]

    model.distance_upper_bound = Constraint(model.clusters,
                                           rule=distance_upper_bound_rule)

    def distance_primary_substation_rule(model, i):
        return model.dist[i] == 0

    model.distance_primary_substation = Constraint(model.substations,
                                                  rule=distance_primary_substation_rule)

    # define loss of load  dependent on distance from connection point
    def lol_calculation_rule(model,i):
        return model.lol_line[i] == model.dist[i]*line_rel*model.energy[i]/8760/proj_lifetime

    model.lol_calculation = Constraint(model.clusters, rule =lol_calculation_rule)


    ####################Define objective function##########################

    # total npc over microgrid lifetime
    def ObjectiveFunctionCost(model):
        return model.obj['cost'] == summation(model.c_microgrids, model.z) \
                + summation(model.c_substations, model.y) + summation(model.c_links, model.x) \
                + sum(model.energy[i] * (1-model.z[i]) for i in model.clusters) * coe

    model.Obj1 = Constraint(rule=ObjectiveFunctionCost)

    # total direct emissions over microgrid lifetime
    def ObjectiveFunctionEmis(model):
        return model.obj['emis'] == summation(model.emission, model.z) +\
               summation(model.em_links, model.x)+\
               sum(model.energy[i] * (1-model.z[i]) for i in model.clusters) * nation_emis

    model.Obj2 = Constraint(rule=ObjectiveFunctionEmis)

    # total energy not supplied [MWh]
    def ObjectiveFunctionRel(model):
        return model.obj['rel'] == \
               sum(model.rel_mg[i]* (model.z[i])*100 for i in model.clusters)  +\
               summation(model.lol_line)+\
               sum(model.energy[i] /8760/proj_lifetime* (1-model.z[i]) for i in model.clusters) * nation_rel

    model.Obj3 = Constraint(rule=ObjectiveFunctionRel)


    # auxiliary variable to allow the activation and deactivation of OF in the for loop
    def AuxiliaryNorm(model, of):
        return model.aux[of] == model.dir[of] * model.obj[of]

    model.AuxNorm = Constraint(model.of, rule=AuxiliaryNorm)

    # aux is null for the OF not optimized in the loop
    def NullAuxiliary(model, of):
        return model.aux[of] == 0

    model.NullAux = Constraint(model.of, rule=NullAuxiliary)

    # objective function for identification of ranges of objective functions (needed for normalization)
    def ObjectiveFunctionNorm(model):
        return sum(model.aux[of] for of in model.of)

    model.ObjNorm = Objective(rule=ObjectiveFunctionNorm, sense=maximize)

    # normalized objective functions
    def DefineNormalizedObj(model, of):
        if model.dir[of] == 1:
            return model.norm_obj[of] == (model.obj[of] - model.min_obj[of]) / (model.max_obj[of] - model.min_obj[of])
        else:
            return model.norm_obj[of] == (model.max_obj[of] - model.obj[of]) / (model.max_obj[of] - model.min_obj[of])

    model.DefNormObj = Constraint(model.of, rule=DefineNormalizedObj)

    # multi-objective optimization through weighted sum approach
    def MultiObjective(model):
        return summation(model.weight, model.norm_obj)

    model.MultiObj = Objective(rule=MultiObjective, sense=maximize)

    #############Solve model##################

    # opt = SolverFactory('cplex',executable=r'C:\Users\silvi\IBM\ILOG\CPLEX_Studio1210\cplex\bin\x64_win64\cplex')
    # opt = SolverFactory('glpk')
    opt = SolverFactory('gurobi')
    opt.options['mipgap'] = 0.01

    instance = model.create_instance(data)
    print('Instance is constructed:', instance.is_constructed())

    obj_list = list(instance.of)  # list of the objective functions
    print(obj_list)
    num_of = len(obj_list)  # number of objective functions
    # payoff_table = np.empty((num_of,num_of))  # table of the ranges of variations of objective functions
    payoff_table = pd.DataFrame(index=obj_list,columns=obj_list)
    payoff_table.index.name = 'optimization'

    # for the first step, ObjNorm is the OF to be used
    instance.MultiObj.deactivate()
    instance.ObjNorm.activate()
    instance.DefNormObj.deactivate()

    print('1) Optimizing one objective function at a time to identify ranges of variations')

    '''
    Options for step 1 (identify ranges for normalizaiton):

    (1) variabile ausiliaria ciao[of] per poter attivare e disattivare i constraints di interesse nel loop
    constraint[of] ciao[of]=dir[of]*obj[of]
    OF sum,of ciao[of]
    for of in of:
        deactivate constraint
        activate constraint in of
        solve

    (2) ciclare sull'indice per definire cosa diventa OF nelle varie iterazioni del loop
    for of in Nof:
        model.ObjNorm = Objective(rule=ObjectiveFunction%of%, sense=maximize)
        (definire le atre come constraint)
        solve
        deactivate objective

    (3) versione esplicitata di (2)
    for of in Nof:
        if of=1
        model.ObjNorm = Objective(rule=ObjectiveFunction1, sense=maximize)
        solve
        deactivate objective

    (4) mettere le OF in un array e indicare con il suo indice quale ottimizzare
    for of in Nof:
        model.ObjNorm = Objective(rule=ObjectiveFunction[1], sense=maximize) # if array of rules of OF
        (definire gli altri elementi dell'array come constraints)
        solve
        deactivate objective
    
    (2),(3),(4) fattibili?
    posso creare l'instance e poi definire cosa è OF e cosa constraint? forse posso solo attivare e disattivare

    '''

    time_i = datetime.now()

    for of in obj_list:
    #for of in instance.of:
        print('Optimize '+of)
        instance.NullAux.activate()
        instance.NullAux[of].deactivate()
        instance.AuxNorm.deactivate()
        instance.AuxNorm[of].activate()
        opt.solve(instance, tee=True)
        payoff_of = []
        for i in obj_list:
            p_of = float(instance.obj.get_values()[i])
            payoff_of.append(p_of)
        payoff_table.loc[of,:] = payoff_of

    print(payoff_table)

    multi_obj=True
    k=0
    print('Find ranges of variation of each objective function:')
    for of in obj_list:
        instance.min_obj[of] = min(payoff_table[of])
        instance.max_obj[of] = max(payoff_table[of])
        print('min' + str(of) + '=' + str(min(payoff_table[of])))
        print('max' + str(of) + '=' + str(max(payoff_table[of])))
        # do not make multiobjective optimization if there is a unique solution
        #that means if all objective functions do not change
        if instance.min_obj[of] == instance.max_obj[of]:
            k=k+1
    if k==num_of:
        multi_obj =False
        print('Multi-obj not needed')


    # for the second step, MultiObj is the OF to be used
    instance.NullAux.deactivate()
    instance.AuxNorm.deactivate()
    instance.ObjNorm.deactivate()
    instance.MultiObj.activate()
    instance.DefNormObj.activate()

    if multi_obj:
        print('2) Multi-objective optimization: Weighted sum approach')
        opt.solve(instance, tee=True)

        for of in obj_list:
            print(str(of)+'='+str(instance.obj.get_values()[of]))

    time_f = datetime.now()
    print('Time required for the two steps is', time_f - time_i)

    ###################Process results#######################
    links = instance.x
    power = instance.P
    microgrids = instance.z
    distance =instance.dist
    lol_line=instance.lol_line
    connections_output = pd.DataFrame(columns=[['id1', 'id2']])
    microgrids_output = pd.DataFrame(columns=['ID'])
    power_output = pd.DataFrame(columns=[['id1', 'id2', 'P']])
    dist_output = pd.DataFrame(columns=[['ID','dist','lol']])
    k = 0
    for index in links:
        if int(round(value(links[index]))) == 1:
            connections_output.loc[k, 'id1'] = index[0]
            connections_output.loc[k, 'id2'] = index[1]
            k = k + 1
    k = 0
    for index in microgrids:
        if int(round(value(microgrids[index]))) == 1:
            microgrids_output.loc[k, 'ID'] = index
            k = k + 1
    k = 0
    for index in power:
        if value(power[index]) != 0:
            power_output.loc[k, 'id1'] = index[0]
            power_output.loc[k, 'id2'] = index[1]
            power_output.loc[k, 'P'] = value(power[index])
            k = k + 1

    k = 0
    for index in distance:
        if value(distance[index]) != 0:
            dist_output.loc[k, 'ID'] = index
            dist_output.loc[k, 'Length'] = value(distance[index])
            dist_output.loc[k, 'lol'] = value(lol_line[index])
            k = k + 1



    connections_output.to_csv('Output/LCOE/MV_connections_output.csv', index=False)
    microgrids_output.to_csv('Output/LCOE/MV_SHS_output.csv', index=False)
    power_output.to_csv('Output/LCOE/MV_power_output.csv', index=False)
    dist_output.to_csv('Output/LCOE/MV_dist_output.csv', index=False)
    return microgrids_output, connections_output
