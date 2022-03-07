''' ALL DIFFERENT MILP MODELS IN ONE FILE
1. MILP WITHOUT MG -> Doesn't consider the Microgrid option and works with only 1 type of cable. However, it considers reliability.
2. MILP2 -> Consider
3. MILP3
4. MILP4
5. MILP5
'''

from __future__ import division

from pyomo.opt import SolverFactory
from pyomo.core import AbstractModel
from pyomo.dataportal.DataPortal import DataPortal
from pyomo.environ import *
import pandas as pd
from datetime import datetime
import os
def MILP_without_MG(gisele_folder,case_study,n_clusters,coe,voltage,resistance,reactance,Pmax,line_cost):
    ############ Create abstract model ###########
    model = AbstractModel()
    data = DataPortal()
    MILP_input_folder = gisele_folder + '/Case studies/' + case_study + '/Intermediate/Optimization/MILP_input'
    MILP_output_folder = gisele_folder + '/Case studies/' + case_study + '/Intermediate/Optimization/MILP_output'
    os.chdir(MILP_input_folder)

    # Define some basic parameter for the per unit conversion and voltage limitation
    Abase = 1
    Vmin = 0.9

    ####################Define sets#####################
    model.N = Set()
    data.load(filename='nodes.csv', set=model.N)  # first row is not read
    model.N_clusters = Set()
    data.load(filename='nodes_clusters.csv', set=model.N_clusters)
    model.N_PS = Set()
    data.load(filename='nodes_PS.csv', set=model.N_PS)
    # Node corresponding to primary substation

    # Allowed connections
    model.links = Set(dimen=2)  # in the csv the values must be delimited by commas
    data.load(filename='links_all.csv', set=model.links)

    model.links_clusters = Set(dimen=2)
    data.load(filename='links_clusters.csv', set=model.links_clusters)

    model.links_decision = Set(dimen=2)
    data.load(filename='links_decision.csv', set=model.links_decision)

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

    #####################Define parameters#####################

    # Electric power in the nodes (injected (-) or absorbed (+))
    model.Psub = Param(model.N_clusters)
    data.load(filename='power_nodes.csv', param=model.Psub)

    model.ps_cost = Param(model.N_PS)
    data.load(filename='PS_costs.csv', param=model.ps_cost)

    model.PSmax = Param(model.N_PS)
    data.load(filename='PS_power_max.csv', param=model.PSmax)

    model.PS_voltage = Param(model.N_PS)
    data.load(filename='PS_voltage.csv', param=model.PS_voltage)

    model.PS_distance = Param(model.N_PS)
    data.load(filename='PS_distance.csv', param=model.PS_distance)

    # Connection distance of all the edges
    model.dist = Param(model.links)
    data.load(filename='distances.csv', param=model.dist)

    model.weights = Param(model.links_decision)
    data.load(filename='weights_decision_lines.csv', param=model.weights)
    # Electrical parameters of all the cables
    model.V_ref = Param(initialize=voltage)
    model.A_ref = Param(initialize=Abase)
    model.E_min = Param(initialize=Vmin)

    model.R_ref = Param(initialize=resistance)
    model.X_ref = Param(initialize=reactance)
    model.P_max = Param(initialize=Pmax)
    model.cf = Param(initialize=line_cost)

    model.Z = Param(initialize=model.R_ref + model.X_ref * 0.5)
    model.Z_ref = Param(initialize=model.V_ref ** 2 / Abase)

    model.n_clusters = Param(initialize=n_clusters)
    model.coe = Param(initialize=coe)

    #####################Define variables#####################

    # binary variable x[i,j]: 1 if the connection i,j is present, 0 otherwise
    model.x = Var(model.links_decision, within=Binary)
    # power[i,j] is the power flow of connection i-j
    model.P = Var(model.links)
    # positive variables E(i) is p.u. voltage at each node
    model.E = Var(model.N, within=NonNegativeReals)

    # binary variable k[i]: 1 if node i is a primary substation, 0 otherwise
    model.k = Var(model.N_PS, within=Binary)
    # Power output of Primary substation
    model.PPS = Var(model.N_PS, within=NonNegativeReals)
    model.positive_p = Var(model.links_clusters, within=Binary)
    model.Distance = Var(model.N, within=Reals)

    model.cable_type = Var(model.links)

    #####################Define constraints###############################

    # Radiality constraint
    def Radiality_rule(model):
        return summation(model.x) == model.n_clusters

    model.Radiality = Constraint(rule=Radiality_rule)

    # Power flow constraints
    def Power_flow_conservation_rule(model, node):
        return (sum(model.P[j, node] for j in model.NodesIn[node]) - sum(
            model.P[node, j] for j in model.NodesOut[node])) == model.Psub[node]

    model.Power_flow_conservation = Constraint(model.N_clusters, rule=Power_flow_conservation_rule)

    def Power_flow_conservation_rule3(model, node):
        return (sum(model.P[j, node] for j in model.NodesIn[node]) - sum(
            model.P[node, j] for j in model.NodesOut[node])) == - model.PPS[node]

    model.Power_flow_conservation3 = Constraint(model.N_PS, rule=Power_flow_conservation_rule3)

    def Power_upper_decision(model, i, j):
        return model.P[i, j] <= model.P_max * model.x[i, j]

    model.Power_upper_decision = Constraint(model.links_decision, rule=Power_upper_decision)

    def Power_lower_decision(model, i, j):
        return model.P[i, j] >= 0  # -model.P_max*model.x[i,j]

    model.Power_lower_decision = Constraint(model.links_decision, rule=Power_lower_decision)

    def Power_upper_clusters(model, i, j):
        return model.P[i, j] <= model.P_max * model.positive_p[i, j]

    model.Power_upper_clusters = Constraint(model.links_clusters, rule=Power_upper_clusters)

    def Power_lower_clusters(model, i, j):
        return model.P[i, j] >= 0  # -model.P_max*model.x[i,j]

    model.Power_lower_clusters = Constraint(model.links_clusters, rule=Power_lower_clusters)

    # Voltage constraints
    def Voltage_balance_rule(model, i, j):
        return (model.E[i] - model.E[j]) + model.x[i, j] - 1 <= model.dist[i, j] / 1000 * model.P[
            i, j] * model.Z / model.Z_ref

    model.Voltage_balance_rule = Constraint(model.links_decision, rule=Voltage_balance_rule)

    def Voltage_balance_rule2(model, i, j):
        return (model.E[i] - model.E[j]) - model.x[i, j] + 1 >= model.dist[i, j] / 1000 * model.P[
            i, j] * model.Z / model.Z_ref

    model.Voltage_balance_rule2 = Constraint(model.links_decision, rule=Voltage_balance_rule2)

    def Voltage_balance_rule3(model, i, j):
        return (model.E[i] - model.E[j]) + model.positive_p[i, j] - 1 <= model.dist[i, j] / 1000 * model.P[
            i, j] * model.Z / model.Z_ref

    model.Voltage_balance_rule3 = Constraint(model.links_clusters, rule=Voltage_balance_rule3)

    def Voltage_balance_rule4(model, i, j):
        return (model.E[i] - model.E[j]) - model.positive_p[i, j] + 1 >= model.dist[i, j] / 1000 * model.P[
            i, j] * model.Z / model.Z_ref

    model.Voltage_balance_rule4 = Constraint(model.links_clusters, rule=Voltage_balance_rule4)

    def Voltage_limit(model, i):
        return model.E[i] >= model.k[i] * (model.PS_voltage[i] - model.E_min) + model.E_min

    model.Voltage_limit = Constraint(model.N_PS, rule=Voltage_limit)

    def Voltage_PS2(model, i):
        return model.E[i] <= model.PS_voltage[i]

    model.Voltage_PS2 = Constraint(model.N_PS, rule=Voltage_PS2)

    def Voltage_limit_clusters2(model, i):
        return model.E[i] >= model.E_min

    model.Voltage_limit_clusters2 = Constraint(model.N_clusters, rule=Voltage_limit_clusters2)

    def PS_power_rule_upper(model, i):
        return model.PPS[i] <= model.PSmax[i] * model.k[i]

    model.PS_power_upper = Constraint(model.N_PS, rule=PS_power_rule_upper)

    def distance_from_PS(model, i):
        return model.Distance[i] <= -model.PS_distance[i]

    model.distance_from_PS = Constraint(model.N_PS, rule=distance_from_PS)

    def distance_from_PS2(model, i):
        return model.Distance[i] >= (model.k[i] - 1) * 200 - model.PS_distance[i] * model.k[i]

    model.distance_from_PS2 = Constraint(model.N_PS, rule=distance_from_PS2)

    def distance_balance_decision(model, i, j):
        return model.Distance[i] - model.Distance[j] + 1000 * (model.x[i, j] - 1) <= model.dist[i, j] / 1000

    model.distance_balance_decision = Constraint(model.links_decision, rule=distance_balance_decision)

    def distance_balance_decision2(model, i, j):
        return (model.Distance[i] - model.Distance[j]) - 1000 * (model.x[i, j] - 1) >= model.dist[i, j] / 1000

    model.distance_balance_decision2 = Constraint(model.links_decision, rule=distance_balance_decision2)

    def distance_balance_clusters(model, i, j):
        return model.Distance[i] - model.Distance[j] + 1000 * (model.positive_p[i, j] - 1) <= model.dist[i, j] / 1000

    model.distance_balance_clusters = Constraint(model.links_clusters, rule=distance_balance_clusters)

    def distance_balance_clusters2(model, i, j):
        return (model.Distance[i] - model.Distance[j]) - 1000 * (model.positive_p[i, j] - 1) >= model.dist[i, j] / 1000

    model.distance_balance_clusters2 = Constraint(model.links_clusters, rule=distance_balance_clusters2)

    def Balance_rule(model):
        return (sum(model.PPS[i] for i in model.N_PS) - sum(model.Psub[i] for i in model.N_clusters)) == 0

    model.Balance = Constraint(rule=Balance_rule)

    def anti_paralel(model, i, j):
        return model.x[i, j] + model.x[j, i] <= 1

    model.anti_paralel = Constraint(model.links_decision, rule=anti_paralel)

    def anti_paralel_clusters(model, i, j):
        return model.positive_p[i, j] + model.positive_p[j, i] == 1

    model.anti_paralel_clusters = Constraint(model.links_clusters, rule=anti_paralel_clusters)
    ####################Define objective function##########################

    ####################Define objective function##########################
    reliability_index = 1000

    def ObjectiveFunction(model):
        return summation(model.weights, model.x) * model.cf / 1000 + summation(model.ps_cost, model.k)

        # return summation(model.weights, model.x) * model.cf / 1000  + summation(model.ps_cost,model.k) - sum(model.Psub[i]*model.Distance[i] for i in model.N_clusters)*reliability_index

    model.Obj = Objective(rule=ObjectiveFunction, sense=minimize)

    #############Solve model##################

    instance = model.create_instance(data)
    print('Instance is constructed:', instance.is_constructed())
    # opt = SolverFactory('cbc',executable=r'C:\Users\Asus\Desktop\POLIMI\Thesis\GISELE\Gisele_MILP\cbc')
    opt = SolverFactory('gurobi')
    opt.options['TimeLimit'] = 300
    # opt.options['numericfocus']=0
    # opt.options['mipgap'] = 0.0002
    # opt.options['presolve']=2
    # opt.options['mipfocus']=2
    # opt = SolverFactory('cbc',executable=r'C:\Users\Asus\Desktop\POLIMI\Thesis\GISELE\New folder\cbc')
    print('Starting optimization process')
    time_i = datetime.now()
    opt.solve(instance, tee=True, symbolic_solver_labels=True)
    time_f = datetime.now()
    print('Time required for optimization is', time_f - time_i)
    links = instance.x
    power = instance.P
    subs = instance.k
    voltage = instance.E
    PS = instance.PPS
    DISTANCE = instance.Distance

    links_clusters = instance.links_clusters
    # voltage_drop=instance.z
    connections_output = pd.DataFrame(columns=[['id1', 'id2', 'power']])
    PrSubstation = pd.DataFrame(columns=[['index', 'power']])
    all_lines = pd.DataFrame(columns=[['id1', 'id2', 'power']])
    Voltages = pd.DataFrame(columns=[['index', 'voltage [p.u]']])
    Links_Clusters = pd.DataFrame(columns=[['id1', 'id2', 'power']])
    distance = pd.DataFrame(columns=[['index', 'length[km]']])
    k = 0
    for index in links:
        if int(round(value(links[index]))) == 1:
            connections_output.loc[k, 'id1'] = index[0]
            connections_output.loc[k, 'id2'] = index[1]
            connections_output.loc[k, 'power'] = value(power[index])
            k = k + 1
    k = 0
    for index in subs:
        if int(round(value(subs[index]))) == 1:
            PrSubstation.loc[k, 'index'] = index
            PrSubstation.loc[k, 'power'] = value(PS[index])
            print(((value(PS[index]))))
            k = k + 1
    k = 0
    for v in voltage:
        Voltages.loc[k, 'index'] = v
        Voltages.loc[k, 'voltage [p.u]'] = value(voltage[v])
        k = k + 1
    k = 0
    for index in power:
        all_lines.loc[k, 'id1'] = index[0]
        all_lines.loc[k, 'id2'] = index[1]
        all_lines.loc[k, 'power'] = value(power[index])
        k = k + 1
    k = 0

    for index in links_clusters:
        Links_Clusters.loc[k, 'id1'] = index[0]
        Links_Clusters.loc[k, 'id2'] = index[1]
        Links_Clusters.loc[k, 'power'] = value(power[index])
        k = k + 1
    k = 0
    for dist in DISTANCE:
        distance.loc[k, 'index'] = dist
        distance.loc[k, 'length[m]'] = value(DISTANCE[dist])
        k = k + 1
    Links_Clusters.to_csv(MILP_output_folder + '/links_clusters.csv', index=False)
    connections_output.to_csv(MILP_output_folder + '/connections_output.csv', index=False)
    PrSubstation.to_csv(MILP_output_folder + '/PrimarySubstations.csv', index=False)
    Voltages.to_csv(MILP_output_folder + '/Voltages.csv', index=False)
    all_lines.to_csv(MILP_output_folder + '/all_lines.csv', index=False)
    distance.to_csv(MILP_output_folder + '/Distances.csv', index=False)

def MILP_MG_reliability(gisele_folder,case_study,n_clusters,coe,voltage,resistance,reactance,Pmax,line_cost):

    model = AbstractModel()
    data = DataPortal()
    MILP_input_folder = gisele_folder + '/Case studies/' + case_study + '/Intermediate/Optimization/MILP_input'
    MILP_output_folder = gisele_folder + '/Case studies/' + case_study + '/Intermediate/Optimization/MILP_output'
    os.chdir(MILP_input_folder)

    # Define some basic parameter for the per unit conversion and voltage limitation
    Abase = 1
    Vmin = 0.9
    # ####################Define sets#####################

    # Name of all the nodes (primary and secondary substations)
    model.N = Set()
    data.load(filename='nodes.csv', set=model.N)  # first row is not read
    model.N_clusters = Set()
    data.load(filename='nodes_clusters.csv', set=model.N_clusters)
    model.N_MG = Set()
    data.load(filename='microgrids_nodes.csv', set=model.N_MG)
    model.N_PS = Set()
    data.load(filename='nodes_PS.csv', set=model.N_PS)
    # Node corresponding to primary substation

    # Allowed connections
    model.links = Set(dimen=2)  # in the csv the values must be delimited by commas
    data.load(filename='links_all.csv', set=model.links)

    model.links_clusters = Set(dimen=2)
    data.load(filename='links_clusters.csv', set=model.links_clusters)

    model.links_decision = Set(dimen=2)
    data.load(filename='links_decision.csv', set=model.links_decision)

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

    #####################Define parameters#####################

    # Electric power in the nodes (injected (-) or absorbed (+))
    model.Psub = Param(model.N_clusters)
    data.load(filename='power_nodes.csv', param=model.Psub)

    # model.PS=Param(model.N)
    # data.load(filename='PS.csv',param=model.PS)

    model.microgrid_power = Param(model.N_MG)
    data.load(filename='microgrids_powers.csv', param=model.microgrid_power)

    model.energy = Param(model.N_MG)
    data.load(filename='energy.csv', param=model.energy)

    model.mg_cost = Param(model.N_MG)
    data.load(filename='microgrids_costs.csv', param=model.mg_cost)

    model.ps_cost = Param(model.N_PS)
    data.load(filename='PS_costs.csv', param=model.ps_cost)

    model.PSmax = Param(model.N_PS)
    data.load(filename='PS_power_max.csv', param=model.PSmax)

    # Power of the primary substation as sum of all the other powers
    # def PPS_init(model):
    #    return sum(model.Psub[i] for i in model.N)
    # model.PPS=Param(model.PS,initialize=PPS_init)

    # Connection distance of all the edges
    model.dist = Param(model.links)
    data.load(filename='distances.csv', param=model.dist)

    model.weights = Param(model.links_decision)
    data.load(filename='weights_decision_lines.csv', param=model.weights)
    # Electrical parameters of all the cables
    model.V_ref = Param()
    model.A_ref = Param()
    model.R_ref = Param()
    model.X_ref = Param()
    model.P_max = Param()
    model.cf = Param()
    model.cPS = Param()
    model.E_min = Param()
    model.PPS_max = Param()
    model.PPS_min = Param()

    model.Z = Param()
    model.Z_ref = Param()

    model.n_clusters = Param()
    model.coe = Param()
    data.load(filename='data2.dat')

    #####################Define variables#####################

    # binary variable x[i,j]: 1 if the connection i,j is present, 0 otherwise
    model.x = Var(model.links_decision, within=Binary)
    # power[i,j] is the power flow of connection i-j
    model.P = Var(model.links)
    # positive variables E(i) is p.u. voltage at each node
    model.E = Var(model.N, within=NonNegativeReals)
    # microgrid
    model.z = Var(model.N_MG, within=Binary)
    # binary variable k[i]: 1 if node i is a primary substation, 0 otherwise
    model.k = Var(model.N_PS, within=Binary)
    # Power output of Primary substation
    model.PPS = Var(model.N_PS, within=NonNegativeReals)
    model.MG_output = Var(model.N_MG)
    model.positive_p = Var(model.links_clusters, within=Binary)
    model.Distance = Var(model.N, within=Reals)

    #####################Define constraints###############################
    # def Make_problem_easy(model,i,j):
    #    return model.x[i,j]+model.weights[i,j]>=1
    # model.easy = Constraint(model.links, rule=Make_problem_easy)

    def Radiality_rule(model):
        # return summation(model.x)==len(model.N)-summation(model.k)
        return summation(model.x) == model.n_clusters

    model.Radiality = Constraint(rule=Radiality_rule)

    def Radiality_rule(model):
        return summation(model.k) + summation(model.z) <= model.n_clusters

    def Power_flow_conservation_rule(model, node):
        return (sum(model.P[j, node] for j in model.NodesIn[node]) - sum(
            model.P[node, j] for j in model.NodesOut[node])) == model.Psub[node]

    model.Power_flow_conservation = Constraint(model.N_clusters, rule=Power_flow_conservation_rule)

    def Power_flow_conservation_rule2(model, node):
        return (sum(model.P[j, node] for j in model.NodesIn[node]) - sum(
            model.P[node, j] for j in model.NodesOut[node])) == - model.MG_output[node]

    model.Power_flow_conservation2 = Constraint(model.N_MG, rule=Power_flow_conservation_rule2)

    def Power_flow_conservation_rule3(model, node):
        return (sum(model.P[j, node] for j in model.NodesIn[node]) - sum(
            model.P[node, j] for j in model.NodesOut[node])) == - model.PPS[node]

    model.Power_flow_conservation3 = Constraint(model.N_PS, rule=Power_flow_conservation_rule3)

    def Power_upper_decision(model, i, j):
        return model.P[i, j] <= model.P_max * model.x[i, j]

    model.Power_upper_decision = Constraint(model.links_decision, rule=Power_upper_decision)

    def Power_lower_decision(model, i, j):
        return model.P[i, j] >= 0  # -model.P_max*model.x[i,j]

    model.Power_lower_decision = Constraint(model.links_decision, rule=Power_lower_decision)

    def Power_upper_clusters(model, i, j):
        return model.P[i, j] <= model.P_max * model.positive_p[i, j]

    model.Power_upper_clusters = Constraint(model.links_clusters, rule=Power_upper_clusters)

    def Power_lower_clusters(model, i, j):
        return model.P[i, j] >= 0  # -model.P_max*model.x[i,j]

    model.Power_lower_clusters = Constraint(model.links_clusters, rule=Power_lower_clusters)

    # Voltage constraints
    def Voltage_balance_rule(model, i, j):
        return (model.E[i] - model.E[j]) + model.x[i, j] - 1 <= model.dist[i, j] / 1000 * model.P[
            i, j] * model.Z / model.Z_ref

    model.Voltage_balance_rule = Constraint(model.links_decision, rule=Voltage_balance_rule)

    def Voltage_balance_rule2(model, i, j):
        return (model.E[i] - model.E[j]) - model.x[i, j] + 1 >= model.dist[i, j] / 1000 * model.P[
            i, j] * model.Z / model.Z_ref

    model.Voltage_balance_rule2 = Constraint(model.links_decision, rule=Voltage_balance_rule2)

    def Voltage_balance_rule3(model, i, j):
        return (model.E[i] - model.E[j]) + model.positive_p[i, j] - 1 <= model.dist[i, j] / 1000 * model.P[
            i, j] * model.Z / model.Z_ref

    model.Voltage_balance_rule3 = Constraint(model.links_clusters, rule=Voltage_balance_rule3)

    def Voltage_balance_rule4(model, i, j):
        return (model.E[i] - model.E[j]) - model.positive_p[i, j] + 1 >= model.dist[i, j] / 1000 * model.P[
            i, j] * model.Z / model.Z_ref

    model.Voltage_balance_rule4 = Constraint(model.links_clusters, rule=Voltage_balance_rule4)

    def Voltage_limit(model, i):
        return model.E[i] >= model.k[i] * (1 - model.E_min) + model.E_min

    model.Voltage_limit = Constraint(model.N_PS, rule=Voltage_limit)

    def Voltage_PS2(model, i):
        return model.E[i] <= 1

    model.Voltage_PS2 = Constraint(model.N_PS, rule=Voltage_PS2)

    def Voltage_limit_MG(model, i):
        return model.E[i] <= 1

    model.Voltage_limit_MG = Constraint(model.N_MG, rule=Voltage_limit_MG)

    def Voltage_limit_MG2(model, i):
        return model.E[i] >= model.z[i] * (1 - model.E_min) + model.E_min

    model.Voltage_limit_MG2 = Constraint(model.N_MG, rule=Voltage_limit_MG2)

    def Voltage_limit_clusters2(model, i):
        return model.E[i] >= model.E_min

    model.Voltage_limit_clusters2 = Constraint(model.N_clusters, rule=Voltage_limit_clusters2)

    def PS_power_rule_upper(model, i):
        return model.PPS[i] <= model.PSmax[i] * model.k[i]

    model.PS_power_upper = Constraint(model.N_PS, rule=PS_power_rule_upper)

    def distance_from_PS(model, i):
        return model.Distance[i] <= 0

    model.distance_from_PS = Constraint(model.N_PS, rule=distance_from_PS)

    def distance_from_PS2(model, i):
        return model.Distance[i] >= (model.k[i] - 1) * 100

    model.distance_from_PS2 = Constraint(model.N_PS, rule=distance_from_PS2)

    def distance_from_MG(model, i):
        return model.Distance[i] <= 0

    model.distance_from_MG = Constraint(model.N_MG, rule=distance_from_MG)

    def distance_from_MG2(model, i):
        return model.Distance[i] >= (model.z[i] - 1) * 100  # length must be <100km in this case

    model.distance_from_MG2 = Constraint(model.N_MG, rule=distance_from_MG2)

    def distance_balance_decision(model, i, j):
        return model.Distance[i] - model.Distance[j] + 1000 * (model.x[i, j] - 1) <= model.dist[i, j] / 1000

    model.distance_balance_decision = Constraint(model.links_decision, rule=distance_balance_decision)

    def distance_balance_decision2(model, i, j):
        return (model.Distance[i] - model.Distance[j]) - 1000 * (model.x[i, j] - 1) >= model.dist[i, j] / 1000

    model.distance_balance_decision2 = Constraint(model.links_decision, rule=distance_balance_decision2)

    def distance_balance_clusters(model, i, j):
        return model.Distance[i] - model.Distance[j] + 1000 * (model.positive_p[i, j] - 1) <= model.dist[i, j] / 1000

    model.distance_balance_clusters = Constraint(model.links_clusters, rule=distance_balance_clusters)

    def distance_balance_clusters2(model, i, j):
        return (model.Distance[i] - model.Distance[j]) - 1000 * (model.positive_p[i, j] - 1) >= model.dist[i, j] / 1000

    model.distance_balance_clusters2 = Constraint(model.links_clusters, rule=distance_balance_clusters2)

    def Balance_rule(model):
        return (sum(model.PPS[i] for i in model.N_PS) + sum(model.MG_output[i] for i in model.N_MG) - sum(
            model.Psub[i] for i in model.N_clusters)) == 0

    model.Balance = Constraint(rule=Balance_rule)

    def MG_power_limit(model, i):
        return model.MG_output[i] == model.z[i] * model.microgrid_power[i]

    model.MG_power_limit = Constraint(model.N_MG, rule=MG_power_limit)

    def anti_paralel(model, i, j):
        return model.x[i, j] + model.x[j, i] <= 1

    model.anti_paralel = Constraint(model.links_decision, rule=anti_paralel)

    def anti_paralel_clusters(model, i, j):
        return model.positive_p[i, j] + model.positive_p[j, i] == 1

    model.anti_paralel_clusters = Constraint(model.links_clusters, rule=anti_paralel_clusters)

    ####################Define objective function##########################

    def ObjectiveFunction(model):
        # return summation(model.weights, model.x) * model.cf / 1000
        return summation(model.weights, model.x) * model.cf / 1000 + summation(model.mg_cost, model.z) * 1000 + sum(
            model.energy[i] * (1 - model.z[i]) for i in model.N_MG) * model.coe
        #      + summation(model.ps_cost,model.k)
        # +sum((model.P[i]/model.A_ref)**2*0.5*1.25*model.R_ref/model.Z_ref*model.dist[i]/1000*24*365*20 for i in model.links)
        # return summation(model.dist,model.x)*model.cf/1000 + summation(model.k) *model.cPS

    model.Obj = Objective(rule=ObjectiveFunction, sense=minimize)

    #############Solve model##################

    instance = model.create_instance(data)
    print('Instance is constructed:', instance.is_constructed())
    # opt = SolverFactory('cbc',executable=r'C:\Users\Asus\Desktop\POLIMI\Thesis\GISELE\Gisele_MILP\cbc')
    opt = SolverFactory('gurobi')
    # opt.options['numericfocus']=0
    # opt.options['mipgap'] = 0.0002
    # opt.options['presolve']=2
    # opt.options['mipfocus']=2
    # opt = SolverFactory('cbc',executable=r'C:\Users\Asus\Desktop\POLIMI\Thesis\GISELE\New folder\cbc')
    print('Starting optimization process')
    time_i = datetime.now()
    opt.solve(instance, tee=True, symbolic_solver_labels=True)
    time_f = datetime.now()
    print('Time required for optimization is', time_f - time_i)
    links = instance.x
    power = instance.P
    subs = instance.k
    voltage = instance.E
    PS = instance.PPS
    mg_output = instance.MG_output
    microGrid = instance.z
    DISTANCE = instance.Distance
    links_clusters = instance.links_clusters
    # voltage_drop=instance.z
    connections_output = pd.DataFrame(columns=[['id1', 'id2', 'power']])
    PrSubstation = pd.DataFrame(columns=[['index', 'power']])
    all_lines = pd.DataFrame(columns=[['id1', 'id2', 'power']])
    Voltages = pd.DataFrame(columns=[['index', 'voltage [p.u]']])
    Microgrid = pd.DataFrame(columns=[['index', 'microgrid', 'power']])
    distance = pd.DataFrame(columns=[['index', 'length[km]']])
    Links_Clusters = pd.DataFrame(columns=[['id1', 'id2', 'power']])
    k = 0
    for index in links:
        if int(round(value(links[index]))) == 1:
            connections_output.loc[k, 'id1'] = index[0]
            connections_output.loc[k, 'id2'] = index[1]
            connections_output.loc[k, 'power'] = value(power[index])
            k = k + 1
    k = 0
    for index in subs:
        if int(round(value(subs[index]))) == 1:
            PrSubstation.loc[k, 'index'] = index
            PrSubstation.loc[k, 'power'] = value(PS[index])
            print(((value(PS[index]))))
            k = k + 1
    k = 0
    for v in voltage:
        Voltages.loc[k, 'index'] = v
        Voltages.loc[k, 'voltage [p.u]'] = value(voltage[v])
        k = k + 1
    k = 0
    for index in power:
        all_lines.loc[k, 'id1'] = index[0]
        all_lines.loc[k, 'id2'] = index[1]
        all_lines.loc[k, 'power'] = value(power[index])
        k = k + 1
    k = 0
    for index in mg_output:
        Microgrid.loc[k, 'index'] = index
        Microgrid.loc[k, 'microgrid'] = value(microGrid[index])
        Microgrid.loc[k, 'power'] = value(mg_output[index])
        k = k + 1
    k = 0
    for dist in DISTANCE:
        distance.loc[k, 'index'] = dist
        distance.loc[k, 'length[m]'] = value(DISTANCE[dist])
        k = k + 1

    k = 0
    for index in links_clusters:
        Links_Clusters.loc[k, 'id1'] = index[0]
        Links_Clusters.loc[k, 'id2'] = index[1]
        Links_Clusters.loc[k, 'power'] = value(power[index])
        k = k + 1
    k = 0
    for dist in DISTANCE:
        distance.loc[k, 'index'] = dist
        distance.loc[k, 'length[m]'] = value(DISTANCE[dist])
        k = k + 1
    Links_Clusters.to_csv(MILP_output_folder + '/links_clusters.csv', index=False)
    connections_output.to_csv(MILP_output_folder + '/connections_output.csv', index=False)
    PrSubstation.to_csv(MILP_output_folder + '/PrimarySubstations.csv', index=False)
    Voltages.to_csv(MILP_output_folder + '/Voltages.csv', index=False)
    all_lines.to_csv(MILP_output_folder + '/all_lines.csv', index=False)
    Microgrid.to_csv(MILP_output_folder + '/Microgrid.csv', index=False)
    distance.to_csv(MILP_output_folder + '/Distances.csv', index=False)

def MILP_multiobjective(p_max_lines, coe, nation_emis, nation_rel, line_rel,input_michele):

    # useful parameters from michele
    proj_lifetime = input_michele['num_years']
    nren = 3

    # Initialize model
    model = AbstractModel()
    data = DataPortal()

    # Define sets

    model.of = Set(initialize=['cost', 'emis', 'rel'])  # Set of objective functions
    model.N = Set()  # Set of all nodes, clusters and substations
    data.load(filename=r'Output/LCOE/set.csv', set=model.N)
    model.clusters = Set()  # Set of clusters
    data.load(filename='Output/LCOE/clusters.csv', set=model.clusters)
    model.renfr = RangeSet(0, nren - 1,1)  # set of microgrids with different ren fractions
    model.mg = Set(dimen=2, within=model.clusters * model.renfr)
    model.substations = Set()  # Set of substations
    data.load(filename='Output/LCOE/subs.csv', set=model.substations)
    model.links = Set(dimen=2,within=model.N * model.N)  # in the csv the values must be delimited by commas
    data.load(filename='Output/LCOE/possible_links_complete.csv',set=model.links)

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
    model.c_microgrids = Param(model.mg)
    # data.load(filename='Output/LCOE/c_npc.csv', param=model.c_microgrids)

    # Total net present cost of substations
    model.c_substations = Param(model.substations)
    data.load(filename='Output/LCOE/sub_npc.csv',
              param=model.c_substations)

    # Connection cost of the possible links
    model.c_links = Param(model.links)
    data.load(filename='Output/LCOE/cost_links_complete.csv',
              param=model.c_links)

    # Energy consumed by each cluster in microgrid lifetime
    model.energy = Param(model.mg)
    # data.load(filename='Output/LCOE/energy.csv', param=model.energy)

    # CO2 emission produced by each cluster in microgrid lifetime
    model.emission = Param(model.mg)
    # data.load(filename='Output/LCOE/emissions.csv', param=model.emission)

    # CO2 emission related to construction of power infrastructure
    model.em_links = Param(model.links)
    data.load(filename='Output/LCOE/em_links.csv', param=model.em_links)

    # lol due to microgrid components_failure
    model.rel_mg = Param(model.mg)
    # data.load(filename='Output/LCOE/mg_rel.csv', param=model.rel_mg)

    # Connection length associated to the nodes
    # todo->put the real length
    model.d_nodes = Param(model.N)
    data.load(filename='Output/LCOE/len_nodes.csv', param=model.d_nodes)

    # Connection length of the possible links
    model.d_links = Param(model.links)
    data.load(filename='Output/LCOE/len_links_complete.csv',
              param=model.d_links)

    # poximum power flowing on lines
    # model.p_max_lines = Param()  # max power flowing on MV lines
    # data.load(filename='Input/data_procedure2.dat')

    # M_max and M_min, values required to linearize the problem
    model.M_max = Param(initialize=10000)
    model.M_min = Param(initialize=-10000)

    data.load(filename='Output/Microgrids/microgrids.csv',
              select=(
              'Cluster', 'Renewable fraction index', 'Total Cost [kEUR]',
              'Energy Demand [MWh]', 'CO2 [kg]', 'Unavailability [MWh/y]'),
              param=(model.c_microgrids, model.energy, model.emission,
                     model.rel_mg), index=model.mg)

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
    model.z = Var(model.mg, within=Binary)
    # power[i,j] is the power flow of connection i-j
    model.P = Var(model.links, within=NonNegativeReals)
    # power[i] is the power provided by substation i
    model.p_substations = Var(model.substations, within=NonNegativeReals)
    # # variables k(i,j) is the variable necessary to linearize
    model.k = Var(model.links)
    # distance of cluster from substation
    model.dist = Var(model.N, within=NonNegativeReals)
    # lol due to MV lines
    model.lol_line = Var(model.clusters, within=NonNegativeReals)

    #####################Define constraints###############################

    def Radiality_rule(model):
        return summation(model.x) == len(model.clusters) - summation(
            model.z)

    model.Radiality = Constraint(
        rule=Radiality_rule)  # all the clusters are either connected to the MV grid or powered by microgrid

    def Power_flow_conservation_rule(model, node):
        return (sum(model.P[j, node] for j in model.NodesIn[node]) - sum(
            model.P[node, j]
            for j in model.NodesOut[node])) == +model.p_clusters[node] * (
                       1 - sum(model.z[node, j] for j in model.renfr))

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
        return model.p_substations[i] <= model.p_max_substations[i] * \
               model.y[i]

    model.Primary_substation_upper_bound = Constraint(model.substations,
                                                      rule=Primary_substation_upper_bound_rule)  # limit to power of PS

    def Number_substations_rule(model, node):
        return model.y[node] <= sum(
            model.x[j, node] for j in model.NodesInSub[node]) + \
               sum(model.x[node, j] for j in model.NodesOutSub[node])

    model.Number_substation = Constraint(model.substations,
                                         rule=Number_substations_rule)

    def Limit_mg_rule(model, i):
        return sum(model.z[i, j] for j in model.renfr) <= 1

    model.Limit_mg = Constraint(model.clusters, rule=Limit_mg_rule)

    #### Distance constraints #####
    def distance_balance_rule(model, i, j):
        return model.k[i, j] == (-model.d_links[i, j] - model.d_nodes[i]) * \
               model.x[i, j]

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
        return model.dist[i] - model.dist[j] - (
                    1 - model.x[i, j]) * model.M_max <= \
               model.k[i, j]

    model.distance_linearization_rule_3 = Constraint(model.links,
                                                     rule=distance_linearization_rule_3)

    def distance_linearization_rule_4(model, i, j):
        return model.dist[i] - model.dist[j] - (
                    1 - model.x[i, j]) * model.M_min >= \
               model.k[i, j]

    model.distance_linearization_rule_4 = Constraint(model.links,
                                                     rule=distance_linearization_rule_4)

    def distance_linearization_rule_5(model, i, j):
        return model.dist[i] - model.dist[j] + (
                    1 - model.x[i, j]) * model.M_max >= \
               model.k[i, j]

    model.distance_linearization_rule_5 = Constraint(model.links,
                                                     rule=distance_linearization_rule_5)

    # if a cluster is electrified with a microgrid, its distance is 0,
    # otherwise it must be less than a max treshold
    def distance_upper_bound_rule(model, i, j):
        return 100 * (1 - model.z[i, j]) >= model.dist[i]

    model.distance_upper_bound = Constraint(model.clusters, model.renfr,
                                            rule=distance_upper_bound_rule)

    def distance_primary_substation_rule(model, i):
        return model.dist[i] == 0

    model.distance_primary_substation = Constraint(model.substations,
                                                   rule=distance_primary_substation_rule)

    # define loss of load  dependent on distance from connection point
    def lol_calculation_rule(model, i):
        return model.lol_line[i] == model.dist[i] * line_rel * \
               model.energy[i, 1] / 8760 / proj_lifetime

    model.lol_calculation = Constraint(model.clusters,
                                       rule=lol_calculation_rule)

    ####################Define objective function##########################

    # total npc over microgrid lifetime
    def ObjectiveFunctionCost(model):
        return model.obj['cost'] == summation(model.c_microgrids, model.z) \
               + summation(model.c_substations, model.y) + summation(
            model.c_links, model.x) \
               + sum(model.energy[i, 1] * (
                    1 - sum(model.z[i, j] for j in model.renfr)) for i in
                     model.clusters) * coe

    model.Obj1 = Constraint(rule=ObjectiveFunctionCost)

    # total direct emissions over microgrid lifetime
    def ObjectiveFunctionEmis(model):
        return model.obj['emis'] == summation(model.emission, model.z) + \
               summation(model.em_links, model.x) + \
               sum(model.energy[i, 1] * (
                           1 - sum(model.z[i, j] for j in model.renfr)) for
                   i in model.clusters) * nation_emis

    model.Obj2 = Constraint(rule=ObjectiveFunctionEmis)

    # minimization of total energy not supplied [MWh]
    def ObjectiveFunctionRel(model):
        return model.obj['rel'] == \
               sum(model.rel_mg[i] * (model.z[i]) for i in model.mg) + \
               summation(model.lol_line) + \
               sum(model.energy[i, 1] / 8760 / proj_lifetime * (
                           1 - sum(model.z[i, j] for j in model.renfr)) for
                   i in model.clusters) * nation_rel

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
            return model.norm_obj[of] == (
                        model.obj[of] - model.min_obj[of]) / (
                               model.max_obj[of] - model.min_obj[of])
        else:
            return model.norm_obj[of] == (
                        model.max_obj[of] - model.obj[of]) / (
                               model.max_obj[of] - model.min_obj[of])

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
    payoff_table = pd.DataFrame(index=obj_list, columns=obj_list)
    payoff_table.index.name = 'optimization'

    # for the first step, ObjNorm is the OF to be used
    instance.MultiObj.deactivate()
    instance.ObjNorm.activate()
    instance.DefNormObj.deactivate()

    print(
        '1) Optimizing one objective function at a time to identify ranges of variations')

    time_i = datetime.now()

    for of in obj_list:
        # for of in instance.of:
        print('Optimize ' + of)
        instance.NullAux.activate()
        instance.NullAux[of].deactivate()
        instance.AuxNorm.deactivate()
        instance.AuxNorm[of].activate()
        opt.solve(instance, tee=True)
        payoff_of = []
        for i in obj_list:
            p_of = float(instance.obj.get_values()[i])
            payoff_of.append(p_of)
        payoff_table.loc[of, :] = payoff_of

    print(payoff_table)

    multi_obj = True
    k = 0
    print('Find ranges of variation of each objective function:')
    for of in obj_list:
        instance.min_obj[of] = min(payoff_table[of])
        instance.max_obj[of] = max(payoff_table[of])
        print('min' + str(of) + '=' + str(min(payoff_table[of])))
        print('max' + str(of) + '=' + str(max(payoff_table[of])))
        # do not make multiobjective optimization if there is a unique solution
        # that means if all objective functions do not change
        if instance.min_obj[of] == instance.max_obj[of]:
            k = k + 1
    if k == num_of:
        multi_obj = False
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
            print(str(of) + '=' + str(instance.obj.get_values()[of]))

    time_f = datetime.now()
    print('Time required for the two steps is', time_f - time_i)

    ###################Process results#######################
    links = instance.x
    power = instance.P
    microgrids = instance.z
    distance = instance.dist
    lol_line = instance.lol_line
    connections_output = pd.DataFrame(columns=[['id1', 'id2']])
    microgrids_output = pd.DataFrame(
        columns=['Cluster', 'Renewable fraction'])
    power_output = pd.DataFrame(columns=[['id1', 'id2', 'P']])
    dist_output = pd.DataFrame(columns=[['ID', 'dist', 'lol']])
    k = 0
    for index in links:
        if int(round(value(links[index]))) == 1:
            connections_output.loc[k, 'id1'] = index[0]
            connections_output.loc[k, 'id2'] = index[1]
            k = k + 1
    k = 0
    for index in microgrids:
        if int(round(value(microgrids[index]))) == 1:
            microgrids_output.loc[k, 'Cluster'] = index[0]
            microgrids_output.loc[k, 'Renewable fraction'] = index[1]

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
            dist_output.loc[k, 'dist'] = value(distance[index])
            dist_output.loc[k, 'lol'] = value(lol_line[index])
            k = k + 1

    connections_output.to_csv('Output/LCOE/MV_connections_output.csv',
                              index=False)
    microgrids_output.to_csv('Output/LCOE/MV_SHS_output.csv', index=False)
    power_output.to_csv('Output/LCOE/MV_power_output.csv', index=False)
    dist_output.to_csv('Output/LCOE/MV_dist_output.csv', index=False)
    return microgrids_output, connections_output

def MILP_MG_noRel(gisele_folder,case_study,n_clusters,coe,voltage,resistance,reactance,Pmax,line_cost):
    model = AbstractModel()
    data = DataPortal()
    MILP_input_folder = gisele_folder + '/Case studies/' + case_study + '/Intermediate/Optimization/MILP_input'
    MILP_output_folder = gisele_folder + '/Case studies/' + case_study + '/Intermediate/Optimization/MILP_output'
    os.chdir(MILP_input_folder)
    # ####################Define sets#####################

    # Define some basic parameter for the per unit conversion and voltage limitation
    Abase = 1
    Vmin = 0.9

    # Name of all the nodes (primary and secondary substations)
    model.N = Set()
    data.load(filename='nodes.csv', set=model.N)  # first row is not read
    model.N_clusters = Set()
    data.load(filename='nodes_clusters.csv', set=model.N_clusters)
    model.N_MG = Set()
    data.load(filename='microgrids_nodes.csv', set=model.N_MG)
    model.N_PS = Set()
    data.load(filename='nodes_PS.csv', set=model.N_PS)
    # Node corresponding to primary substation

    # Allowed connections
    model.links = Set(dimen=2)  # in the csv the values must be delimited by commas
    data.load(filename='links_all.csv', set=model.links)

    model.links_clusters = Set(dimen=2)
    data.load(filename='links_clusters.csv', set=model.links_clusters)

    model.links_decision = Set(dimen=2)
    data.load(filename='links_decision.csv', set=model.links_decision)

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

    #####################Define parameters#####################

    # Electric power in the nodes (injected (-) or absorbed (+))
    model.Psub = Param(model.N_clusters)
    data.load(filename='power_nodes.csv', param=model.Psub)

    # model.PS=Param(model.N)
    # data.load(filename='PS.csv',param=model.PS)

    model.microgrid_power = Param(model.N_MG)
    data.load(filename='microgrids_powers.csv', param=model.microgrid_power)

    model.energy = Param(model.N_MG)
    data.load(filename='energy.csv', param=model.energy)
    # TODO also calculate npv
    model.mg_cost = Param(model.N_MG)
    data.load(filename='microgrids_costs.csv', param=model.mg_cost)
    # TODO also calculate npv
    model.ps_cost = Param(model.N_PS)
    data.load(filename='PS_costs.csv', param=model.ps_cost)

    model.PSmax = Param(model.N_PS)
    data.load(filename='PS_power_max.csv', param=model.PSmax)

    model.PS_voltage = Param(model.N_PS)
    data.load(filename='PS_voltage.csv', param=model.PS_voltage)

    # Power of the primary substation as sum of all the other powers
    # def PPS_init(model):
    #    return sum(model.Psub[i] for i in model.N)
    # model.PPS=Param(model.PS,initialize=PPS_init)

    # Connection distance of all the edges
    model.dist = Param(model.links)
    data.load(filename='distances.csv', param=model.dist)
    #TODO use the npv cost of the lines
    model.weights = Param(model.links_decision)
    data.load(filename='weights_decision_lines.csv', param=model.weights)
    #data.load(filename='weights_decision_lines_npv.csv', param=model.weights)

    # Electrical parameters of all the cables
    model.V_ref = Param(initialize=voltage)
    model.A_ref = Param(initialize=Abase)
    model.E_min = Param(initialize=Vmin)

    model.R_ref = Param(initialize=resistance)
    model.X_ref = Param(initialize=reactance)
    model.P_max = Param(initialize=Pmax)
    model.cf = Param(initialize=line_cost)

    model.Z = Param(initialize=model.R_ref + model.X_ref * 0.5)
    model.Z_ref = Param(initialize=model.V_ref ** 2 / Abase)

    model.n_clusters = Param(initialize=n_clusters)
    model.coe = Param(initialize=coe)


    #####################Define variables#####################

    # binary variable x[i,j]: 1 if the connection i,j is present, 0 otherwise
    model.x = Var(model.links_decision, within=Binary)
    # power[i,j] is the power flow of connection i-j
    model.P = Var(model.links)
    # positive variables E(i) is p.u. voltage at each node
    model.E = Var(model.N, within=NonNegativeReals)
    # microgrid
    model.z = Var(model.N_MG, within=Binary)
    # binary variable k[i]: 1 if node i is a primary substation, 0 otherwise
    model.k = Var(model.N_PS, within=Binary)
    # Power output of Primary substation
    model.PPS = Var(model.N_PS, within=NonNegativeReals)
    model.MG_output = Var(model.N_MG)

    #####################Define constraints###############################
    # def Make_problem_easy(model,i,j):
    #    return model.x[i,j]+model.weights[i,j]>=1
    # model.easy = Constraint(model.links, rule=Make_problem_easy)

    def Radiality_rule(model):
        # return summation(model.x)==len(model.N)-summation(model.k)
        return summation(model.x) == model.n_clusters

    model.Radiality = Constraint(rule=Radiality_rule)

    def Power_flow_conservation_rule(model, node):
        return (sum(model.P[j, node] for j in model.NodesIn[node]) - sum(
            model.P[node, j] for j in model.NodesOut[node])) == model.Psub[node]

    model.Power_flow_conservation = Constraint(model.N_clusters, rule=Power_flow_conservation_rule)

    def Power_flow_conservation_rule2(model, node):
        return (sum(model.P[j, node] for j in model.NodesIn[node]) - sum(
            model.P[node, j] for j in model.NodesOut[node])) == - model.MG_output[node]

    model.Power_flow_conservation2 = Constraint(model.N_MG, rule=Power_flow_conservation_rule2)

    def Power_flow_conservation_rule3(model, node):
        return (sum(model.P[j, node] for j in model.NodesIn[node]) - sum(
            model.P[node, j] for j in model.NodesOut[node])) == - model.PPS[node]

    model.Power_flow_conservation3 = Constraint(model.N_PS, rule=Power_flow_conservation_rule3)

    def Power_upper_decision(model, i, j):
        return model.P[i, j] <= model.P_max * model.x[i, j]

    model.Power_upper_decision = Constraint(model.links_decision, rule=Power_upper_decision)

    def Power_lower_decision(model, i, j):
        return model.P[i, j] >= -model.P_max * model.x[i, j]

    model.Power_lower_decision = Constraint(model.links_decision, rule=Power_lower_decision)

    def Power_upper_clusters(model, i, j):
        return model.P[i, j] <= model.P_max

    model.Power_upper_clusters = Constraint(model.links_clusters, rule=Power_upper_clusters)

    def Power_lower_clusters(model, i, j):
        return model.P[i, j] >= -model.P_max

    model.Power_lower_clusters = Constraint(model.links_clusters, rule=Power_lower_clusters)

    # Voltage constraints
    def Voltage_balance_rule(model, i, j):
        return (model.E[i] - model.E[j]) + model.x[i, j] - 1 <= model.dist[i, j] / 1000 * model.P[
            i, j] * model.Z / model.Z_ref

    model.Voltage_balance_rule = Constraint(model.links_decision, rule=Voltage_balance_rule)

    def Voltage_balance_rule2(model, i, j):
        return (model.E[i] - model.E[j]) - model.x[i, j] + 1 >= model.dist[i, j] / 1000 * model.P[
            i, j] * model.Z / model.Z_ref

    model.Voltage_balance_rule2 = Constraint(model.links_decision, rule=Voltage_balance_rule2)

    def Voltage_balance_rule3(model, i, j):
        return (model.E[i] - model.E[j]) <= model.dist[i, j] / 1000 * model.P[i, j] * model.Z / model.Z_ref

    model.Voltage_balance_rule3 = Constraint(model.links_clusters, rule=Voltage_balance_rule3)

    def Voltage_balance_rule4(model, i, j):
        return (model.E[i] - model.E[j]) >= model.dist[i, j] / 1000 * model.P[i, j] * model.Z / model.Z_ref

    model.Voltage_balance_rule4 = Constraint(model.links_clusters, rule=Voltage_balance_rule4)

    def Voltage_limit(model, i):
        return model.E[i] >= model.k[i] * (model.PS_voltage[i] - model.E_min) + model.E_min

    model.Voltage_limit = Constraint(model.N_PS, rule=Voltage_limit)

    def Voltage_PS2(model, i):
        return model.E[i] <= model.PS_voltage[i]

    model.Voltage_PS2 = Constraint(model.N_PS, rule=Voltage_PS2)

    def Voltage_limit_MG(model, i):
        return model.E[i] <= 1

    model.Voltage_limit_MG = Constraint(model.N_MG, rule=Voltage_limit_MG)

    def Voltage_limit_MG2(model, i):
        return model.E[i] >= model.z[i] * (1 - model.E_min) + model.E_min

    model.Voltage_limit_MG2 = Constraint(model.N_MG, rule=Voltage_limit_MG2)


    def Voltage_limit_clusters2(model, i):
        return model.E[i] >= model.E_min

    model.Voltage_limit_clusters2 = Constraint(model.N_clusters, rule=Voltage_limit_clusters2)

    def PS_power_rule_upper(model, i):
        return model.PPS[i] <= model.PSmax[i] * model.k[i]

    model.PS_power_upper = Constraint(model.N_PS, rule=PS_power_rule_upper)

    def Balance_rule(model):
        return (sum(model.PPS[i] for i in model.N_PS) + sum(model.MG_output[i] for i in model.N_MG) - sum(
            model.Psub[i] for i in model.N_clusters)) == 0

    model.Balance = Constraint(rule=Balance_rule)

    def MG_power_limit(model, i):
        return model.MG_output[i] == model.z[i] * model.microgrid_power[i]

    model.MG_power_limit = Constraint(model.N_MG, rule=MG_power_limit)

    ####################Define objective function##########################
    print(coe)
    def ObjectiveFunction(model):
        # model.weights is in euro, model.coe is euro/MWh,
        return summation(model.weights, model.x)  + summation(model.mg_cost, model.z) * 1000 + \
               sum(model.energy[i] * (1 - model.z[i]) for i in model.N_MG) * model.coe + summation(model.ps_cost,
                                                                                                   model.k)
        # +sum((model.P[i]/model.A_ref)**2*0.5*1.25*model.R_ref/model.Z_ref*model.dist[i]/1000*24*365*20 for i in model.links)
        # return summation(model.dist,model.x)*model.cf/1000 + summation(model.k) *model.cPS

    model.Obj = Objective(rule=ObjectiveFunction, sense=minimize)

    #############Solve model##################

    instance = model.create_instance(data)
    print('Instance is constructed:', instance.is_constructed())
    # opt = SolverFactory('cbc',executable=r'C:\Users\Asus\Desktop\POLIMI\Thesis\GISELE\Gisele_MILP\cbc')
    opt = SolverFactory('gurobi')
    # opt.options['numericfocus']=0
    opt.options['mipgap'] = 0.02
    opt.options['presolve']=2
    # opt.options['mipfocus'] = 3
    print('Starting optimization process')
    time_i = datetime.now()
    opt.solve(instance, tee=True, symbolic_solver_labels=True)
    time_f = datetime.now()
    print('Time required for optimization is', time_f - time_i)
    links = instance.x
    power = instance.P
    subs = instance.k
    voltage = instance.E
    PS = instance.PPS
    mg_output = instance.MG_output
    microGrid = instance.z

    links_clusters = instance.links_clusters
    # voltage_drop=instance.z
    connections_output = pd.DataFrame(columns=[['id1', 'id2', 'power']])
    PrSubstation = pd.DataFrame(columns=[['index', 'power']])
    all_lines = pd.DataFrame(columns=[['id1', 'id2', 'power']])
    Voltages = pd.DataFrame(columns=[['index', 'voltage [p.u]']])
    Microgrid = pd.DataFrame(columns=[['index', 'microgrid', 'power']])
    Links_Clusters = pd.DataFrame(columns=[['id1', 'id2', 'power']])
    k = 0
    for index in links:
        if int(round(value(links[index]))) == 1:
            connections_output.loc[k, 'id1'] = index[0]
            connections_output.loc[k, 'id2'] = index[1]
            connections_output.loc[k, 'power'] = value(power[index])
            k = k + 1
    k = 0
    for index in subs:
        if int(round(value(subs[index]))) == 1:
            PrSubstation.loc[k, 'index'] = index
            PrSubstation.loc[k, 'power'] = value(PS[index])
            print(((value(PS[index]))))
            k = k + 1
    k = 0
    for v in voltage:
        Voltages.loc[k, 'index'] = v
        Voltages.loc[k, 'voltage [p.u]'] = value(voltage[v])
        k = k + 1
    k = 0
    for index in power:
        all_lines.loc[k, 'id1'] = index[0]
        all_lines.loc[k, 'id2'] = index[1]
        all_lines.loc[k, 'power'] = value(power[index])
        k = k + 1
    k = 0
    for index in mg_output:
        Microgrid.loc[k, 'index'] = index
        Microgrid.loc[k, 'microgrid'] = value(microGrid[index])
        Microgrid.loc[k, 'power'] = value(mg_output[index])
        k = k + 1
    k = 0
    for index in links_clusters:
        Links_Clusters.loc[k, 'id1'] = index[0]
        Links_Clusters.loc[k, 'id2'] = index[1]
        Links_Clusters.loc[k, 'power'] = value(power[index])
        k = k + 1

    Links_Clusters.to_csv(MILP_output_folder + '/links_clusters.csv', index=False)
    connections_output.to_csv(MILP_output_folder + '/connections_output.csv', index=False)
    PrSubstation.to_csv(MILP_output_folder + '/PrimarySubstations.csv', index=False)
    Voltages.to_csv(MILP_output_folder + '/Voltages.csv', index=False)
    all_lines.to_csv(MILP_output_folder + '/all_lines.csv', index=False)
    Microgrid.to_csv(MILP_output_folder + '/Microgrid.csv', index=False)

    # k=0
    # for index in links:
    #     if int(round(value(links[index])))==1:
    #         Voltage_drop.loc[k,'id1']=index[0]
    #         Voltage_drop.loc[k,'id2']=index[1]
    #         Voltage_drop.loc[k,'v_drop']=value(voltage_drop[index])
    #         k=k+1
    # Voltage_drop.to_csv('MILP_results/solution_new/Voltage_drops.csv',index=False)

def MILP_MG_2cables(gisele_folder,case_study,n_clusters,coe,voltage,resistance,reactance,Pmax,line_cost,
                    resistance2,reactance2,Pmax2,line_cost2):


    ############ Create abstract model ###########
    model = AbstractModel()
    data = DataPortal()
    MILP_input_folder = gisele_folder + '/Case studies/' + case_study + '/Intermediate/Optimization/MILP_input'
    MILP_output_folder = gisele_folder + '/Case studies/' + case_study + '/Intermediate/Optimization/MILP_output'
    os.chdir(MILP_input_folder)
    # ####################Define sets#####################

    # Define some basic parameter for the per unit conversion and voltage limitation
    Abase = 1
    Vmin = 0.9

    # Name of all the nodes (primary and secondary substations)
    model.N = Set()
    data.load(filename='nodes.csv', set=model.N)  # first row is not read
    model.N_clusters = Set()
    data.load(filename='nodes_clusters.csv', set=model.N_clusters)
    model.N_MG = Set()
    data.load(filename='microgrids_nodes.csv', set=model.N_MG)
    model.N_PS = Set()
    data.load(filename='nodes_PS.csv', set=model.N_PS)
    # Node corresponding to primary substation

    # Allowed connections
    model.links = Set(dimen=2)  # in the csv the values must be delimited by commas
    data.load(filename='links_all.csv', set=model.links)

    model.links_clusters = Set(dimen=2)
    data.load(filename='links_clusters.csv', set=model.links_clusters)

    model.links_decision = Set(dimen=2)
    data.load(filename='links_decision.csv', set=model.links_decision)

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

    #####################Define parameters#####################

    # Electric power in the nodes (injected (-) or absorbed (+))
    model.Psub = Param(model.N_clusters)
    data.load(filename='power_nodes.csv', param=model.Psub)

    # model.PS=Param(model.N)
    # data.load(filename='PS.csv',param=model.PS)

    model.microgrid_power = Param(model.N_MG)
    data.load(filename='microgrids_powers.csv', param=model.microgrid_power)

    model.energy = Param(model.N_MG)
    data.load(filename='energy.csv', param=model.energy)

    model.mg_cost = Param(model.N_MG)
    data.load(filename='microgrids_costs.csv', param=model.mg_cost)

    model.ps_cost = Param(model.N_PS)
    data.load(filename='PS_costs.csv', param=model.ps_cost)

    model.PSmax = Param(model.N_PS)
    data.load(filename='PS_power_max.csv', param=model.PSmax)

    model.PS_voltage = Param(model.N_PS)
    data.load(filename='PS_voltage.csv', param=model.PS_voltage)

    # Power of the primary substation as sum of all the other powers
    # def PPS_init(model):
    #    return sum(model.Psub[i] for i in model.N)
    # model.PPS=Param(model.PS,initialize=PPS_init)

    # Connection distance of all the edges
    model.dist = Param(model.links)
    data.load(filename='distances.csv', param=model.dist)

    model.weights = Param(model.links_decision)
    data.load(filename='weights_decision_lines.csv', param=model.weights)
    # Electrical parameters of all the cables
    model.V_ref = Param(initialize=voltage)
    model.A_ref = Param(initialize=Abase)
    model.E_min = Param(initialize=Vmin)

    model.R_ref = Param(initialize=resistance)
    model.R_ref2 = Param(initialize=resistance2)
    model.X_ref = Param(initialize=reactance)
    model.X_ref2 = Param(initialize = reactance2)
    model.P_max = Param(initialize=Pmax)
    model.P_max2 = Param(initialize = Pmax2)
    model.cf = Param(initialize=line_cost)
    model.cf2 = Param(initialize = line_cost2)

    model.Z = Param(initialize=model.R_ref + model.X_ref * 0.5)
    model.Z2 = Param(initialize = model.Ref2+model.X_ref2 * 0.5)
    model.Z_ref = Param(initialize=model.V_ref ** 2 / Abase)

    model.n_clusters = Param(initialize=n_clusters)
    model.coe = Param(initialize=coe)


    #####################Define variables#####################

    # binary variable x[i,j]: 1 if the connection i,j is present, 0 otherwise
    model.x = Var(model.links_decision, within=Binary)
    model.x1 = Var(model.links_decision, within=Binary)
    model.y = Var(model.links_clusters, within=Binary)
    # x is if decision link is of cable type 1, x1 if it is of cable type 2. On the other hand, y is simply 1 decision variable
    # since we know that it is chosen, we just want to know if it's the better cable or not. So, we consider that the base cable is the smaller one,
    # and we just add the additional power availability, voltage drop and cost.
    # power[i,j] is the power flow of connection i-j
    model.P = Var(model.links)
    # positive variables E(i) is p.u. voltage at each node
    model.E = Var(model.N, within=NonNegativeReals)
    # microgrid
    model.z = Var(model.N_MG, within=Binary)
    # binary variable k[i]: 1 if node i is a primary substation, 0 otherwise
    model.k = Var(model.N_PS, within=Binary)
    # Power output of Primary substation
    model.PPS = Var(model.N_PS, within=NonNegativeReals)
    model.MG_output = Var(model.N_MG)

    model.delta_Pmax = model.P_max - model.P_max2
    model.delta_Z = model.Z - model.Z2
    model.delta_cf = model.cf - model.cf2

    #####################Define constraints###############################
    # def Make_problem_easy(model,i,j):
    #    return model.x[i,j]+model.weights[i,j]>=1
    # model.easy = Constraint(model.links, rule=Make_problem_easy)

    def Radiality_rule(model):
        # return summation(model.x)==len(model.N)-summation(model.k)
        return summation(model.x) + summation(model.x1) == model.n_clusters

    model.Radiality = Constraint(rule=Radiality_rule)

    def Radiality_rule(model):
        return summation(model.k) + summation(model.z) <= model.n_clusters

    def Power_flow_conservation_rule(model, node):
        return (sum(model.P[j, node] for j in model.NodesIn[node]) - sum(
            model.P[node, j] for j in model.NodesOut[node])) == model.Psub[node]

    model.Power_flow_conservation = Constraint(model.N_clusters, rule=Power_flow_conservation_rule)

    def Power_flow_conservation_rule2(model, node):
        return (sum(model.P[j, node] for j in model.NodesIn[node]) - sum(
            model.P[node, j] for j in model.NodesOut[node])) == - model.MG_output[node]

    model.Power_flow_conservation2 = Constraint(model.N_MG, rule=Power_flow_conservation_rule2)

    def Power_flow_conservation_rule3(model, node):
        return (sum(model.P[j, node] for j in model.NodesIn[node]) - sum(
            model.P[node, j] for j in model.NodesOut[node])) == - model.PPS[node]

    model.Power_flow_conservation3 = Constraint(model.N_PS, rule=Power_flow_conservation_rule3)

    def Power_upper_decision(model, i, j):
        return model.P[i, j] <= model.P_max * model.x[i, j] + model.P_max2 * model.x1[i, j]

    model.Power_upper_decision = Constraint(model.links_decision, rule=Power_upper_decision)

    def Power_lower_decision(model, i, j):
        return model.P[i, j] >= -model.P_max * model.x[i, j] - model.P_max2 * model.x1[i, j]

    model.Power_lower_decision = Constraint(model.links_decision, rule=Power_lower_decision)

    def Power_upper_clusters(model, i, j):
        return model.P[i, j] <= model.P_max2 + model.y[i, j] * model.delta_Pmax

    model.Power_upper_clusters = Constraint(model.links_clusters, rule=Power_upper_clusters)

    def Power_lower_clusters(model, i, j):
        return model.P[i, j] >= -model.P_max2 - model.y[i, j] * model.delta_Pmax

    model.Power_lower_clusters = Constraint(model.links_clusters, rule=Power_lower_clusters)

    # Voltage constraints
    def Voltage_balance_rule(model, i, j):
        return (model.E[i] - model.E[j]) + model.x[i, j] - 1 <= model.dist[i, j] / 1000 * model.P[
            i, j] * model.Z / model.Z_ref

    model.Voltage_balance_rule = Constraint(model.links_decision, rule=Voltage_balance_rule)

    def Voltage_balance_rule11(model, i, j):
        return (model.E[i] - model.E[j]) + model.x1[i, j] - 1 <= model.dist[i, j] / 1000 * model.P[
            i, j] * model.Z2 / model.Z_ref

    model.Voltage_balance_rule11 = Constraint(model.links_decision, rule=Voltage_balance_rule11)

    def Voltage_balance_rule2(model, i, j):
        return (model.E[i] - model.E[j]) - model.x[i, j] + 1 >= model.dist[i, j] / 1000 * model.P[
            i, j] * model.Z / model.Z_ref

    model.Voltage_balance_rule2 = Constraint(model.links_decision, rule=Voltage_balance_rule2)

    def Voltage_balance_rule22(model, i, j):
        return (model.E[i] - model.E[j]) - model.x1[i, j] + 1 >= model.dist[i, j] / 1000 * model.P[
            i, j] * model.Z2 / model.Z_ref

    model.Voltage_balance_rule22 = Constraint(model.links_decision, rule=Voltage_balance_rule22)

    def Voltage_balance_rule3(model, i, j):
        return (model.E[i] - model.E[j]) + model.y[i, j] - 1 <= model.dist[i, j] / 1000 * model.P[
            i, j] * model.Z / model.Z_ref

    model.Voltage_balance_rule3 = Constraint(model.links_clusters, rule=Voltage_balance_rule3)

    def Voltage_balance_rule33(model, i, j):
        return (model.E[i] - model.E[j]) - model.y[i, j] <= model.dist[i, j] / 1000 * model.P[
            i, j] * model.Z2 / model.Z_ref

    model.Voltage_balance_rule33 = Constraint(model.links_clusters, rule=Voltage_balance_rule33)

    def Voltage_balance_rule4(model, i, j):
        return (model.E[i] - model.E[j]) - model.y[i, j] + 1 >= model.dist[i, j] / 1000 * model.P[
            i, j] * model.Z / model.Z_ref

    model.Voltage_balance_rule4 = Constraint(model.links_clusters, rule=Voltage_balance_rule4)

    def Voltage_balance_rule44(model, i, j):
        return (model.E[i] - model.E[j]) + model.y[i, j] >= model.dist[i, j] / 1000 * model.P[
            i, j] * model.Z2 / model.Z_ref

    model.Voltage_balance_rule44 = Constraint(model.links_clusters, rule=Voltage_balance_rule44)

    def Voltage_limit(model, i):
        return model.E[i] >= model.k[i] * (model.PS_voltage[i] - model.E_min) + model.E_min

    model.Voltage_limit = Constraint(model.N_PS, rule=Voltage_limit)

    def Voltage_PS2(model, i):
        return model.E[i] <= model.PS_voltage[i]

    model.Voltage_PS2 = Constraint(model.N_PS, rule=Voltage_PS2)

    def Voltage_limit_MG(model, i):
        return model.E[i] <= 1

    model.Voltage_limit_MG = Constraint(model.N_MG, rule=Voltage_limit_MG)

    def Voltage_limit_MG2(model, i):
        return model.E[i] >= model.z[i] * (1 - model.E_min) + model.E_min

    model.Voltage_limit_MG2 = Constraint(model.N_MG, rule=Voltage_limit_MG2)


    def Voltage_limit_clusters2(model, i):
        return model.E[i] >= model.E_min

    model.Voltage_limit_clusters2 = Constraint(model.N_clusters, rule=Voltage_limit_clusters2)

    def PS_power_rule_upper(model, i):
        return model.PPS[i] <= model.PSmax[i] * model.k[i]

    model.PS_power_upper = Constraint(model.N_PS, rule=PS_power_rule_upper)

    def Balance_rule(model):
        return (sum(model.PPS[i] for i in model.N_PS) + sum(model.MG_output[i] for i in model.N_MG) - sum(
            model.Psub[i] for i in model.N_clusters)) == 0

    model.Balance = Constraint(rule=Balance_rule)

    def anti_paralel(model, i, j):
        return model.x[i, j] + model.x1[i, j] <= 1

    model.anti_paralel = Constraint(model.links_decision, rule=anti_paralel)

    def MG_power_limit(model, i):
        return model.MG_output[i] == model.z[i] * model.microgrid_power[i]

    model.MG_power_limit = Constraint(model.N_MG, rule=MG_power_limit)

    ####################Define objective function##########################

    def ObjectiveFunction(model):
        # return summation(model.weights, model.x) * model.cf / 1000 + summation(model.weights, model.x1) * model.cf2 / 10000
        return summation(model.weights, model.x) * model.cf / 10000 + summation(model.mg_cost, model.z) * 1000 + sum(
            model.energy[i] * (1 - model.z[i]) for i in model.N_MG) * model.coe \
               + summation(model.ps_cost, model.k) + summation(model.weights, model.x1) * model.cf2 / 10000 + sum(
            model.y[i] for i in model.links_clusters) * model.delta_cf

    model.Obj = Objective(rule=ObjectiveFunction, sense=minimize)

    #############Solve model##################

    instance = model.create_instance(data)
    print('Instance is constructed:', instance.is_constructed())
    # opt = SolverFactory('cbc',executable=r'C:\Users\Asus\Desktop\POLIMI\Thesis\GISELE\Gisele_MILP\cbc')
    opt = SolverFactory('gurobi')
    # opt.options['numericfocus']=0
    opt.options['mipgap'] = 0.01
    opt.options['presolve'] = 2
    # opt.options['mipfocus']=2
    # opt = SolverFactory('cbc',executable=r'C:\Users\Asus\Desktop\POLIMI\Thesis\GISELE\New folder\cbc')
    print('Starting optimization process')
    time_i = datetime.now()
    opt.solve(instance, tee=True, symbolic_solver_labels=True)
    time_f = datetime.now()
    print('Time required for optimization is', time_f - time_i)
    links = instance.x
    links_small = instance.x1
    power = instance.P
    subs = instance.k
    voltage = instance.E
    PS = instance.PPS
    mg_output = instance.MG_output
    microGrid = instance.z
    links_clusters_type = instance.y
    links_clusters = instance.links_clusters
    # voltage_drop=instance.z
    connections_output = pd.DataFrame(columns=[['id1', 'id2', 'power', 'Type']])
    PrSubstation = pd.DataFrame(columns=[['index', 'power']])
    all_lines = pd.DataFrame(columns=[['id1', 'id2', 'power']])
    Voltages = pd.DataFrame(columns=[['index', 'voltage [p.u]']])
    Microgrid = pd.DataFrame(columns=[['index', 'microgrid', 'power']])
    Links_Clusters = pd.DataFrame(columns=[['id1', 'id2', 'power', 'Type']])
    k = 0
    for index in links:
        if int(round(value(links[index]))) == 1:
            connections_output.loc[k, 'id1'] = index[0]
            connections_output.loc[k, 'id2'] = index[1]
            connections_output.loc[k, 'power'] = value(power[index])
            connections_output.loc[k, 'Type'] = 'Large'
            k = k + 1
        elif int(round(value(links_small[index]))) == 1:
            connections_output.loc[k, 'id1'] = index[0]
            connections_output.loc[k, 'id2'] = index[1]
            connections_output.loc[k, 'power'] = value(power[index])
            connections_output.loc[k, 'Type'] = 'Small'
            k = k + 1
    k = 0
    for index in subs:
        if int(round(value(subs[index]))) == 1:
            PrSubstation.loc[k, 'index'] = index
            PrSubstation.loc[k, 'power'] = value(PS[index])
            print(((value(PS[index]))))
            k = k + 1
    k = 0
    for v in voltage:
        Voltages.loc[k, 'index'] = v
        Voltages.loc[k, 'voltage [p.u]'] = value(voltage[v])
        k = k + 1
    k = 0
    for index in power:
        all_lines.loc[k, 'id1'] = index[0]
        all_lines.loc[k, 'id2'] = index[1]
        all_lines.loc[k, 'power'] = value(power[index])
        k = k + 1
    k = 0
    for index in mg_output:
        Microgrid.loc[k, 'index'] = index
        Microgrid.loc[k, 'microgrid'] = value(microGrid[index])
        Microgrid.loc[k, 'power'] = value(mg_output[index])
        k = k + 1
    k = 0
    for index in links_clusters:
        print(index)
        Links_Clusters.loc[k, 'id1'] = index[0]
        Links_Clusters.loc[k, 'id2'] = index[1]
        Links_Clusters.loc[k, 'power'] = value(power[index])
        if int(round(value(links_clusters_type[index]))) == 1:
            Links_Clusters.loc[k, 'Type'] = 'Large'
        else:
            Links_Clusters.loc[k, 'Type'] = 'Small'
        k = k + 1

    Links_Clusters.to_csv(MILP_output_folder + '/links_clusters.csv', index=False)
    connections_output.to_csv(MILP_output_folder + '/connections_output.csv', index=False)
    PrSubstation.to_csv(MILP_output_folder + '/PrimarySubstations.csv', index=False)
    Voltages.to_csv(MILP_output_folder + '/Voltages.csv', index=False)
    all_lines.to_csv(MILP_output_folder + '/all_lines.csv', index=False)
    Microgrid.to_csv(MILP_output_folder + '/Microgrid.csv', index=False)

    # k=0
    # for index in links:
    #     if int(round(value(links[index])))==1:
    #         Voltage_drop.loc[k,'id1']=index[0]
    #         Voltage_drop.loc[k,'id2']=index[1]
    #         Voltage_drop.loc[k,'v_drop']=value(voltage_drop[index])
    #         k=k+1
    # Voltage_drop.to_csv('MILP_results/solution_new/Voltage_drops.csv',index=False)

def MILP_base_no_voltage(gisele_folder,case_study,n_clusters,coe,voltage,resistance,reactance,Pmax,line_cost):
    ############ Create abstract model ###########
    model = AbstractModel()
    data = DataPortal()
    MILP_input_folder = gisele_folder + '/Case studies/' + case_study + '/Intermediate/Optimization/MILP_input'
    MILP_output_folder = gisele_folder + '/Case studies/' + case_study + '/Intermediate/Optimization/MILP_output'
    os.chdir(MILP_input_folder)

    # Define some basic parameter for the per unit conversion and voltage limitation
    Abase = 1
    Vmin = 0.9

    ####################Define sets#####################
    model.N = Set()
    data.load(filename='nodes.csv', set=model.N)  # first row is not read
    model.N_clusters = Set()
    data.load(filename='nodes_clusters.csv', set=model.N_clusters)
    model.N_PS = Set()
    data.load(filename='nodes_PS.csv', set=model.N_PS)
    # Node corresponding to primary substation

    # Allowed connections
    model.links = Set(dimen=2)  # in the csv the values must be delimited by commas
    data.load(filename='links_all.csv', set=model.links)

    model.links_clusters = Set(dimen=2)
    data.load(filename='links_clusters.csv', set=model.links_clusters)

    model.links_decision = Set(dimen=2)
    data.load(filename='links_decision.csv', set=model.links_decision)
    # Connection distance of all the edges
    model.dist = Param(model.links)
    data.load(filename='distances.csv', param=model.dist)
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

    #####################Define parameters#####################

    # Electric power in the nodes (injected (-) or absorbed (+))
    model.Psub = Param(model.N_clusters)
    data.load(filename='power_nodes.csv', param=model.Psub)

    model.ps_cost = Param(model.N_PS)
    data.load(filename='PS_costs.csv', param=model.ps_cost)

    model.PSmax = Param(model.N_PS)
    data.load(filename='PS_power_max.csv', param=model.PSmax)


    model.weights = Param(model.links_decision)
    data.load(filename='weights_decision_lines.csv', param=model.weights)
    # Electrical parameters of all the cables
    model.V_ref = Param(initialize=voltage)
    model.A_ref = Param(initialize=Abase)
    model.E_min = Param(initialize=Vmin)

    model.R_ref = Param(initialize=resistance)
    model.X_ref = Param(initialize=reactance)
    model.P_max = Param(initialize=Pmax)
    model.cf = Param(initialize=line_cost)

    model.Z = Param(initialize=model.R_ref + model.X_ref * 0.5)
    model.Z_ref = Param(initialize=model.V_ref ** 2 / Abase)

    model.n_clusters = Param(initialize=n_clusters)
    model.coe = Param(initialize=coe)

    #####################Define variables#####################

    # binary variable x[i,j]: 1 if the connection i,j is present, 0 otherwise
    model.x = Var(model.links_decision, within=Binary)
    # power[i,j] is the power flow of connection i-j
    model.P = Var(model.links)

    # binary variable k[i]: 1 if node i is a primary substation, 0 otherwise
    model.k = Var(model.N_PS, within=Binary)
    # Power output of Primary substation
    model.PPS = Var(model.N_PS, within=NonNegativeReals)
    model.cable_type = Var(model.links)

    #####################Define constraints###############################

    # Radiality constraint
    def Radiality_rule(model):
        return summation(model.x) == model.n_clusters

    model.Radiality = Constraint(rule=Radiality_rule)

    # Power flow constraints
    def Power_flow_conservation_rule(model, node):
        return (sum(model.P[j, node] for j in model.NodesIn[node]) - sum(
            model.P[node, j] for j in model.NodesOut[node])) == model.Psub[node]

    model.Power_flow_conservation = Constraint(model.N_clusters, rule=Power_flow_conservation_rule)

    def Power_flow_conservation_rule3(model, node):
        return (sum(model.P[j, node] for j in model.NodesIn[node]) - sum(
            model.P[node, j] for j in model.NodesOut[node])) == - model.PPS[node]

    model.Power_flow_conservation3 = Constraint(model.N_PS, rule=Power_flow_conservation_rule3)

    def Power_upper_decision(model, i, j):
        return model.P[i, j] <= model.P_max * model.x[i, j]

    model.Power_upper_decision = Constraint(model.links_decision, rule=Power_upper_decision)

    def Power_lower_decision(model, i, j):
        return model.P[i, j] >= -model.P_max * model.x[i, j]

    model.Power_lower_decision = Constraint(model.links_decision, rule=Power_lower_decision)
    def Power_upper_cluster(model, i, j):
        return model.P[i, j] <= model.P_max

    model.Power_upper_cluster = Constraint(model.links_clusters, rule=Power_upper_cluster)

    def Power_lower_cluster(model, i, j):
        return model.P[i, j] >= -model.P_max

    model.Power_lower_cluster= Constraint(model.links_clusters, rule=Power_lower_cluster)


    def PS_power_rule_upper(model, i):
        return model.PPS[i] <= model.PSmax[i] * model.k[i]

    model.PS_power_upper = Constraint(model.N_PS, rule=PS_power_rule_upper)

    def Balance_rule(model):
        return (sum(model.PPS[i] for i in model.N_PS) - sum(model.Psub[i] for i in model.N_clusters)) == 0

    model.Balance = Constraint(rule=Balance_rule)

    ####################Define objective function##########################

    ####################Define objective function##########################
    reliability_index = 1000

    def ObjectiveFunction(model):
        return summation(model.weights, model.x) * model.cf / 1000 + summation(model.ps_cost, model.k)

        # return summation(model.weights, model.x) * model.cf / 1000  + summation(model.ps_cost,model.k) - sum(model.Psub[i]*model.Distance[i] for i in model.N_clusters)*reliability_index

    model.Obj = Objective(rule=ObjectiveFunction, sense=minimize)

    #############Solve model##################

    instance = model.create_instance(data)
    print('Instance is constructed:', instance.is_constructed())
    # opt = SolverFactory('cbc',executable=r'C:\Users\Asus\Desktop\POLIMI\Thesis\GISELE\Gisele_MILP\cbc')
    opt = SolverFactory('gurobi')
    opt.options['TimeLimit'] = 300
    # opt.options['numericfocus']=0
    # opt.options['mipgap'] = 0.0002
    # opt.options['presolve']=2
    # opt.options['mipfocus']=2
    # opt = SolverFactory('cbc',executable=r'C:\Users\Asus\Desktop\POLIMI\Thesis\GISELE\New folder\cbc')
    print('Starting optimization process')
    time_i = datetime.now()
    opt.solve(instance, tee=True, symbolic_solver_labels=True)
    time_f = datetime.now()
    print('Time required for optimization is', time_f - time_i)
    links = instance.x
    power = instance.P
    subs = instance.k
    PS = instance.PPS


    links_clusters = instance.links_clusters
    # voltage_drop=instance.z
    connections_output = pd.DataFrame(columns=[['id1', 'id2', 'power']])
    PrSubstation = pd.DataFrame(columns=[['index', 'power']])
    all_lines = pd.DataFrame(columns=[['id1', 'id2', 'power']])
    Links_Clusters = pd.DataFrame(columns=[['id1', 'id2', 'power']])

    k = 0
    for index in links:
        if int(round(value(links[index]))) == 1:
            connections_output.loc[k, 'id1'] = index[0]
            connections_output.loc[k, 'id2'] = index[1]
            connections_output.loc[k, 'power'] = value(power[index])
            k = k + 1
    k = 0
    for index in subs:
        if int(round(value(subs[index]))) == 1:
            PrSubstation.loc[k, 'index'] = index
            PrSubstation.loc[k, 'power'] = value(PS[index])
            print(((value(PS[index]))))
            k = k + 1

    for index in power:
        all_lines.loc[k, 'id1'] = index[0]
        all_lines.loc[k, 'id2'] = index[1]
        all_lines.loc[k, 'power'] = value(power[index])
        k = k + 1
    k = 0

    for index in links_clusters:
        Links_Clusters.loc[k, 'id1'] = index[0]
        Links_Clusters.loc[k, 'id2'] = index[1]
        Links_Clusters.loc[k, 'power'] = value(power[index])
        k = k + 1

    Links_Clusters.to_csv(MILP_output_folder + '/links_clusters.csv', index=False)
    connections_output.to_csv(MILP_output_folder + '/connections_output.csv', index=False)
    PrSubstation.to_csv(MILP_output_folder + '/PrimarySubstations.csv', index=False)
    all_lines.to_csv(MILP_output_folder + '/all_lines.csv', index=False)

def MILP_base(gisele_folder,case_study,n_clusters,coe,voltage,resistance,reactance,Pmax,line_cost):
    ############ Create abstract model ###########
    model = AbstractModel()
    data = DataPortal()
    MILP_input_folder = gisele_folder + '/Case studies/' + case_study + '/Intermediate/Optimization/MILP_input'
    MILP_output_folder = gisele_folder + '/Case studies/' + case_study + '/Intermediate/Optimization/MILP_output'
    os.chdir(MILP_input_folder)

    # Define some basic parameter for the per unit conversion and voltage limitation
    Abase = 1
    Vmin = 0.9

    ####################Define sets#####################
    model.N = Set()
    data.load(filename='nodes.csv', set=model.N)  # first row is not read
    model.N_clusters = Set()
    data.load(filename='nodes_clusters.csv', set=model.N_clusters)
    model.N_PS = Set()
    data.load(filename='nodes_PS.csv', set=model.N_PS)
    # Node corresponding to primary substation

    # Allowed connections
    model.links = Set(dimen=2)  # in the csv the values must be delimited by commas
    data.load(filename='links_all.csv', set=model.links)

    model.links_clusters = Set(dimen=2)
    data.load(filename='links_clusters.csv', set=model.links_clusters)

    model.links_decision = Set(dimen=2)
    data.load(filename='links_decision.csv', set=model.links_decision)
    # Connection distance of all the edges
    model.dist = Param(model.links)
    data.load(filename='distances.csv', param=model.dist)
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

    #####################Define parameters#####################

    # Electric power in the nodes (injected (-) or absorbed (+))
    model.Psub = Param(model.N_clusters)
    data.load(filename='power_nodes.csv', param=model.Psub)

    model.ps_cost = Param(model.N_PS)
    data.load(filename='PS_costs.csv', param=model.ps_cost)

    model.PSmax = Param(model.N_PS)
    data.load(filename='PS_power_max.csv', param=model.PSmax)

    model.PS_voltage = Param(model.N_PS)
    data.load(filename='PS_voltage.csv', param=model.PS_voltage)


    model.weights = Param(model.links_decision)
    data.load(filename='weights_decision_lines.csv', param=model.weights)
    # Electrical parameters of all the cables
    model.V_ref = Param(initialize=voltage)
    model.A_ref = Param(initialize=Abase)
    model.E_min = Param(initialize=Vmin)

    model.R_ref = Param(initialize=resistance)
    model.X_ref = Param(initialize=reactance)
    model.P_max = Param(initialize=Pmax)
    model.cf = Param(initialize=line_cost)

    model.Z = Param(initialize=model.R_ref + model.X_ref * 0.5)
    model.Z_ref = Param(initialize=model.V_ref ** 2 / Abase)

    model.n_clusters = Param(initialize=n_clusters)
    model.coe = Param(initialize=coe)

    #####################Define variables#####################

    # binary variable x[i,j]: 1 if the connection i,j is present, 0 otherwise
    model.x = Var(model.links_decision, within=Binary)
    # power[i,j] is the power flow of connection i-j
    model.P = Var(model.links)
    # positive variables E(i) is p.u. voltage at each node
    model.E = Var(model.N, within=NonNegativeReals)

    # binary variable k[i]: 1 if node i is a primary substation, 0 otherwise
    model.k = Var(model.N_PS, within=Binary)
    # Power output of Primary substation
    model.PPS = Var(model.N_PS, within=NonNegativeReals)
    model.cable_type = Var(model.links)

    #####################Define constraints###############################

    # Radiality constraint
    def Radiality_rule(model):
        return summation(model.x) == model.n_clusters

    model.Radiality = Constraint(rule=Radiality_rule)

    # Power flow constraints
    def Power_flow_conservation_rule(model, node):
        return (sum(model.P[j, node] for j in model.NodesIn[node]) - sum(
            model.P[node, j] for j in model.NodesOut[node])) == model.Psub[node]

    model.Power_flow_conservation = Constraint(model.N_clusters, rule=Power_flow_conservation_rule)

    def Power_flow_conservation_rule3(model, node):
        return (sum(model.P[j, node] for j in model.NodesIn[node]) - sum(
            model.P[node, j] for j in model.NodesOut[node])) == - model.PPS[node]

    model.Power_flow_conservation3 = Constraint(model.N_PS, rule=Power_flow_conservation_rule3)

    def Power_upper_decision(model, i, j):
        return model.P[i, j] <= model.P_max * model.x[i, j]

    model.Power_upper_decision = Constraint(model.links_decision, rule=Power_upper_decision)

    def Power_lower_decision(model, i, j):
        return model.P[i, j] >= -model.P_max * model.x[i, j]

    model.Power_lower_decision = Constraint(model.links_decision, rule=Power_lower_decision)
    def Power_upper_cluster(model, i, j):
        return model.P[i, j] <= model.P_max

    model.Power_upper_cluster = Constraint(model.links_clusters, rule=Power_upper_cluster)

    def Power_lower_cluster(model, i, j):
        return model.P[i, j] >= -model.P_max

    model.Power_lower_cluster= Constraint(model.links_clusters, rule=Power_lower_cluster)

    # Voltage constraints
    def Voltage_balance_rule(model, i, j):
        return (model.E[i] - model.E[j]) + model.x[i, j] - 1 <= model.dist[i, j] / 1000 * model.P[
            i, j] * model.Z / model.Z_ref

    model.Voltage_balance_rule = Constraint(model.links_decision, rule=Voltage_balance_rule)

    def Voltage_balance_rule2(model, i, j):
        return (model.E[i] - model.E[j]) - model.x[i, j] + 1 >= model.dist[i, j] / 1000 * model.P[
            i, j] * model.Z / model.Z_ref

    model.Voltage_balance_rule2 = Constraint(model.links_decision, rule=Voltage_balance_rule2)

    def Voltage_balance_rule3(model, i, j):
        return (model.E[i] - model.E[j]) <= model.dist[i, j] / 1000 * model.P[i, j] * model.Z / model.Z_ref

    model.Voltage_balance_rule3 = Constraint(model.links_clusters, rule=Voltage_balance_rule3)

    def Voltage_balance_rule4(model, i, j):
        return (model.E[i] - model.E[j]) >= model.dist[i, j] / 1000 * model.P[i, j] * model.Z / model.Z_ref

    model.Voltage_balance_rule4 = Constraint(model.links_clusters, rule=Voltage_balance_rule4)



    def Voltage_limit(model, i):
        return model.E[i] >= model.k[i] * (model.PS_voltage[i] - model.E_min) + model.E_min

    model.Voltage_limit = Constraint(model.N_PS, rule=Voltage_limit)

    def Voltage_PS2(model, i):
        return model.E[i] <= model.PS_voltage[i]

    model.Voltage_PS2 = Constraint(model.N_PS, rule=Voltage_PS2)

    def Voltage_limit_clusters2(model, i):
        return model.E[i] >= model.E_min

    model.Voltage_limit_clusters2 = Constraint(model.N_clusters, rule=Voltage_limit_clusters2)

    def PS_power_rule_upper(model, i):
        return model.PPS[i] <= model.PSmax[i] * model.k[i]

    model.PS_power_upper = Constraint(model.N_PS, rule=PS_power_rule_upper)

    def Balance_rule(model):
        return (sum(model.PPS[i] for i in model.N_PS) - sum(model.Psub[i] for i in model.N_clusters)) == 0

    model.Balance = Constraint(rule=Balance_rule)

    ####################Define objective function##########################

    ####################Define objective function##########################
    reliability_index = 1000

    def ObjectiveFunction(model):
        return summation(model.weights, model.x) * model.cf / 1000 + summation(model.ps_cost, model.k)
        # return summation(model.cost_asset,model.flex_positive) + return summation(model.cost_asset, model.flex_negative)

    model.Obj = Objective(rule=ObjectiveFunction, sense=minimize)

    #############Solve model##################

    instance = model.create_instance(data)
    print('Instance is constructed:', instance.is_constructed())
    # opt = SolverFactory('cbc',executable=r'C:\Users\Asus\Desktop\POLIMI\Thesis\GISELE\Gisele_MILP\cbc')
    opt = SolverFactory('gurobi')
    opt.options['TimeLimit'] = 600
    # opt.options['numericfocus']=0
    # opt.options['mipgap'] = 0.0002
    opt.options['presolve']=0
    opt.options['PreCrush']=1
    # opt.options['mipfocus']=2
    # opt = SolverFactory('cbc',executable=r'C:\Users\Asus\Desktop\POLIMI\Thesis\GISELE\New folder\cbc')
    print('Starting optimization process')
    time_i = datetime.now()
    opt.solve(instance, tee=True, symbolic_solver_labels=True)
    time_f = datetime.now()
    print('Time required for optimization is', time_f - time_i)
    links = instance.x
    power = instance.P
    subs = instance.k
    voltage = instance.E
    PS = instance.PPS


    links_clusters = instance.links_clusters
    # voltage_drop=instance.z
    connections_output = pd.DataFrame(columns=[['id1', 'id2', 'power']])
    PrSubstation = pd.DataFrame(columns=[['index', 'power']])
    all_lines = pd.DataFrame(columns=[['id1', 'id2', 'power']])
    Voltages = pd.DataFrame(columns=[['index', 'voltage [p.u]']])
    Links_Clusters = pd.DataFrame(columns=[['id1', 'id2', 'power']])
    distance = pd.DataFrame(columns=[['index', 'length[km]']])
    k = 0
    for index in links:
        if int(round(value(links[index]))) == 1:
            connections_output.loc[k, 'id1'] = index[0]
            connections_output.loc[k, 'id2'] = index[1]
            connections_output.loc[k, 'power'] = value(power[index])
            k = k + 1
    k = 0
    for index in subs:
        if int(round(value(subs[index]))) == 1:
            PrSubstation.loc[k, 'index'] = index
            PrSubstation.loc[k, 'power'] = value(PS[index])
            print(((value(PS[index]))))
            k = k + 1
    k = 0
    for v in voltage:
        Voltages.loc[k, 'index'] = v
        Voltages.loc[k, 'voltage [p.u]'] = value(voltage[v])
        k = k + 1
    k = 0
    for index in power:
        all_lines.loc[k, 'id1'] = index[0]
        all_lines.loc[k, 'id2'] = index[1]
        all_lines.loc[k, 'power'] = value(power[index])
        k = k + 1
    k = 0

    for index in links_clusters:
        Links_Clusters.loc[k, 'id1'] = index[0]
        Links_Clusters.loc[k, 'id2'] = index[1]
        Links_Clusters.loc[k, 'power'] = value(power[index])
        k = k + 1

    Links_Clusters.to_csv(MILP_output_folder + '/links_clusters.csv', index=False)
    connections_output.to_csv(MILP_output_folder + '/connections_output.csv', index=False)
    PrSubstation.to_csv(MILP_output_folder + '/PrimarySubstations.csv', index=False)
    Voltages.to_csv(MILP_output_folder + '/Voltages.csv', index=False)
    all_lines.to_csv(MILP_output_folder + '/all_lines.csv', index=False)