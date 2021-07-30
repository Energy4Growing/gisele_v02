from pyomo.environ import Param, RangeSet, Set, Var, NonNegativeReals, NonNegativeIntegers, Binary, Constraint, Objective, minimize, AbstractModel
from pyomo.opt import SolverFactory
from pyomo.dataportal.DataPortal import DataPortal
from math import sqrt
from pandas import read_excel, read_csv

model= AbstractModel()
data = DataPortal()


# Parameters related to set definition
model.project_duration = Param()
model.pv_types = Param()
model.wt_types = Param()
model.bess_types = Param()
model.dg_types = Param()

# SETS
model.hours = RangeSet(1, model.project_duration)
model.pv = RangeSet(1, model.pv_types)
model.wt = RangeSet(1, model.wt_types)
model.bess = RangeSet(1, model.bess_types)
model.dg = RangeSet(1, model.dg_types)

# PARAMETERS


# Parameters of the PV
model.pv_nominal_capacity = Param(model.pv)  # Nominal capacity of the PV in kW/unit
model.pv_investment_cost = Param(model.pv)  # Cost of solar panel in €/kW
model.pv_OM_cost = Param(model.pv)  # Cost of O&M solar panel in €/kW/y
model.pv_max_units = Param(model.pv)  # Maximum number of installed units [-]
model.pv_life = Param(model.pv)  # Lifetime of panels [y]

# Parameters of the Wind Turbine
model.wt_nominal_capacity = Param(model.wt)  # Nominal capacity of the WT in kW/unit
model.wt_investment_cost = Param(model.wt)  # Cost of WT in €/kW
model.wt_OM_cost = Param(model.wt)  # Cost of O&M WT in €/kW/y
model.wt_max_units = Param(model.wt)  # Maximum number of installed units [-]
model.wt_life = Param(model.wt)  # Lifetime of WT [y]

# Parameters of the Storage System
model.bess_nominal_capacity = Param(model.bess)  # Nominal capacity of the BESS in kWh/unit
model.bess_investment_cost = Param(model.bess)  # Cost of BESS in €/kWh
model.bess_OM_cost = Param(model.bess)  # Cost of O&M BESS in €/kWh/y
model.bess_max_units = Param(model.bess)  # Maximum number of installed units [-]
model.bess_life = Param(model.bess)  # Lifetime of BESS [y]
model.bess_cycle_life = Param(model.bess)  # cycling lifetime of BESS [kwh]
model.bess_depth_of_discharge = Param(model.bess)  # Depth of discharge of the battery (DOD) [0-1]
model.bess_P_ratio_max = Param(model.bess)  # Maximum Ratio of Power/Energy [kW/kWh]
model.bess_ch_efficiency_nom = Param(model.bess)  # Efficiency of the charge of the battery [0-1]
model.bess_dis_efficiency_nom = Param(model.bess)  # Efficiency of the discharge of the battery [0-1]
model.bess_initial_SOC = Param(model.bess)  # initiali state of charge [0-1]

# Parameters of the Diesel Generator
model.dg_nominal_capacity = Param(model.dg)  # Nominal capacity of the DG in kW/unit
model.dg_investment_cost = Param(model.dg)  # Cost of DG in €/kW
model.dg_OM_cost = Param(model.dg)  # Cost of O&M DG in €/kW/y
model.dg_max_units = Param(model.dg)  # Maximum number of installed units [-]
model.dg_life = Param(model.dg)  # Lifetime of DG [functioning hours]
model.dg_P_min = Param(model.dg)  # Minimum power [0-1]
model.dg_cost_coeff_A = Param(model.dg)  # linear Fuel consumption curve coefficient [l/h]
model.dg_cost_coeff_B = Param(model.dg)  # linear Fuel consumption curve coefficient [l/h/kW]
model.dg_SU_cost = Param(model.dg)  # Start up cost [€]
model.dg_SD_cost = Param(model.dg)  # Shut down cost [€]
model.dg_UT = Param(model.dg)  # minimum up time [h]
model.dg_DT = Param(model.dg)  # minimum down time [h]
model.dg_RU = Param(model.dg)  # ramp up limit [kW/h]
model.dg_RD = Param(model.dg)  # ramp down limit [kW/h]

# Scalars
model.inflation_rate = Param()  # inflation rate [0-1]
model.nominal_interest_rate = Param()  # nominal interest rate [0-1]
model.lost_load_max = Param()  # maximum admitted loss of load [0-1]
model.lost_load_value = Param()
model.fuel_cost = Param()  # cost of diesel [€/l]
model.inverter_cost = Param()  # investment cost of inverter [€/kW]
model.inverter_life = Param()  # lifetime of inverter [y]
model.load_forecast_error = Param()  # [0-1]
model.pv_forecast_error = Param()  # error on power produced from pv panels[0-1]
model.wt_forecast_error = Param()  # [0-1]
model.demand_growth = Param()  # yearly demand growth [0-1]
model.M = Param()  # big number
model.epsilon = Param()  # small number

# initialize parameters

def Initialize_wear_cost(model,b):
    return model.bess_investment_cost[b] / (model.bess_cycle_life[b] * sqrt(model.bess_ch_efficiency_nom[b] * model.bess_dis_efficiency_nom[b]))

model.bess_wear_cost = Param(model.bess, initialize=Initialize_wear_cost)  # cycling degradation cost [€/kwh]

def Initialize_dg_repl_cost(model,g):
    return model.dg_investment_cost[g] / model.dg_life[g]

model.dg_repl_cost = Param(model.dg, initialize=Initialize_dg_repl_cost)  # unitary replacement dg cost [€/h ON]

def Initialize_ir(model):
    return (model.nominal_interest_rate - model.inflation_rate) / (1 + model.inflation_rate)

model.ir = Param(initialize=Initialize_ir)  # real interest rate [0-1]

def Initialize_Discount_Rate(model,h):
    return 1 / (1+model.ir)**(h//8760)

model.discount_rate = Param(model.hours, initialize=Initialize_Discount_Rate)  # discount rate [0-1]

'''
def Initialize_hours_last(model):
    if model.hours % 8760 == 0:
        return model.hours

model.hours_last = Set(model.hours, initialize=Initialize_hours_last)  # set of last hour of each year


def Initialize_end_of_project(model, h):
    if h == model.project_duration:
        return h

model.end_of_project = Set(model.hours, initialize=Initialize_end_of_project)  # set of last hour of each year
'''

# input profiles
model.input_load = Param(model.hours)  # hourly load profile [kWh]
model.input_pv_prod = Param(model.hours, model.pv)  # hourly PV production [kWh]
model.input_wt_prod = Param(model.hours, model.wt)  # hourly WT production [kWh]

#model.input_load = read_csv('load_Soroti_20.csv')
#model.input_pv_prod = read_csv('solarPV_Soroti_20.csv')
#model.input_wt_prod = read_csv('windPower_Soroti_20.csv')

data.load(filename='load_Soroti.csv', param=model.input_load)
data.load(filename='solarPV_Soroti_10.csv', param=model.input_pv_prod)
data.load(filename='windPower_Soroti_10.csv', param=model.input_wt_prod)

data.load(filename='data.dat', param=model.project_duration)
data.load(filename='data_scalars.dat')


# VARIABLES

# Variables associated to the project
model.NPC = Var()  # Net Present Cost [k€]
model.initial_investment = Var(within=NonNegativeReals)
model.OM_cost = Var(within=NonNegativeReals)
model.replacement_cost = Var(within=NonNegativeReals)
model.salvage_value = Var(within=NonNegativeReals)

# Variables associated to RES
model.pv_units = Var(model.pv, within=NonNegativeReals)  # Number of units of solar panels
model.wt_units = Var(model.wt, within=NonNegativeIntegers)  # Number of units of wind turbines
model.total_power_res = Var(model.hours, within=NonNegativeReals)  # Power generated from the PV and WT [kW]

# Variables associated to the battery bank
model.bess_units = Var(model.bess, within=NonNegativeReals)  # Number of units of batteries
model.bess_dis_power = Var(model.hours, model.bess, within=NonNegativeReals)  # Battery discharging power in kW
model.bess_ch_power = Var(model.hours, model.bess, within=NonNegativeReals)  # Battery charging power in kW
model.bess_total_energy = Var(model.hours, model.bess, within=NonNegativeReals)  # Battery charge level at h [kWh]
model.bess_power_max = Var(model.bess,
                           within=NonNegativeReals)  # maximum power withdrawn or injected by the batteries [kW]
model.bess_bin = Var(model.hours, model.bess, within=Binary)  # Binary variable, 1 if charging mode

# Variables associated to the diesel generator
model.dg_units = Var(model.dg, within=NonNegativeIntegers)  # Number of units of diesel generators
model.dg_power = Var(model.hours, model.dg, within=NonNegativeReals)  # Power level the Diesel generator [kWh]
model.dg_fuel_consumption = Var(model.hours, model.dg, within=NonNegativeReals)  # diesel consumption [L]
model.dg_units_on = Var(model.hours, model.dg, within=NonNegativeIntegers)  # number of active DG in h

# Variables associated to the energy balance
model.load_hourly = Var(model.hours, within=NonNegativeReals)
model.lost_load = Var(model.hours, within=NonNegativeReals)  # Power not supplied by the system [kW]
model.load_total = Var(within=NonNegativeReals)  # Total energy requirement of the project [kWh]

# Variables associated to reserve needs
model.reserve = Var(model.hours, within=NonNegativeReals)  # total reserve needed per hour [kW]
model.reserve_dg = Var(model.hours, model.dg, within=NonNegativeReals)  # reserve provided by DG [kW]
model.reserve_bess = Var(model.hours, model.bess, within=NonNegativeReals)  # reserve provided by BESS [kW]



#--------------------------------OBJECTIVE FUNTION

def total_net_present_cost(model):
    return model.NPC == model.initial_investment + model.OM_cost + model.replacement_cost
           #+ model.salvage_value


def total_initial_investment(model):
    return model.initial_investment == 0.001 * (
                +sum(model.dg_units[g] * model.dg_investment_cost[g] for g in model.dg)
                + sum(model.pv_units[p] * model.pv_investment_cost[p] for p in model.pv)
                + sum(model.wt_units[w] * model.wt_investment_cost[w] for w in model.wt)
                + sum(model.bess_units[b] * model.bess_investment_cost[b] for b in model.bess)
        # +sum(model.bess_power_max[b]*model.inverter_cost for b in model.bess) \
    )


def total_replacement_cost(model):
    return model.replacement_cost == 0.001 * (
                +sum(sum(model.discount_rate[h] * model.bess_wear_cost[b] * model.bess_dis_power[h, b] for b in model.bess) for h in model.hours)
                + sum(sum(model.discount_rate[h] * model.dg_units_on[h, g] * model.dg_repl_cost[g] for g in model.dg) for h in model.hours)
        )


def total_OM_cost(model):
    return model.OM_cost == 0.001 * (
                + sum(sum(model.discount_rate[h] * model.fuel_cost * model.dg_fuel_consumption[h, g] for g in model.dg)
                     for h in model.hours)
                + sum(sum(model.discount_rate[h] * model.dg_repl_cost[g] * model.dg_units_on[h, g] for g in model.dg) for h in
                     model.hours)
#                + sum(sum(model.discount_rate[h] * model.pv_units[p] * model.pv_OM_cost[p] for p in model.pv) for h in
#                     model.hours_last)
#                + sum(sum(model.discount_rate[h] * model.wt_units[w] * model.wt_OM_cost[w] for w in model.wt) for h in
#                     model.hours_last)
#                + sum(sum(model.discount_rate[h] * model.bess_units[b] * model.bess_OM_cost[b] for b in model.bess) for h in
#                     model.hours_last)
                + sum(model.discount_rate[h]*model.lost_load[h]*model.lost_load_value for h in model.hours) \
    )

'''
def total_salvage_value(model):
    return model.salvage_value == 0.001 * (
            + sum(sum(model.discount_rate[h] * model.pv_units[p] * model.pv_investment_cost[p] * (
                    model.pv_life[p] - model.project_duration // 8760) / model.pv_life[p] for p in model.pv)for h in model.end_of_project)
            + sum(sum(model.discount_rate[h] * model.wt_units[w] * model.wt_investment_cost[w] * (
                    model.wt_life[w] - model.project_duration // 8760) / model.wt_life[w] for w in model.wt)for h in model.end_of_project)
            + sum(sum(model.discount_rate[h] * model.bess_power_max[b] * model.inverter_cost * (
                    model.inverter_life - model.project_duration // 8760) / model.inverter_life for b in model.bess)for h in model.end_of_project)
            )
'''


# ------------------------------CONSTRAINTS

# this constraint defines the hourly demand, subject to a linear growth rate starting from the second year
def total_load(model, h):
    if h <= 8760:
        return model.load_hourly[h] == model.input_load[h]
    else:
        return model.load_hourly[h] == model.load_hourly[h - 8760] * (1 + model.demand_growth)


# this group of constraints limits the number of units installed
def pv_installed(model, p):
    return model.pv_units[p] <= model.pv_max_units[p]


def wt_installed(model, w):
    return model.wt_units[w] <= model.wt_max_units[w]


def bess_installed(model, b):
    return model.bess_units[b] <= model.bess_max_units[b]


def dg_installed(model, g):
    return model.dg_units[g] <= model.dg_max_units[g]


# this constraint defines the maximum power produced by renewables
def res_energy(model, h):
    return model.total_power_res[h] <= sum(model.pv_units[p] * model.input_pv_prod[h, p] for p in model.pv) + sum(
        model.wt_units[w] * model.input_wt_prod[h, w] for w in model.wt)


# this constraints expresses the balance of the system
def system_balance(model, h):
    return model.load_hourly[h] == model.total_power_res[h] + sum(model.dg_power[h, g] for g in model.dg) + model.lost_load[h] \
           + sum(model.bess_dis_power[h, b] * model.bess_dis_efficiency_nom[b] - model.bess_ch_power[h, b] / model.bess_ch_efficiency_nom[b]
        for b in model.bess)


# these constraints define the total energy requirement along the project and the maximum allowable unmet portion
def total_energy_req(model):
    return model.load_total == sum(model.load_hourly[h] for h in model.hours)


def total_lost_load(model):
    return sum(model.lost_load[h] for h in model.hours) <= model.load_total * model.lost_load_max


# these constraints set the reserve requirement and its allocation among DG and BESS
def total_reserve_req(model, h):
    return model.reserve[h] == model.load_forecast_error * model.load_hourly[h] + model.pv_forecast_error * sum(
        model.pv_units[p] * model.input_pv_prod[h, p] for p in model.pv) + model.wt_forecast_error * sum(
        model.wt_units[w] * model.input_wt_prod[h, w] for w in model.wt)


def reserve_allocation(model, h):
    return sum(model.reserve_dg[h, g] for g in model.dg) + sum(model.reserve_bess[h, b] for b in model.bess) >= \
           model.reserve[h]


#########################################
# constraints related to diesel generators

def fuel_consumption_curve(model, h, g):  # linear characteristic
    return model.dg_fuel_consumption[h, g] == model.dg_cost_coeff_A[g] * model.dg_units_on[h, g] + model.dg_cost_coeff_B[g] * \
           model.dg_power[h, g]


def dg_power_max(model, h, g):
    # with reserve
    return model.dg_power[h, g] + model.reserve_dg[h, g] <= model.dg_nominal_capacity[g] * model.dg_units_on[h, g]
    # without reserve
    # return model.dg_power[h,g] <= model.dg_nominal_capacity[g]*model.dg_units_on[h,g]


def dg_power_min(model, h, g):
    return model.dg_power[h, g] >= model.dg_nominal_capacity[g] * model.dg_P_min[g] * model.dg_units_on[h, g]


def dg_online(model, h, g):
    return model.dg_units_on[h, g] <= model.dg_units[g]


###########################################
# constraints related to batteries
def battery_power_max(model, h, b):  # maximum power flowing through batteries, to size converters
    return model.bess_dis_power[h, b] * model.bess_dis_efficiency_nom[b] + model.bess_ch_power[h, b] / \
           model.bess_ch_efficiency_nom[b] <= model.bess_power_max[b]


# following two constraints to avoid charging and discharging at the same time
def bess_condition1(model, h, b):
    return model.bess_ch_power[h, b] <= model.bess_bin[h, b] * model.M


def bess_condition2(model, h, b):
    return model.bess_dis_power[h, b] <= (1 - model.bess_bin[h, b]) * model.M


def bess_charging_level(model, h, b):
    if h == 1:
        return model.bess_total_energy[h, b] == model.bess_units[b] * model.bess_nominal_capacity[b] * \
               model.bess_initial_SOC[b] + model.bess_ch_power[h, b] - model.bess_dis_power[h, b]
    elif h < model.project_duration:
        return model.bess_total_energy[h, b] == model.bess_total_energy[h - 1, b] + model.bess_ch_power[h, b] - \
               model.bess_dis_power[h, b]
    else:
        return model.bess_total_energy[h, b] == 0.5 * model.bess_units[b] * model.bess_nominal_capacity[b] * \
               model.bess_initial_SOC[b]


def bess_charging_level_min(model, h, b):
    return model.bess_total_energy[h, b] >= model.bess_units[b] * model.bess_nominal_capacity[b] * (
                1 - model.bess.depth_of_discharge[b]) + model.reserve_bess[h, b]


def bess_charging_level_max(model, h, b):
    return model.bess_total_energy[h, b] <= model.bess_units[b] * model.bess_nominal_capacity[b]


# maximum charging and discharging power depending on the P-ratio
def bess_ch_power_max(model, h, b):
    return model.bess_ch_power[h, b] <= model.bess_units[b] * model.bess_nominal_capacity[b] * model.bess_P_ratio_max[b]


def bess_dis_power_max(model, h, b):
    return model.bess_dis_power[h, b] <= model.bess_units[b] * model.bess_nominal_capacity[b] * model.bess_P_ratio_max[b]



# OBJETIVE FUNTION:
model.ObjectiveFuntion = Objective(rule=total_net_present_cost, sense=minimize)

# CONSTRAINTS
# to compute OF
model.TotalInitialInvestment = Constraint(rule=total_initial_investment)
model.TotalReplacementCost = Constraint(rule=total_replacement_cost)
model.TotalOMCost = Constraint(rule=total_OM_cost)
#model.TotalSalvageValue = Constraint(rule=total_salvage_value)
# to design the system
model.TotalLoad = Constraint(model.hours, rule=total_load)
model.PvInstalled = Constraint(model.pv, rule=pv_installed)
model.WtInstalled = Constraint(model.wt, rule=wt_installed)
model.BessInstalled = Constraint(model.bess, rule=bess_installed)
model.DgInstalled = Constraint(model.dg, rule=dg_installed)
model.ResEnergy = Constraint(model.hours, rule=res_energy)
model.SystemBalance = Constraint(model.hours, rule=system_balance)
model.TotalEnergyReq = Constraint(rule=total_energy_req)
model.TotalLostLoad = Constraint(rule=total_lost_load)
model.TotalReserveReq = Constraint(model.hours, rule=total_reserve_req)
model.ReserveAllocation = Constraint(model.hours, rule=reserve_allocation)
# constraints related to diesel generators
model.FuelConsumptionCurve = Constraint(model.hours,model.dg, rule=fuel_consumption_curve)
model.DgPowerMax = Constraint(model.hours, model.dg, rule=dg_power_max)
model.DgPowerMin = Constraint(model.hours, model.dg, rule=dg_power_min)
model.DgOnline = Constraint(model.hours, model.dg, rule=dg_online)
# constraints related to batteries
model.BatteryPowerMax = Constraint(model.hours,model.bess, rule=battery_power_max)
model.BessCondition1 = Constraint(model.hours, model.hours, rule=bess_condition1)
model.BessCondition2 = Constraint(model.hours, model.hours, rule=bess_condition2)
model.BessChargingLevel = Constraint(model.hours,model.bess, rule=bess_charging_level)
model.BessChargingLevelMin = Constraint(model.hours,model.bess, rule=bess_charging_level_min)
model.BessChargingLevelMax = Constraint(model.hours,model.bess, rule=bess_charging_level_max)
model.BessChPowerMax = Constraint(model.hours,model.bess, rule=bess_ch_power_max)
model.BessDisPowerMax = Constraint(model.hours, model.bess, rule=bess_dis_power_max)

instance = model.create_instance(data) # load parameters
opt = SolverFactory('cplex', executable='C:\Program Files\cplex_studio1210.win-x86-64')  # Solver use during the optimization
#model.parameters.mip.tolerances.mipgap(float(0.1))
results = opt.solve(instance, options="threads=10", tee=True)  # Solving a model instance
instance.solutions.load_from(results)  # Loading solution into instance
