from pyomo.environ import Objective, minimize, Constraint
from pyomo.opt import SolverFactory
from pyomo.dataportal.DataPortal import DataPortal
import pandas as pd


#def Model_Instance(model):

def Model_Resolution(model):

    instance = model.create_instance(r'gisele\michele\Inputs\data.json')  # load parameters
    # opt = SolverFactory('cplex',executable=r'C:\Program Files\IBM\ILOG\CPLEX_Studio1210\cplex\bin\x64_win64\cplex')
    opt = SolverFactory('gurobi')  # Solver use during the optimization
    opt.options['mipgap'] = 0.05

    print('Begin Optimization')

    n_ren=3
    print('Solving for '+str(n_ren)+' different values of minimum renewable fraction, ranging from 0% to 100%')

    renfr=-1  # initialize at a value lower than 0

    #instances={}

    for i in range(n_ren):
        min_renfr=i*1/(n_ren-1)
        instance.ren_fraction=min_renfr
        print('Solving for minimum renewable fraction of '+str(min_renfr*100)+'%')
        if renfr>=min_renfr:
            print('Skipping solve as minimum renewable fraction already fulfilled')
            #instances[i]=instances[i-1]
        else:
            results = opt.solve(instance, tee=True)  # Solving a model instance
            instance.solutions.load_from(results)  # Loading solution into instance
            # Compute renewable fraction
            project_duration = int(instance.project_duration.extract_values()[None])
            h_weight = int(instance.h_weight.extract_values()[None])
            dg_power = pd.DataFrame.from_dict(instance.dg_power.get_values(), orient='index')
            dg_tot = 0.001 * h_weight * sum(dg_power.iloc[h, 0] for h in range(project_duration))
            load = pd.DataFrame.from_dict(instance.Load.get_values(), orient='index')
            load_tot = 0.001 * h_weight * sum(load.iloc[h, 0] for h in range(project_duration))
            renfr=1-dg_tot/load_tot
            print('Renewable fraction='+str(renfr*100)+'%')
            #instances[i]=instance
    return instance



