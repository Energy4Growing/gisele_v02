from pyomo.environ import AbstractModel, Var, value

from gisele.michele.components_creation import Model_Creation
#from gisele.michele.model_solve import Model_Instance
from gisele.michele.model_solve import Model_Resolution
from gisele.michele.results import Load_results
from gisele.michele.components_initialization import importing
from pyomo.opt import SolverFactory
import pandas as pd

def start(load_profile, pv_avg, wt_avg,input_michele, ht_avg,n_mg):
    input_load, wt_prod, pv_prod = importing(load_profile, pv_avg, wt_avg, ht_avg)

    model = AbstractModel()  # define type of optimization problem

    # Optimization model
    print('Starting model creation')
    Model_Creation(model, input_load, wt_prod, pv_prod,input_michele)  # Creation of the Sets, parameters and variables.
    print('Starting model resolution')

    instance = model.create_instance(
        r'gisele\michele\Inputs\data.json')  # load parameters

    # opt = SolverFactory('cplex',executable=r'C:\Program Files\IBM\ILOG\CPLEX_Studio1210\cplex\bin\x64_win64\cplex')
    opt = SolverFactory('gurobi')  # Solver use during the optimization
    opt.options['mipgap'] = 0.05
    print('Begin Optimization')
    result_michele = {}

    if n_mg==1:

        print('Solving for minimum renewable fraction of '+str(instance.ren_fraction.extract_values()[None]*100)+'%')
        results = opt.solve(instance, tee=True)  # Solving a model instance
        instance.solutions.load_from(results)  # Loading solution into instance
        result_michele['0'] = Load_results(instance)

    else:

        print('Solving for '+str(n_mg)+' different values of minimum renewable fraction, ranging from 0% to 100%')

        renfr=-1  # initialize at a value lower than 0

        for i in range(n_mg):
            min_renfr = i*1/(n_mg-1)
            instance.ren_fraction=min_renfr
            print('Solving for minimum renewable fraction of '+str(min_renfr*100)+'%')

            if renfr >= min_renfr:
                print('Skipping solve as minimum renewable fraction already fulfilled')

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
                renfr = 1-dg_tot/load_tot
                print('Renewable fraction='+str(renfr*100)+'%')

            print('Show results')
            result_michele[str(i)]= Load_results(instance)

    return result_michele



