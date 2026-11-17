#The same libraries we had for SAA
import xtrack as xt
import xobjects as xo
import numpy as np
import config as cfg
import matplotlib.pyplot as plt
from lattice_design import ring

#NSGA2 using pymoo library
from pymoo.core.problem import ElementwiseProblem
from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.optimize import minimize
from pymoo.termination import get_termination

#For tunes Qx, Qy
from SAA_optimisation import third_int

class LatticeOptimisationProblem(ElementwiseProblem):
    def __init__(self, set_vars, get_metrics, lo, hi, ref, tol_qx=1e-3, tol_qy=1e-3):
        super().__init__(n_var=4, n_obj=4, n_constr=2, xl=lo, xu=hi)
        self.set_vars = set_vars
        self.get_metrics = get_metrics
        self.ref = ref
        self.tol_qx = tol_qx
        self.tol_qy = tol_qy

    def _evaluate(self, x, out, *args, **kwargs):
        self.set_vars(x)
        metrics = self.get_metrics()
        try:
            m=self.get_metrics()
        except Exception:
            m=dict(stable=False)

        if not m.get('stable', False):
            out["F"] = [1e6, 1e6, 1e6, 1e6]
            out["G"] = [1e6, 1e6]
            return

        #Objectives
        phi = (
            m['betax_max'] / self.ref['betax_max'],
            m['bety_max'] / self.ref['bety_max'],
            m['dx_max'] / self.ref['dx_max'],
            m['dx_start_zero'] / self.ref['dx_start_zero'],
        )

        out["F"] = [phi[0], phi[1], phi[2], phi[3]]

        #Constraints
        c = {
                'qx': max(0.0, third_int(m['qx']) / self.tol_qx - 1.0),
                'qy': max(0.0, third_int(m['qy']) / self.tol_qy - 1.0)
            }
        out["G"] = [c['qx'], c['qy']]

#Main program
if __name__ == '__main__':
    quad_names = ['qf1', 'qd1', 'qf2', 'qd2']
    
    k_init = np.array([ring[name].k1 for name in quad_names])
    lo = np.array([ 0.1, -2.0,  0.1, -2.0])
    hi = np.array([ 2.0, -0.1,  2.0, -0.1])
    ref=dict(betax_max=30.0, bety_max=12.0, dx_max=13.0, dx_start_zero=1e-3, qx=1.672, qy=1.72)

#Same as in SAA
    def set_vars(x):
        for name, val in zip(quad_names, x):
            ring[name].k1 = float(val)

    def get_metrics():
        tw = ring.twiss(method='4d')
        return dict(
            stable=True,
            betax_max=float(np.max(tw.betx)),
            bety_max=float(np.max(tw.bety)),
            dx_max=float(np.max(np.abs(tw.dx))),
            dx_start_zero=float(np.abs(tw.dx[0])),
            qx=float(tw.qx),
            qy=float(tw.qy),
        )
    problem = LatticeOptimisationProblem(set_vars, get_metrics, lo, hi, ref, tol_qx=1e-3, tol_qy=1e-3)

    algorithm = NSGA2(pop_size=40)
    termination = get_termination("n_eval", 800)

    res = minimize(problem,
                   algorithm,
                   termination,
                   seed=0,
                   verbose=True)

    print(res.X)
    print(res.F)

