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
from pymoo.visualization.scatter import Scatter
from pymoo.indicators.hv import HV
from pymoo.visualization.pcp import PCP

#For tunes Qx, Qy
from SAA_optimisation import third_int

#from SAA_optimisation.py
data = np.load("sa_results.npz")
sa_obj = data["sa_objectives"]

class LatticeOptimisationProblem(ElementwiseProblem):
    def __init__(self, set_vars, get_metrics, lo, hi, ref, tol_qx=0.01, tol_qy=0.01):
        super().__init__(n_var=4, n_obj=4, n_constr=2, xl=lo, xu=hi)
        self.set_vars = set_vars
        self.get_metrics = get_metrics
        self.ref = ref
        self.tol_qx = tol_qx
        self.tol_qy = tol_qy

    def _evaluate(self, x, out, *args, **kwargs):
        self.set_vars(x)
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
            m['betx_max'] / self.ref['betx_max'],
            m['bety_max'] / self.ref['bety_max'],
            m['dx_max'] / self.ref['dx_max'],
            m['dx_start_zero'] / self.ref['dx_start_zero'],
        )

        out["F"] = [phi[0], phi[1], phi[2], phi[3]]

        #Constraints
        c = {
                'qx': self.tol_qx - third_int(m['qx']),
                'qy': self.tol_qy - third_int(m['qy'])
            }
        out["G"] = [c['qx'], c['qy']]

#Main program
if __name__ == '__main__':
    quad_names = ['qf1', 'qd1', 'qf2', 'qd2']
    
    lo = np.array([ 0.1, -2.0,  0.1, -2.0])
    hi = np.array([ 2.0, -0.1,  2.0, -0.1])
    ref=dict(betx_max=30.0, bety_max=12.0, dx_max=13.0, dx_start_zero=1e-3, qx=1.672, qy=1.72)

#Same as in SAA
    def set_vars(x):
        for name, val in zip(quad_names, x):
            ring[name].k1 = float(val)

    def get_metrics():
        tw = ring.twiss(method='4d')
        return dict(
            stable=True,
            betx_max=float(np.max(tw.betx)),
            bety_max=float(np.max(tw.bety)),
            dx_max=float(np.max(np.abs(tw.dx))),
            dx_start_zero=float(np.abs(tw.dx[0])),
            qx=float(tw.qx),
            qy=float(tw.qy),
        )
    
    problem = LatticeOptimisationProblem(set_vars, get_metrics, lo, hi, ref, tol_qx=0.01, tol_qy=0.01)

    algorithm = NSGA2(pop_size=40)
    termination = get_termination("n_eval", 3000)

    res = minimize(problem,
                   algorithm,
                   termination,
                   seed=0,
                   save_history=True,
                   verbose=True)

    print(res.X)
    print(res.F)

    n_evals = []
    hv_vals = []
    best_J_nsga = []

    ref_point = np.array([10.0, 10.0, 10.0, 10.0])
    ind_hv = HV(ref_point=ref_point)

    weights_o = np.array([0.20, 0.50, 0.15, 0.15])
    weights_c = np.array([0.5, 0.5])

    for entry in res.history:
        n_evals.append(entry.evaluator.n_eval)

        F_pop = entry.pop.get("F")
        G_pop = entry.pop.get("G")

        cons = np.sum(weights_c * np.maximum(0.0, G_pop), axis=1)
        objs = np.sum(weights_o * F_pop, axis=1)
        J_pop = objs + cons
        best_J_nsga.append(np.min(J_pop))

        feas_mask = np.all(G_pop <= 0, axis=1)
        if np.any(feas_mask):
            F_feas = F_pop[feas_mask]
            valid = np.all(F_feas < ref_point, axis=1)
            hv_vals.append(ind_hv(F_feas[valid]) if np.any(valid) else 0.0)
        else:
            hv_vals.append(0.0)

    plt.figure(figsize=(8, 4.5))
    plt.plot(n_evals, hv_vals, marker='o', color='tab:purple', linewidth=2)
    plt.xlabel("Evaluations")
    plt.ylabel("Hypervolume (HV)")
    plt.title("NSGA-II Convergence: Hypervolume over Evaluations")
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.tight_layout()
    plt.show()

    feas_final = np.all(res.G <= 0, axis=1) if res.G is not None else np.ones(len(res.F), dtype=bool)
    F_pareto = res.F[feas_final]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    ax1.scatter(F_pareto[:, 0], F_pareto[:, 1], color='tab:blue', alpha=0.8, label='NSGA-II Pareto Front')
    ax1.scatter(sa_obj[0], sa_obj[1], color='red', marker='*', s=200, zorder=5, label='SA Best Solution')
    ax1.axvline(1.0, color='gray', linestyle=':', alpha=0.7)
    ax1.axhline(1.0, color='gray', linestyle=':', alpha=0.7)
    ax1.set_xlabel(r'$\beta_{x,\max} / \mathrm{ref}$')
    ax1.set_ylabel(r'$\beta_{y,\max} / \mathrm{ref}$')
    ax1.set_title(r'Trade-off: $\beta_{y,\max}$ vs $\beta_{x,\max}$')
    ax1.legend()
    ax1.grid(True, linestyle='--', alpha=0.6)

    ax2.scatter(F_pareto[:, 2], F_pareto[:, 3], color='tab:green', alpha=0.8, label='NSGA-II Pareto Front')
    ax2.scatter(sa_obj[2], sa_obj[3], color='red', marker='*', s=200, zorder=5, label='SA Best Solution')    
    ax2.axvline(1.0, color='gray', linestyle=':', alpha=0.7)
    ax2.axhline(1.0, color='gray', linestyle=':', alpha=0.7)
    ax2.set_xlabel(r'$D_{x,\max} / \mathrm{ref}$')
    ax2.set_ylabel(r'$D_{x,0} / \mathrm{ref}$')
    ax2.set_title(r'Trade-off: $D_{x,0}$ vs $D_{x,\max}$')
    ax2.legend()
    ax2.grid(True, linestyle='--', alpha=0.6)

    plt.tight_layout()
    plt.show()

    if len(F_pareto) > 0:
        pcp = PCP(
            labels=[r"$\beta_{x,\max}$", r"$\beta_{y,\max}$", r"$D_{x,\max}$", r"$D_{x,0}$"],
            normalize_each_axis=True
        )
        pcp.set_axis_style(color="grey", alpha=0.5)
        pcp.add(F_pareto, color="tab:blue", alpha=0.5)
        best_idx = np.argmin(np.sum(weights_o * F_pareto, axis=1))
        pcp.add(F_pareto[best_idx], color="red", linewidth=2.5, label="Min Weighted J in Pareto")
        pcp.do()
        pcp.fig.suptitle("Parallel Coordinates of Pareto-Optimal Solutions", y=0.98, fontsize=13)
        pcp.fig.subplots_adjust(top=0.82)
        pcp.show()

    plt.figure(figsize=(8.5, 5))
    plt.plot(n_evals, best_J_nsga, marker='s', markersize=4, linewidth=2, color='black', label='NSGA-II (Seed 0, best J in pop)')

    plt.yscale("log")
    plt.xlabel("Function Evaluations")
    plt.ylabel("Best-so-far Equivalent Cost J")
    plt.title("Benchmarking: Simulated Annealing vs NSGA-II (3000 Evaluations)")
    plt.legend()
    plt.grid(True, which="both", linestyle='--', alpha=0.6)
    plt.tight_layout()
    plt.show()