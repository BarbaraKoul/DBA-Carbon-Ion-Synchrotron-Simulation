import xtrack as xt
import xobjects as xo
import numpy as np
import config as cfg
import math
import random
from lattice_design import ring

J_UNSTABLE=1e3

#computing the distance of Q from the nearest n+-1/3
def third_int(q):
    frac=q-np.floor(q)
    return min(abs(frac-1.0/3.0),abs(frac-2.0/3.0))

#encoder function
class Encoder:
    def __init__(self, ref, weights=(0.0, 0.0, 0.5, 0.5), w_con=None, power=2.0, mu=100.0, tol_q=0.01, extra=None):
        self.ref = ref
        self.weights = weights
        self.w_con = w_con or {}
        self.power = power
        self.mu = mu
        self.tol_q = tol_q
        self.extra = extra or {}

    def __call__(self,m):
        if not m.get('stable', False):
            return J_UNSTABLE, dict(unstable=True)

        phi= (
            m['circ'] / self.ref['circ'],
            m['n_magnets'] / self.ref['n_magnets'],
            m['bety_max'] / self.ref['bety_max'],
            m['dx_max'] / self.ref['dx_max'],
        )

        obj = sum(w * p  for w, p in zip(self.weights, phi))

        v = {'c4_tune': max(0.0, third_int(m['qx']) / self.tol_q - 1.0)}
        for name, (key, limit) in self.extra.items():
            v[name] = max(0.0, m[key] / limit - 1.0)

        pen = sum(self.w_con.get(k, self.mu) * (x ** self.power) for k, x in v.items())
        return obj + pen, dict(phi=phi, objective=obj, violations=v, penalty=pen)

class ScalarEvaluator:
    def __init__(self, cost, set_vars, get_metrics, lo, hi, budget):
        self.cost, self.set_vars, self.get_metrics = cost, set_vars, get_metrics
        self.lo, self.hi, self.budget = np.asarray(lo, float), np.asarray(hi, float), budget
        self.n, self.best_J, self.best_u, self.best_parts, self.history = 0, np.inf, None, None, []

#Count how many times the evaluator has been called, and stop when the budget is reached
    def __call__(self, u):
        if self.n >= self.budget:
            raise StopIteration
        self.n += 1
        u = np.clip(np.asarray(u, float), 0.0, 1.0)
        self.set_vars(self.lo + u * (self.hi - self.lo))
        try:
            m = self.get_metrics()
        except Exception:
            m = dict(stable=False)
        J, parts = self.cost(m)
        if J < self.best_J:
            self.best_J, self.best_u, self.best_parts = J, u.copy(), dict(parts, metrics=m)
        self.history.append(self.best_J)
        return J

#Function that keeps the variables within the range [0,1]
def _reflect(u):
    u = np.abs(u)
    u = np.where(u > 1.0, 2.0 - u, u)
    return np.clip(u, 0.0, 1.0)


def anneal(ev, u0, rng, T0=0.3, T_end=0.002, sigma0=0.2, sigma_min=1e-3, sigma_max=0.5,
           adapt_every=20, target_acc=(0.2, 0.5)):
    n = len(u0)
    u = np.array(u0, float)
    E = np.log(ev(u))
    adapt, n_prop, n_acc = 1.0, 0, 0
    try:
        while True:
            f = min(ev.n / ev.budget, 1.0)
            T = T0 * (T_end / T0) ** f
            sigma = float(np.clip(sigma0 * (T / T0) ** 0.5 * adapt, sigma_min, sigma_max))
            un = _reflect(u + sigma * rng.standard_normal(n))
            En = np.log(ev(un))
            n_prop += 1
            if En <= E or rng.random() < np.exp(-(En - E) / T):
                u, E = un, En
                n_acc += 1
            if n_prop == adapt_every:
                rate = n_acc / n_prop
                adapt = adapt * 1.3 if rate > target_acc[1] else (adapt / 1.3 if rate < target_acc[0] else adapt)
                n_prop = n_acc = 0
    except StopIteration:
        pass
    return ev.best_u, ev.best_J

if __name__ == '__main__':

    quad_names = ['qf1', 'qd1', 'qf2', 'qd2']

    k_init = np.array([ring[name].k1 for name in quad_names])
    lo = k_init - 1.5
    hi = k_init + 1.5

    def set_vars(x):
        for name, val in zip(quad_names, x):
            ring[name].k1 = float(val)

    def get_metrics():
        tw = ring.twiss(method='4d')
        return dict(
            stable=True,
            circ=float(ring.get_length()),
            n_magnets=4,
            bety_max=float(np.max(tw.bety)),
            dx_max=float(np.max(np.abs(tw.dx))),
            dx_start=float(np.abs(tw.dx[0])),
            qx=float(tw.qx),
        )

    #Cost
    cost = Encoder(
        ref=dict(circ=ring.get_length(), n_magnets=4, bety_max=15.0, dx_max=1.0),
        weights=(0.0, 0.0, 2, 1),
        w_con={'c4_tune': 100.0, 'dx_zero_start': 100.0},
        power=1.0,
        tol_q=0.01,
        extra={'dx_zero_start': ('dx_start', 1e-3)}
    )

    #4 seeds
    best_overall_x = None
    best_overall_J = np.inf
    best_overall_parts = None

    for seed in range(4):
        rng = np.random.default_rng(seed)
        ev = ScalarEvaluator(cost, set_vars, get_metrics, lo=lo, hi=hi, budget=800)
        u0 = rng.random(4)
        u_opt, J_opt = anneal(ev, u0, rng)
        x_opt = ev.lo + u_opt * (ev.hi - ev.lo)
        print(f"Seed {seed}: J = {J_opt:.4f} | k1 = {np.round(x_opt, 4)}")

        if J_opt < best_overall_J:
            best_overall_J = J_opt
            best_overall_x = x_opt
            best_overall_parts = ev.best_parts

    #Optimum solution
    set_vars(best_overall_x)
    print("\n Optimum solution:")
    print("k1 [m^-2]:", np.round(best_overall_x, 5))
    print("Cost J :", round(best_overall_J, 5))
    print("Analysis  :", best_overall_parts)