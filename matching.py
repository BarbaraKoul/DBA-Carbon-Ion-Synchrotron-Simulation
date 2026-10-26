import xtrack as xt
import xobjects as xo
import numpy as np
import config as cfg

ctx = xo.ContextCpu()
env = xt.Environment()

def gamma_t_from(tw):
    try:
        a = tw.momentum_compaction_factor
        if a is not None and np.isfinite(a) and a != 0:
            return 1 / np.sqrt(abs(a))
    except Exception:
        pass
    L = cfg.L_DIPOLE_M
    th = np.deg2rad(cfg.BEND_ANGLE_PER_DIPOLE_DEG)
    rho = L / th
    I = 0.0
    for i, n in enumerate(tw.name):
        if str(n).startswith('mb'):
            D0, Dp0 = tw.dx[i], tw.dpx[i]
            I += (D0 * rho * np.sin(th) + rho**2 * Dp0 * (1 - np.cos(th))
                  + rho * (L - rho * np.sin(th))) / rho
    return 1 / np.sqrt(I / tw.s[-1])

# Carbon-ion reference particle
ion = xt.Particles(_context=ctx, mass0=cfg.ION_MASS_EV, q0=cfg.ION_Z, p0c=cfg.P0C)
env.set_particle_ref(ion)
BETA0 = float(ion.beta0[0])

# Drift lengths
env.vars['l1'] = 1.978
env.vars['l2'] = 0.403
env.vars['l3'] = 0.382
env.vars['l4'] = 0.367
env.vars['l5'] = 1.477
env.vars['l6'] = 0.412
env.vars['l7'] = 0.404
env.vars['kqf1']= cfg.K1_QF1
env.vars['kqf2']= cfg.K1_QF2
env.vars['kqd1']= cfg.K1_QD1
env.vars['kqd2']= cfg.K1_QD2

# Drifts
env.new('d_str', xt.Drift, length='l1')
env.new('d1', xt.Drift, length='l2')
env.new('d2', xt.Drift, length='l3')
env.new('d3', xt.Drift, length='l4')
env.new('d4', xt.Drift, length='l5')
env.new('d5', xt.Drift, length='l6')
env.new('d6', xt.Drift, length='l7')

# Magnets
env.new('mb', xt.Bend, length=cfg.L_DIPOLE_M, angle=np.deg2rad(cfg.BEND_ANGLE_PER_DIPOLE_DEG), k0_from_h=True)
env.new('mq', xt.Quadrupole, length=cfg.L_QUAD_M)

# Quadrupole strengths
env.new('qf1', 'mq', k1='kqf1')
env.new('qf2', 'mq', k1='kqf2')
env.new('qd1', 'mq', k1='kqd1')
env.new('qd2', 'mq', k1='kqd2')


full_dba = [
    'd_str', 'qf1', 'd1', 'qd1', 'd2', 'mb', 'd3', 'mb', 'd4', 'qd2', 'd5', 'mb', 'd6', 'qf2',
    'd6', 'mb', 'd5', 'qd2', 'd4', 'mb', 'd3', 'mb', 'd2', 'qd1', 'd1', 'qf1', 'd_str'
]

ring = env.new_line(components=full_dba * 2)
# Matching in DBA
opt = ring.match(chrom=False,
    method='4d',
    vary=[
        xt.Vary('l1', step=1e-4, limits=(1.8, 2.2)),
        xt.Vary('l2', step=1e-4, limits=(0.30, 0.50)),
        xt.Vary('l3', step=1e-4, limits=(0.30, 0.50)),
        xt.Vary('l4', step=1e-4, limits=(0.30, 0.50)),
        xt.Vary('l5', step=1e-4, limits=(1.30, 1.70)),
        xt.Vary('l6', step=1e-4, limits=(0.30, 0.50)),
        xt.Vary('l7', step=1e-4, limits=(0.30, 0.50)),
    ],
    targets=[
        #Tunes for DBA
        xt.Target('qx', cfg.QX_TARGET, tol=2e-3),
        xt.Target('qy', cfg.QY_TARGET, tol=2e-3),
        xt.Target(lambda tw: tw.s[-1], cfg.CIRCUMFERENCE_M, tol=3e-2),
        xt.Target(lambda tw: tw.betx.max(), cfg.BETA_X_MAX_M, tol=5e-2),
        xt.Target(lambda tw: tw.bety.max(), cfg.BETA_Y_MAX_M, tol=5e-2),
        xt.Target(lambda tw: tw.dx.max() / BETA0, cfg.DISPERSION_X_MAX_M, tol=5e-2),
        xt.Target(gamma_t_from, cfg.TRANSITION_GAMMA, tol=3e-3),
    ]

)

opt.show()
opt.target_status()
print(env.vars.get_table())