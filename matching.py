import xtrack as xt
import xobjects as xo
import numpy as np
import config as cfg

ctx = xo.ContextCpu()
env = xt.Environment()

# Carbon-ion reference particle
ion = xt.Particles(_context=ctx, mass0=cfg.ION_MASS_EV, q0=cfg.ION_Z, p0c=cfg.P0C)
env.set_particle_ref(ion)

# Drift lengths
env.vars['l1'] = 1.6113
env.vars['l2'] = 0.4265
env.vars['l3'] = 1.0508
env.vars['l4'] = 0.6859
env.vars['l5'] = 0.5551
env.vars['l6'] = 0.6388
env.vars['l7'] = 0.4515
env.vars['kqf1']= cfg.K1_QF1
env.vars['kqf2']= cfg.K1_QF2
env.vars['kqd1']= cfg.K1_QD1
env.vars['kqd2']= cfg.K1_QD2

# Drifts
env.new('d_str', xt.Drift, length=env.vars['l1'])
env.new('d1', xt.Drift, length=env.vars['l2'])
env.new('d2', xt.Drift, length=env.vars['l3'])
env.new('d3', xt.Drift, length=env.vars['l4'])
env.new('d4', xt.Drift, length=env.vars['l5'])
env.new('d5', xt.Drift, length=env.vars['l6'])
env.new('d6', xt.Drift, length=env.vars['l7'])

# Magnets
env.new('mb', xt.Bend, length=cfg.L_DIPOLE_M, angle=np.deg2rad(cfg.BEND_ANGLE_PER_DIPOLE_DEG))
env.new('mq', xt.Quadrupole, length=cfg.L_QUAD_M)

# Quadrupole strengths
env.new('qf1', 'mq', k1=env.vars['kqf1'])
env.new('qf2', 'mq', k1=env.vars['kqf2'])
env.new('qd1', 'mq', k1=env.vars['kqd1'])
env.new('qd2', 'mq', k1=env.vars['kqd2'])

half_arc = [
    xt.Drift(length=1.6113),
    'qf1',
    xt.Drift(length=0.4265),
    'qd1',
    xt.Drift(length=1.0508),
    'mb',
    xt.Drift(length=0.6859),
    'mb',
    xt.Drift(length=0.5551),
    'qd2',
    xt.Drift(length=0.6388),
    'mb',
    xt.Drift(length=0.4515),
    'qf2'
]

full_dba = [
    'd_str', 'qf1', 'd1', 'qd1', 'd2', 'mb', 'd3', 'mb', 'd4', 'qd2', 'd5', 'mb', 'd6', 'qf2',
    'd6', 'mb', 'd5', 'qd2', 'd4', 'mb', 'd3', 'mb', 'd2', 'qd1', 'd1', 'qf1', 'd_str'
]
line_dba = env.new_line(components=full_dba)
ring = env.new_line(components=[line_dba, line_dba])
# Matching in DBA
opt = ring.match(
    method='4d',
    vary=[
        xt.Vary('kqf1', limits=(0.8, 1.3)),
        xt.Vary('kqf2', limits=(0.6, 1.1)),
        xt.Vary('kqd1', limits=(-1.6, -1.0)),
        xt.Vary('kqd2', limits=(-0.9, -0.4)),
    ],
    targets=[
        #Tunes for DBA
        xt.Target('qx', value=cfg.QX_TARGET, tol=1e-5),
        xt.Target('qy', value=cfg.QY_TARGET, tol=1e-5),

        xt.Target('dx',  value=0.0, at=xt.START, tol=1e-3),
        xt.Target('dpx', value=0.0, at=xt.START, tol=1e-3),
    ]

)

opt.show()
print(env.vars.get_table())