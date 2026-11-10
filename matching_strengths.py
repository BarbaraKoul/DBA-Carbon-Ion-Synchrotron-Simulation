import xtrack as xt
import xobjects as xo
import numpy as np
import config as cfg

ctx = xo.ContextCpu()
env = xt.Environment()

def frac(k, f=0.4):
    return tuple(sorted((k * (1 - f), k * (1 + f))))


# Carbon-ion reference particle
ion = xt.Particles(_context=ctx, mass0=cfg.ION_MASS_EV, q0=cfg.ION_Z, p0c=cfg.P0C)
env.set_particle_ref(ion)

env.vars['kqf1']= 1.1
env.vars['kqf2']= 0.8
env.vars['kqd1']= -1.4
env.vars['kqd2']= -0.5

#magnets
env.new('mb', xt.Bend, length=cfg.L_DIPOLE_M, angle=np.deg2rad(cfg.BEND_ANGLE_PER_DIPOLE_DEG),k0_from_h=True)
env.new('mq', xt.Quadrupole, length=cfg.L_QUAD_M)

#quadruple strengths
env.new('qf1', 'mq', k1='kqf1')
env.new('qf2', 'mq', k1='kqf2')
env.new('qd1', 'mq', k1='kqd1')
env.new('qd2', 'mq', k1='kqd2')

dba=env.new_line(length=27.501, components=[
    env.place('qf1', anchor='start', at=1.98296),
    env.place('qd1', anchor='start', at=0.403088, from_='qf1@end'),
    env.place('mb', anchor='start', at=0.399686, from_='qd1@end'),
    env.place('mb', anchor='start', at=0.339825, from_='mb@end'),
    env.place('qd2', anchor='start', at=1.48352, from_='mb@end'),
    env.place('mb', anchor='start', at=0.422807, from_='qd2@end'),
    env.place('qf2', anchor='start', at=0.390371, from_='mb@end'),
    env.place('mb', anchor='start', at=0.390371, from_='qf2@end'),
    env.place('qd2', anchor='start', at=0.422807, from_='mb@end'),
    env.place('mb', anchor='start', at=1.48352, from_='qd2@end'),
    env.place('mb', anchor='start', at=0.339825, from_='mb@end'),
    env.place('qd1', anchor='start', at=0.399686, from_='mb@end'),
    env.place('qf1', anchor='start', at=0.403088, from_='qd1@end'),
])

ring=2*dba

opt = ring.match(chrom=False,
    method='4d',
    vary=[xt.Vary('kqf1', step=1e-5, limits=frac(cfg.K1_QF1)),
          xt.Vary('kqf2', step=1e-5, limits=frac(cfg.K1_QF2)),
          xt.Vary('kqd1', step=1e-5, limits=frac(cfg.K1_QD1)),
          xt.Vary('kqd2', step=1e-5, limits=frac(cfg.K1_QD2))], 
    targets=[
        xt.Target('qx', cfg.QX_TARGET, tol=2e-3),
        xt.Target('qy', cfg.QY_TARGET, tol=2e-3),
        xt.Target('dx', 0.0, at=xt.START, tol=1e-3),
        xt.Target(lambda tw: tw.betx.max(), cfg.BETA_X_MAX_M, tol=2e-3),])

opt.show()
opt.target_status()
print(env.vars.get_table())
tw = ring.twiss(method='4d', chrom=False)
print('bety_max', tw.bety.max(), '(11.839)')