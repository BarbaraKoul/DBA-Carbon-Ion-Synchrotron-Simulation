import xtrack as xt
import xobjects as xo
import numpy as np
import matplotlib.pyplot as plt
import config as cfg

ctx= xo.ContextCpu()

env = xt.Environment()

#setting the carbon-ion particle as reference particle
ion = xt.Particles(_context=ctx, mass0=cfg.ION_MASS_EV, q0=cfg.ION_Z, p0c=cfg.P0C)
env.set_particle_ref(ion)

#magnets
env.new('mb', xt.Bend, length=cfg.L_DIPOLE_M, angle=np.deg2rad(cfg.BEND_ANGLE_PER_DIPOLE_DEG))
env.new('mq', xt.Quadrupole, length=cfg.L_QUAD_M)

#quadruple strengths
env.new('qf1', 'mq', k1=cfg.K1_QF1)
env.new('qf2', 'mq', k1=cfg.K1_QF2)
env.new('qd1', 'mq', k1=cfg.K1_QD1)
env.new('qd2', 'mq', k1=cfg.K1_QD2)

tt_elem = env.elements.get_table()
print(tt_elem)