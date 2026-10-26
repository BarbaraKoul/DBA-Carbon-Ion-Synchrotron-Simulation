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
env.new('ms', xt.Sextupole, length=cfg.L_SEXTUPOLE_M)

#quadruple strengths
env.new('qf1', 'mq', k1=cfg.K1_QF1)
env.new('qf2', 'mq', k1=cfg.K1_QF2)
env.new('qd1', 'mq', k1=cfg.K1_QD1)
env.new('qd2', 'mq', k1=cfg.K1_QD2)
env.new('sxr1', 'ms', k2=cfg.K2_SXR1)


tt_elem = env.elements.get_table()
print(tt_elem)

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

ring.insert([
    env.place('sxr1', anchor='start', at=-14.729, from_='qf2@end')
])
plt.close('all')
sv=ring.survey()
sv.plot()

tw=ring.twiss4d()
tw.plot()

chrom_x=tw.dqx
chrom_y=tw.dqy
print("Horizontal chromaticity:", chrom_x)
print("Vertical chromaticity:", chrom_y)
print("betx max: ", tw.betx.max())

print("bety max: ", tw.bety.max())
print("dx max: ", tw.dx.max())
print("dpx max: ", tw.dpx.max())
print("qx: ", tw.qx)
print("qy: ", tw.qy)

#tt = ring.get_table()
#tt.show(cols=['s_start', 's_center', 's_end'])

plt.show()