import xtrack as xt
import xpart as xp
import xobjects as xo
import numpy as np
import matplotlib.pyplot as plt
import configRF as cfg

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

#rf cavity
env.new('rf1', xt.Cavity, length=cfg.L_CAVITY_RF, voltage=cfg.RF_VOLT, harmonic=cfg.RF_HARM)


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
    env.place('sxr1', anchor='start', at=-13.129, from_='qf2@end'),
    env.place('rf1', anchor='start', at=-15.003, from_='qf2@end'),
])

plt.close('all')
sv = ring.survey()
sv.plot()

tw = ring.twiss6d()
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
print("tune: ",tw.qs)

#tt = ring.get_table()
#tt.show(cols=['s_start', 's_center', 's_end'])

plt.show()

monitor = xt.BeamStatsMonitor(
    start_at_turn=0,
    stop_at_turn=300,
    stats=["num_particles", "sigma_zeta", "mean_zeta", "betx", "bety", "dx", "dy"],
)

ring.insert("beam_monitor", monitor, at=0)

np.random.seed(42)
ion = ring.xpart.generate_matched_gaussian_bunch(
    num_particles=100,
    total_intensity_particles=3,
    nemitt_x=11e-5,
    nemitt_y=29e-5,
    sigma_z=1e-3,
)

ring.track(ion, num_turns=300,
           with_progress=1
)

#monitor.stats  
#monitor.available_levels 
#monitor.default_level  
#monitor.turns  
#monitor.mean_x.shape
#monitor.num_particles.shape

plt.close('all')
fig, axes = plt.subplots(2, 1, sharex=True)
axes[0].plot(monitor.turns, monitor.sigma_zeta * 1000)
axes[0].set_ylabel("beam size zeta [mm]")
axes[1].plot(monitor.turns, monitor.mean_zeta * 1000)
axes[1].set_ylabel("mean_zeta [mm]")
axes[1].set_xlabel("turn")
plt.show()

plt.close('all')
fig, axes = plt.subplots(2, 1, sharex=True)
axes[0].plot(monitor.turns, monitor.betx)
axes[0].set_ylabel("betx [m]")
axes[1].plot(monitor.turns, monitor.bety * 10)
axes[1].set_ylabel("bety [m]")
axes[1].set_xlabel("turn")
plt.show()

plt.close('all')
fig, axes = plt.subplots(2, 1, sharex=True)
axes[0].plot(monitor.turns, monitor.dx)
axes[0].set_ylabel("dispersionX")
axes[1].plot(monitor.turns, monitor.dy)
axes[1].set_ylabel("dispersionY")
axes[1].set_xlabel("turn")
plt.show()
