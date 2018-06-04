# Import the rebound module
import rebound
import os


def replace_snapshot(sim, filename):
    if os.path.isfile(filename):
        os.remove(filename)
    sim.simulationarchive_snapshot(filename)


div = 100

filename = "simulationarchive.bin"
tmax = 1 * 3.1415 * 1e5

try:
    sim = rebound.Simulation.from_archive(filename)
    sim.automateSimulationArchive(filename, interval=tmax / div)
    print("Restarting from simulation archive. Last snapshot found at t=%.1f" % sim.t)
except:
    print("Cannot load SimulationArchive. Creating new simulation.")

    sim = rebound.Simulation()
    sim.integrator = 'whfast'
    sim.G = 1
    sim.ri_whfast.safe_mode = 0

    sim.collision = 'line'
    # sim.collision_resolve = collision

    sim.add(m=1)  # star
    r = 0.0001
    sim.add(m=1e-3, P=1, e=0.1, r=r)  # planet 1
    sim.add(m=1e-3, P=1.3, e=0.1, r=r)  # planet 2
    sim.dt = 3.1415 * 2. * 6. / 365.25  # 6 days in units where G=1
    sim.move_to_com()
    # sim.automateSimulationArchive(filename, interval=tmax/1000, deletefile=True)
    sim.automateSimulationArchive(filename, interval=tmax / div)


try:
    sim.integrate(tmax)
except rebound.Collision:
    replace_snapshot(sim, "final_collide.bin")

sa = rebound.SimulationArchive(filename)
print(sa[-2].t)
print(sa[-1].t)
print(sa[-1].t - sa[-2].t)
print(tmax / div)
print(sim.dt)
print(len(sa))

replace_snapshot(sim, "final.bin")
