#The purpose of this script is to run the jobs generated by "make_4planet.py". Actual 5-body simulations
import pandas as pd
import numpy as np
import os
import sys
import rebound
import time
import math

#if a collision occurs, end the simulation
def collision(reb_sim, col):
    reb_sim.contents._status = 5 # causes simulation to stop running and have flag for whether sim stopped due to collision
    return 0

#arguments
system = sys.argv[1]
id_ = int(sys.argv[2])
maxorbs = float(sys.argv[3])
Nplanets = int(sys.argv[4])
#name = sys.argv[5] #replaces because of name restrictions imposed by ics-aci
#name = "%s_%d_%d"%(system,int(np.log10(maxorbs)),id_)
name = "%s_%d"%(system,id_)

#load data
data = pd.read_csv('systems/%s_data.csv'%system)
d = data.iloc[id_]

Ms=d["Ms"]

#find minimum mutual hill radius (where collisions will be detected)
a=np.zeros(Nplanets)
for i in range(Nplanets):
    a[i]=((d["P%d"%(i+1)]/365*(2*np.pi))**2 * Ms)**(1./3.)
solar_mass_2_earth_mass = 0.000003003

#mut_hill=np.zeros(Nplanets-1)
#for i in range(Nplanets-1):
#    mut_hill[i]=np.mean([a[i],a[i+1]])*((d["m%d"%(i+1)]+d["m%d"%(i+2)])*solar_mass_2_earth_mass/Ms/3.)**(1./3.)
#minradius=min(mut_hill)

radii=np.zeros(Nplanets)
for i in range(Nplanets):
    radii[i]=a[i]*(d["m%d"%(i+1)]*solar_mass_2_earth_mass/Ms/3.)**(1./3.)/2 #half hill radius

#set up simulation
sim = rebound.Simulation()
sim.integrator = 'whfast'
sim.G = 1
sim.ri_whfast.safe_mode = 0
sim.collision = 'direct'
sim.collision_resolve = collision

#add star
sim.add(m=Ms)

#add planets
for i in range(1,Nplanets+1):
    m = d["m%d"%i]*solar_mass_2_earth_mass
    P = d["P%d"%i]
    e = d["e%d"%i]
    w = d["w%d"%i]
    M = d["MA%d"%i]
    sim.add(m=m, P=P*2*np.pi/365., e=e, omega=w, M=M, r=radii[i-1]) #G=1 units!
sim.move_to_com()

#timestep
#To ensure a low energy error (dE/E(0) < 1e-9) the timestep needs to be a small fraction of the innermost orbital period (in this case roughly 3%).
#Studies have shown that if the timestep is an exact fraction of the innermost orbital period then numerical artifacts can be introduced 
dt = 2.*math.sqrt(3)/100. #0.0346410162
P1 = sim.particles[1].P
sim.dt = dt*P1 # ~3% orbital period
tmax = maxorbs*P1

#save simulation archive
#dir_path = os.path.dirname(os.path.realpath(__file__)) #directory of this program
#out_dir = dir_path+"/output/"
#os.system('mkdir %s'%out_dir)
sim.initSimulationArchive('output/%s_SA.bin'%name, interval=tmax/1000.)     #save checkpoints.

#simulate
E0 = sim.calculate_energy()
t0 = time.time()
print("starting simulation")

try:
    sim.integrate(tmax) #will stop if collision occurs
except:
    pass
print("finished simulation")

Ef = sim.calculate_energy()
Eerr = abs((Ef-E0)/E0)

#need to store the result somewhere
f = open('systems/%s_Nbodyresults.csv'%system, "a")
f.write('%s, %d, %e, %e, %e, %e, %e \n'%(name,id_,maxorbs,P1,sim.t,Eerr,time.time()-t0))

