# Script Description
# generate_compact_system - ---
# 	---
# 
# Inputs:
# 
# 	---
# 
# Outputs:
# 
# 	--- 
# 
# Other files required: ---
# 
# See also: ---
# 
# Author: Christian Gilbertson, BS., Physics and Engineering Science and Mechanics
#     Graduate student at The Pennsylvania State University
# Email address: chrisgil@psu.edu  
# Website: https://sites.psu.edu/chrisgil/
# February 2018; Last revision: 13-February-2018

# ------------- BEGIN CODE --------------
import numpy as np
from scipy import optimize

#This function returns the combined mass for two planets needed to be at a certain hill radius limit. 
#r is the period ratio between the two planets
#mass_limit scales the result. Effectively used to subtract the value of one planet to find a limit on the other
def f(r,hill_limit,mass_limit):
	Mstar=Ms/earth
	return ((24*Mstar*(r**(2/3)-1)**3)/(hill_limit**3*(r**(2/3)+1)**3)-mass_limit)

#Draws a valid period that is not close to a resonance and within the passed limits
def next_period(p,index,lower_ratio_limit,upper_ratio_limit):
	closest2resonance=0
	j=index+1

	#if the current period is too close to a resonance, draw another
	while closest2resonance<resonance_tolerance:

		#draw the period within the given limits
		p[j-1]=p[j-2]*10**(np.random.uniform(low=np.log10(lower_ratio_limit), high=np.log10(upper_ratio_limit)))
		
		#period ratios of planets
		ratios=np.zeros((j,j))
		for x in range(j):
			for y in range(j):
				ratios[x,y]=p[x]/p[y]		

		#check for closest period ratio to resonance
		close2resonance=np.zeros((j,j))
		for x in range(j):
			for y in range(j):
				if y<x:
					close2resonance[x,y]=min(abs(ratios[x,y]-resonances))
		
		#find the smallest nonzero value in the close to resonance matrix
		closest2resonance=min(i for i in close2resonance.flatten() if i > 0)

	#return the period ratio (for convenient printing later)
	return ratios

#Setting some constants
earth = 0.000003003 #Earth mass/Sun mass
Ms=1 #Mass of star for the test system in solar masses

#Number of planets we want in the system
Np=4

#Initializing period and mass lists
P=np.zeros(Np)
P[0]=1
m=np.zeros(Np)

#Setting hard limits
lowerplanetmass=1/10 #Earth masses
upperplanetmass=100 #Earth masses
lowerhill=8 #mutual hill radii
upperhill=12 #mutual hill radii
secondaryhill=30 #mutual hill radii
resonance_tolerance=0.05 #how close to a resonance a period ratio can be

#resonances to check for
orders=4
far=10 #ends at (far+order)/far
resonances=np.zeros(far*orders+1)
resonances[-1]=1
for x in range(orders):
	for y in range(far):
		resonances[far*x+y]=(y+2+x)/(y+1)
resonances=np.unique(resonances)
#print(resonances)

#find the smallest period ratio allowed by the closeness to resonance criterion
i=0
min_resonance_period_ratio=1
while i<len(resonances):

	#If the current minimum period ratio is too close to the next highest resonance, increase it
	if min_resonance_period_ratio>=resonances[i]-resonance_tolerance:
		min_resonance_period_ratio=resonances[i]+resonance_tolerance
		i+=1

	#if it is far enough from the next highest resonance, leave it where it is and end the loop
	else:
		endloop=len(resonances)

		#take out the values from the resonance list that won't be used
		resonances=resonances[i:]
		
		#break the loop
		i=endloop
#print(min_resonance_period_ratio)
#print(resonances)

#draw periods and masses interatively until you get a working result
i=0
while i<Np:

	#for the first draw, you just need to get a mass within the hard limits
	if i==0:

		m[i]=10**np.random.uniform(low=np.log10(lowerplanetmass), high=np.log10(upperplanetmass))
		i+=1

	#for other draws, get period and mass
	elif i>0:

		#get limits on period ratios from constraints on total mass and period ratios
		lower_ratio_limit = max(min_resonance_period_ratio,optimize.brenth(f, 1, 5, args=(lowerhill,upperplanetmass+m[i-1])))
		upper_ratio_limit = optimize.brenth(f, 1, 5, args=(upperhill,upperplanetmass+m[i-1]))

		#if the constraints are impossible, start fresh
		if lower_ratio_limit>upper_ratio_limit:
			i=0
			print("reset b/c of period constraints")

		#if a period can be drawn, draw it and get a mass
		else:
			ratios=next_period(P,i,lower_ratio_limit,upper_ratio_limit)

			#after the second planet, look at secondary hill radii separations as well
			if i>1:
				lower_mass_limit=max(lowerplanetmass, f(P[i]/P[i-1],upperhill,m[i-1]), f(P[i]/P[i-2],secondaryhill,m[i-2]))
			else:
				lower_mass_limit=max(lowerplanetmass, f(P[i]/P[i-1],upperhill,m[i-1]))
			upper_mass_limit=min(upperplanetmass, f(P[i]/P[i-1],lowerhill,m[i-1]))

			#if the constraints are impossible, start fresh
			if lower_mass_limit>upper_mass_limit:
				i=0
				print("reset b/c of mass constraints")

			#if a mass can be drawn, draw it
			else:
				m[i]=10**np.random.uniform(low=np.log10(lower_mass_limit), high=np.log10(upper_mass_limit))
				
				#move on to the next planet
				i+=1

scalefactor=365.25/10 #the smallest orbital period is 1/10 of a year (in days)
Preal=P*scalefactor #I'm pretty sure this doesn't matter

#get semi-major axes
a=np.zeros(Np)
for i in range(Np):
    a[i]=((Preal[i]/365*(2*np.pi))**2 * Ms)**(1./3.)

#find mutual hill radii separations for printing
mut_hill=np.zeros((Np,Np))
separation=np.zeros((Np,Np))
for x in range(Np):
	for y in range(Np):
		if y<x:
			mut_hill[x,y]=np.mean([a[x],a[y]])*((m[x]+m[y])*earth/Ms/3.)**(1./3.)
			separation[x,y]=abs(a[x]-a[y])/mut_hill[x,y]

print("\nPeriods (scaled to first period)\n")
print(P)
print("\nPeriod Ratios\n")
print(ratios)
print("\nPlanet Masses (Earth masses)\n")
print(m)
print("\nHill Separations\n")
print(separation)

# # ------------- END OF CODE --------------