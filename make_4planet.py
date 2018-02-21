import numpy as np
import pandas as pd
import mr_forecast as mr
import matplotlib.pyplot as plt
import math
import sys
import os
import warnings

#This function draws a distribution from the values listed in the planets.csv file from NASA's exoplanet archive and also returns whether it found anything
def returnvalues(sample,parname,nsamples,ind=0):
	
	#get the mean
	mean=sample[parname][ind:ind+1].astype(np.float)
	
	#if a mean isn't found
	if math.isnan(mean):
		found=False

		#return a zero (used later)
		values=0
	
	#if a mean is found
	else:
		found=True
		
		#get the positive error
		err1=sample[parname+"err1"][ind:ind+1].astype(np.float)
		
		#if no error is found
		if math.isnan(err1):

			#just return a list of the mean repeated
			values=np.repeat(mean,nsamples)

			limitflag=sample[parname+"lim"][ind:ind+1].astype(np.int)
			#print a warning
			basewarning="No "+parname+" error information found for planet"+str(ind+1)+". The mean value is being repeated. "
			if limitflag==0:
				warnings.warn(basewarning+"The value does not appear to be a limit")
			elif limitflag==1:
				warnings.warn(basewarning+"The value appears to be a limit")
			else:
				warnings.warn(basewarning+"The value does not appear to be a limit, but an improperly accounted for special case inserted by the programmer modifying the orginal data file")
			warnings.resetwarnings()
			
		#if an error is found
		else:

			#draw values from an assymetric normal distribution
			err1=abs(err1)
			err2=abs(sample[parname+"err2"][ind:ind+1].astype(np.float))
			rand=np.random.normal(0,1,nsamples)
			error=np.zeros(nsamples)
			for x in range(0, nsamples):
				
				#if the drawn value will cause a negative final value (or a value greater than 1 for eccentricities), throw it out
				if parname=="pl_orbeccen":
					while float(mean+rand[x]*err2)<0 or float(mean+rand[x]*err1)>1:
						rand[x]=np.random.normal(0,1,1)
				elif parname=="pl_ecos_omegabar" or parname=="pl_esin_omegabar":
					while float(mean+rand[x]*err2)<-1 or float(mean+rand[x]*err1)>1:
						rand[x]=np.random.normal(0,1,1)
				else:
					while float(mean+rand[x]*err2)<0:
						rand[x]=np.random.normal(0,1,1)
				
				#multiply the error by its respective value
				if rand[x]>0:
					error[x]=rand[x]*err1
				else:
					error[x]=rand[x]*err2

			values=error+np.repeat(mean,nsamples)
			
			#this could be used if you didn't care about unequal errors
			#err=np.mean(np.array(abs(sample[parname+"err1"][ind:ind+1].astype(np.float)),abs(sample[parname+"err2"][ind:ind+1].astype(np.float))))	
			#values=np.random.normal(mean,err,nsamples)

	#values[values<0]=0 #should be accounted for by earlier while loop
	return np.array(values), found

#Assigns probabilistic mass and periods, and randomly draws arguments of periapsis, inclinations, and mean anomaly.
def planetparameters(sample,index, nsamples):

	#planet mass
	m, found=returnvalues(sample,"pl_bmassj",nsamples,ind=index)
	
	#if no mass is given, use the radius information to predict the mass distribution
	if not found:
		R, found=returnvalues(sample,"pl_radj",nsamples,ind=index)
		#use Rpost2M because it actually generates values instead of converting statistics
		m=mr.Rpost2M(R, unit='Jupiter')
	m=m*317.8 #Converting to Earth masses

	#period
	P, found=returnvalues(sample,"pl_orbper",nsamples,ind=index)

	#inclination
	inc=np.zeros(nsamples)

	#longitude of ascending node not necessary because inclination is 0

	#argument of periapsis
	w=2*np.pi*np.random.sample(nsamples)

	#mean anomaly
	MA=2*np.pi*np.random.sample(nsamples)

	return m, P, inc, w, MA

#Draws eccentricities for those planets that dont have any
def eccentricities(sample,index,nsamples,ms,p1,w,p0=np.zeros(2),p2=np.zeros(2),m0=np.zeros(2),m2=np.zeros(2)):

	limitflag=sample["pl_orbeccenlim"][index:index+1].astype(np.int)
	limitflagint=int(limitflag)

	#if a no limit flag or a standard one is found, procede normally
	if math.isnan(limitflag) or limitflagint<2:
		e, found=returnvalues(sample,"pl_orbeccen",nsamples,ind=index)

	#this should only happen for KOI-94_mas right now
	elif limitflagint==2:

		ecos, found=returnvalues(sample,"pl_ecos_omegabar",nsamples,ind=index)
		esin, found=returnvalues(sample,"pl_esin_omegabar",nsamples,ind=index)

		#use the values from the Masuda et al. paper to contruct new e and w distributions
		e=np.sqrt(np.multiply(ecos,ecos)+np.multiply(esin,esin))
		w=np.add(np.arctan2(esin,ecos),2*np.pi)%(2*np.pi)
		found=True

	#if a totally weird limit flag is found
	else:
		found=False
	
	if not found:

		#cannot reference nsamples in function definition
		if not p0[0]:
			p0=np.zeros(nsamples)
		if not p2[0]:
			p2=np.zeros(nsamples)
		if not m0[0]:
			m0=np.zeros(nsamples)
		if not m2[0]:
			m2=np.zeros(nsamples)

		e=np.zeros(nsamples)
		earth = 0.000003003
		
		#get semi-major axes
		for x in range(0, nsamples):
			a0 = ((p0[x]/365)**2 * ms[x])**(1./3.)
			a1 = ((p1[x]/365)**2 * ms[x])**(1./3.)
			a2 = ((p2[x]/365)**2 * ms[x])**(1./3.)
			
			#place where orbits would cross if the inner orbit (a0) had 0 eccentricity (if no p0 is passed, this equals 1)
			ecrit10 = (a1-a0)/a1

			#place where orbits would cross if the outer orbit (a2) had 0 eccentricity
			if not p2[0]:
				ecrit12 = 1
			else:
				ecrit12 = (a2-a1)/a1

			#maximum eccentricty bound by inner and outer orbits (maybe include hill radii in this calculation later?)
			logemax = np.log10(min(ecrit10, ecrit12))

			#minimum eccentricity based on perturbations by the inner or outer planet (derived by Dan Tamayo)
			logemin = np.log10(max(m0[x]*earth/ecrit10**2, m2[x]*earth/ecrit12**2))

			#seeding evenly (in log space) eccentricities between the upper and lower bounds
			e[x]=10**np.random.uniform(logemin, max(logemax,0))

	#else:
		#e[e<0]=0 #should be taken care of by while loop in return values
		#e[e>1]=1 # keeps eccentricities less than 1 (should be properly corrected in returnvalues)

	return e, w

#function that writes the job files
def write_jobs(indices, system, jobs_dir, norbits, Np):

	#generate jobs
	#currentplace=list(data.index)[0]
	os.system('mkdir %s'%jobs_dir)
	for i in range(len(indices)):
		id_ = indices[i]             #id number of sample
		job_name = "%s_%d_%d"%(system,int(np.log10(norbits)),id_)
		sh_script_name = "%s%s"%(jobs_dir,job_name)
		with open(sh_script_name, 'w') as f:
			f_head = open('job_header_icsaci','r')
			f.write(f_head.read())
			f_head.close()
			f.write('#PBS -N %s \n'%job_name) #Names the job
			f.write('cd $PBS_O_WORKDIR\n')      #This will be the directory the jobs are submitted from
			f.write('source activate stability \n')
			#f.write('python run_4planet.py %s %d %f %d %d %s >& batch.output\n'%(system,id_,Ms[id_-currentplace],norbits,Np,job_name))
			f.write('python run_Nplanet.py %s %d %d %d %s >& batch.output\n'%(system,id_,norbits,Np,job_name))
		f.close()

	return 1

#collects generated parameters
def collect_parameters(sample, n_sims):

	#number of planets based on passed sample
	Np = sample.shape[0]

	#stellar mass
	Ms, found=returnvalues(sample,"st_mass",n_sims)

	#arrays to hold all planet parameters
	m = np.zeros((Np,n_sims))
	P = np.zeros((Np,n_sims))
	inc = np.zeros((Np,n_sims))
	w = np.zeros((Np,n_sims))
	MA = np.zeros((Np,n_sims))
	e = np.zeros((Np,n_sims))

	#get all planet parameters
	for i in range(Np):
		m[i,:], P[i,:], inc[i,:], w[i,:], MA[i,:] = planetparameters(sample, i, n_sims)

	#get all eccentricities (first and last have to be done separately as they dont fit the syntax pattern of planets that are surrounded by other planets)
	e[0,:], w[0,:] = eccentricities(sample,0,n_sims,Ms,P[0,:],w[0,:],p2=P[1,:],m2=m[1,:])
	e[Np-1,:], w[Np-1,:] = eccentricities(sample,Np-1,n_sims,Ms,P[Np-1,:],w[Np-1,:],p0=P[Np-2,:],m0=m[Np-2,:])
	for i in range(1,Np-1):
		e[i,:], w[i,:] = eccentricities(sample,i,n_sims,Ms,P[i,:],w[i,:],p0=P[i-1,:],p2=P[i+1,:],m0=m[i-1,:],m2=m[i+1,:])
	return Ms, m, P, inc, w, MA, e, Np

#saves data to csv
def save_data(Ms, m, P, w, MA, e, Np, n_sims, dat_dir, system):

	#generate file header based on amount of planets
	orb_elements = ["Ms"]
	for i in range(1,Np+1):
		orb_elements.append("m%d"%i)
		orb_elements.append("MA%d"%i)
		orb_elements.append("P%d"%i)
		orb_elements.append("e%d"%i)
		orb_elements.append("w%d"%i)

	#create properly formatted data frame (this could probably be done more efficiently)
	data = []
	for i in range(n_sims):
		intermediate=np.zeros(5*Np+1)
		intermediate[0]=Ms[i]
		for x in range(Np):
			intermediate[5*x+1]=m[x,i]
			intermediate[5*x+2]=MA[x,i]
			intermediate[5*x+3]=P[x,i]
			intermediate[5*x+4]=e[x,i]
			intermediate[5*x+5]=w[x,i]
		data.append(intermediate)
		#data.append([Ms[i],m[0,i],MA[0,i],P[0,i],e[0,i],w[0,i],m2[i],MA2[i],P2[i],e2[i],w2[i],m3[i],MA3[i],P3[i],e3[i],w3[i],m4[i],MA4[i],P4[i],e4[i],w4[i]])
	data = pd.DataFrame(np.asarray(data),columns=orb_elements)

	#save data to csv
	#os.system('mkdir %s'%dat_dir)
	incl_header=True
	data_file = "%s%s_data.csv"%(dat_dir,system)
	if os.path.isfile(data_file) == True:
		data.index += pd.read_csv(data_file).index[-1] + 1   #start current index number at previous entry number
		incl_header=False
	data.to_csv(data_file, mode="a", header=incl_header)
	
	indices = list(data.index)

	return indices

#wrapper function to generate the job files
def generate_jobs(sample,system,dat_dir,jobs_dir,dir_path,n_sims,norbits,permute=0,plotstuff=0):

	#create and save parameters and create job files for original system
	Ms, m, P, inc, w, MA, e, Np = collect_parameters(sample, n_sims)
	indices = save_data(Ms, m, P, w, MA, e, Np, n_sims, dat_dir, system)
	out = write_jobs(indices, system, jobs_dir, norbits, Np)

	if plotstuff:
		#plotting
		for i in range(sample.shape[0]):
			plt.hist(e[i,:], bins='auto', facecolor='red', alpha=0.75)
			fig = plt.gcf()
			fig.savefig("%s planet %d eccentricities"%(samplename,i)+".png")   # save the figure to file
			plt.close(fig)

	#if you want to create jobs of all equivalent 3-planet systems
	if permute:

		#for each planet
		for i in range(Np):
			
			#name the systems based on which planet is removed
			new_name=system+"_r%d"%(i+1)
			new_jobs_dir = dir_path+"/jobs/"+new_name+"/"     #output directory for jobs
			
			#save the data and write the jobs for the new systems
			new_indices = save_data(Ms, np.delete(m,i,0), np.delete(P,i,0), np.delete(w,i,0), np.delete(MA,i,0), np.delete(e,i,0), Np-1, n_sims, dat_dir, new_name)
			out = write_jobs(new_indices, new_name, new_jobs_dir, norbits, Np-1)
		
		print("Generated %d (%d permuted) simulations for %s"%(samps, samps*Np, samplename))
	
	else:
		
		print("Generated %d simulations for %s"%(samps,samplename))

	return 1

#reading in NASA exoplanets archive information for 4 planets as a panda array
data = pd.read_csv('planets_mod.csv', header=40)

if __name__ == '__main__': #do this if the file is called directly from the console (not by another function)

	names=list(set(data["pl_hostname"])) #list of unique names
	names.sort() #sort the list
	#print(names)
	#['GJ 3293','GJ 676 A','GJ 876','HD 141399','HD 160691','HD 20794','HR 8799','K2-72','KOI-94','Kepler-106','Kepler-107','Kepler-132',
	#'Kepler-1388','Kepler-1542','Kepler-167','Kepler-172','Kepler-176','Kepler-197','Kepler-208','Kepler-215','Kepler-220','Kepler-221',
	#'Kepler-223','Kepler-224','Kepler-235','Kepler-24','Kepler-245','Kepler-251','Kepler-256','Kepler-26','Kepler-265','Kepler-282',
	#'Kepler-286','Kepler-299','Kepler-304','Kepler-306','Kepler-338','Kepler-341','Kepler-342','Kepler-37','Kepler-402','Kepler-48',
	#'Kepler-49','Kepler-758','Kepler-79','Kepler-82','Kepler-85','WASP-47','tau Cet']
	
	#system you want to generate jobs for. NEEDS TO BE CHANGED EVERYTIME!
	samplename= "KOI-94_mas"

	#finding and isolating the system data
	sampleindex=np.where(data["pl_hostname"]==samplename)[0]
	sample=data[sampleindex[0]:sampleindex[-1]+1]

	sample=sample.sort_values(by="pl_orbper") #sort the samples in order of their periods
	samplename=samplename.replace(" ","_") #take out spaces in the sample name (for creating directories)

	dir_path = os.path.dirname(os.path.realpath(__file__)) #directory of this program
	jobs_dir = dir_path+"/jobs/"+samplename+"/"     #output directory for jobs
	dat_dir = dir_path+"/systems/"    #output directory for storing _data.csv files
	
	#number of sims created. NEEDS TO BE CHANGED EVERYTIME!
	samps = 200
	
	#number of orbits of innermost planet. NEEDS TO BE CHANGED EVERYTIME!
	norbits = 1e9
	
	#do you want to create systems where individual planets are removed? NEEDS TO BE CHANGED EVERYTIME!
	perm=1

	#do you want to plot the distribution of eccentricities
	plots=1

	#generalize to multiple systems later?
	# systems = [samplename]
	# for system in systems:

	#generate jobs for original system
	out = generate_jobs(sample,samplename,dat_dir,jobs_dir,dir_path,samps,norbits,permute=perm,plotstuff=plots)
