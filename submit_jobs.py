#This script submits jobs to the ICS-ACI cluster
import os
import os.path
import glob
import numpy as np

#takes all jobs in folder and submits them from the directory that this program is in
def submit_job(f, job_name):
    os.system('mv %s %s'%(f, job_name))
    os.system('qsub %s'%job_name)
    os.system('mv %s %s'%(job_name,f))

#retrieves list of unsubmitted or incomplete job files in the specified jobs directory
def find_unfinished_jobs(jobs_dir):
    unsub_jobs = []
    N_jobs_submitted = 0
    jobs = glob.glob('%s/*'%jobs_dir)
    for j in jobs:
        basename = os.path.basename(j)
        #if a final simulation archive exits, the job with that name has already been submitted
        if os.path.isfile('output/'+samplename+'/%s_final.bin'%basename) == False:
            unsub_jobs.append(j)
        else:
            N_jobs_submitted += 1
    print('N_jobs_submitted=%d'%N_jobs_submitted)
    print('found %d jobs'%len(unsub_jobs))
    return unsub_jobs

###############################
#['GJ 3293','GJ 676 A','GJ 876','HD 141399','HD 160691','HD 20794','HR 8799','K2-72','KOI-94','Kepler-106','Kepler-107','Kepler-132',
#'Kepler-1388','Kepler-1542','Kepler-167','Kepler-172','Kepler-176','Kepler-197','Kepler-208','Kepler-215','Kepler-220','Kepler-221',
#'Kepler-223','Kepler-224','Kepler-235','Kepler-24','Kepler-245','Kepler-251','Kepler-256','Kepler-26','Kepler-265','Kepler-282',
#'Kepler-286','Kepler-299','Kepler-304','Kepler-306','Kepler-338','Kepler-341','Kepler-342','Kepler-37','Kepler-402','Kepler-48',
#'Kepler-49','Kepler-758','Kepler-79','Kepler-82','Kepler-85','WASP-47','tau Cet']

permute=0

#what systems to use. NEEDS TO BE CHANGED EVERYTIME!
#systems = ["Ari Fake 10 0.05","Ari Fake 20 0.05","Ari Fake 30 0.05","Ari Fake 40 0.05","Ari Fake 50 0.05", "Ari Fake 10 0.1","Ari Fake 20 0.1","Ari Fake 30 0.1","Ari Fake 40 0.1","Ari Fake 50 0.1"]
#systems = ["Rasio Fake 2 Jup 1"]
systems=["bully 0.045", "super earth 0.06", "neptune 0.1", "Ari Fake 15 0.09", "analog 0.04", "3.5 planet 0.02"]
# for i in ["bully ", "super earth ", "neptune ", "Ari Fake 15 ", "analog ", "3.5 planet "]:
# # for i in ["bully "]:
#     for j in ["0.01", "0.02", "0.03", "0.04", "0.05", "0.075", "0.1"]:
#     # for j in ["0.1"]:
#         systems+=[i+j]

duplicate=list(systems)
if permute:
    for samplename in duplicate:
        for i in range(4):
            systems += [samplename+" r%d"%(i+1)]

#print(systems)
for samplename in systems:

    samplename=samplename.replace(" ","_") #take out spaces in the sample name (for creating directories)
    jobs_dir="jobs/"+samplename+"/"

    files = find_unfinished_jobs(jobs_dir)
    Njobs_counter = 0

    # submit every found job
    for f in files:
        job_name = f.split(jobs_dir)[1]
        submit_job(f, job_name)
        Njobs_counter += 1 #seems unnecessary to keep this

    print('submitted %d jobs'%Njobs_counter)