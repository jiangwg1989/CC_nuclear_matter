####################################################
#initialize the jobs for nuclear matter
####################################################
import numpy as np
import os



def output_ccm_in_file(file_path,cD,cE,LEC_ci,c1s0,c3s1,cnlo,particle_num,matter_type,density,nmax):
    with open(file_path,'w') as f_1:
        f_1.write('!Chiral order for Deltas(LO = 0,NLO=2,NNLO=3,N3LO=4) and cutoff'+'\n')
        f_1.write('3, 450\n')
        f_1.write('! cE and cD 3nf parameters:'+'\n' )
        f_1.write('%.2f, %.2f\n' % (cE,cD))
        f_1.write('! LEC ci \n')
        f_1.write('%.12f, %.12f, %.12f, %.12f \n' % (LEC_ci[0],LEC_ci[1],LEC_ci[2],LEC_ci[3]))
        f_1.write('!c1s0 & c3s1 \n')
        f_1.write('%.12f, %.12f, %.12f, %.12f, %.12f, %.12f \n' % (c1s0[0],c1s0[1],c1s0[2],c3s1[0],c3s1[1],c3s1[2]))
        f_1.write('! cnlo(7) \n')
        f_1.write('%.12f, %.12f, %.12f, %.12f, %.12f, %.12f, %.12f \n' % (cnlo[0],cnlo[1],cnlo[2],cnlo[3],cnlo[4],cnlo[5],cnlo[6]))
        f_1.write('! number of particles'+'\n')
        f_1.write('%d\n' % (particle_num) )
        f_1.write('! specify: pnm/snm, input type: density/kfermi'+'\n')
        f_1.write(matter_type+', density'+'\n')
        f_1.write('! specify boundary conditions (PBC/TABC/TABCsp)'+'\n')
        f_1.write('PBC'+'\n')
        f_1.write('! dens/kf, ntwist,  nmax'+'\n')
        f_1.write('%.2f, 1, %d\n' % (density, nmax))
        f_1.write('! specify cluster approximation: CCD, CCDT'+'\n')
        f_1.write('CCD(T)'+'\n')
        f_1.write('! tnf switch (T/F) and specify 3nf approximation: 0=tnf0b, 1=tnf1b, 2=tnf2b'+'\n')
        f_1.write('T, 3'+'\n')
        f_1.write('! 3nf cutoff(MeV),non-local reg. exp'+'\n')
        f_1.write('450, 3'+'\n')


#def output_job_file(file_path,cD,cE,particle_num,matter_type,density,nmax):
#    ccm_in_file_path = './output/ccm_in_'+matter_type+'_%d_cD_%.2f_cE_%.2f_rho_%.2f' % (particle_num,cD,cE,density)
#    ccm_out_file_path = './'+matter_type+'_%d_cD_%.2f_cE_%.2f_rho_%.2f.out &' % (particle_num,cD,cE,density)
#    with open(file_path,'w') as f_1:
#        f_1.write('export OMP_NUM_THREADS=4'+'\n\n')
#        f_1.write('/home/g1u/sw/intel_openmpi/bin/mpiexec -np 4 ./prog_ccm.exe '+ccm_in_file_path+' > '+ccm_out_file_path)
#        f_1.write('wait'+'\n')



####################################################
#  set up all the parameters
####################################################
nmax = 2
particle_num = 28
neutron_num = 14
cE_min = -1
cE_max = -1
cE_gap = 1 
cE_count = int( (cE_max - cE_min) / cE_gap + 1 )
cD_min = -3
cD_max = -3
cD_gap = 1
cD_count = int( (cD_max - cD_min) / cD_gap + 1 )
density_min = 0.12
density_max = 0.20
density_gap = 0.02
density_count = int( (density_max - density_min) / density_gap +1 )
CADES_run_path = '/lustre/or-hydra/cades-nucphys/w01'
nodes_num = 16
threads_num = 16
walltime = '00:05:00'
LEC_ci = np.zeros(4)
c1s0   = np.zeros(3)
c3s1   = np.zeros(3)
cnlo   = np.zeros(7)
LEC_ci[0] = -0.74 
LEC_ci[1] = -0.49
LEC_ci[2] = -0.65
LEC_ci[3] =  0.96
c1s0[0]   =  -0.337137 #c1s0pp 
c1s0[1]   =  -0.338139 #c1s0np
c1s0[2]   =  -0.338023 #c1s0nn
c3s1[0]   =  -0.229310 #c3s1
c3s1[1]   =  -0.229310
c3s1[2]   =  -0.229310
cnlo[0]   =   2.476589  #1s1
cnlo[1]   =   0.645550  #3p0
cnlo[2]   =  -0.028541  #1p1
cnlo[3]   =  -1.022359  #3p1
cnlo[4]   =   0.695953  #3s1
cnlo[5]   =   0.358330  #3s13d1
cnlo[6]   =  -0.870203  #3p2  


os.system('mkdir /lustre/or-hydra/cades-nucphys/w01/output')
os.system('mkdir /lustre/or-hydra/cades-nucphys/w01/output/output')
os.system('cp /home/w01/CC_nuclear_matter/CCM_kspace_deltafull/prog_ccm.exe /lustre/or-hydra/cades-nucphys/w01/prog_ccm.exe')
####################################################
#  set up all the ccm_in files for snm 
####################################################
for loop1 in range(cE_count):
    for loop2 in range(cD_count):
        cE = cE_min + cE_gap * loop1
        cD = cD_min + cD_gap * loop2
        for loop3 in range(density_count): 
            density = density_min + density_gap * loop3
            file_path = CADES_run_path+'/output/ccm_in_snm_%d_cD_%.2f_cE_%.2f_rho_%.2f' % (particle_num,cD,cE,density)
            output_ccm_in_file(file_path,cD,cE,LEC_ci,c1s0,c3s1,cnlo,particle_num,'snm',density,nmax)


####################################################
#  set up all the ccm_in files for pnm /home/w01/CC_nuclear_matter/output
####################################################
for loop1 in range(cE_count):
    for loop2 in range(cD_count):
        cE = cE_min + cE_gap * loop1
        cD = cD_min + cD_gap * loop2
        for loop3 in range(density_count):
            density = density_min + density_gap * loop3
            file_path = CADES_run_path+'/output/ccm_in_pnm_%d_cD_%.2f_cE_%.2f_rho_%.2f' % (neutron_num,cD,cE,density)
            output_ccm_in_file(file_path,cD,cE,LEC_ci,c1s0,c3s1,cnlo,neutron_num,'pnm',density,nmax)


####################################################
#  set up job script for snm 
####################################################
matter_type = 'snm'
file_path = './output/snm_job.script'
with open(file_path,'w') as f_1: 
    f_1.write('#PBS -N nuclearmatter'+'\n')
    f_1.write('#PBS -j eo'+'\n')
    f_1.write('#PBS -l walltime='+walltime+',nodes='+str(int(nodes_num))+'\n')
    f_1.write('export OMP_NUM_THREADS='+str(int(threads_num))+'\n\n')
    f_1.write('cd '+CADES_run_path+'\n')
    for loop1 in range(cE_count):
        for loop2 in range(cD_count):
            cE = cE_min + cE_gap * loop1
            cD = cD_min + cD_gap * loop2
            for loop3 in range(density_count):
                density = density_min + density_gap * loop3
                ccm_in_file_path = './ccm_in_'+matter_type+'_%d_cD_%.2f_cE_%.2f_rho_%.2f' % (particle_num,cD,cE,density)
                ccm_out_file_path = './output/'+matter_type+'_%d_cD_%.2f_cE_%.2f_rho_%.2f.out &' % (particle_num,cD,cE,density)
                f_1.write('aprun  -n'+str(int(nodes_num))+' -d'+str(int(threads_num))+' ./prog_ccm.exe '+ccm_in_file_path+' > '+ccm_out_file_path+'\n')
                f_1.write('wait'+'\n')





####################################################
#  set up all the ccm_in files for pnm 
####################################################
matter_type = 'pnm'
file_path = './output/pnm_job.script'
with open(file_path,'w') as f_1:
    f_1.write('#!/bin/bash\n')
    f_1.write('#PBS -l nodes=%d:ppn=%d\n' % (nodes_num,threads_num))
    #f_1.write('#PBS -l mem=12gb\n')
    f_1.write('#PBS -l walltime='+walltime+'\n')
    f_1.write('cd '+CADES_run_path+'\n')
    f_1.write('module purge\n')
    f_1.write('module load PE-gnu R\n')
    for loop1 in range(cE_count):
        for loop2 in range(cD_count):
            cE = cE_min + cE_gap * loop1
            cD = cD_min + cD_gap * loop2
            density = 0.16
            ccm_in_file_path = './output/ccm_in_'+matter_type+'_%d_cD_%.2f_cE_%.2f_rho_%.2f' % (neutron_num,cD,cE,density)
            ccm_out_file_path = './output/output/'+matter_type+'_%d_cD_%.2f_cE_%.2f_rho_%.2f.out &' % (neutron_num,cD,cE,density)
            #f_1.write('aprun  -n'+str(int(nodes_num))+' -d'+str(int(threads_num))+' ./prog_ccm.exe '+ccm_in_file_path+' > '+ccm_out_file_path+'\n')
            f_1.write('Rscript ./prog_ccm.exe '+ccm_in_file_path+' > '+ccm_out_file_path+'\n')
            #f_1.write('wait'+'\n')


