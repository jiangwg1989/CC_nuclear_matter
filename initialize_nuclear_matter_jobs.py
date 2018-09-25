####################################################
#initialize the jobs for nuclear matter
####################################################
import numpy as np
import os



def output_ccm_in_file(file_path,cD,cE,particle_num,matter_type,density,nmax):
    with open(file_path,'w') as f_1:
        f_1.write('!Chiral order for Deltas(LO = 0,NLO=2,NNLO=3,N3LO=4) and cutoff'+'\n')
        f_1.write('3, 450\n')
        f_1.write('! cE and cD 3nf parameters:'+'\n' )
        f_1.write('%.2f, %.2f\n' % (cE,cD))
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
nmax = 4
particle_num = 132
neutron_num = 66
cE_min = 0
cE_max = 0
cE_gap = 1
cE_count = int( (cE_max - cE_min) / cE_gap + 1 )
cD_min = 0
cD_max = 0
cD_gap = 1
cD_count = int( (cD_max - cD_min) / cD_gap + 1 )
density_min = 0.12
density_max = 0.20
density_gap = 0.02
density_count = int( (density_max - density_min) / density_gap +1 )
titan_run_path = '/lustre/atlas/scratch/w01/nph123/runs/nucmat_test'
nodes_num = 300
threads_num = 16
walltime = '03:00:00'

os.system('mkdir output')


####################################################
#  set up all the ccm_in files for snm 
####################################################
for loop1 in range(cE_count):
    for loop2 in range(cD_count):
        cE = cE_min + cE_gap * loop1
        cD = cD_min + cD_gap * loop2
        for loop3 in range(density_count):
            density = density_min + density_gap * loop3
            file_path = titan_run_path+'/output/ccm_in_snm_%d_cD_%.2f_cE_%.2f_rho_%.2f' % (particle_num,cD,cE,density)
            output_ccm_in_file(file_path,cD,cE,particle_num,'snm',density,nmax)


####################################################
#  set up all the ccm_in files for pnm 
####################################################
for loop1 in range(cE_count):
    for loop2 in range(cD_count):
        cE = cE_min + cE_gap * loop1
        cD = cD_min + cD_gap * loop2
        for loop3 in range(density_count):
            density = density_min + density_gap * loop3
            file_path = titan_run_path+'/output/ccm_in_pnm_%d_cD_%.2f_cE_%.2f_rho_%.2f' % (neutron_num,cD,cE,density)
            output_ccm_in_file(file_path,cD,cE,neutron_num,'pnm',density,nmax)



####################################################
#  set up job script for snm 
####################################################
matter_type = 'snm'
file_path = './output/snm_job.script'
with open(file_path,'w') as f_1:
    f_1.write('#PBS -N nuclearmatter'+'\n')
    f_1.write('#PBS -j eo'+'\n')
    f_1.write('#PBS -q batch'+'\n')
    f_1.write('#PBS -l walltime='+walltime+',nodes='+str(int(nodes_num))+'\n')
    f_1.write('#PBS -A NPH123'+'\n\n')
    f_1.write('export OMP_NUM_THREADS='+str(int(threads_num))+'\n\n')
    f_1.write('cd '+titan_run_path+'\n')
    for loop1 in range(cE_count):
        for loop2 in range(cD_count):
            cE = cE_min + cE_gap * loop1
            cD = cD_min + cD_gap * loop2
            for loop3 in range(density_count):
                density = density_min + density_gap * loop3
                ccm_in_file_path = './output/ccm_in_'+matter_type+'_%d_cD_%.2f_cE_%.2f_rho_%.2f' % (particle_num,cD,cE,density)
                ccm_out_file_path = './output/output/'+matter_type+'_%d_cD_%.2f_cE_%.2f_rho_%.2f.out &' % (particle_num,cD,cE,density)
                f_1.write('aprun  -n'+str(int(nodes_num))+' -d'+str(int(threads_num))+' ./prog_ccm.exe '+ccm_in_file_path+' > '+ccm_out_file_path+'\n')
                f_1.write('wait'+'\n')





####################################################
#  set up all the ccm_in files for pnm 
####################################################
matter_type = 'pnm'
file_path = './output/pnm_job.script'
with open(file_path,'w') as f_1:
    f_1.write('#PBS -N nuclearmatter'+'\n')
    f_1.write('#PBS -j eo'+'\n')
    f_1.write('#PBS -q batch'+'\n')
    f_1.write('#PBS -l walltime='+walltime+',nodes='+str(int(nodes_num))+'\n')
    f_1.write('#PBS -A NPH123'+'\n\n')
    f_1.write('export OMP_NUM_THREADS='+str(int(threads_num))+'\n\n')
    f_1.write('cd '+titan_run_path+'\n')
    for loop1 in range(cE_count):
        for loop2 in range(cD_count):
            cE = cE_min + cE_gap * loop1
            cD = cD_min + cD_gap * loop2
            for loop3 in range(density_count):
                density = density_min + density_gap * loop3
                ccm_in_file_path = './output/ccm_in_'+matter_type+'_%d_cD_%.2f_cE_%.2f_rho_%.2f' % (neutron_num,cD,cE,density)
                ccm_out_file_path = './output/output/'+matter_type+'_%d_cD_%.2f_cE_%.2f_rho_%.2f.out &' % (neutron_num,cD,cE,density)
                f_1.write('aprun  -n'+str(int(nodes_num))+' -d'+str(int(threads_num))+' ./prog_ccm.exe '+ccm_in_file_path+' > '+ccm_out_file_path+'\n')
                f_1.write('wait'+'\n')


