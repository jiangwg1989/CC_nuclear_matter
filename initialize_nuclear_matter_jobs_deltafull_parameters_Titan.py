####################################################
#initialize the jobs for nuclear matter
####################################################
import numpy as np
import os
import re

def start_from_break_point(file_path,line_num):
    with open(file_path,'r') as f_1:
        data = f_1.readlines()
        wtf = re.match('#', 'abc',flags=0)
        temp_1 = re.findall(r"[-+]?\d+\.?\d*",data[line_num-1])
        DNNLO450_input[0:6] = temp_1[0:6]
        temp_1 = re.findall(r"[-+]?\d+\.?\d*",data[line_num])
        DNNLO450_input[6:12] = temp_1[0:6]
        temp_1 = re.findall(r"[-+]?\d+\.?\d*",data[line_num+1])
        DNNLO450_input[12:17] = temp_1[0:5]




def output_ccm_in_file(file_path,cD,cE,LEC_ci,c1s0,c3s1,cnlo,particle_num,matter_type,density,nmax):
    with open(file_path,'w') as f_1:
        f_1.write('!Chiral order for Deltas(LO = 0,NLO=2,NNLO=3,N3LO=4) and cutoff'+'\n')
        f_1.write('3, 450\n')
        f_1.write('! cE and cD 3nf parameters:'+'\n' )
        f_1.write('%.8f, %.8f\n' % (cE,cD))
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
nmax = 4
particle_num = 132
neutron_num = 66
cE_min = 0.03000471
cE_max = 0.03000471
cE_gap = 1 
cE_count = int( (cE_max - cE_min) / cE_gap + 1 )
cD_min = 0.68668995
cD_max = 0.68668995
cD_gap = 1
cD_count = int( (cD_max - cD_min) / cD_gap + 1 )
density_min = 0.12
density_max = 0.20
density_gap = 0.02
density_count = int( (density_max - density_min) / density_gap +1 )
titan_run_path = '/lustre/atlas/scratch/w01/nph123/runs/nucmat_test'
nodes_num = 300
threads_num = 16
pnm_walltime = '00:30:00'
snm_walltime = '03:00:00'
LEC_ci = np.zeros(4)
c1s0   = np.zeros(3)
c3s1   = np.zeros(3)
cnlo   = np.zeros(7)
vec_input = np.zeros(15)


DNNLO450_input = np.array([-0.74, #ci
            -0.49,
            -0.65,
             0.96,
            -0.33813946528363, #Ct_1S0np
            -0.33802308849055, #Ct_1S0nn
            -0.33713665635790, #Ct_1S0pp 
            -0.22931007944035, #Ct_3S1(pp,nn,np)
             2.47658908242147, #C_1S0
             0.64555041107198, #C_3P0
            -1.02235931835913, #C_3P1
            -0.87020321739728, #C_3P2
            -0.02854100307153, #C_1P1
             0.69595320984261, #C_3S1
             0.35832984489387, #C_3S1-3D1
             0.79,             #cD
             0.017])           #cE

start_from_break_point('./pounders.out',1)

vec_input[0:15] = DNNLO450_input[0:15]
cE_min = DNNLO450_input[16]
cE_max = DNNLO450_input[16]
cD_min = DNNLO450_input[15]
cD_max = DNNLO450_input[15]

print ('vec_input='+str(vec_input))

#vec_input[0]  = -0.74   
#vec_input[1]  = -0.49
#vec_input[2]  = -0.65
#vec_input[3]  = 0.96
#vec_input[4]  = -0.33982408 
#vec_input[5]  = -0.33847639
#vec_input[6]  = -0.3387161
#vec_input[7]  = -0.2860546
#vec_input[8]  = 2.52024031
#vec_input[9]  = 0.64477023
#vec_input[10] = -0.90129773
#vec_input[11] = -0.88829717
#vec_input[12] = 0.18364315
#vec_input[13] = 1.31945937
#vec_input[14] = 0.72331122


LEC_ci[0] = vec_input[0]
LEC_ci[1] = vec_input[1]
LEC_ci[2] = vec_input[2]
LEC_ci[3] = vec_input[3]
c1s0[0]   = vec_input[6]
c1s0[1]   = vec_input[4]
c1s0[2]   = vec_input[5]
c3s1[0]   = vec_input[7]
c3s1[1]   = vec_input[7]
c3s1[2]   = vec_input[7]
cnlo[0]   = vec_input[8]
cnlo[1]   = vec_input[9]
cnlo[2]   = vec_input[12]
cnlo[3]   = vec_input[10]
cnlo[4]   = vec_input[13]
cnlo[5]   = vec_input[14]
cnlo[6]   = vec_input[11]


#LEC_ci[0] = -0.74 
#LEC_ci[1] = -0.49
#LEC_ci[2] = -0.65
#LEC_ci[3] =  0.96
#c1s0[0]   =  -0.337137 #c1s0pp 
#c1s0[1]   =  -0.338139 #c1s0np
#c1s0[2]   =  -0.338023 #c1s0nn
#c3s1[0]   =  -0.229310 #c3s1
#c3s1[1]   =  -0.229310
#c3s1[2]   =  -0.229310
#cnlo[0]   =   2.476589  #1s1
#cnlo[1]   =   0.645550  #3p0
#cnlo[2]   =  -0.028541  #1p1
#cnlo[3]   =  -1.022359  #3p1
#cnlo[4]   =   0.695953  #3s1
#cnlo[5]   =   0.358330  #3s13d1
#cnlo[6]   =  -0.870203  #3p2  


#LEC_ci[0] =-0.728365308352
#LEC_ci[1] =-0.487307350938
#LEC_ci[2] =-0.650986995096
#LEC_ci[3] = 0.974417059712
#c1s0[0]   =-0.336931595688
#c1s0[1]   =-0.338011916329
#c1s0[2]   =-0.337290791852
#c3s1[0]   =-0.228842794936
#c3s1[1]   =-0.228842794936
#c3s1[2]   =-0.228842794936
#cnlo[0]   =2.468442159468
#cnlo[1]   =0.731389531683
#cnlo[2]   =-0.047052011951
#cnlo[3]   =-0.940503433140
#cnlo[4]   =0.679254714232
#cnlo[5]   =0.369756349283
#cnlo[6]   =-0.851354838008


os.system('mkdir output')
os.system('cp /ccs/home/w01/programs_titian/CCM_kspace_deltafull/prog_ccm.exe /ccs/home/w01/run/runs/nucmat_test/')
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
            file_path = titan_run_path+'/output/ccm_in_pnm_%d_cD_%.2f_cE_%.2f_rho_%.2f' % (neutron_num,cD,cE,density)
            output_ccm_in_file(file_path,cD,cE,LEC_ci,c1s0,c3s1,cnlo,neutron_num,'pnm',density,nmax)


####################################################
#  set up job script for snm 
####################################################
matter_type = 'snm'
file_path = './output/snm_job.script'
with open(file_path,'w') as f_1:
    f_1.write('#PBS -N nuclearmatter'+'\n')
    f_1.write('#PBS -j eo'+'\n')
    f_1.write('#PBS -q batch'+'\n')
    f_1.write('#PBS -l walltime='+snm_walltime+',nodes='+str(int(nodes_num))+'\n')
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
    f_1.write('#PBS -l walltime='+pnm_walltime+',nodes='+str(int(nodes_num))+'\n')
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


