######################################################
# processing the .out data 
######################################################
import numpy as np
import re
import matplotlib
matplotlib.use('PS')
import matplotlib.pyplot as plt
from scipy import interpolate


######################################################
######################################################
#####      read .out data 
######################################################
######################################################
def input_file(file_path):
    with open(file_path,'r') as f_1:
        converge_flag = int(1)
        count = len(open(file_path,'rU').readlines())
        data =  f_1.readlines()
        wtf = re.match('#', 'abc',flags=0)
        for loop1 in range(0,count):
            if ( re.search('E/N', data[loop1],flags=0) != wtf):
                temp_1 = re.findall(r"[-+]?\d+\.?\d*",data[loop1+1])
                energy_per_nucleon = float(temp_1[2])
                return energy_per_nucleon,converge_flag
        print ('No "E/A" found in the file:'+file_path)


######################################################
######################################################
#####     setup all the parameters
######################################################
######################################################
particle_num = 28
neutron_num = 14
cE_min = -1
cE_max = 1
cE_gap = 1
cE_count = int( (cE_max - cE_min) / cE_gap + 1 )
cD_min = -3
cD_max = 3
cD_gap = 2
cD_count = int( (cD_max - cD_min) / cD_gap + 1 )
density_min = 0.12
density_max = 0.20
density_gap = 0.02
density_count = int( (density_max - density_min) / density_gap +1 )
CADES_run_path = '/lustre/or-hydra/cades-nucphys/w01'

data_num = cE_count*cD_count*density_count
raw_data = np.zeros((data_num,7),dtype = np.float)  # cD,cE,density,snm,pnm


######################################################
######################################################
#####     read all the .out file for pnm 
######################################################
######################################################
#matter_type = 'pnm'
#loop4 = 0
#for loop1 in range(cE_count):
#    for loop2 in range(cD_count):
#        cE = cE_min + cE_gap * loop1
#        cD = cD_min + cD_gap * loop2
#        for loop3 in range(density_count):
#            density = density_min + density_gap * loop3
#            ccm_out_file_path = CADES_run_path + '/output/output/'+matter_type+'_%d_cD_%.2f_cE_%.2f_rho_%.2f.out' % (neutron_num,cD,cE,density)
#            raw_data[loop4][4],raw_data[loop4][6] = input_file(ccm_out_file_path) # 4 is the pnm
#            loop4 = loop4 + 1
#print ('total_count='+str(loop4))



######################################################
######################################################
#####      read all the .out file for snm 
######################################################
######################################################
matter_type = 'snm'
loop4 = 0
for loop1 in range(cE_count):
    for loop2 in range(cD_count):
        cE = cE_min + cE_gap * loop1
        cD = cD_min + cD_gap * loop2
        for loop3 in range(density_count):
            density = density_min + density_gap * loop3
            ccm_out_file_path = CADES_run_path + '/output/output/'+matter_type+'_%d_cD_%.2f_cE_%.2f_rho_%.2f.out' % (particle_num,cD,cE,density)
            raw_data[loop4][0] = cD # 0 is the cD
            raw_data[loop4][1] = cE # 0 is the cE
            raw_data[loop4][2] = density # 3 is the density
            raw_data[loop4][3],raw_data[loop4][5] = input_file(ccm_out_file_path) # 3 is the snm
            loop4 = loop4 + 1
if (data_num != loop4):
    print ('total_count error!!!')


######################################################
######################################################
#####      write raw data file 
######################################################
######################################################
file_path = './cD%.2f-%.2f_cE%.2f-%.2f.dat' % (cD_min,cD_max,cE_min,cE_max)
loop1 = 0
#with open (file_path, 'w') as f:
#    f.write('#  cD    cE     density       E/A_snm       E/A_pnm \n')
#    for loop1 in range(data_num):
#        f.write('%.2f   %.2f    %.2f    %.8f    %.8f \n' % (raw_data[loop1][0],raw_data[loop1][1],raw_data[loop1][2],raw_data[loop1][3],raw_data[loop1][4]))

with open (file_path, 'w') as f:
    f.write('#  cD     cE       density     E/A_snm       E/A_pnm     snm_converge_flag \n')
    for loop1 in range(data_num):
        if (raw_data[loop1][5] == 1):
            f.write('% .2f   % .2f    % .2f    % .8f    % .8f         %d\n' % (raw_data[loop1][0],raw_data[loop1][1],raw_data[loop1][2],raw_data[loop1][3],raw_data[loop1][4],raw_data[loop1][5]))
        elif (raw_data[loop1][5] == 0):
            f.write('#% .2f   % .2f    % .2f    % .8f    % .8f         %d\n' % (raw_data[loop1][0],raw_data[loop1][1],raw_data[loop1][2],raw_data[loop1][3],raw_data[loop1][4],raw_data[loop1][5]))


interpol_count = 1000
saturation_point = np.zeros((cE_count*cD_count,2))
kind = "quadratic"

for loop1 in range(cE_count*cD_count):
    x = raw_data[loop1*density_count:loop1*density_count+density_count,2]
    y = raw_data[loop1*density_count:loop1*density_count+density_count,3]
    spldens = np.linspace(density_min,density_max,num=interpol_count) 
    f       = interpolate.interp1d(x,y,kind=kind)
    y_new   = np.zeros(interpol_count)
    y_new   = f(spldens)
    saturation_point[loop1,1] = np.min(y_new)
    temp = spldens[np.where(y_new == np.min(y_new))]
    saturation_point[loop1,0] = temp[0]

x_list = saturation_point[:,0]
y_list = saturation_point[:,1]
fig_saturation_point = plt.figure('fig_saturation_point')

l1 = plt.scatter(x_list,y_list,color = 'b', s = 10, marker='x')

plot_path = 'saturation_point_cD_%.2f_%.2f_cE_%.2f_%.2f.eps' % (cD_min,cD_max,cE_min,cE_max)
plt.savefig(plot_path)
plt.show()


