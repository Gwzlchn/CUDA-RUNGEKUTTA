import numpy as np
import os



open_file_name = 'K1K2K3K4_add.dat'
save_file_name = 'wzl_dx_1_'

os.chdir("./wzl_data")  


os.listdir()
der_list = ['x1','px1','x2','px2','y1','py1',
            'y2','py2','z1','pz1','z2','pz2']
raw_dx_1 = np.loadtxt(open_file_name)


for i ,value in enumerate(der_list):
    now_dat_name = save_file_name + value + ".dat"
    np.savetxt(now_dat_name,raw_dx_1[i:-1:12,:],
    fmt='%20.12lf',newline='\n')


