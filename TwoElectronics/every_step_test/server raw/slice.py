import numpy as np



der_list = ['px1','fx1','px2','fx2','py1','fy1',
            'py2','fy2','pz1','fz1','pz2','fz2']
raw_dx_1 = np.loadtxt('K1K2K3K4_dx_10.dat')


for i ,value in enumerate(der_list):
    now_dat_name = "wzl_dx_10_" + value + ".dat"
    np.savetxt(now_dat_name,raw_dx_1[i:-1:12,:],
    fmt='%20.12lf',newline='\n')


