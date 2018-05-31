import numpy as np



der_list = ['g1','f1','g2','f2','g3','f3',
            'g4','f4','g5','f5','g6','f6']
raw_dx_1 = np.loadtxt('K1K2K3K4.dat')


for i ,value in enumerate(der_list):
    now_dat_name = "wzl_dx_1_" + value + ".dat"
    np.savetxt(now_dat_name,raw_dx_1[i:-1:12,:],
    fmt='%20.12lf',newline='\n')


