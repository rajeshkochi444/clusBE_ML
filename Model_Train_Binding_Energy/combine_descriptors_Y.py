import numpy as np

metal_list = ['Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu' ]


for metal in metal_list:
    fname = 'y_all_adsorb_1eV_random_'+ metal + '.npy'
    if metal == 'Sc':
        a = np.load(fname)
        print(metal, a.shape)
    else:
        b = np.load(fname)
        print(metal, b.shape)
        a = np.concatenate((a,b))
        print("stacked a", a.shape)



np.save('y_all_adsorb_firstrowTM.npy', a)


np.load('y_all_adsorb_firstrowTM.npy').shape

