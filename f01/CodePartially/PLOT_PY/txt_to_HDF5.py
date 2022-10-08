import h5py
import numpy as np
import os
import sys
from glob import glob

cfields = [ 'cooper0',
            'cooper1',
            'cooper3',
            'delta',
            'delta-00',
            'delta-01',
            'delta-02',
            'delta-03',
            'density',
            'jydown',
            'jytot',
            'jyup',
            'jyxspin',
            'jyyspin',
            'jyzspin',
            'jzdown',
            'jztot',
            'jzup',
            'jzxspin',
            'jzyspin',
            'jzzspin',
            'mx',
            'my',
            'mz',
            'ndown',
            'nup']

def read_1D_file(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
        while(len(lines[-1]) < 3):
            lines.pop()

        npx = len(lines)

        width = len(lines[0].strip().split())

        data = np.ndarray((npx, width), dtype=np.float32)

        for i in range(npx):
            data[i,:] = [float(x) for x in lines[i].split()]

        return data

def read_2D_file(npy, npz, filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
        while(len(lines[-1]) < 3):
            lines.pop()

        width = len(lines[0].strip().split())

        data = np.ndarray((npy, npz, width), dtype = np.float32)

        for y in range(npy):
            for z in range(npz):
                idx = y*npz + z
                data[y,z,:] = [float(x) for x in lines[idx].split()]

        return data

def read_single_line(filename):
    with open(filename, 'r') as f:
        data = f.readline()
        values = np.array([float(x) for x in data.split()], np.float32)
        return values

def read_eigen_values(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
        while(len(lines[-1]) < 3):
            lines.pop()

    tokens = lines[-1].strip().split()
    nkx = int(tokens[0])
    N = int(tokens[1])
    data = np.ndarray((nkx, N), dtype = np.float32)
    for kx in range(nkx):
        for col in range(N):
            idx = kx*N + col + 1
            tokens = lines[idx].strip().split()
            data[kx,col] = float(tokens[2])

    return data

if __name__ == '__main__':
    try:
        os.mkdir('hdf5_results')
    except FileExistsError:
        pass

    h5 = h5py.File('hdf5_results/spin_flip.h5', 'w')

    # make all of the Complex Field groups
    for cf in cfields:
        h5.create_group(cf)

    # get npy, npz, y_coord, and z_coord (mm_y and mm_z in the process)
    # mm_y
    data = read_1D_file('mm_y.txt')
    y_coord = data[:,0]
    npy = len(y_coord)
    tmp = np.zeros((2,npy), dtype=np.float32)
    tmp[0,:] = data[:,1]
    h5['mx']['y_scan'] = tmp
    tmp[0,:] = data[:,2]
    h5['my']['y_scan'] = tmp
    tmp[0,:] = data[:,3]
    h5['mz']['y_scan'] = tmp

    # mm_z
    data = read_1D_file('mm_z.txt')
    z_coord = data[:,0]
    npz = len(z_coord)
    tmp = np.zeros((2,npz), dtype=np.float32)

    tmp[0,:] = data[:,1]
    h5['mx']['z_scan'] = tmp
    tmp[0,:] = data[:,2]
    h5['my']['z_scan'] = tmp
    tmp[0,:] = data[:,3]
    h5['mz']['z_scan'] = tmp

    h5['y_coord'] = y_coord
    h5['z_coord'] = z_coord 

    # avgf{i}
    for i in ['0', '1', '3']:
        tmp = read_single_line(f'avgf{i}.txt')
        values = np.array((tmp[3], tmp[0], tmp[1], tmp[2], tmp[4]),
                                dtype=np.float32)
        h5[f'cooper{i}']['averages'] = values

    # avgjcharge
    values = read_single_line('avgjcharge.txt')
    h5['avg_j_charge'] = values

    # delta-??
    deltas = glob('delta-*')
    for d in deltas:
        data = read_2D_file(npy, npz, d)
        data = np.transpose(data,(2,0,1))
        h5[d[0:-4]]['values'] = data

    # density, nup, ndown
    data = read_2D_file(npy, npz, 'density.txt')
    tmp = np.zeros((2, npy, npz), dtype=np.float32)
    tmp[0,:,:] = data[:,:,2]
    h5['density']['values'] = tmp
    tmp[0,:,:] = data[:,:,3]
    h5['nup']['values'] = tmp
    tmp[0,:,:] = data[:,:,4]
    h5['ndown']['values'] = tmp

    # density_y
    data = read_1D_file('density_y.txt')
    tmp = np.zeros((2, npy), dtype = np.float32)
    tmp[0,:] = data[:,1]
    h5['density']['y_scan'] = tmp

    # density_z
    data = read_1D_file('density_z.txt')
    tmp = np.zeros((2, npz), dtype = np.float32)
    tmp[0,:] = data[:,1]
    h5['density']['z_scan'] = tmp

    # divjc
    values = read_single_line('divjc.txt')
    values = np.reshape(values, (2,3))
    h5['jytot']['c_line_integrals'] = values

    # djspin
    values = read_2D_file(npy-1, npz-1, 'djspin.txt')
    tmp = np.ndarray((npy-1, npz-1), dtype=np.float32)
    tmp[:,:] = values[:,:,2]
    h5['divx'] = tmp
    tmp[:,:] = values[:,:,3]
    h5['divy'] = tmp
    tmp[:,:] = values[:,:,4]
    h5['divz'] = tmp

    # dtavg_tot
    values = read_single_line('dtavg_tot.txt')
    h5['dtavg'] = values

    # free0
    values = read_single_line('free0.txt')
    h5['free0'] = values

    # g/g_coord
    values = read_2D_file(npy, npz, 'g.txt')
    values = np.transpose(values,(2,0,1))
    h5['delta']['values'] = values

    # iter
    with open('iter.txt', 'r') as f:
        lines = f.readlines()
        tmp = -1 * np.ones((100,2), dtype=np.float32)
        i = 0
        for l in lines:
            if(len(l) > 3):
                tmp[i,:] = [float(x) for x in l.split()]
            i += 1
        h5['iter'] = tmp


    # W_final
    data = read_eigen_values('W_final.txt')
    h5['final_eigen_values'] = data

    # W_initial_min
    data = read_1D_file('W_initial_min.txt')
    tmp = np.zeros(len(data), dtype = np.float32)
    tmp[:] = data[:,0]
    h5['initial_eigen_values-min'] = tmp

    # W_initial
    data = read_eigen_values('W_initial.txt')
    h5['initial_eigen_values'] = data

