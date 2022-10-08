#!/usr/bin/env python3

import math
import os

max_elements = 2**31

class Model_Params(object):
    '''A class that holds all of the information that is set by the user'''
    def __init__(   self,
                    restart_flag=0,
                    kfl1=5,
                    kfl2=5,
                    kfl3=5,
                    kfds=150,
                    kfw=50,
                    kfdx=5.0,
                    min_kx_factor=1.0,
                    H0s=[1.0, 1.0, 1.0],
                    omega=0.1,
                    nkx=16,
                    itmax=8,
                    minerror=0.001,
                    time=100.0,
                    coher_len=50.0,
                    rel_temp=0.09,
                    deg_vphi1=0.0,
                    deg_vphi2=0.0,
                    grid_int=4,
                    deg_theta=[0.0, 90.0, 90.0],
                    deg_phi=[0.0, 90.0, 0.0],
                    wtxt_range=[-1.0, 1.0],
                    DOS_range=[-4,4],
                    DOS_inter=1.0,
                    DOS_energy_out = [0, 0, 0],
                    DOS_loc_out_y = [0.0, 0.0, 0.0],
                    DOS_loc_out_z = [0.0, 0.0, 0.0],
                    scatt_str = 0.0,
                    imp_conc = [0.0, 0.0, 0.0, 0.0, 0.0],
                    rand_opt = 0,
                    rand_seed = 12345
                    ):

        self.restart_flag=0
        self.kfl1 = kfl1
        self.kfl2 = kfl2
        self.kfl3 = kfl3
        self.kfds = kfds
        self.kfw = kfw
        self.kfdx = kfdx
        self.min_kx_factor = min_kx_factor
        self.H0s = H0s
        self.omega = omega

        self.kfl = kfl1 + kfl2 + kfl3 + (2*kfds)
        self.max_H0 = max(H0s)
        self.max_Kx = 1.0 + omega + self.max_H0
        self.kx_factor = min_kx_factor # starting point

        self.nkx = nkx
        self.itmax = itmax
        self.minerror = minerror
        self.time = time

        self.coher_len = coher_len
        self.rel_temp = rel_temp

        self.deg_vphi1 = deg_vphi1
        self.deg_vphi2 = deg_vphi2

        self.grid_int = grid_int

        self.deg_theta = deg_theta
        self.deg_phi = deg_phi

        self.wtxt_range=wtxt_range

        self.scatt_str = scatt_str
        self.imp_conc = imp_conc
        self.rand_opt = rand_opt
        self.rand_seed = rand_seed
        
        self.DOS_range = DOS_range
        self.DOS_inter = DOS_inter
        self.DOS_energy_out = DOS_energy_out
        self.DOS_loc_out_y = DOS_loc_out_y
        self.DOS_loc_out_z = DOS_loc_out_z
            
    def calculate_N(self):
        self.Np = math.sqrt(self.max_Kx) * self.kx_factor * self.kfl / math.pi
        self.Np = int(self.Np + 0.5)

        self.Nq = math.sqrt(self.max_Kx) * self.kx_factor * self.kfw / math.pi
        self.Nq = int(self.Nq + 0.5)

    def size_check(self, n_proc_rows, n_proc_cols, clustersize, twodims):
        self.n_proc_rows = n_proc_rows
        self.n_proc_cols = n_proc_cols

        self.calculate_N()
        self.grid_size = n_proc_rows * n_proc_cols
        self.nmaxnum = self.Np * self.Nq
        self.tot_N = 4 * self.nmaxnum
        self.blocking = self.Np
        self.NN = max(self.tot_N, self.blocking, 2)
        self.local_rows = self.tot_N / n_proc_rows
        self.local_cols = self.tot_N / n_proc_cols
        if (twodims) and (self.kx_factor > 2.0):
            clustersize = int(self.blocking * self.blocking)
        self.lrwork =   4*self.tot_N + \
                        max(5*self.NN, self.local_rows * self.local_cols) + \
                        math.ceil(self.tot_N / self.grid_size) * self.NN + \
                        self.tot_N*clustersize

        self.ram_per_rank = (self.lrwork + 4*(self.local_rows * self.local_cols))
        self.ram_per_rank *= 8.0
        self.ram_per_rank /= (1024 * 1024 * 1024)

        return (max_elements > self.lrwork,
                self.nmaxnum % (n_proc_rows * self.blocking) == 0,
                self.nmaxnum % (n_proc_cols * self.blocking) == 0)

    def read_env(self):
        try:
            restart = os.environ['spin_restart']
            if ((restart == 'true') or (restart == 'True') or (restart == 'TRUE')):
                self.restart_flag = 1
        except:
            pass

        try:
            self.kfl1 = float(os.environ['spin_kfl1'])
        except:
            pass

        try:
            self.kfl2 = float(os.environ['spin_kfl2'])
        except:
            pass

        try:
            self.kfl3 = float(os.environ['spin_kfl3'])
        except:
            pass

        try:
            self.kfds = float(os.environ['spin_kfds'])
        except:
            pass

        try:
            self.kfw = float(os.environ['spin_kfw'])
        except:
            pass

        try:
            self.kfdx = float(os.environ['spin_kfdx'])
        except:
            pass

        try:
            self.min_kx_factor = float(os.environ['spin_min_kx_factor'])
            self.kx_factor = self.min_kx_factor
        except:
            pass

        try:
            H0s = os.environ['spin_H0']
            H0s = H0s.replace('[', '')
            H0s = H0s.replace(']', '')
            self.H0s = [float(x) for x in H0s.split(',')]
        except:
            raise
            # pass

        try:
            self.omega = float(os.environ['spin_omega'])
        except:
            pass

        self.max_H0 = max(self.H0s)
        self.max_Kx = 1.0 + float(self.omega) + self.max_H0

        try:
            self.nkx = int(os.environ['spin_nkx'])
        except:
            pass

        try:
            self.itmax = int(os.environ['spin_itmax'])
        except:
            pass

        try:
            self.minerror = float(os.environ['spin_minerror'])
        except:
            pass

        try:
            self.time = float(os.environ['spin_time'])
        except:
            pass

        try:
            self.coher_len = float(os.environ['spin_cl'])
        except:
            pass

        try:
            self.rel_temp = float(os.environ['spin_rel_temp'])
        except:
            pass

        try:
            self.deg_vphi1 = float(os.environ['spin_deg_vphi1'])
        except:
            pass

        try:
            self.deg_vphi2 = float(os.environ['spin_deg_vphi2'])
        except:
            pass

        try:
            self.deg_vphi1 = float(os.environ['spin_vphi1']) * 180.0 / math.pi
        except:
            pass

        try:
            self.deg_vphi2 = float(os.environ['spin_vphi2']) * 180.0 / math.pi
        except:
            pass

        try:
            self.grid_int = int(os.environ['spin_grid_int'])
        except:
            pass

        try:
            theta = os.environ['spin_deg_theta']
            theta = theta.replace('[','')
            theta = theta.replace(']','')
            self.deg_theta = [float(x) for x in theta.split(',')]
        except:
            pass

        try:
            theta = os.environ['spin_theta']
            theta = theta.replace('[','')
            theta = theta.replace(']','')
            self.deg_theta = [(180.0 * float(x) / math.pi) for x in theta.split(',')]
        except:
            pass

        try:
            phi = os.environ['spin_deg_phi']
            phi = phi.replace('[','')
            phi = phi.replace(']','')
            self.deg_phi = [float(x) for x in phi.split(',')]
        except:
            pass

        try:
            phi = os.environ['spin_phi']
            phi = phi.replace('[','')
            phi = phi.replace(']','')
            self.deg_phi = [(180.0 * float(x) / math.pi) for x in phi.split(',')]
        except:
            pass

        try:
            w_range = os.environ['spin_Wtxt_range']
            w_range = w_range.replace('[','')
            w_range = w_range.replace(']','')
            self.wtxt_range = [float(x) for x in w_range.split(',')]
        except:
            pass


        try:
            self.scatt_str = float(os.environ['spin_scatt_str'])
        except:
            pass

        try:
            self.rand_opt = int(os.environ['spin_rand_opt'])
        except:
            pass

        if self.rand_opt == 1:
            try:
                self.rand_seed = int(os.environ['spin_rand_seed'])
            except:
                pass

        try:
            conc = os.environ['spin_imp_conc']
            conc = conc.replace('[','')
            conc = conc.replace(']','')
            self.imp_conc = [float(x) for x in conc.split(',')]
        except:
            pass

        try:
            DOS_r = os.environ['spin_DOS_range']
            DOS_r = DOS_r.replace('[','')
            DOS_r = DOS_r.replace(']','')
            self.DOS_range = [int(x) for x in DOS_r.split(',')]
        except:
            pass

        try:
            self.DOS_inter = float(os.environ['spin_DOS_inter'])
        except:
            pass

        try:
            DOS_e = os.environ['spin_DOS_energy']
            DOS_e = DOS_e.replace('[','')
            DOS_e = DOS_e.replace(']','')
            self.DOS_energy_out = [float(x) for x in DOS_e.split(',')]
        except:
            pass
        try:
            DOS_y = os.environ['spin_DOS_loc_y']
            DOS_y = DOS_y.replace('[','')
            DOS_y = DOS_y.replace(']','')
            self.DOS_loc_out_y = [float(x) for x in DOS_y.split(',')]
        except:
            pass

        try:
            DOS_z = os.environ['spin_DOS_loc_z']
            DOS_z = DOS_z.replace('[','')
            DOS_z = DOS_z.replace(']','')
            self.DOS_loc_out_z = [float(x) for x in DOS_z.split(',')]
        except:
            pass

    def write_input(self, filename):
        with open(filename, 'w') as f:
            f.write(f'restart   = {self.restart_flag}\n')
            
            f.write(f'nkx       = {self.nkx}\n')
            f.write(f'itmax     = {self.itmax}\n')
            f.write(f'minerror  = {self.minerror}\n')

            f.write(f'time      = {self.time}\n')
            f.write(f'omega     = {self.omega}\n')
            f.write(f'rel_temp  = {self.rel_temp}\n')
            f.write(f'cl        = {self.coher_len}\n')

            f.write(f'deg_vphi1 = {self.deg_vphi1}\n')
            f.write(f'deg_vphi2 = {self.deg_vphi2}\n')

            f.write(f'kfl1      = {self.kfl1}\n')
            f.write(f'kfl2      = {self.kfl2}\n')
            f.write(f'kfl3      = {self.kfl3}\n')
            f.write(f'kfds      = {self.kfds}\n')
            f.write(f'kfw       = {self.kfw}\n')
            f.write(f'kfdx       = {self.kfdx}\n')
            f.write(f'kx_factor = {self.kx_factor:.3f}\n')

            f.write(f'grid_int  = {self.grid_int}\n')

            f.write(f'n_proc_rows  = {self.n_proc_rows}\n')
            f.write(f'n_proc_cols  = {self.n_proc_cols}\n')

            f.write('deg_theta = [ ' + \
                        ', '.join([str(x) for x in self.deg_theta]) + ' ]\n')
            f.write('deg_phi   = [ ' + \
                        ', '.join([str(x) for x in self.deg_phi]) + ' ]\n')
            f.write('h0        = [ ' + ', '.join([str(x) for x in self.H0s]) + ' ]\n')

            f.write('wtxt_range   = [ ' + ', '.join([str(x) for x in self.wtxt_range]) + ' ]\n')

            f.write(f'scatt_str = {self.scatt_str}\n')
            f.write('imp_conc   = [ ' + ', '.join([str(x) for x in self.imp_conc]) + ' ]\n')
            f.write(f'rand_opt  = {self.rand_opt}\n')
            if self.rand_opt == 1:
                f.write(f'rand_seed = {self.rand_seed}\n')

            f.write('dos_range    = [ ' + ', '.join([str(x) for x in self.DOS_range]) + ' ]\n')
            f.write(f'dos_inter       = {self.DOS_inter}\n')

            f.write('dos_energies   = [ ' + ', '.join([str(x) for x in self.DOS_energy_out]) + ' ]\n')

            f.write('dos_locs_y  = [ ' + ', '.join([str(x) for x in self.DOS_loc_out_y]) + ' ]\n')
            f.write('dos_locs_z  = [ ' + ', '.join([str(x) for x in self.DOS_loc_out_z]) + ' ]\n')            

if __name__ == '__main__':
    
    params = Model_Params(5, 5, 5, 150, 50, 1.2, (1.0, 0.1, 0.1), 0.1)
    print(params.kfl, params.max_H0)

    while(False in params.size_check(2, 16, 0)):
        params.kx_factor += 0.001

    print(f'kx_factor = {params.kx_factor:.3f}')
    print(f'ram_per_rank = {params.ram_per_rank:.2f}')

    params.read_env()

    params.write_input('test_input.txt')
