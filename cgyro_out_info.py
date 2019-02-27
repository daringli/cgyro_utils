from __future__ import division
import numpy as np





class Cgyro_out_info(object):
    def read_file(self):
        self.z = []
        self.n = []
        self.T = []
        self.m = []
        self.a_Ln = []
        self.a_Lt = []
        self.nu = []
        self.s_a_Ln = []
        self.s_a_Lt = []
                
        
        def read_default(line):
            ls = line.split()
            next_parser = read_default
            """Read line and determine how to read the next line."""
            if ls[0] == "n_theta":
                next_parser = read_resolution
            elif ls[0] == "nc_loc":
                next_parser = read_parallel
            elif ls[0] == "n":
                next_parser = read_rhokx
            elif ls[0] == "D-theta:":
                self.D_theta,self.D_radial,self.D_alpha = [int(x) for x in ls[1:6:2]]
                next_parser = read_default
            elif ls[0] == "up_theta:":
                self.up_theta,self.up_radial,self.up_alpha = [float(x) for x in ls[1:6:2]]
                next_parser = read_default
            elif ls[0] == "C(theta):":
                self.C_theta = float(ls[1])
                next_parser = read_default
            elif ls[0] == "r/a:":
                self.r_a = float(ls[1])
                next_parser = read_default
            elif ls[0] == "R/a:":
                self.R_a,self.shift,self.betae = [float(x) for x in ls[1:6:2]]
                next_parser = read_default
            elif ls[0] == "q:":
                self.q,self.s,self.beta_star = [float(x) for x in ls[1:6:2]]
                next_parser = read_default
            elif ls[0] == "kappa:":
                self.kappa,self.s_kappa,self.lamb_star = [float(x) for x in ls[1:6:2]]
                next_parser = read_default
            elif ls[0] == "delta:":
                self.delta,self.s_delta,self.gamma_e = [float(x) for x in ls[1:6:2]]
                next_parser = read_default
            elif ls[0] == "zeta:":
                self.zeta,self.s_zeta,self.gamma_p = [float(x) for x in ls[1:6:2]]
                next_parser = read_default
            elif ls[0] == "zmag:":
                self.zmag,self.dzmag,self.mach = [float(x) for x in ls[1:6:2]]
            elif ls[0] == "[rho/a]:":
                self.rho_a,self.z_eff,self.w_E_dt = [float(x) for x in ls[1:6:2]]
                next_parser = read_default
            elif ls[0] in [str(i) for i in range(1,11)]:
                i = int(ls[0])
                z,n,T,m,a_Ln,a_Lt,nu,s_a_Ln,s_a_Lt = [float(x) for x in ls[1:]]
                self.z.append(z)
                self.n.append(n)
                self.T.append(T)
                self.m.append(m)
                self.a_Ln.append(a_Ln)
                self.a_Lt.append(a_Lt)
                self.nu.append(nu)
                self.s_a_Ln.append(s_a_Ln)
                self.s_a_Lt.append(s_a_Lt)
                next_parser = read_default
            elif ls[0] == "a[m]:":
                keys = ls[0::2]
                for i,key in enumerate(keys):
                    if key=="a[m]:":
                        self.a = float(ls[1+2*i])
                    elif key=="b_unit[T]:":
                        self.b_unit = float(ls[1+2*i])
                    elif key=="rhos/a:":
                        self.rhos_a = float(ls[1+2*i])
                    else:
                        raise ValueError("Unrecognized key: " + key)
                next_parser = read_default
            elif ls[0] == "n_norm[e19/m^3]:":
                self.n_norm,self.v_norm,self.T_norm = [float(x) for x in ls[1:6:2]]
                next_parser = read_default
            else:
                pass
            return next_parser
                
        def read_resolution(line):
            self.n_theta,self.n_species,self.n_energy,self.n_xi = [int(x) for x in line.split()]
            next_parser = read_default
            return next_parser
        def read_parallel(line):
            self.nc_loc,self.nv_loc,self.nsplit,self.n_MPI,self.n_OMP = [int(x) for x in line.split()]
            next_parser = read_default
            return next_parser
        def read_rhokx(line):
            ls = line.split()
            try:
                self.n_rhokx = int(ls[1])
            except ValueError:
                self.n_rhokx = np.nan
            try:
                self.Delta_rhokx = float(ls[2])
            except:
                self.Delta_rhokx = np.nan
            try:
                self.max_rhokx = float(ls[3])
            except ValueError:
                self.max_rhokx = np.nan
            try:
                self.L_rhox = float(ls[4])
            except ValueError:
                self.L_rhox= np.nan
            next_parser = read_rhoky
            return next_parser
        def read_rhoky(line):
            ls = line.split()
            try:
                self.n_rhoky = int(ls[1])
            except ValueError:
                self.n_rhoky = np.nan
            try:
                self.Delta_rhoky = float(ls[2])
            except:
                self.Delta_rhoky = np.nan
            try:
                self.max_rhoky = float(ls[3])
            except ValueError:
                self.max_rhoky = np.nan
            try:
                self.L_rhoy = float(ls[4])
            except ValueError:
                self.L_rhoy= np.nan
            next_parser = read_default
            return next_parser
        
        parser = read_default
        with open(self.filename,'r') as f:
            for l in f.readlines():
                l = l.strip()
                if len(l)==0:
                    continue
                parser = parser(l)
        self.z = np.array(self.z)
        self.n = np.array(self.n)
        self.T = np.array(self.T)
        self.m = np.array(self.m)
        self.a_Ln = np.array(self.a_Ln)
        self.a_Lt = np.array(self.a_Lt)
        self.nu = np.array(self.nu)
        self.s_a_Ln = np.array(self.s_a_Ln)
        self.s_a_Lt = np.array(self.s_a_Lt)
                
    
    def __init__(self,filename="out.cgyro.info"):
        self.filename = filename
        self.read_file()

if __name__ == "__main__":
    cgoi = Cgyro_out_info("../../../../../../cgyroSimuls/AUG/30701/H_scan/KY100.00_RMIN0.90_base/out.cgyro.info")
    print cgoi.a_Ln
    print cgoi.Delta_rhokx
    print cgoi.n_rhoky
