from __future__ import division
import numpy as np
from cgyro_input import Cgyro_input

from shutil import copytree
import os
import subprocess

from get_input_profiles_extra import get_dlnnedr,get_dlntedr,get_dlnnidr,get_dlntidr, get_ni_new, get_dlnnidr_new, get_rhos
from get_input_profiles import get_r, get_psi, get_Ti, get_ni, get_Te, get_ne, get_IPCCW
from diff_matrix import diff_matrix

from cgyro_out_info import Cgyro_out_info

from EXPRO_2_perfect_geometry import EXPRO_geometry

import matplotlib.pyplot as plt

elecharge = 1.60217662e-19
pi = np.pi
eps0 = 8.854187817e-12
u=1.660539040e-27
mD=2.013553212745*u # deuterium mass


class Cgyro_simul(object):

    def read_output1(self,filename,col):
        """Reads a column 'col' of a certain type of CGYRO output used for:
        out.cgyro.time
        out.cgyro.freq
        ...
        Where col is the 0-indexed column"""

        ret = []
        try:
            with open(filename,'r') as f:
                for l in f.readlines():
                    ret.append(float(l.split()[col]))
        except IOError:
            print "File '" + filename + "' does not exists for '" + self.dirname + "'!"
        return np.array(ret)

    def copy(self,destination):
        copytree(self.dirname,destination)
        return Cgyro_simul(destination)

    def clean(self):
        this_dir = os.getcwd()
        os.chdir(self.dirname)
        subprocess.call("cgyro -clean", shell=True)
        os.chdir(this_dir)

    def run(self,n=4,):
        this_dir = os.getcwd()
        os.chdir(self.dirname)
        platform = os.environ['GACODE_PLATFORM']
        if platform == 'KEBNEKAISE':
            runstr= "gacode_qsub -code cgyro -n" + str(n) + " -nomp 2 -numa 2 -s -e"
        else:
            runstr= "cgyro -n " + str(n) + " -e ."
        print ">" + runstr
        subprocess.call(runstr, shell=True)
        os.chdir(this_dir)

        
    @property
    def t(self):
        time_filename = "out.cgyro.time"
        return self.read_output1(self.fn(time_filename),col=0)

    @property
    def converged(self):
        ret = False
        out_filename = "out.cgyro.info"
        try:
            with open(self.fn(out_filename),'r') as f:
                for l in f.readlines():
                    if l[:4] == "EXIT":
                        if l.find("converged")>= 0:
                            ret= True
        except IOError:
            print "File 'out.cgyro.info' does not exists for '" + self.dirname + "'!"
        return ret

    
    @property
    def error(self):
        ret = False
        out_filename = "out.cgyro.info"
        try:
            with open(self.fn(out_filename),'r') as f:
                for l in f.readlines():
                    if l[:5] == "ERROR":
                        ret = True
        except IOError:
            print "File 'out.cgyro.info' does not exists for '" + self.dirname + "'!"
        return ret

    
    @property
    def omega(self):
        freq_filename = "out.cgyro.freq"
        return self.read_output1(self.fn(freq_filename),col=0)

    @property
    def gamma(self):
        freq_filename = "out.cgyro.freq"
        return self.read_output1(self.fn(freq_filename),col=1)

    def avg_omega(self,i=100):
        return np.average(self.omega[i:])
        
    def avg_gamma(self,i=100):
        return np.average(self.gamma[i:])

    @property
    def last_omega(self):
        try:
            return self.omega[-1]
        except IndexError:
            print "File 'out.cgyro.freq' is empty for '" + self.dirname + "'!"
        
    @property
    def last_gamma(self):
        try:
            return self.gamma[-1]
        except IndexError:
            print "File 'out.cgyro.freq' is empty for '" + self.dirname + "'!"
        
        
    @property
    def KY(self):
        return self.input.KY

    @property
    def RMIN(self):
        return self.input.RMIN

    @property
    def r(self):
        return self.RMIN * self.a

    
    @property
    def N_RADIAL(self):
        return self.input.N_RADIAL

    @property
    def N_THETA(self):
        return self.input.N_THETA
    
    @property
    def DELTA_T(self):
        return self.input.DELTA_T

    
    def __init__(self,dirname):
        self.dirname = dirname
        self.input = Cgyro_input(self.fn("input.cgyro"))
        self.has_read_out_info = False
        self.has_read_geometry = False

    
    def read_out_info(self):
        self.out_info = Cgyro_out_info(self.fn("out.cgyro.info"))
        self.has_read_out_info = True
        
    def read_geometry(self):
        # this is seperated from __init__ since it is rather slow
        # and only used for a few quantities at the moment
        # it is called automatically when it is needed for the first time
        self.expro_geometry = EXPRO_geometry(self.dirname)
        self.has_read_geometry = True
        
    def fn(self,fn):
        """Turns a filename into a filename in the right dir"""
        try:
            ret = self.dirname + "/" + fn
        except TypeError as e:
            raise TypeError(str(e) + "dirname : " +str(self.dirname))
        return ret
    
    @property
    def a_Ln(self):
        if not self.has_read_out_info:
            self.read_out_info()
        return self.out_info.a_Ln
        
    @property
    def a_LT(self):
        if not self.has_read_out_info:
            self.read_out_info()
        return self.out_info.a_Lt

    @property
    def r_profile(self):
        return get_r(self.fn("input.profiles"))

    @property
    def psi_profile(self):
        return get_psi(self.fn("input.profiles"))

    
    @property
    def dpsidr_profile(self):
        psi = self.psi_profile
        r = self.r_profile
        #D=diff_matrix(r[1],r[-2],len(r)-2,order=4)
        #dpsidr = np.dot(D,psi[1:-1])
        #return np.concatenate(([dpsidr[0]],dpsidr,[dpsidr[-1]]))
        drdpsi = np.diff(psi[1:])/np.diff(r[1:])
        return np.concatenate(([drdpsi[0],drdpsi[0]],drdpsi))
        
        
    
    @property
    def T_profile(self):
        T = []
        for i in range(1,self.Nspecies):
            T.append(get_Ti(self.fn("input.profiles"),i))
        T.append(get_Te(self.fn("input.profiles")))
        return np.array(T).transpose()

    @property
    def Te_profile(self):
        return self.T_profile[:,self.electron_index]

    @property
    def n_profile(self):
        n = []
        for i in range(1,self.Nspecies):
            n.append(get_ni(self.fn("input.profiles"),i))
        n.append(get_ne(self.fn("input.profiles")))
        return np.array(n).transpose()

    @property
    def IPCCW(self):
        return get_IPCCW(self.fn("input.profiles"))
    
    @property
    def RBp_profile(self):
        if self.has_read_geometry == False:
            self.read_geometry()
        
        RBp2= (self.expro_geometry.R*self.expro_geometry.B)**2-self.expro_geometry.I[:,np.newaxis]**2
        # the sign of Bp is given by -IPCCW
        # see: http://gafusion.github.io/doc/geometry.html
        Ntheta = RBp2.shape[1]
        #plt.plot(np.sqrt(RBp2)[50])
        #plt.show()
        return -self.IPCCW * np.sqrt(RBp2)[:,Ntheta//2+1]

    @property
    def dlnndr_profile(self):
        dlnndr = []
        dlnndr.append(get_dlnnidr_new(self.fn("input.profiles.extra")))
        for i in range(1,self.Nspecies-1):
            dlnndr.append(get_dlnnidr(self.fn("input.profiles.extra"))[i])
        dlnndr.append(get_dlnnedr(self.fn("input.profiles.extra")))
        # minus sign since dlnnidr from input.profiles.extra in fact is -dln(n)/dr
        return -np.array(dlnndr).transpose()

    @property
    def dlnTdr_profile(self):
        dlnTdr = []
        for i in range(0,self.Nspecies-1):
            dlnTdr.append(get_dlntidr(self.fn("input.profiles.extra"))[i])
        dlnTdr.append(get_dlntedr(self.fn("input.profiles.extra")))
        # minus sign since dlnnidr from input.profiles.extra in fact is -dln(T)/dr
        return -np.array(dlnTdr).transpose()

    @property
    def a_Ln_profile(self):
        return -self.a*self.dlnndr_profile
    
    @property
    def a_LT_profile(self):
        return -self.a*self.dlnTdr_profile

    @property
    def Er_estimate_profile(self):
        # E_r ~ T_i \nabla ln p_i/e
        i=0 #assumes bulk ions are the first species
        return self.T_profile[:,i] * (self.dlnTdr_profile[:,i] + self.dlnndr_profile[:,i])/(self.Z[i] * elecharge)
    
    @property
    def omega0_estimate_profile(self):
        return self.Er_estimate_profile/(self.dpsidr_profile)

    @property
    def gamma_E_profile(self):
        r = self.r_profile
        # exclude first and last grid point, since those are potentially artificial
        D=diff_matrix(r[1],r[-2],len(r)-2,order=4)
        omega0 = self.omega0_estimate_profile[1:-1]
        ddr_omega0 = np.dot(D,omega0)
        ddr_omega0 = np.concatenate(([ddr_omega0[0]],ddr_omega0,[ddr_omega0[-1]]))
        return -r*ddr_omega0/(self.q)
    
    @property
    def gamma_E(self):
        return  np.interp(self.r,self.r_profile,self.gamma_E_profile)

    @property
    def gamma_E_norm_profile(self):
        return self.gamma_E_profile*self.a/self.c_s_profile

    @property
    def gamma_E_norm(self):
        return self.gamma_E*self.a/self.c_s
                
    @property
    def T(self):
        if not self.has_read_out_info:
            self.read_out_info()
        return self.out_info.T * self.out_info.T_norm *1e3 * elecharge

    @property
    def n(self):
        if not self.has_read_out_info:
            self.read_out_info()
        return self.out_info.n * self.out_info.n_norm * 1e19

    @property
    def electron_index(self):
        return np.argwhere(self.Z==-1)[0][0]
    
    @property
    def ne(self):
        return self.n[self.electron_index]

    @property
    def Te(self):
        return self.T[self.electron_index]

    @property
    def Z(self):
        if not self.has_read_out_info:
            self.read_out_info()
        return self.out_info.z

    @property
    def m(self):
        if not self.has_read_out_info:
            self.read_out_info()
        # placeholder for actual mnorm
        # this will be correct for PROFILE_MODE=2
        mnorm = mD
        return self.out_info.m * mnorm

    @property
    def me(self):
        return self.m[self.electron_index]

    @property
    def nu_norm(self):
        nu_normed = 0 #nu/(c_s/a)
        # nu defined as \sum_{b} 1/tau_{ab} with tau_{ab} in Appendix of Helander+Sigmar
        # i.e. nu_ab = 4/(3sqrt(pi)) \hat{nu}_ab

        for b in range(self.Nspecies):
            nu_normed = nu_normed + self.Z[b]**2 * self.n[b]
        nu_normed = (self.Z**2/(self.T**(3/2)*(self.m)**(1/2))) * (3* self.a *np.sqrt(mD) *elecharge**4*self.logLambda/(np.sqrt(self.Te *2 *pi) *2*eps0**2)) * nu_normed
        return nu_normed

    @property
    def nuee_norm(self):
        # defined as in the CGYRO manual
        # but in SI units e_G^4 = e_SI^4/(4*pi eps0)**2
        nu = self.ne*elecharge**4*self.logLambda/((2*self.Te)**(3/2) * 4 * pi * eps0**2 * np.sqrt(self.me))
        return nu/(self.c_s/self.a)
        
    @property
    def Nspecies(self):
        return len(self.Z)
    
    @property
    def a(self):
        if not self.has_read_out_info:
            self.read_out_info()
        return self.out_info.a

    @property
    def c_s(self):
        return np.sqrt(self.Te/mD)

    @property
    def c_s_profile(self):
        return np.sqrt(self.Te_profile/mD)

    @property
    def q(self):
        if not self.has_read_out_info:
            self.read_out_info()
        return self.out_info.q

    @property
    def s_q(self):
        if not self.has_read_out_info:
            self.read_out_info()
        return self.out_info.s

    @property
    def R(self):
        if not self.has_read_out_info:
            self.read_out_info()
        return self.out_info.R_a*self.a

    @property
    def shift(self):
        if not self.has_read_out_info:
            self.read_out_info()
        return self.out_info.shift

    @property
    def kappa(self):
        if not self.has_read_out_info:
            self.read_out_info()
        return self.out_info.kappa

    @property
    def s_kappa(self):
        if not self.has_read_out_info:
            self.read_out_info()
        return self.out_info.s_kappa

    @property
    def delta(self):
        if not self.has_read_out_info:
            self.read_out_info()
        return self.out_info.delta

    @property
    def s_delta(self):
        if not self.has_read_out_info:
            self.read_out_info()
        return self.out_info.s_delta

    @property
    def zeta(self):
        if not self.has_read_out_info:
            self.read_out_info()
        return self.out_info.zeta

    @property
    def s_zeta(self):
        if not self.has_read_out_info:
            self.read_out_info()
        return self.out_info.s_zeta

    @property
    def zmag(self):
        if not self.has_read_out_info:
            self.read_out_info()
        return self.out_info.zmag

    @property
    def dzmag(self):
        if not self.has_read_out_info:
            self.read_out_info()
        return self.out_info.dzmag
    

    @property
    def Er_estimate(self):
        # E_r ~ T_i \nabla ln p_i/e
        i = 0
        self.T[i] * (a_LT[i] + a_Ln[i])/(self.a * self.Z[i] * elecharge)

    @property
    def logLambda(self):
        # from EXPRO_compute_derived.f90
        #loglam(:) = 24.0 - log(sqrt(EXPRO_ne(:)*1e13)/(EXPRO_te(:)*1e3))
        # which is the second electron-ion formula in NRL
        loglam = 24.0 - np.log(np.sqrt(self.ne/1e6)/(self.Te/elecharge))
        return loglam
    
if __name__=="__main__":
    import matplotlib.pyplot as plt
    cs1= Cgyro_simul('../AUG/30701/r_0.9_K_Y_0.2/N1')

    css = [cs1]
    colors = ['b']
    
    fig,ax = plt.subplots(nrows=2,sharex=True)

    istart = 50
    for i,cs in enumerate(css):
        ax[0].plot(cs.t[istart:],cs.omega[istart:],color=colors[i])
        ax[1].plot(cs.t[istart:],cs.gamma[istart:],color=colors[i])

        if False:
            istartavg = 50
            ax[0].plot([cs.t[istartavg],cs.t[-1]],[cs.avg_omega(istartavg)]*2,color=colors[i],linestyle='dashed')
            ax[1].plot([cs.t[istartavg],cs.t[-1]],[cs.avg_gamma(istartavg)]*2,color=colors[i],linestyle='dashed')
        
    
    plt.xlabel(r"$t/[a/c_s]$")
    ax[0].set_ylabel(r"$\omega/[c_s/a]$")
    ax[1].set_ylabel(r"$\gamma/[c_s/a]$")
    plt.savefig("30701_r0.9_ky0.2_resolution.pdf")
    
