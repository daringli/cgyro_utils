from __future__ import division
import numpy as np
from cgyro_input import Cgyro_input

from shutil import copytree
import os
import subprocess

def read_output1(filename,col):
    """Reads a column 'col' of a certain type of CGYRO output used for:
    out.cgyro.time
    out.cgyro.freq
    ...
    Where col is the 0-indexed column
    """

    ret = []
    with open(filename,'r') as f:
        for l in f.readlines():
            ret.append(float(l.split()[col]))
    return np.array(ret)


class Cgyro_simul(object):

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
            pass
        else:
            runstr= "cgyro -n " + str(n) + " -e ."
            print ">" + runstr
            subprocess.call(runstr, shell=True)
        os.chdir(this_dir)

        
    @property
    def t(self):
        time_filename = "out.cgyro.time"
        return read_output1(self.fn(time_filename),col=0)

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
        return read_output1(self.fn(freq_filename),col=0)

    @property
    def gamma(self):
        freq_filename = "out.cgyro.freq"
        return read_output1(self.fn(freq_filename),col=1)

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

    def fn(self,fn):
        """Turns a filename into a filename in the right dir"""
        return self.dirname + "/" + fn
        
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
    
