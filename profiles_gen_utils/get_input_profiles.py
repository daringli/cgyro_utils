from __future__ import division
import numpy as np
from scipy.constants import e

class Input_profiles(object):


    def read_header_data(self):        
        findstring1 = "RADIAL GRIDPOINTS"
        findstring2 = "N_EXP"
        findstring3 = "IPCCW"
        findstring4 = "BTCCW"
        N1 = None
        N2 = None
        IPCCW = None
        BTCCW = None
        f1 = False
        f2 = False
        with  open(self.filename,'r') as f:
            # first pass, we read comment lines for meta data
            for l in f.readlines():
                if l.find(findstring1) != -1:
                    try:
                        N1 = int(l[1:].strip().split()[3])
                    except ValueError:
                        pass
                    else:
                        f1= True
                elif l.find(findstring2) != -1:
                    try:
                        N2 = int(l[1:].strip().split('=')[1])
                    except ValueError:
                        pass
                    else:
                        f2 = True
                elif  l.find(findstring3) != -1:
                    IPCCW = int(l[1:].strip().split()[2])
                elif  l.find(findstring4) != -1:
                    BTCCW = int(l[1:].strip().split()[2])
                else:
                    pass
            else:
                if (N1 is None) and (N2 is None):
                    raise ValueError("File does not contain information on the number of gridpoints")
                elif N1 is not None:
                    N2 = N1
                elif N2 is not None:
                    N1 = N2
                else:
                   raise IllegalStateException("This should not happen") 
                
        self.N_radial_gridpoints = N1
        self.IPCCW = IPCCW
        self.BTCCW = BTCCW
        
    
    def read_data(self):
        metadata = []
        data = {}
        header_end_string = "rho(-)"
        past_header = False
        read_data = False
        with  open(self.filename,'r') as f:
            for l in f.readlines():
                if not past_header:
                    # first we search for the string signifying
                    # that we are past header lines
                    if l.find(header_end_string) != -1:
                        # header found!
                        # but we are not reading data yet
                        past_header = True
                if past_header:
                    # we are past header, will only read data and metadata now
                    if read_data:
                        tmp = l.strip().split()
                        if len(tmp) != len(metadata):
                            read_data = False
                            continue
                        for i,d in enumerate(tmp):
                            data[metadata[i]][iradial] = float(d)
                        iradial = iradial + 1
                    else:
                        metadata = l[1:].split()
                        for md in metadata:
                            data[md] = np.zeros(self.N_radial_gridpoints)
                        iradial = 0
                        read_data = True
        self.data = data
                    

    def __init__(self,filename):
        self.filename = filename
        self.read_header_data()
        self.read_data()

def get_rhoN(filename):
    return Input_profiles(filename).data["rho(-)"]

def get_q(filename):
    return Input_profiles(filename).data["q(-)"]

def get_rmaj(filename):
    return Input_profiles(filename).data["rmaj(m)"]

def get_zmag(filename):
    return Input_profiles(filename).data["zmag(m)"]

def get_kappa(filename):
    return Input_profiles(filename).data["kappa(-)"]

def get_delta(filename):
    return Input_profiles(filename).data["delta(-)"]

def get_zeta(filename):
    return Input_profiles(filename).data["zeta(-)"]

def get_r(filename):
    return Input_profiles(filename).data["rmin(m)"]

def get_ne(filename):
    return Input_profiles(filename).data["ne(10^19/m^3)"] * 1e19

def get_ni(filename,number):
    return Input_profiles(filename).data["ni_"+str(number)+"(10^19/m^3)"] * 1e19

def get_Te(filename):
    return Input_profiles(filename).data["Te(keV)"] * 1e3 * e

def get_Ti(filename,number):
    return Input_profiles(filename).data["Ti_"+str(number)+"(keV)"] * 1e3 * e

def get_psi(filename):
    return Input_profiles(filename).data["polflux(Wb/rad)"]

def get_IPCCW(filename):
    return Input_profiles(filename).IPCCW
    
if __name__=="__main__":
    IP= Input_profiles("../JET89454_2018-08-16/input.profiles")
    print IP.N_radial_gridpoints
    print IP.data
