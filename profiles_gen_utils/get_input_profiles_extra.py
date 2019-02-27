from __future__ import division
import numpy as np

class Input_profiles_extra(object):

    
    def get_metadata(self):
        with  open(self.filename,'r') as f:
            # first pass, we read comments line for meta data
            metadata = []
            for l in f.readlines():
                if l[0] == '#':
                    _l = l[1:-1].strip()
                    if _l != '':
                        metadata.append(_l)
        return metadata

    def get_Nentries(self):
        try:
            Nentries= int(self.metadata[-1].split()[-1])
        except ValueError as e:
            raise ValueError(str(e) + ": " + self.filename)
        return Nentries
        
    @staticmethod
    def get_index_range(intervalstring):
        l=intervalstring.split('-')
        if len(l) == 1:
            return [int(l[0]),int(l[0])]
        elif len(l) == 2:
            return [int(l[0]),int(l[1])]
        else:
            raise ValueError("indexstring should be an integer 'X' or an interval on the form 'X-Y'")

        
    def get_quantity_dict(self):
        quantity_dict = {}
        for md in self.metadata[1:-1]:
            mds = md.split()
            quantity_dict[mds[1]] = Input_profiles_extra.get_index_range(mds[0])
        return quantity_dict
    
    def __init__(self,filename):
        self.filename = filename
        self.metadata = self.get_metadata()
        self.array_length = self.get_Nentries()
        self.quantity_dict = self.get_quantity_dict()
        

    def get_quantity(self,quantity_name):
        indices = self.quantity_dict[quantity_name]
        start_location = (indices[0]-1)*self.array_length
        stop_location = (indices[1]) *self.array_length
        data = np.zeros(self.array_length*(indices[1]-(indices[0]-1)))

        with  open(self.filename,'r') as f:
            # second pass, we ignore comments
            i = 0
            for l in f.readlines():
                if l[0] != '#':
                    if (i>= start_location) and (i< stop_location):
                        data[i-start_location] = float(l.strip())
                    i = i + 1
        return data


def get_dVdr(filename):
    """Function to get volp (dV/dr) from input.profiles.extra"""
    IPE = Input_profiles_extra(filename)
    return IPE.get_quantity("EXPRO_volp(:)")

def get_Bt0(filename):
    """Function to get B_tor at theta=0 (T) from input.profiles.extra"""
    IPE = Input_profiles_extra(filename)
    return IPE.get_quantity("EXPRO_bt0(:)")

def get_Bp0(filename):
    """Function to get B_pol at theta=0 (T) from input.profiles.extra"""
    IPE = Input_profiles_extra(filename)
    return IPE.get_quantity("EXPRO_bp0(:)")

def get_dlnnedr(filename):
    """Function to get logarithmic electron density gradient from input.profiles.extra"""
    IPE = Input_profiles_extra(filename)
    ret = IPE.get_quantity("EXPRO_dlnnedr(:)")
    return ret


def get_dlntedr(filename):
    """Function to get logarithmic electron temperature gradient from input.profiles.extra"""
    IPE = Input_profiles_extra(filename)
    ret = IPE.get_quantity("EXPRO_dlntedr(:)")
    return ret

def get_dlnnidr(filename):
    """Function to get logarithmic ion density gradient from input.profiles.extra"""
    IPE = Input_profiles_extra(filename)
    ret = IPE.get_quantity("EXPRO_dlnnidr(1:10,:)")
    Nspecies = 10
    Nr = int(len(ret)/Nspecies)
    return ret.reshape((Nspecies,Nr))

def get_dlntidr(filename):
    """Function to get logarithmic ion temperature gradient from input.profiles.extra"""
    IPE = Input_profiles_extra(filename)
    ret = IPE.get_quantity("EXPRO_dlntidr(1:10,:)")
    Nspecies = 10
    Nr = int(len(ret)/Nspecies)
    return ret.reshape((Nspecies,Nr))

def get_ni_new(filename):
    """Function to get ion density of the first ion when QN has been enforced, from input.profiles.extra"""
    IPE = Input_profiles_extra(filename)
    ret = IPE.get_quantity("EXPRO_ni_new(:)")
    return ret

def get_dlnnidr_new(filename):
    """Function to get logarithmic ion density gradient of the first ion when QN has been enforced, from input.profiles.extra """
    IPE = Input_profiles_extra(filename)
    ret = IPE.get_quantity("EXPRO_dlnnidr_new(:)")
    return ret

def get_rhos(filename):
    IPE = Input_profiles_extra(filename)
    ret = IPE.get_quantity("EXPRO_rhos(:)")
    return ret
        
if __name__=="__main__":
    IPE = Input_profiles_extra("../JET89454_2018-08-16/input.profiles.extra")
    print IPE.metadata
    print IPE.array_length
    print IPE.quantity_dict
    q= IPE.get_quantity("EXPRO_nuee(:)")
    print q[0]
    print q[-1]
    print repr(get_dVdr("../JET89454_2018-08-16/input.profiles.extra"))

    t=get_dlntidr("../JET89454_2018-08-16/input.profiles.extra")
    print t[0]
    #print repr(t[0])
    #print repr(t[1])
