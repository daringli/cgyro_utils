from __future__ import division

from cgyro_simul import Cgyro_simul

import os

class Cgyro_collection(object):
    def __getattribute__(self, attr):
        # NOTE: can only get attributes directly exposed on Cgyro_simul,
        # i.e. not attributes on Cgyro_simul.input
        return [getattr(s,attr) for s in object.__getattribute__(self,"simuls")]

    def __setattr__(self, attr, value):
        for s in object.__getattribute__(self,"simuls"):
            setattr(s,attr,value)

    def __len__(self):
        return len(object.__getattribute__(self,"simuls"))
    
    def __init__(self,dirlist):
        object.__setattr__(self,"dirlist",dirlist)
        object.__setattr__(self,"simuls",[])
        for d in dirlist:
            object.__getattribute__(self,"simuls").append(Cgyro_simul(d))


def Cgyro_collection_from_subdirs(superdir):
    subdirs = next(os.walk(superdir))[1]
    subdirs.sort()
    subdirs = [superdir + "/" + d for d in subdirs]
    return Cgyro_collection(subdirs)

            
if __name__=="__main__":
    dirlist =['../AUG/30701/r_0.9_K_Y_0.2/N1','../AUG/30701/r_0.9_K_Y_1.0/N1','../AUG/30701/r_0.9_K_Y_10.0/N1','../AUG/30701/r_0.9_K_Y_100.0/N1']
    
    cc=Cgyro_collection(dirlist)
    print cc.ky
    print cc.last_gamma
    print cc.last_omega
