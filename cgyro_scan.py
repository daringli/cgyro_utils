from cgyro_simul import Cgyro_simul
from cgyro_collection import Cgyro_collection_from_subdirs

import os

class Cgyro_scan(object):
    def __init__(self,basedir,location):
        self.basedir = basedir
        self.location = location
        self.basesimul = Cgyro_simul(self.basedir)
        self.modified_simuls = []
        os.mkdir(self.location)

    def create_modified(self,name,delta_dict):
        modified_simul = self.basesimul.copy(self.location +"/" + name)
        modified_simul.clean()
        for key in delta_dict:
            setattr(modified_simul.input,key,delta_dict[key])
        return modified_simul

    def add_simulation(self,name,delta_dict):
        self.modified_simuls.append(self.create_modified(name,delta_dict))

    def add_simulations(self,names,delta_dicts):
        for name,delta_dict in zip(names,delta_dicts):
            self.add_simulation(name,delta_dict)

    def run_simulations(self,n=4):
        for simul in self.modified_simuls:
            print "Running in directory: " + simul.dirname
            simul.run(n)

if __name__=="__main__":
    scan_basedir = "test_scan"
    basedir ='../AUG/30701/r_0.9_K_Y_0.2/test'
    names = ["1","2","3"]
    delta_dicts = [{"N_THETA":51,"N_RADIAL":10},{"N_THETA":61,"N_XI":32},{"N_THETA":71,"N_RADIAL":16}]
    cs= Cgyro_scan(basedir,scan_basedir)
    cs.add_simulations(names,delta_dicts)

    cs.run_simulations()
