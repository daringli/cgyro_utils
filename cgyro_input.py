from __future__ import division
import numpy as np
import subprocess

from warnings import warn

class Cgyro_input(object):

    def copy_file(self,new_filename):
        subprocess.call("cp " + self.input_filename + " " + new_filename, shell=True)
        return Cgyro_input(new_filename)
    
    def change_attr(self,attr,value):
        # Taken from perfectInputFile and modified to ignore fortran namelist features like
        # groups, arrays, strings
        # Warning: this command will fail silently if the pattern is not found in the file. Sorry about that.
        # Warning: case insensitive
        subprocess.call("sed -i '{s/^"+attr+"=.*/"+attr+"="+str(value)+"/I }' "+self.input_filename, shell=True)

    #def __getattribute__(self, attr):
    #    return object.__getattribute__(self, attr)

    def __setattr__(self, attr, value):
        if attr in self.file_attributes:
            self.change_attr(attr,value)
        if attr not in self.__dict__:
            warn("Cgyro_input is setting an attribute that does not already exist. Is this really what you wanted? Note that file attributes must already exist in the underyling cgyro.input file for Cgyro_input to be able to set them.")
        object.__setattr__(self, attr, value)
    
    def read_inputfile(self):
        def trim_line(line):
            line = line.split('#',1)[0]
            line = line.strip()
            return line

        object.__setattr__(self, "file_attributes", [])
        
        with open(self.input_filename,'r') as input_file:
            for line in input_file.readlines():
                line=trim_line(line)
                if len(line)>0:
                    try:
                        attr,value = line.split('=',1)
                    except ValueError:
                        raise ValueError("Line is not on form: X=Y, line='" + line + "'")
                    #TODO cast attr to upper or lower if input.cgyro is not case sensitive
                    attr = attr.strip()
                    try:
                        value =  float(value.strip())
                    except ValueError:
                        raise ValueError("Attribute +'" +attr+"' has value that could not be cast to float. Value: '" + value + "'")
                    object.__setattr__(self,attr,value)
                    self.file_attributes.append((attr))

    def __init__(self,input_filename):
        # use object's setattr since we override our own
        object.__setattr__(self, "file_attributes", [])
        object.__setattr__(self,"input_filename",input_filename)
        self.read_inputfile()

if __name__=="__main__":
    ci = Cgyro_input("input.cgyro")
    print ci.__dict__.keys()
    ci.Z_EFF_METHOD=2
    ci.Zeff_method="bjs"
    print ci.Zeff_method
    ci2 = ci.copy_file("input.cgyro2")
    print ci2.input_filename
    print ci2.Z_EFF_METHOD
