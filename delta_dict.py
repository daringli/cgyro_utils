# This file specifies a functions for creating "delta_dicts", which are dictionaries that specify parameters (keys) and parameter values (values).
# The "delta" in the name comes from the fact that these key-values descrube how to modify parameters in a cgyro simulation.
# List of delta_dicts are used for creating scans, these are called dictlists.

def add_dicts(dictionary,dict2add):
    """Merge delta_dicts into one dictionary, provided that there are no duplicate keys"""
    dictionary = dict(dictionary)
    for key in dict2add:
        if not key in dictionary.keys():
            dictionary[key]=dict2add[key]
        else:
            print "Warning: diction already contains key from dict2add, key: '" + str(key) + "' thus cannot be added!"
    return dictionary

def mesh2d_dictlists(dictlist1,dictlist2):
    """Takes 2 dictlists and merges every entry in dictlist 1 with every entery in dictlist 2."""
    ret = []
    for d1 in dictlist1:
        for d2 in dictlist2:
            ret.append(add_dicts(d1,d2))
    return ret
            
def mesh_dictlists(delta_dict_lists):
    """Takes n dictlists and merges every entry in every dictlist."""
    ret = [{}]
    for dl in delta_dict_lists:
        ret=mesh2d_dictlists(ret,dl)
    return ret
        
def linear_dictlist(params,values):
    """Generates a linear dictlist from a key and a list of values"""
    if not type(params) is tuple:
        params = (params,)
    if not type(values) is tuple:
        values = (values,)
        
    assert(len(params) == len(values))
    Nvalues = len(values[0])
    assert(all([len(v) == Nvalues for v in values]))
    
    linear_dictlist = []
    for i in range(Nvalues):
        delta_dict = {}
        for param,value in zip(params,values):
            delta_dict[param] = value[i]
        linear_dictlist.append(delta_dict)
    return linear_dictlist
        
if __name__=="__main__":
    params = "N_THETA"
    values = [1,51,101]

    dictlist1= linear_dictlist(params,values)
    print dictlist1
    
    dictlist11 = add_dicts(dictlist1[0],{"gurka":2})
    print dictlist11

    dictlist12 = mesh2d_dictlists(dictlist1,[{"gurka":2},{"tomat":3}])
    print dictlist12

    print "13"
    dictlist13 = mesh2d_dictlists(dictlist1,[{}])
    print dictlist13

    print "14"
    dictlist14 = mesh_dictlists((dictlist1,[{"N_XI":1,"N_X":1},{"N_XI":2}],[{"N_R":0, "N_Y":1},{"N_R":2}]))
    print dictlist14

    
    params = ("N_THETA","N_XI")
    values = ([1,51,101],[16,32,96])
    dictlist2= linear_dictlist(params,values)
    print dictlist2
    
