# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

""" Generates a text table with mineral properties. Run 'python table.py latex' to write a tex version of the table to mytable.tex """


import os, sys, numpy as np, matplotlib.pyplot as plt
#hack to allow scripts to be placed in subdirectories next to burnman:
if not os.path.exists('burnman') and os.path.exists('../burnman'):
    sys.path.insert(1,os.path.abspath('..')) 

import burnman
import inspect

from burnman import minerals
from burnman import tools

if __name__ == "__main__":    

    names = set()
    phasenames = []

    libs = dir(minerals)
    for l in libs:
        mineralgroup = getattr(minerals, l)
        
        if mineralgroup.__class__.__name__ == "module":
            for m in dir(mineralgroup):
                mineral = getattr(mineralgroup, m)
                if inspect.isclass(mineral) and mineral!=burnman.Mineral and issubclass(mineral, burnman.Mineral):
                    #print mineral.__module__ + mineral.__name__
                    name1 = str(mineralgroup.__name__) +"." + str(m)
                    name = name1.replace("burnman.minerals.","")
                    #name = mineral.__module__.replace("minlib_","").replace("burnman.","").replace("minerals.","") + "." + mineral.__name__
                    #print mineral, name, name1
                    if not name in names:
                        names.add(name)
                        try:
                            x=mineral()
                            phasenames.append((name,x))
                        except:
                            print "Could not create '%s'" % name

    
    params = ['V_0','K_0','Kprime_0','G_0','Gprime_0','molar_mass','n','Debye_0','grueneisen_0','q_0','eta_s_0']
    
    
    def create_list(name, mineral):
        ownname = mineral.to_string().replace("'","").replace("burnman.minerals.","")
        if name!=ownname:
            name=name+" ("+ownname+")"
        row = [ name ]
        for param in params:
            
            row.append(str(p.params[param] if param in p.params else ""))
        return row
    
    
    
    table = [ ['Name'] + params ]
    tablel = []

    sortedlist = sorted(phasenames, key=lambda x: x[0])

    for (name,p) in sortedlist:
        p.set_state(1e9,300)
        row = create_list(name,p)
        table.append(row)
        tablel.append(row)    
 
    if (len(sys.argv)==1):
        tools.pretty_print_table(table, False)
    elif sys.argv[1]=="tab":
        tools.pretty_print_table(table, True)
