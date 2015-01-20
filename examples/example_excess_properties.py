# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

"""

excess_properties
----------------------
    
This script includes the code required to recreate figures for the paper on
the modification to the asymmetric regular solid solution model of
Holland and Powell (2003)

"""

import os, sys, numpy as np, matplotlib.pyplot as plt
#hack to allow scripts to be placed in subdirectories next to burnman:
if not os.path.exists('burnman') and os.path.exists('../burnman'):
    sys.path.insert(1,os.path.abspath('..'))

import burnman
from burnman import minerals
from burnman.processchemistry import *

if __name__ == "__main__":
    '''
    First, let's pick a starting pressure and temperature
    '''
    P=1.e5
    T=298.15

    class mg_fe_ca_garnet(burnman.SolidSolution):
        def __init__(self):
            self.name='Asymmetric pyrope-almandine-grossular garnet'
            endmembers = [[minerals.HP_2011_ds62.py(), '[Mg]3[Al]2Si3O12'],[minerals.HP_2011_ds62.alm(), '[Fe]3[Al]2Si3O12'],[minerals.HP_2011_ds62.gr(), '[Ca]3[Al]2Si3O12']]
            alphas=[1.0, 1.0, 2.7]
            enthalpy_interaction=[[2.5e3,30.1e3],[1.0e3]]
            volume_interaction=[[0.,0.4e-5],[0.122e-5]]
            burnman.SolidSolution.__init__(self, endmembers, \
                          burnman.solutionmodel.AsymmetricRegularSolution(endmembers, alphas, enthalpy_interaction, volume_interaction) )


    g=mg_fe_ca_garnet()
    comp=np.linspace(0.0,1.0,101)
    g_K_T = np.empty_like(comp)
    g_K_T_KV = np.empty_like(comp)

    g.set_composition([1.0, 0.0, 0.0])
    g.set_state(P,T)
    KV0=g.K_T*g.V
    g.set_composition([0.0, 0.0, 1.0])
    g.set_state(P,T)
    KV1=g.K_T*g.V

    for i,c in enumerate(comp):
        molar_fraction=[1.0-c, 0., c]
        g.set_composition(molar_fraction)
        g.set_state(P,T)
        g_K_T[i] = g.K_T

        KV=KV0*(1.0-c) + KV1*c
        g_K_T_KV[i] = KV/g.V

    # see Du et al., 2015
    for endmember in g.endmembers:
        print endmember[0].params['name'], endmember[0].params['K_0'], endmember[0].params['Kprime_0'] 

    plt.plot( comp, g_K_T/1.e9, 'r-', linewidth=1., label='Excess volume')
    plt.plot( comp, g_K_T_KV/1.e9, 'b-', linewidth=1., label='KV=constant')
    plt.plot( [comp[0], comp[100]], [g_K_T[0]/1.e9, g_K_T[100]/1.e9], 'g-',  linewidth=1., label='linear')
    plt.title("Pyrope-grossular join (asymmetric model)")
    plt.ylabel("Bulk modulus (calculated, GPa)")
    plt.xlabel("Grossular fraction")
    plt.legend(loc='lower left')
    plt.show()




    # Volume data from Huckenholz et al., 1974, plot in Woodland and Ross, 1994.
    # Skew towards grossular
    # Bulk modulus data from Babuska et al., 1978, -ve, skew towards andradite
    # Also Fan et al., 2011 (+ve, ~ symmetric)
    # Also Lacivita et al, 2014 (slightly -ve Kex for slight +ve Vex)
    class gr_andr(burnman.SolidSolution):
        def __init__(self):
            self.name='Symmetric grossular-andradite garnet'
            endmembers = [[minerals.HP_2011_ds62.gr(), '[Ca]3[Al]2Si3O12'],[minerals.HP_2011_ds62.andr(), '[Ca]3[Fe]2Si3O12']]
            enthalpy_interaction=[[0.0e3]]
            volume_interaction=[[-0.2e-5]]
            burnman.SolidSolution.__init__(self, endmembers, \
                          burnman.solutionmodel.SymmetricRegularSolution(endmembers, enthalpy_interaction, volume_interaction) )



    g=gr_andr()
    comp=np.linspace(0.0,1.0,101)
    g_Vex = np.empty_like(comp)

    for i,c in enumerate(comp):
        molar_fraction=[1.0-c, c]
        g.set_composition(molar_fraction)
        g.set_state(P,T)
        g_Vex[i] = g.excess_volume

    plt.plot( comp, g_Vex*1.e6, 'r-', linewidth=1., label='Excess volume (cm^3)')
    plt.title("Grossular-andradite join (asymmetric model)")
    plt.ylabel("Excess volume")
    plt.xlabel("Andradite fraction")
    plt.legend(loc='lower left')
    plt.show()

    P=1.e10
    g.set_composition([1.0, 0.0])
    g.set_state(P,T)
    KV0=g.K_T*g.V

    g.set_composition([0.0, 1.0])
    g.set_state(P,T)
    KV1=g.K_T*g.V

    g_K_T_KV = np.empty_like(comp)
    g_K_T_test= np.empty_like(comp)

    for i,c in enumerate(comp):

        molar_fraction=[1.0-c, c]
        g.set_composition(molar_fraction)
        g.set_state(P,T)
        g_K_T[i] = g.K_T

        KV=KV0*(1.0-c) + KV1*c
        g_K_T_KV[i] = KV/g.V


        # Test
        a=sum([g.endmembers[e][0].V*molar_fraction[e] for e in range(2)])
        b=sum([g.endmembers[e][0].V*molar_fraction[e]/g.endmembers[e][0].K_T for e in range(2)])

        Kex=000000.e9
        v0=g.excess_volume
        v1=1.0e-11*v0
        v2=-2.e-22*v0
        Vex=v0 + v1*P
        Vprimeex=v1
        g_K_T_test[i]=(a+Vex)/(b+Vprimeex)

    #print g_K_T_test
    plt.plot( comp, g_K_T/1.e9, 'r-', linewidth=1., label='Regular solution model')
    plt.plot( comp, g_K_T_KV/1.e9, 'b-', linewidth=1., label='KV=constant')
    plt.plot( comp, g_K_T_test/1.e9, 'b--', linewidth=1., label='test')
    plt.plot( [comp[0], comp[100]], [g_K_T[0]/1.e9, g_K_T[100]/1.e9], 'g-',  linewidth=1., label='linear')
    plt.title("Grossular-andradite join (asymmetric model)")
    plt.ylabel("Bulk modulus (GPa)")
    plt.xlabel("Andradite fraction")
    plt.legend(loc='lower left')
    plt.show()


    atomic_masses=read_masses()

    class grandite (burnman.Mineral):
        def __init__(self):
            formula='Ca3.0Fe1.0Al1.0Si3.0O12.0'
            formula = dictionarize_formula(formula)
            self.params = {
                'name': 'grandite',
                'formula': formula,
                'equation_of_state': 'hp_tmt',
                'H_0': -5769100.0 ,
                'S_0': 316.4 ,
                'V_0': (0.00013204+0.00012535)/2. -0.05e-5 ,
                'Cp': [638.6, 0.0, -4955100.0, -3989.2] ,
                'a_0': (2.86e-05+2.2e-05)/2. ,
                'K_0': (1.588e+11+1.72e+11)/2. ,
                'Kprime_0': (5.68+5.53)/2. ,
                'Kdprime_0': (-3.6e-11+-3.2e-11)/2. ,
                'n': sum(formula.values()),
                'molar_mass': formula_mass(formula, atomic_masses)}
            burnman.Mineral.__init__(self)


    grand=grandite()
    g.set_composition([0.5, 0.5])

    pressures=np.linspace(1.e5, 1.e11,101)
    g_V=np.linspace(1.e5, 6.e10,101)
    grand_V=np.linspace(1.e5, 6.e10,101)
    diff_V=np.linspace(1.e5, 6.e10,101)
    diff_K=np.linspace(1.e5, 6.e10,101)

    P=1.e5
    grand.set_state(P, T)
    g.set_state(P,T)

    dVdP=grand.V/grand.K_T - 0.5*(g.endmembers[0][0].V/g.endmembers[0][0].K_T + g.endmembers[1][0].V/g.endmembers[1][0].K_T)
    d2VdP2=grand.V*(1.+grand.params['Kprime_0'])/(grand.K_T*grand.K_T*grand.params['Kprime_0']) - 0.5*(g.endmembers[0][0].V*(1.+g.endmembers[0][0].params['Kprime_0'])/(g.endmembers[0][0].K_T*g.endmembers[0][0].K_T*g.endmembers[0][0].params['Kprime_0']) + g.endmembers[1][0].V*(1.+g.endmembers[1][0].params['Kprime_0'])/(g.endmembers[1][0].K_T*g.endmembers[1][0].K_T*g.endmembers[1][0].params['Kprime_0']))

    for i,P in enumerate(pressures):
        grand.set_state(P, T)
        g.set_state(P,T)
        g_V[i] = g.V
        diff_V[i] = (grand.V-g.V+dVdP*P-d2VdP2*P*P)/g.V*100.
        diff_K[i] = (grand.K_T-g.K_T)/g.K_T*100.

    plt.plot( pressures/1.e9, diff_V, marker='.', linewidth=1., label='V, Grandite endmember - 50% solid solution')
    plt.plot( pressures/1.e9, diff_K, marker='.', linewidth=1., label='K, Grandite endmember - 50% solid solution')
    plt.title("Grossular-andradite join (asymmetric model)")
    plt.ylabel("% difference")
    plt.xlabel("pressure (GPa)")
    plt.legend(loc='upper left')
    plt.show()
