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
    diff_orig_V=np.linspace(1.e5, 6.e10,101)
    diff_K=np.linspace(1.e5, 6.e10,101)

    P=1.e5
    grand.set_state(P, T)
    g.set_state(P,T)

    import burnman.eos.modified_tait as mt
    a, b, c = mt.tait_constants(grand.params)
    a0, b0, c0 = mt.tait_constants(g.endmembers[0][0].params)
    a1, b1, c1 = mt.tait_constants(g.endmembers[1][0].params)


    pressures=np.linspace(1.e5, 5.e10, 101)
    V_series=np.linspace(1.e5, 1.e10, 101)
    V_series2=np.linspace(1.e5, 1.e10, 101)
    V_actual=np.linspace(1.e5, 1.e10, 101)

    for i, P in enumerate(pressures):
        v=grand.params['V_0']
        V_series[i]=v-a*b*c*v*P #+ (1./2.)*a*b*b*c*(c+1)*v*P*P - (1./6.)*a*b*b*b*c*(c+1.)*(c+2.)*v*P*P*P + (1./24.)*a*b*b*b*b*c*(c+1.)*(c+2.)*(c+3.)*v*P*P*P*P
        V_series2[i]=v-a*b*c*v*P + (1./2.)*a*b*b*c*(c+1)*v*P*P - (1./6.)*a*b*b*b*c*(c+1.)*(c+2.)*v*P*P*P + (1./24.)*a*b*b*b*b*c*(c+1.)*(c+2.)*(c+3.)*v*P*P*P*P - (1./120.)*a*b*b*b*b*b*c*(c+1.)*(c+2.)*(c+3.)*(c+4.)*v*P*P*P*P*P
        V_actual[i]=v*(1.-a*(1.-pow(1.+b*P,-c)))

    plt.plot( pressures/1.e9, V_series, marker='.', linewidth=1., label='V series')
    plt.plot( pressures/1.e9, V_series2, marker='.', linewidth=1., label='V series 2')
    plt.plot( pressures/1.e9, V_actual, marker='.', linewidth=1., label='V full')
    plt.title("Grossular-andradite join (asymmetric model)")
    plt.ylabel("Volume")
    plt.xlabel("pressure (GPa)")
    plt.legend(loc='upper left')
    plt.show()


    VoK =grand.params['V_0']/grand.params['K_0']
    VoK0=g.endmembers[0][0].params['V_0']/g.endmembers[0][0].params['K_0']
    VoK1=g.endmembers[1][0].params['V_0']/g.endmembers[1][0].params['K_0']

    #b*(c+1)
    dVdP=VoK - 0.5*(VoK0 + VoK1)
    d2VdP2=0.#VoK*b*(c+1.)/grand.params['Kprime_0'] - 0.5*(VoK0*b0*(c0+1.)/g.endmembers[0][0].params['Kprime_0'] + VoK1*b1*(c1+1.)/g.endmembers[1][0].params['Kprime_0'])
    d2VdP2=0.#VoK*b*(c+1.) - 0.5*(VoK0*b0*(c0+1.) + VoK1*b1*(c1+1.))

    d3VdP3=0.#(1./2.)*VoK*b*b*(c+1.)*(c+2.) - (0.5/2.)*(VoK0*b0*b0*(c0+1.)*(c0+2.) + VoK1*b1*b1*(c1+1.)*(c1+2.)) # need to check this equation

    d4VdP4=0.#(1./6.)*VoK*b*b*b*(c+1.)*(c+2.)*(c+3.) - (0.5/6.)*(VoK0*b0*b0*b0*(c0+1.)*(c0+2.)*(c0+3.) + VoK1*b1*b1*b1*(c1+1.)*(c1+2.)*(c1+3.)) # need to check this equation

    d5VdP5=0.#(1./24.)*VoK*b*b*b*b*(c+1.)*(c+2.)*(c+3.)*(c+4.) - (0.5/24.)*(VoK0*b0*b0*b0*b0*(c0+1.)*(c0+2.)*(c0+3.)*(c0+4.) + VoK1*b1*b1*b1*b1*(c1+1.)*(c1+2.)*(c1+3.)*(c1+4.))

    d6VdP6=0.#(1./120.)*VoK*b*b*b*b*b*(c+1.)*(c+2.)*(c+3.)*(c+4.)*(c+5.) - (0.5/24.)*(VoK0*b0*b0*b0*b0*b0*(c0+1.)*(c0+2.)*(c0+3.)*(c0+4.)*(c0+5.) + VoK1*b1*b1*b1*b1*b1*(c1+1.)*(c1+2.)*(c1+3.)*(c1+4.)*(c1+5.))

    #d2VdP2=VoK*(1.+grand.params['Kprime_0'])/(grand.K_T*grand.params['Kprime_0']) - 0.5*(g.endmembers[0][0].V*(1.+g.endmembers[0][0].params['Kprime_0'])/(g.endmembers[0][0].K_T*g.endmembers[0][0].K_T*g.endmembers[0][0].params['Kprime_0']) + g.endmembers[1][0].V*(1.+g.endmembers[1][0].params['Kprime_0'])/(g.endmembers[1][0].K_T*g.endmembers[1][0].K_T*g.endmembers[1][0].params['Kprime_0']))

    for i,P in enumerate(pressures):
        grand.set_state(P, T)
        g.set_state(P,T)
        g_V[i] = g.V
        diff_V[i] = (grand.V-g.V + dVdP*P - d2VdP2*P*P + d3VdP3*P*P*P -  d4VdP4*P*P*P*P + d5VdP5*P*P*P*P*P - d6VdP6*P*P*P*P*P*P)/g.V*100.
        diff_orig_V[i] =  (grand.V-g.V)/g.V*100.
        diff_K[i] = (grand.K_T-g.K_T - grand.V/VoK + g.V/(0.5*(VoK0 + VoK1))  )/g.K_T*100.

    plt.plot( pressures/1.e9, diff_V, marker='.', linewidth=1., label='V diff, Grandite endmember - 50% solid solution')
    plt.plot( pressures/1.e9, diff_orig_V, marker='.', linewidth=1., label='V diff_original, Grandite endmember - 50% solid solution')
    plt.plot( pressures/1.e9, diff_K, marker='.', linewidth=1., label='K, Grandite endmember - 50% solid solution')
    plt.title("Grossular-andradite join (asymmetric model)")
    plt.ylabel("% difference")
    plt.xlabel("pressure (GPa)")
    plt.ylim(-1, 1)
    plt.legend(loc='upper left')
    plt.show()



# The other way ... 
    class gr_andr(burnman.SolidSolution):
        def __init__(self):
            self.name='Symmetric grossular-andradite garnet'
            endmembers = [[minerals.HP_2011_ds62.gr(), '[Ca]3[Al]2Si3O12'],[minerals.HP_2011_ds62.andr(), '[Ca]3[Fe]2Si3O12']]
            enthalpy_interaction=[[0.0e3]]
            volume_interaction=[[-0.0e-5]]
            burnman.SolidSolution.__init__(self, endmembers, \
                          burnman.solutionmodel.SymmetricRegularSolution(endmembers, enthalpy_interaction, volume_interaction) )

    g=gr_andr()

    V_excess=-0.2e-5
    c_average=np.average([g.endmembers[i][0].params['V_0']*g.endmembers[i][0].params['K_0'] for i in range(g.n_endmembers)])
    K_excess=c_average / grand.params['V_0'] - np.average([g.endmembers[i][0].params['K_0'] for i in range(g.n_endmembers)])

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
                'V_0': np.average([g.endmembers[i][0].params['V_0'] for i in range(g.n_endmembers)]) + V_excess/4. ,
                'Cp': [np.average([g.endmembers[i][0].params['Cp'][j] for i in range(g.n_endmembers)]) for j in range(4)] ,
                'a_0': np.average([g.endmembers[i][0].params['a_0'] for i in range(g.n_endmembers)]) ,
                'K_0': np.average([g.endmembers[i][0].params['K_0'] for i in range(g.n_endmembers)]) + K_excess ,
                'Kprime_0': np.average([g.endmembers[i][0].params['Kprime_0'] for i in range(g.n_endmembers)]) ,
                'Kdprime_0': np.average([g.endmembers[i][0].params['Kdprime_0'] for i in range(g.n_endmembers)]) ,
                'n': sum(formula.values()),
                'molar_mass': formula_mass(formula, atomic_masses)}
            burnman.Mineral.__init__(self)


    grand=grandite()

    pressures=np.linspace(1.e5, 1.2e11, 101)
    V_ex=np.linspace(1.e5, 1.e10, 101)
    K_diff=np.linspace(1.e5, 1.e10, 101)
    for i, P in enumerate(pressures):
        grand.set_state(P, T)
        g.set_composition([0.5, 0.5])
        g.set_state(P,T)
        V_excess = grand.V - g.V
        dPdV_excess = 0.0 - grand.K_T/grand.V + g.K_T/g.V
        g_K_T = (g.V + V_excess)*((g.K_T/(g.V)-dPdV_excess))
        V_ex[i]=V_excess
        K_diff[i]=grand.K_T - g_K_T # should be zero


    plt.plot( pressures/1.e9, V_ex, marker='.', linewidth=1., label='V excess')
    plt.plot( pressures/1.e9, K_diff, marker='.', linewidth=1., label='V excess')
    plt.title("Grossular-andradite join (asymmetric model)")
    plt.ylabel("Volume")
    plt.xlabel("pressure (GPa)")
    plt.legend(loc='upper left')
    plt.show()
