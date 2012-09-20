#/usr/bin/env python
"""
Tests for nuclide_data

"""

import numpy as np
import uncertainties as unc

import nuclide_data


def test_isotopes():
    """Do we have the correct isotopes for select elements?"""

    zs = [1, 8, 56, 95]
    isos = [range(1, 8), 
            range(12, 29),
            range(112, 154),
            range(230, 250),]
    for z, ref_iso in zip(zs, isos):
        iso = nuclide_data.isotopes[z]
        assert ref_iso == iso


def test_isomers():
    """Do we have the correct isomers for select nuclides?"""

    nuclides = [ (4,10), (7,14), (52,115), (81,185), (115,290) ]
    isomer_states = [ [0., 3.3680, 5.9584, 6.1793],
                      [0., 8.49, 8.964, 9.129],
                      [0., 0.02, 0.2801],
                      [0., 0.4548], 
                      [0.,] ]

    for n, ref_isomers in zip(nuclides, isomer_states):
        isomers = nuclide_data.isomers(*n)
        assert ref_isomers == isomers

def ufloat_equiv(a, b):
    try:
        return ( np.allclose([a.nominal_value], [b.nominal_value]) 
             and np.allclose([a.std_dev()], [b.std_dev()]) )
    except AttributeError:
        return a == b 
        #return np.allclose([a], [b]) 

def test_data():
    """
    Do we have the correct weight, abundance, and half-life for select nuclides?
    """
    nuclides = [ (1,1), (19, 49), (63, 148), (96, 240) ]
    weights  = [ unc.ufloat("1.00782503207(10)"),
                 unc.ufloat("48.967450(80)"),
                 unc.ufloat("147.918086(11)"),
                 unc.ufloat("240.0555295(25)"), ]
    abundances  = [ unc.ufloat("0.999885(70)"), 0., 0., 0. ]
    half_lifes = [ 0., 1.26E+00, 4.71E+06, 2.33E+06]
    stable = [True, False, False, False]

    for i in range(len(nuclides)):
        d = nuclide_data.nuc(*nuclides[i])

        assert ufloat_equiv(d['weight'], weights[i])
        assert ufloat_equiv(d['abundance'], abundances[i])
        assert d['half-life'] == half_lifes[i]
        assert d['stable'] == stable[i]


def test_weight():
    """
    Does weight() function work as expected?
    """
    symbols = [ 'H', 'K', 'eu', 'CM' ]
    nuclides = [ (1,1), (19, 49), (63, 148), (96, 240) ]
    weights  = [ 1.00782503207, 48.967450, 147.918086, 240.0555295, ]

    for i in range(len(nuclides)):

        assert weights[i] == nuclide_data.weight(*nuclides[i])
        assert weights[i] == nuclide_data.weight(symbols[i], nuclides[i][1])
        assert weights[i] == nuclide_data.weight(
                              "{0}-{1}".format(symbols[i], nuclides[i][1]))


def test_zaids():
    """Does zaid conversion work correctly?"""

    zaids = [92235, 3006, "54135", "08016"]
    zas = [(92,235), (3,6), (54, 135), (8, 16)]
    for zaid, ref_za in zip(zaids, zas):
        za = nuclide_data.zaid2za(zaid)
        assert ref_za == za

def test_Nuclide_class_init():
    """Does Nuclide class correctly identify nuclide for a variety of inputs?"""

    class Foo:
       pass

    nuc_obj = Foo()
    nuc_obj.Z = 92
    nuc_obj.A = 235

    nuc_ids =   [ 'U235', 'U-235', '235U', '235-U',
                  'u235', 'u-235', '235u', '235-u',
                   92235, "92235",
                   (92,235), [92, 235],
                   {'Z':92, 'A':235},
                   nuc_obj
                ]

    for nuc_id in nuc_ids:
        nuclide = nuclide_data.Nuclide(nuc_id)

        assert nuclide.Z == 92
        assert nuclide.A == 235
        assert nuclide.element == 'U'

def test_Nuclide_class_init_2():
    """Does Nuclide class correctly identify nuclide for a variety of inputs?"""

    class Foo:
       pass

    nuc_obj = Foo()
    nuc_obj.Z = 3
    nuc_obj.A = 6

    nuc_ids =   [ 'Li6', 'LI-6', '6LI', '6-Li',
                  'li6', 'lI-6', '6li', '6-li',
                   3006, "03006",
                   (3,6), [3, 6],
                   {'Z': 3, 'A':6},
                   nuc_obj
                ]

    for nuc_id in nuc_ids:
        nuclide = nuclide_data.Nuclide(nuc_id)

        assert nuclide.Z == 3
        assert nuclide.A == 6
        assert nuclide.element == 'Li'

def test_Nuclide_class_init_with_E():
    """Does Nuclide class correctly get isomeric energy?"""
    E_ref = 0.5

    class Foo:
       pass

    nuc_obj = Foo()
    nuc_obj.Z = 3
    nuc_obj.A = 6
    nuc_obj.E = E_ref

    nuc_ids =   [ ('Li6', E_ref),  (3006, E_ref),
                   ((3,6, E_ref),), ([3, 6, E_ref],),
                   ({'Z': 3, 'A':6, 'E': E_ref},),
                   (nuc_obj,)
                ]

    for nuc_id in nuc_ids:
        nuclide = nuclide_data.Nuclide(*nuc_id)

        # Primary check
        assert np.allclose([nuclide.E,], [E_ref])

        # Secondary check
        assert nuclide.Z == 3
        assert nuclide.A == 6
        assert nuclide.element == 'Li'


