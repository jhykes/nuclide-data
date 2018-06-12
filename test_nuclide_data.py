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
             and np.allclose([a.std_dev], [b.std_dev]) )
    except AttributeError:
        return a == b 
        #return np.allclose([a], [b]) 

def test_data():
    """
    Do we have the correct weight, abundance, and half-life for select nuclides?
    """
    nuclides = [ (1,1), (19, 49), (63, 148), (96, 240) ]
    weights  = [ unc.ufloat_fromstr("1.00782503207(10)"),
                 unc.ufloat_fromstr("48.967450(80)"),
                 unc.ufloat_fromstr("147.918086(11)"),
                 unc.ufloat_fromstr("240.0555295(25)"), ]
    abundances  = [ unc.ufloat_fromstr("0.999885(70)"), 0., 0., 0. ]
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
        assert nuclide.metastable == False
        assert nuclide.E == 0.

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
        assert nuclide.metastable == False
        assert nuclide.E == 0.


def test_Nuclide_class_init_with_E():
    """Does Nuclide class correctly get isomeric energy?"""
    E_ref = 0.2283

    class Foo:
       pass

    nuc_obj = Foo()
    nuc_obj.Z = 13
    nuc_obj.A = 26
    nuc_obj.E = E_ref

    nuc_ids =   [ ('Al26', E_ref),  (13026, E_ref),
                   ((13, 26, E_ref),), ([13, 26, E_ref],),
                   ({'Z': 13, 'A': 26, 'E': E_ref},),
                   (nuc_obj,)
                ]

    for nuc_id in nuc_ids:
        nuclide = nuclide_data.Nuclide(*nuc_id)

        # Primary check
        assert np.allclose([nuclide.E,], [E_ref])

        # Secondary check
        assert nuclide.Z == 13
        assert nuclide.A == 26
        assert nuclide.element == 'Al'
        assert nuclide.metastable == True


def test_Nuclide_class_init_with_metastable():
    """Does Nuclide class correctly understand metastable notation?"""
    E_ref = 0.5

    nuclide = nuclide_data.Nuclide('Li6m', E_ref)

    # Primary check
    assert np.allclose([nuclide.E,], [E_ref])
    assert nuclide.Z == 3
    assert nuclide.A == 6
    assert nuclide.element == 'Li'
    assert nuclide.metastable == True


    nuclide = nuclide_data.Nuclide('LI-6M')

    # Primary check
    assert nuclide.E is np.inf
    assert nuclide.Z == 3
    assert nuclide.A == 6
    assert nuclide.element == 'Li'
    assert nuclide.metastable == True

    nuclide = nuclide_data.Nuclide((3,6), metastable=True)

    # Primary check
    assert nuclide.E is np.inf
    assert nuclide.Z == 3
    assert nuclide.A == 6
    assert nuclide.element == 'Li'
    assert nuclide.metastable == True



def test_Nuclide_class_MAT():
    """Does Nuclide class correctly set MAT?"""

    nuc_ids = {
        'H -  1 ' : 125,
        'He-  4 ' : 228,
        'N - 15 ' : 728,
        'O - 16 ' : 825,
        'Si- 29 ' : 1428,
        'Ca- 44 ' : 2037,
        'Sc- 45 ' : 2125,
        'Fe- 54 ' : 2625,
        'Co- 58 ' : 2722,
        'Co- 58M' : 2723,
        'Se- 82 ' : 3449,
        'Sr- 87 ' : 3834,
        'Ag-109 ' : 4731,
        'Ag-110M' : 4735,
        'Cd-115M' : 4853,
        'Sn-112 ' : 5025,
        'Gd-152 ' : 6425,
        'Ra-226 ' : 8834,
        'U -235 ' : 9228,
        'Am-242 ' : 9546,
        'Am-242M' : 9547,
        'Am-243 ' : 9549,
        'Am-244 ' : 9552,
        'Am-244M' : 9553,
        'Cm-240 ' : 9625,
        'Cm-241 ' : 9628,
        'Bk-246 ' : 9743,
        'Cf-254 ' : 9867,
        'Es-253 ' : 9913,
        'Es-254 ' : 9914,
        'Es-254M' : 9915,
        'Es-255 ' : 9916,
        'Fm-255 ' : 9936,
       }

    for nuc_id in nuc_ids:
        assert nuclide_data.Nuclide(nuc_id).mat == nuc_ids[nuc_id]

