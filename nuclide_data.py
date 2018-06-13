#!/usr/bin/env python
"""
Access to nuclide atomic weights, abundances, and decay constants.

Data from NIST and NNDC
 * http://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=&ascii=ascii2&isotype=all
 * http://www.nndc.bnl.gov/wallet/

"""

import os.path
import warnings
import re
import gzip
import copy
import string
from functools import total_ordering

import numpy as np

import uncertainties as unc

# Absolute path to this module
basepath = os.path.dirname(__file__)

# conversion factor from MeV/c^2 to amu
mev_per_c_2_amu = 1. / 931.494061




    

# NIST data -------------------------------------------------------------
def split_line(line):
    return map(str.strip, line.split('='))

def parse_one_chunk(chunk):
    d = {}
    for line in chunk:
        k,v = split_line(line)
        if v.find('.') >= 0:
            if v.endswith('#'):
                v = v[:-1]
            d[k] = unc.ufloat_fromstr(v)
        else:
            try:
                d[k] = int(v)
            except ValueError:
                d[k] = v

    return d

# NIST data file
data_file = os.path.join(basepath, "nist-nuclide-data.txt")

# chunk file into nuclides
nist_nuclide_raw_list = []
current_nuclide = []
for line in open(data_file):
    if line == '\n':
        nist_nuclide_raw_list.append(current_nuclide)
        current_nuclide = []
    else:
        current_nuclide.append(line.rstrip())


nist_nuclide_processed_list = []
for raw_chunk in nist_nuclide_raw_list:
    nist_nuclide_processed_list.append( parse_one_chunk(raw_chunk) )

nist_per_element = {}
for nuclide in nist_nuclide_processed_list:
    Z = nuclide['Atomic Number']
    try:
        nist_per_element[Z].append(nuclide)
    except KeyError:
        nist_per_element[Z] = []
        nist_per_element[Z].append(nuclide)

nist_nuclides = {}
for nuclide in nist_nuclide_processed_list:
    Z = nuclide['Atomic Number']
    A = nuclide['Mass Number']

    nist_nuclides[(Z,A)] = nuclide

z2sym = dict(
          [ (Z, nist_per_element[Z][0]['Atomic Symbol']) for Z in range(1,119) ]
            )

sym2z = dict( [ (z2sym[k], k) for k in z2sym ] )

atomic_weights = {}
for Z in nist_per_element:
    w = nist_per_element[Z][0]['Standard Atomic Weight']
    if type(w) is not str:
        atomic_weights[Z] = w

# Nuclear wallet cards data ------------------------------

def nndc_unc(string, delimiter):
    a, b = map(float, string.split(delimiter))
    return unc.ufloat(a,b)

def nndc_abun(string, delimiter):
    a, b = string.split(delimiter)
    unc_string = "{0}({1})".format(float(a), int(b))
    return unc.ufloat_fromstr(unc_string)

def do_if_present(string, func, default=None):
    string = string.strip()
    if string:
        return func(string)
    else:
        return default
    
def process_branch(s):
    try:
        return float(s) / 100.
    except ValueError:
        return None

def process_abundance(s):
    if s.startswith('100'):
        return 1.
    else:
        return nndc_abun(s,'%') / 100.

def parse_one_wallet_line(line):
    d = {}
    d['A'] = int(line[1:4])
    d['Z'] = int(line[6:9])
    d['symbol'] = line[10:12].strip().title()

    d['mass excess'] = unc.ufloat(*map(float, (line[97:105], line[105:113]))) # in MeV
    d['systematics mass'] = (line[114] == 'S')
    d['abundance'] = do_if_present(line[81:96], process_abundance, default=0.)

    d['Jpi'] = line[16:26].strip()

    d['isomeric'] = (line[4] == 'M')
    d['excitation energy'] = do_if_present(line[42:49], lambda x: float(x), default=0.)

    d['decay mode'] = do_if_present(line[30:34], str)
    d['branch fraction'] = do_if_present(line[35:41], process_branch)
    d['Q-value'] = do_if_present(line[49:56], float)   # in MeV

    d['half-life string'] = line[63:80].strip()
    d['stable'] = (d['half-life string'] == 'STABLE')
    d['half-life'] = float(line[124:133])   # in seconds

    if d['stable']: d['half-life'] = np.inf

    return d


wallet_filename = os.path.join(basepath, 'nuclear-wallet-cards.txt.gz')
wallet_file = gzip.open(wallet_filename, 'rb')
wallet_content = wallet_file.read()
wallet_file.close()

wallet_lines = wallet_content.split('\n')[:-1]

wallet_nuclide_processed_list = []
for line in wallet_lines:
    wallet_nuclide_processed_list.append( parse_one_wallet_line(line) )
    d = parse_one_wallet_line(line)


isomer_keys = ['symbol', 'mass excess', 'abundance', 'isomeric',
               'Jpi', 'stable', 'half-life string', 'half-life']

decay_keys = ['branch fraction', 'Q-value']

# -------------------------------------------------------------------------
# Build master dictionary
nuclides = {}
for el in wallet_nuclide_processed_list:
    Z, A, E = [el[i] for i in ['Z', 'A', 'excitation energy']]

    # Pick the nuclide (Z,A) (or create new entry)
    if not ((Z,A) in nuclides):
        nuclides[(Z,A)] = {}

    isomers = nuclides[(Z,A)]

    # Pick the isomer [(Z,A)][E] (or create new entry)
    if not (E in isomers):
        isomers[E] = {}
        isomers[E]['decay modes'] = {}

        isomer = isomers[E]

        # nuclide data not associated with decay
        for k in isomer_keys:
            isomer[k] = el[k]

        if isomer['half-life'] == 0.:
            isomer['lambda'] = np.inf
        else:
            isomer['lambda'] = np.log(2.) / isomer['half-life']

        if (Z,A) in nist_nuclides:
            isomer['weight'] = nist_nuclides[(Z,A)]['Relative Atomic Mass']

    else:
        isomer = isomers[E]

    # decay data
    isomer['decay modes'][el['decay mode']] = {}
    for k in decay_keys:
        isomer['decay modes'][el['decay mode']][k] = el[k]


default_isomer_E = {}
meta_suffixes = 'mnopqrs'
for n in nuclides:
    Es = nuclides[n].keys()

    if n[0] == 0: continue

    nuc_string = '{}-{}'.format(z2sym[n[0]], n[1])

    if len(Es) > 1:
        for i,E in enumerate(Es[1:]):
            default_isomer_E[nuc_string+meta_suffixes[i]] = E


def return_nominal_value(Z_or_symbol, A, E, attribute):
    """
    Input
    -----
     * Z_or_symbol : can be either 'U', 92, or 'U-235'. If it is 'U-235',
       then the A argument is ignored.
     * A : atomic mass number
     * E : excitation energy of isomer
     * attribute : a valid nuclide data dictionary key
    """
    try:
        # Is Z_or_symbol a fully specified nuclide, e.g., 'U-235'
        if Z_or_symbol.find('-') > -1:
            symbol, A = Z_or_symbol.split('-')
            A = int(A)
            Z = sym2z[symbol.title()]
        else:
            Z = sym2z[Z_or_symbol.title()]
    except AttributeError:
        Z = Z_or_symbol
   
    # testing for no A, then return elemental value
    if A is None:
        return atomic_weights[Z].nominal_value
        
    try:
        return nuclides[(Z,A)][E][attribute].nominal_value
    except (ValueError, AttributeError):
        return nuclides[(Z,A)][E][attribute]


# ENDF-6 MAT data -------------------------------------------------------
#  mats is dictionary with
#    key : (Z, A, metastable), Z, A are int, metastable is bool
#    value : MAT nuclide id, integer, from ENDF-6 list
mat_file = os.path.join(basepath, "n-ENDF-B-VII.1.endf.list")
mats = {}
for line in open(mat_file):
    # Skip comment line
    if line.startswith('#'): continue

    # Grab Z, A, and MAT
    Z = int(line[6:9])
    A = int(line[13:16])
    mat = int(line[72:76])

    # Is it metastable?
    metastable = (line[16] == 'M')

    key = (Z, A, metastable)

    mats[key] = int(mat)



# ---------------------------------------------------------------------------- #
# means intended for public access of data

# list_of_As = isotopes[Z]
isotopes = {}
for (Z,A) in nuclides:

    if not (Z in isotopes):
        isotopes[Z] = []
    
    isotopes[Z].append(copy.copy(A))

    isotopes[Z].sort()


def zaid2za(zaid):
    """
    Convert ZZAAA to (Z,A) tuple.

    """
    # Ignores decimal and stuff after decimal.
    zaid = str(int(zaid))

    Z = int(zaid[:-3])
    A = int(zaid[-3:])
    return (Z, A)



def nuc(Z, A, E=0.):
    """
    Return nuclide data for Z, A, and (optionally) E of isomeric state.
    """
    return nuclides[(Z,A)][E]


def isomers(Z, A):
    """
    Return energy levels of isomeric states for particular Z & A.

    Energies in MeV.
    """
    isom = nuclides[(Z,A)].keys()
    isom.sort()
    return isom


def weight(Z_or_symbol, A=None, E=0.):
    """
    Return atomic weight for Z, A, and (optionally) E of isomeric state.
    """
    return return_nominal_value(Z_or_symbol, A, E, 'weight')


@total_ordering
class Nuclide:
    """
    Provide a convenient interface to various nuclide names and data.

    Input nuc_id can be:
       * Alphanumeric: 'U235', 'U-235', '235U', '235-U'
           -- letters may be lower or uppercase
       * ZAID: 92235, "92235"
       * Tuple/list: (92, 235), [92, 235] 
       * Tuple/list with energy: (92, 235, 0.5), [92, 235, 0.5] 
       * Dictionary: {'Z':92, 'A':235}
       * Object x with x.Z and x.A integer attributes
       * Metastable, only as "Am242m" or "AM-242M"
            that is "Sym[-]AAAm", case insensitive

    Energy level E in MeV.

    For a select set of nuclides (see default_isomer_E.keys())
      a default value for E is set if no E is supplied.

    If metastable and no E provided, it is set to np.inf.

    If E is provides as part of nuc_id, then it overrides the value
    provided to the keyword.
    """


    def __init__(self, nuc_id, E=0., metastable=False):

        try:
            # Object with attributes
            self.Z, self.A = nuc_id.Z, nuc_id.A
            try:
                E = nuc_id.E
            except AttributeError:
                E = 0.
        except AttributeError:

            try:
                # Dictionary
                self.Z, self.A = nuc_id['Z'], nuc_id['A']
                try:
                    E = nuc_id['E']
                except (KeyError, TypeError):
                    E = 0.
            except (KeyError, TypeError):

                # Integer ZAID
                if type(nuc_id) is int:
                    self.Z, self.A = zaid2za(nuc_id)

                # List or tuple
                if type(nuc_id) in [list, tuple]:
                    if len(nuc_id) == 2:
                        self.Z, self.A = nuc_id
                    if len(nuc_id) == 3:
                        self.Z, self.A, E = nuc_id

                # String
                if type(nuc_id) is str:

                    # Alphanumeric
                    if re.search('[a-zA-Z]', nuc_id):

                        # upper case to make comparison easier
                        nuc_id = nuc_id.upper()

                        # Metastable
                        if ( nuc_id[0] in string.ascii_letters 
                                          and nuc_id[-1] == 'M'):
                            metastable = True
                            nuc_id = nuc_id[:-1]


                        # Hyphenated
                        if re.search('-', nuc_id):
                            s1, s2 = nuc_id.split('-')
                            s1 = s1.strip()
                            s2 = s2.strip()
                        else:
                            s1 = filter(lambda x: x in string.ascii_letters, nuc_id)
                            s2 = filter(lambda x: not (x in string.ascii_letters), nuc_id).strip()
                            

                        # Not sure of the order of s1 & s2,
                        #  so try one, then the other.
                        try:
                            self.Z = sym2z[s1.title()]
                            self.A = int(s2)
                        except:
                            self.Z = sym2z[s2.title()]
                            self.A = int(s1)
                                
                        
                    else: # assume it is a ZAID string
                        self.Z, self.A = zaid2za(nuc_id)
        
        # Metastable can be specified by either E or metastable flag.
        #  If flag is given but E is not, then set E to inf as
        #  an indication that it is not stable, but that the exact
        #  E value isn't given.
        if metastable and E==0.:
            E = np.inf

        self.metastable = E > 0.

        # Assign E, unless it has already been set
        self.E = E

        self.element = z2sym[self.Z]

        # Assign E for list of metastable nuclides if E wasn't provided
        if (self.E is np.inf and 
               self.__repr__() in default_isomer_E.keys()):
            self.E = default_isomer_E[self.__repr__()]
            
            
        try:
            self.weight = return_nominal_value(self.Z, self.A, self.E, 'weight')
        except:
            warnings.warn("nuclide weight not available for {}".format(self))

        # Set MAT for ENDF6
        try:
            self.mat = mats[(self.Z, self.A, self.metastable)]
        except:
            warnings.warn("nuclide {} not on ENDFB-VII.1 neutron library".format(self))

    def zaid(self):
        return self.Z*1000 + self.A

    def decay_const(self):
        return return_nominal_value(self.Z, self.A, self.E, 'lambda')

    def half_life(self):
        return return_nominal_value(self.Z, self.A, self.E, 'half-life')

    def __repr__(self):
        if self.E==0.:
            return "{x.element}-{x.A}".format(x=self)
        else:
            return "{x.element}-{x.A}m".format(x=self)

    def __key__(self):
        return (self.Z, self.A, self.E)


    def __hash__(self):
        return hash(self.__key__())


    def __eq__(self, other):
        E_test = np.allclose([self.E,], [other.E])
        return ( (self.Z, self.A) == (other.Z, other.A) ) and E_test


    def __lt__(self, other):
        return ( (self.Z, self.A, self.E) < (other.Z, other.A, other.E) )


