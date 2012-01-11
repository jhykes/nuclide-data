#!/usr/bin/env python
"""
Access to isotope atomic weights and abundances.

Data from NIST.
http://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=&ascii=ascii2&isotype=all

"""

import numpy as np
import uncertainties as unc

# NIST data file
data_file = "isotope_data.txt"

# chunk file into isotopes
isotope_raw_list = []
current_isotope = []
for line in open(data_file):
    if line == '\n':
        isotope_raw_list.append(current_isotope)
        current_isotope = []
    else:
        current_isotope.append(line.rstrip())

def split_line(line):
    return map(str.strip, line.split('='))

def parse_one_chunk(chunk):
    d = {}
    for line in chunk:
        k,v = split_line(line)
        if v.find('.') >= 0:
            if v.endswith('#'):
                v = v[:-1]
            d[k] = unc.ufloat(v)
        else:
            try:
                d[k] = int(v)
            except ValueError:
                d[k] = v

    return d

isotope_processed_list = []
for raw_chunk in isotope_raw_list:
    isotope_processed_list.append( parse_one_chunk(raw_chunk) )


per_element = {}
for isotope in isotope_processed_list:
    Z = isotope['Atomic Number']
    try:
        per_element[Z].append(isotope)
    except KeyError:
        per_element[Z] = []
        per_element[Z].append(isotope)

