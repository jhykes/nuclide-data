#!/usr/bin/env python
"""
Access to isotope atomic weights and abundances.

Data from NIST.
http://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=&ascii=ascii2&isotype=all

"""

import gzip

import numpy as np

import uncertainties as unc

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
            d[k] = unc.ufloat(v)
        else:
            try:
                d[k] = int(v)
            except ValueError:
                d[k] = v

    return d

# NIST data file
data_file = "isotope_data.txt"

# chunk file into isotopes
nist_isotope_raw_list = []
current_isotope = []
for line in open(data_file):
    if line == '\n':
        nist_isotope_raw_list.append(current_isotope)
        current_isotope = []
    else:
        current_isotope.append(line.rstrip())


nist_isotope_processed_list = []
for raw_chunk in nist_isotope_raw_list:
    nist_isotope_processed_list.append( parse_one_chunk(raw_chunk) )


nist_per_element = {}
for isotope in nist_isotope_processed_list:
    Z = isotope['Atomic Number']
    try:
        nist_per_element[Z].append(isotope)
    except KeyError:
        nist_per_element[Z] = []
        nist_per_element[Z].append(isotope)


# Nuclear wallet cards data ------------------------------

def nndc_unc(string, delimiter):
    a, b = map(float, string.split(delimiter))
    return unc.ufloat((a,b))

def do_if_present(string, func):
    string = string.strip()
    if string:
        return func(string)
    else:
        return None
    

def parse_one_wallet_line(line):
    d = {}
    d['A'] = int(line[1:4])
    d['isomeric'] = (line[4] == 'M')
    d['Z'] = int(line[6:9])
    d['Atomic Symbol'] = line[10:12].title()
    d['J_pi'] = line[16:26]
    d['decay mode'] = do_if_present(line[30:34], str)

    d['branch fraction'] = do_if_present(line[34:41], lambda x: float(x) / 100.)

    d['excitation energy'] = do_if_present(line[42:49], lambda x: float(x))

    d['Q-value'] = float(line[49:56])                  # in MeV
    d['half-life string'] = line[63:80]        

    d['abundance'] = do_if_present(line[81:96], lambda x: nndc_unc(x,'%')/100.)
   
    d['atomic mass'] = unc.ufloat((line[97:105], line[105:113])) # in MeV
    d['systematics mass'] = (line[114] == 'S')
    d['half-life'] = float(line[124:133])   # in seconds


wallet_filename = 'nuclear-wallet-cards.txt.gz'
wallet_file = gzip.open(wallet_filename, 'rb')
wallet_content = wallet_file.read()
wallet_file.close()

wallet_lines = wallet_content.split('\n')[:-1]

wallet_isotope_processed_list = []
for line in wallet_lines:
    wallet_isotope_processed_list.append( parse_one_wallet_line(line) )

