isotopic-data -- a Python interface to isotopic data
----------------------------------------------------

This small module reads NIST and NNDC data files and makes them accessible
via Python functions, lists, and dictionaries.

The data for each nuclide is contained in a Python dictionary with
the following keys:

  * 'A' : atomic mass number (int)
  * 'Z' : atomic number or number of protons (int)
  * 'symbol' : element symbol for this Z (string)
  * 'isomeric' : is this nuclide in an excited state (Boolean)
  * 'Jpi' : spin and parity (string)
  * 'decay mode' : the type of radioactive decay, if applicable (string)
  * 'branch fraction' : the branching fraction for a particular decay
                         mode (float in (0, 1])
  * 'excitation energy' : excitation energy of isomer in MeV (float)
  * 'Q' : Q-value of decay mode in MeV (float)
  * 'half-life' : half life in seconds (float)
  * 'lambda' : decay constant in 1/seconds (float)
  * 'abundance' : isotopic abundance (float in [0, 1])
  * 'weight' : atomic weight, dimensionless (float)
  * 'weight mev' : atomic weight in MeV/c^2 (float)


Contents
--------

 * isotope_data.py -- Python script to read data files and create interface 
                      to data.
 * isotope_data.txt -- NIST file with atomic weights and abundances
 * nuclear-wallet-cards.txt.gz -- Nuclear Wallet Card ASCII file
 * WC-format -- explanation of Nuclear Wallet Card ASCII format

References
----------

 1. Jagdish K. Tuli, **Nuclear Wallet Cards**,
    National Nuclear Data Center, April 2005. http://www.nndc.bnl.gov/wallet
 2. J. S. Coursey, D. J. Schwab, J. J. Tsai, and R. A. Dragoset,
    **Atomic Weights and Isotopic Compositions with Relative Atomic
    Masses**, NIST Physical Measurement Laboratory,
    updated July 2010. http://www.nist.gov/pml/data/comp.cfm
