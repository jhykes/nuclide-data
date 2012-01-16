nuclide-data -- a Python interface to nuclide data
--------------------------------------------------

This script reads nuclide abundances, weights, and decay constants from
NIST and NNDC data files and makes them accessible
via Python functions, lists, and dictionaries.


The recommended ways to access the data are:

 * ``list_of_As = isotopes[Z]`` dictionary whose keys are Z, and values are list of isotopes' A
 * ``nuclide_dict = nuc(Z, A, E=0.)`` function to return nuclide data for Z & A, 
   and (optionally) E, determining isomeric state.
 * ``list_of_isomer_energies = isomers(Z, A)`` Return energy levels (in MeV) of 
   isomeric states for particular Z & A.
 * ``w = weights(Z_or_symbol, A=None, E=0.)`` Return atomic or isotopic weight
   for an element or nuclide. The input to this function is flexible:

      * ``element_weight = weights('U')`` element weight of U
      * ``element_weight = weights(92)`` element weight of U
      * ``nuclide_weight = weights('U-235')`` isotopic weight of U-235
      * ``nuclide_weight = weights('U', 235)`` isotopic weight of U-235
      * ``nuclide_weight = weights(92, 235)`` isotopic weight of U-235

These commands will return errors if the data is unavailable.


The data for each nuclide is contained in a Python dictionary with
the following keys:

  * 'symbol' : element symbol for this Z (string)
  * 'isomeric' : is this nuclide in an excited state (Boolean)
  * 'Jpi' : spin and parity (string)
  * 'half-life' : half life in seconds (float)
  * 'stable' : is nuclide stable (Boolean)
  * 'lambda' : decay constant in 1/seconds (float)
  * 'abundance' : isotopic abundance (float in [0, 1])
  * 'weight' : atomic weight, dimensionless (ufloat)
  * 'mass excess' : mass excess in MeV (ufloat)
  * 'decay modes' : a dictionary with the type of decay as the keys. Each
    dictionary contains
      * 'Q-value' : Q-value of decay reaction in MeV (float)
      * 'branch fraction' : the branching fraction for a particular decay
        mode (float in (0, 1])


Dependencies
------------
 * numpy
 * uncertainties

Contents
--------

 * nuclide_data.py -- Python script to read data files and create interface 
   to data.
 * test_nuclide_data.py -- unit tests to verify implementation of 
   nuclide_data.py. These tests can be useful as examples.
 * nist-nuclide-data.txt -- NIST file with atomic weights and abundances
 * nuclear-wallet-cards.txt.gz -- Nuclear Wallet Card ASCII file
 * WC-format.pdf -- explanation of Nuclear Wallet Card ASCII format

References
----------

 1. Jagdish K. Tuli, **Nuclear Wallet Cards**,
    National Nuclear Data Center, April 2005. http://www.nndc.bnl.gov/wallet
 2. J. S. Coursey, D. J. Schwab, J. J. Tsai, and R. A. Dragoset,
    **Atomic Weights and Isotopic Compositions with Relative Atomic
    Masses**, NIST Physical Measurement Laboratory,
    updated July 2010. http://www.nist.gov/pml/data/comp.cfm
