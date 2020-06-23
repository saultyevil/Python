# CNO Processed Abundances

This directory contains two atomic data masterfiles which assume CNO processed
abundances instead of solar abundances. For CNO processed abundances, the abundance
of He, C, N and O are modified.

In the following data files, the abundances represent a 2 solar mass star at the
end of its main sequence. The following abundances are used:

* He ~ 2 x Solar
* C ~ 0.5 x Solar
* N ~ 7 x Solar
* O ~ 1 x Solar

These values were obtained from Gallegos-Garcia et al. 2018 (https://ui.adsabs.harvard.edu/abs/2018ApJ...857..109G/abstract).

## Masterfiles

* h20_hetop_standard80_cno_abun.dat : macro atom data with 20 level H and He atoms 
* standard80_cno_abun.dat : the standard data files from the py81c release

## Parameter Files

Included in this directory is an example parameter file, `tde_cno.pf`, for a
model of a TDE accretion disc wind using the macro atom CNO processed abundances 
atomic data. It is recommended that you run this example with the "-p 2" option.
