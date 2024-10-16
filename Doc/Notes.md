
Uses Hollerith constants and some work would be required 
to replace them with character data.

Uses integer arrays as FORMAT statements. These are passed
through the code. Should fix at the same time as 1. above.

The package requires a machine dependent set of routines. The ones
provided in the Drivers directory in the file epc_mpupk.f are specific for
the EPC Fortran 90 compiler and need to be changed for other compilers.
If the name is changed from epc_mpupk.f the makefile will also need
changing.
