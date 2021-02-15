# B-Anomalies

These are adapted bender scripts from https://gitlab.cern.ch/lhcb-imperial-public/examples/bender with added calculations for distance of closest approach (DOCA), impact parameter (IP), and flight distance (FD) of the tau in the Bs2->KKtautau decay. 

## Changes

Changes were made to TupleUtilities.py and Bs2Phi3TauTau.py. TupleUtilities.py now also contains wrappers for the functions get_IP, get_FD, get_DOCA, and get_vertex. The first three utilize [IDistanceCalculator](https://gitlab.cern.ch/lhcb/Phys/blob/2dd7d41d845627ea82e5b5262014bce22c3a2597/Phys/DaVinciInterfaces/Kernel/IDistanceCalculator.h), the last one the [IVertexFit](https://gitlab.cern.ch/lhcb/Phys/blob/2dd7d41d845627ea82e5b5262014bce22c3a2597/Phys/DaVinciInterfaces/Kernel/IVertexFit.h).

The three defined functions are then used in Bs2Phi3TauTau.py to find the required quantities and simply appended to the resultant tuple.
