# B-Anomalies

These are adapted bender scripts from https://gitlab.cern.ch/lhcb-imperial-public/examples/bender with added calculations for distance of closest approach (DOCA), impact parameter (IP), and flight distance (FD) of the tau in the Bs2->KKtautau decay. 

## Changes

Changes were made to TupleUtilities.py and Bs2Phi3TauTau.py. TupleUtilities.py now also contains wrappers for the functions get_IP, get_FD, get_DOCA, and get_vertex. The first three utilize the [IDistanceCalculator](https://gitlab.cern.ch/lhcb/Phys/blob/2dd7d41d845627ea82e5b5262014bce22c3a2597/Phys/DaVinciInterfaces/Kernel/IDistanceCalculator.h) tool, the last one the [IVertexFit](https://gitlab.cern.ch/lhcb/Phys/blob/2dd7d41d845627ea82e5b5262014bce22c3a2597/Phys/DaVinciInterfaces/Kernel/IVertexFit.h) tool.

The three defined functions are then used in Bs2Phi3TauTau.py to find the required quantities and simply appended to the resultant tuple.

The functions get_IP and get_DOCA are directly acted on the particles/tracks/vertices we want distances from. For the FD, a vertex between the muon and phi3 is first found using get_vertex and then the tau FD is taken to be the distance between this newly found vertex and phi3.endVertex(). 

All changes are noted in the respective .py files.

## Running

lb-run -c best Bender/v35r6 bash
python BenderAnalysis.py /vols/lhcb/masmith/gangadir_bs2kktautau/workspace/mesmith/LocalXML/59/0/output/Brunel.dst

## Output

The output tuple is the same name as the .dst Brunel input file, in our case this was tested on /vols/lhcb/masmith/gangadir_bs2kktautau/workspace/mesmith/LocalXML/59/0/output/Brunel.dst. An example output .root file is the Bender.root file.

