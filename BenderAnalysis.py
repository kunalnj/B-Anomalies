###############################################################################
#                                                                             #
# This file contains the boiler plate for running over data in Bender. Unless #
# you want to fiddle with the data options, you probably don't need to do     #
# anything in here.                                                           #
#                                                                             #
###############################################################################

#import everything from bender 
from Bender.Main import *

__author__  = "Michael McCann"
__date__    = "2020/10/18" 
__version__ = "0.01"

from Bender.Logger import getLogger 
if '__main__' == __name__ : logger = getLogger ( 'Something' )
else                      : logger = getLogger ( __name__ )


import os

#This is the work horse of the code, you probably want to be looking there
from Bs2Phi3TauTau import Bs2Phi3TauTau

# =============================================================================

# =============================================================================
## The configuration of the job
def configure ( inputdata        ,    ## the list of input files  
                catalogs = []    ,    ## xml-catalogs (filled by GRID)
                castor   = False ,    ## use the direct access to castor/EOS ? 
                params   = {}    ) :
    

    ## import DaVinci 
    from Configurables import DaVinci, MessageSvc
    froot=os.path.splitext(os.path.basename(inputdata[0]))[0].replace('Brunel', '')
    dv = DaVinci ( DataType   = '2018' ,
                   InputType  = 'DST'  ,
                   TupleFile = 'Bender{}.root'.format(froot),
                   HistogramFile = 'Bender{}-histos.root'.format(froot),
                   Simulation = True   )
    DaVinci().DDDBtag   = "dddb-20170721-3"
    DaVinci().CondDBtag = "sim-20190430-vc-mu100"

    ## define the input data
    setData(inputdata, catalogs, castor, useDBtags = False)
    from GaudiConf import IOHelper

    ## get/create application manager
    gaudi = appMgr() 
    
    #
    ## modify/update the configuration:
    #
    ## (1) create the algorithm
    alg = Bs2Phi3TauTau( 'Bs2phitautau' )
    
    ## (2) replace the list of top level algorithm by
    #     new list, which contains only *THIS* algorithm
    gaudi.setAlgorithms( [ alg ] )
             
    return SUCCESS 
# =============================================================================



# =============================================================================
## Job steering 
import sys
if __name__ == '__main__' :

  logger.info ( 80*'*'  ) 
  logger.info ( __doc__ ) 
  logger.info ( ' Author  : %s ' %  __author__  ) 
  logger.info ( ' Version : %s ' %  __version__ ) 
  logger.info ( ' Date    : %s ' %  __date__    ) 
  logger.info ( 80*'*'  ) 

  ## job configuration
  inputdata = [arg for arg in sys.argv if arg.endswith('dst')]
  if not inputdata:
    print 'You need to specify some input data, none found.'
    sys.exit(0)

  print 'Algorithm to run over', inputdata

  configure( inputdata , castor = False )
    
  ## event loop 
  nEvents = -1  #note this shuold be -1 to run over all events or the approprate number of events
  try:
    nEvents = int(sys.argv[-1])
  except:
    nEvents = -1
  run(nEvents)

# =============================================================================
# The END
# =============================================================================


