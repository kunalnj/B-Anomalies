###############################################################################
#                                                                             #
# This file contains a number of utility functions used to help with writing  #
# things to the tuples.                                                       #
#                                                                             #
###############################################################################


from Bender.Main import *

import array
import numpy as np
from GaudiPython.Bindings import gbl


#value to write to tuple on error
BadFloat = -9999.0


#List of trigger lines to include in the tuple
TriggerLines = [
'L0HadronDecision',
'L0MuonDecision',
'L0DiMuonDecision',
'L0ElectronDecision',
'L0PhotonDecision',
'L0LocalPi0Decision',
'L0GlobalPi0Decision',
'L0MuonNoGlobDecision',
'L0B1gasDecision',
'L0B2gasDecision',
'L0CALODecision',
'L0DiMuon,lowMultDecision',
'L0HCALDecision',
'L0MUON,minbiasDecision',
'L0Muon,lowMultDecision',
'L0MuonHighDecision',
'L0PU20Decision',
'L0PUDecision',
'L0SPD40Decision',
'L0SPDDecision',

'Hlt1TrackAllL0Decision',
'Hlt1TrackMuonDecision',
'Hlt1TrackMVADecision',
#'Hlt1GlobalDecision',

'Hlt2Topo2BodyBBDTDecision',
'Hlt2Topo3BodyBBDTDecision',
'Hlt2Topo4BodyBBDTDecision',
'Hlt2TopoMu2BodyBBDTDecision',
'Hlt2TopoMu3BodyBBDTDecision',
'Hlt2TopoMu4BodyBBDTDecision',
'Hlt2SingleMuonDecision',
'Hlt2DiMuonDecision',
'Hlt2DiMuonDetachedDecision',
#'Hlt2GlobalDecision'

]

def getIp(particle, pv, disttool):
  '''
  Compute the IP & chi2
  '''
  ip = array.array('d', [0])
  chi2 = array.array('d', [0])
  disttool.distance(particle, pv, ip,chi2, False)
  return ip[0], chi2[0]

def getFd(endv, pv, disttool):
  '''
  Compute the FD & chi2
  '''
  fd = array.array('d', [0])
  chi2 = array.array('d', [0])
  disttool.distance(endv, pv, fd, chi2)
  return fd[0], chi2[0]

def getdoca(disttool, track1, track2):
  '''
  Compute DOCA & chi2
  '''
  doca = array.array("d", [0])
  chi2 = array.array("d", [0])
  disttool.distance(track1, track2, doca, chi2)
  return doca[0], chi2[0]

def getvertex(vfittool, particle1, particle2):
  '''
  Compute best vertex between two particle LOFs
  '''
  vertex = gbl.LHCb.Vertex()
  chi2 = array.array("d", [0])
  fit = vfittool.fit(vertex, particle1, particle2)
  return vertex, vertex.chi2(), fit


#Make a dictionary containing all the momentum variables
def getMomentumVariables(momentum=None,**options):
  if not momentum:
    return ['M', 'PX', 'PY', 'PZ', 'E', 'P', 'PT', 'ETA', 'PHI']
  vars = {}
  vars['M']   = momentum.M()
  vars['PX']  = momentum.Px()
  vars['PY']  = momentum.Py()
  vars['PZ']  = momentum.Pz()
  vars['E']   = momentum.E()
  vars['P']   = momentum.P()
  vars['PT']  = momentum.Pt()
  vars['ETA'] = momentum.Eta()
  vars['PHI'] = momentum.Phi()
  return vars

#Make a dictionary containing all the trigger variables
def getTriggerVariables(particle=None,**options):
  if not particle:
    names = []
    for triggerline in TriggerLines + ['L0GlobalDecision', 'Hlt1GlobalDecision', 'Hlt2GlobalDecision']:
      names.append(triggerline + '_Dec')
      names.append(triggerline + '_TOS')
      names.append(triggerline + '_TIS')
    return names
  vars = {}
  if particle.proto():
    particle = particle.proto()
  if 'triggertool' in options:
    options['triggertool'].setOfflineInput(particle)
    def fillLine(triggerline, outname = None):
      if not outname: outname = triggerline
      options['triggertool'].setTriggerInput(triggerline)
      trigger = options['triggertool'].triggerTisTos()
      if trigger.decision():
        print outname, ' fired', particle
      vars[outname + '_Dec'] = trigger.decision()
      vars[outname + '_TIS'] = trigger.tis()
      vars[outname + '_TOS'] = trigger.tos()


    for triggerline in TriggerLines:
      fillLine(triggerline)
    fillLine('L0.*Decision',   'L0GlobalDecision')
    fillLine('Hlt1.*Decision', 'Hlt1GlobalDecision')
    fillLine('Hlt1.*Decision', 'Hlt2GlobalDecision')
  return vars

#Make a dictionary containing all the track variables
def getTrackVariables(track = None,**options):
  if not track:
    return ['X', 'Y', 'Z', 'PX', 'PY', 'PZ', 'TX', 'TY', 'P', 'PT', 'ETA', 'PHI', 'CHARGE', 'nStates', 'nMeasurements', 'CHI2', 'DOF', 'TYPE']
  vars = {}
  vars['X'] = track.position().x()
  vars['Y'] = track.position().y()
  vars['Z'] = track.position().z()
  vars['PX'] = track.momentum().x()
  vars['PY'] = track.momentum().y()
  vars['PZ'] = track.momentum().z()
  vars['TX'] = track.slopes().x()
  vars['TX'] = track.slopes().y()
  vars['P'] = track.p()
  vars['PT'] = track.pt()
  vars['ETA'] = track.pseudoRapidity()
  vars['PHI'] = track.phi()
  vars['CHARGE'] = track.charge()
  vars['nStates'] = track.nStates()
  vars['nMeasurements'] = track.nMeasurements()
  vars['CHI2'] = track.chi2()
  vars['DOF'] = track.nDoF()
  vars['TYPE'] = track.type()
  return vars

#Make a dictionary containing all the vertex variables
def getVertexVariables(vertex=None,**options):
  if not vertex:
    return ['VX', 'VY', 'VZ', 'VCHI2', 'VCHI2NDOF', 'VNDOF']
  vars = {}
  try:
    vars['VX']        = vertex.position().X()
    vars['VY']        = vertex.position().Y()
    vars['VZ']        = vertex.position().Z()
    vars['VCHI2']     = vertex.chi2()
    vars['VNDOF']     = vertex.nDoF()
    vars['VCHI2NDOF'] = vertex.chi2PerDoF()
  except:
    vars['VX']        = BadFloat
    vars['VY']        = BadFloat
    vars['VZ']        = BadFloat
    vars['VCHI2']     = BadFloat
    vars['VNDOF']     = BadFloat
    vars['VCHI2NDOF'] = BadFloat
  return vars

#Make a dictionary containing all the particle ID variables
def getPIDVariables(proto=None,**options):
  if not proto:
    return ['RichDLLe','RichDLLmu','RichDLLpi', 'RichDLLk', 'RichDLLp', 'RichDLLd', 'RichDLLbt', 'ProbNNe', 'ProbNNmu', 'ProbNNpi', 'ProbNNk', 'ProbNNp', 'ProbNNd', 'ProbNNghost', 'IsMuon', 'IsMuonLoose', 'IsMuonTight']
  vars = {}
  try:
    vars['RichDLLe']    = proto.info(proto.RichDLLe,    BadFloat)
    vars['RichDLLmu']   = proto.info(proto.RichDLLmu,   BadFloat)
    vars['RichDLLpi']   = proto.info(proto.RichDLLpi,   BadFloat)
    vars['RichDLLk']    = proto.info(proto.RichDLLk,    BadFloat)
    vars['RichDLLp']    = proto.info(proto.RichDLLp,    BadFloat)
    vars['RichDLLd']    = proto.info(proto.RichDLLd,    BadFloat)
    vars['RichDLLbt']   = proto.info(proto.RichDLLbt,   BadFloat)
    vars['ProbNNe']     = proto.info(proto.ProbNNe,     BadFloat)
    vars['ProbNNmu']    = proto.info(proto.ProbNNmu,    BadFloat)
    vars['ProbNNpi']    = proto.info(proto.ProbNNpi,    BadFloat)
    vars['ProbNNk']     = proto.info(proto.ProbNNk,     BadFloat)
    vars['ProbNNp']     = proto.info(proto.ProbNNp,     BadFloat)
    vars['ProbNNd']     = proto.info(proto.ProbNNd,     BadFloat)
    vars['ProbNNghost'] = proto.info(proto.ProbNNghost, BadFloat)
  except:
    vars['RichDLLe']    = BadFloat
    vars['RichDLLmu']   = BadFloat
    vars['RichDLLpi']   = BadFloat
    vars['RichDLLk']    = BadFloat
    vars['RichDLLp']    = BadFloat
    vars['RichDLLd']    = BadFloat
    vars['RichDLLbt']   = BadFloat
    vars['ProbNNe']     = BadFloat
    vars['ProbNNmu']    = BadFloat
    vars['ProbNNpi']    = BadFloat
    vars['ProbNNk']     = BadFloat
    vars['ProbNNp']     = BadFloat
    vars['ProbNNd']     = BadFloat
    vars['ProbNNghost'] = BadFloat

  try:
    vars['IsMuon']      = proto.muonPID().IsMuon()
    vars['IsMuonLoose'] = proto.muonPID().IsMuonLoose()
    vars['IsMuonTight'] = proto.muonPID().IsMuonTight()
  except:
    vars['IsMuon']      = BadFloat
    vars['IsMuonLoose'] = BadFloat
    vars['IsMuonTight'] = BadFloat
  return vars

#Make a dictionary containing all the best primary vertex variables variables
def getBestPVVariables(particle=None, **options):
  if not particle:
    return ['PV', 'FD', 'FDCHI2']  #added 'FD'
  vars = {}
  try:
    table = options['pvtool'].relatedPVs(particle, options['pvs'])
    relPVs = table.relations(particle)
    if not len(relPVs):
      raise
    bestVertex = relPVs.back().to()
    vars['PV'] = options['pvs'].index(bestVertex)
    _fd_chi2 = getFd(particle.endVertex(), bestVertex, options['disttool'])[1]
    _fd = getFd(particle.endVertex(), bestVertex, options['disttool'])[0]
    vars['FDCHI2'] = _fd_chi2
    vars['FD'] = _fd

  except:
    vars['PV']     = -1
    vars['FDCHI2'] = BadFloat
    vars['FD'] = BadFloat
  return vars

#Make a dictionary containing all the 'Standard' variables
def getStandardVariables(particle=None,**options):
  if not particle:
    return getMomentumVariables() + getVertexVariables() + getPIDVariables() + ['ID', 'IP_0', 'IPCHI2_0', 'MINIP', 'MINIPCHI2', 'TRACKCHI2', 'TRACKCHI2NDOF', 'TRACKNDOF', 'CHARGE', 'ghostProb', 'RefPoint_X', 'RefPoint_Y', 'RefPoint_Z'] + getBestPVVariables() + getTriggerVariables() #added 'IP'
  vars = getMomentumVariables(particle.momentum())
  vars.update(getVertexVariables(particle.endVertex() if particle.endVertex() else -1))
  vars.update(getPIDVariables(particle.proto() if particle.proto() else -1))
  vars.update(getBestPVVariables(particle, **options))
  vars.update(getTriggerVariables(particle, **options))
  vars['ID'] = particle.particleID().pid()
  if 'pvs' in options and len(options['pvs']) > 0 and 'disttool' in options:
    _ip = [getIp(particle,testpv,options['disttool'])[0] for testpv in options['pvs']] #Added
    _ip_chi2 = [getIp(particle,testpv,options['disttool'])[1] for testpv in options['pvs']] #Added
    vars['MINIPCHI2'] = min(_ip_chi2)
    vars['IP_0'] = _ip[0]
    vars['IPCHI2_0'] = _ip_chi2[0]
    if np.isnan(_ip_chi2).any():
      vars['MINIP'] = np.nan
    else:
      vars['MINIP'] = np.array(_ip)[np.where(np.array(_ip_chi2) == min(np.array(_ip_chi2)))[0][0]]
  else:
    vars['MINIPCHI2'] = BadFloat
    vars['MINIP'] = BadFloat
  if particle.proto():
    vars['TRACKCHI2']     = particle.proto().track().chi2()
    vars['TRACKCHI2NDOF'] = particle.proto().track().chi2PerDoF()
    vars['TRACKCHI2NDOF'] = particle.proto().track().nDoF()
    vars['CHARGE']        = particle.proto().track().charge()
    vars['ghostProb']     = particle.proto().track().ghostProbability()
  else:
    vars['TRACKCHI2']     = BadFloat
    vars['TRACKCHI2NDOF'] = BadFloat
    vars['TRACKCHI2NDOF'] = BadFloat
    vars['CHARGE']        = BadFloat
    vars['ghostProb']     = BadFloat
  if particle.referencePoint():
    vars['RefPoint_X']     = particle.referencePoint().X()
    vars['RefPoint_Y']     = particle.referencePoint().Y()
    vars['RefPoint_Z']     = particle.referencePoint().Z()
  else:
    vars['RefPoint_X']     = BadFloat
    vars['RefPoint_Y']     = BadFloat
    vars['RefPoint_Z']     = BadFloat
    
  return vars

#Make a dictionary containing all the 'Standard' variables for an MC particle
def getStandardMCVariables(particle=None,**options):
  if not particle:
    return getMomentumVariables() + ['VX', 'VY', 'VZ', 'OVX', 'OVY', 'OVZ', 'ID']
  vars = getMomentumVariables(particle.momentum())
  endvertex = None
  for vtx in particle.endVertices():
    if vtx.isDecay():
      endvertex = vtx  
  if endvertex:
    vars['VX']        = endvertex.position().X()
    vars['VY']        = endvertex.position().Y()
    vars['VZ']        = endvertex.position().Z()
  else:
    vars['VX']        = BadFloat
    vars['VY']        = BadFloat
    vars['VZ']        = BadFloat
  vars['OVX']        = particle.originVertex().position().X()
  vars['OVY']        = particle.originVertex().position().Y()
  vars['OVZ']        = particle.originVertex().position().Z()
  vars['ID'] = particle.particleID().pid()
  return vars

###############################################################################
# What follows is a bunch of functions and classes to write various lists to  #
# the tuples                                                                  #
###############################################################################
def addGenericVariables(t, name, p, generator):
  vars = generator(p)
  for var in generator():
    t.column(name+'_'+var, vars[var])


def addStandardVariables(t, name, particle):
  vars = getStandardVariables(particle)
  for var in getStandardVariables():
    t.column(name+'_'+var, vars[var])

def addStandardMCVariables(t, name, particle):
  vars = getStandardMCVariables(particle)
  for var in getStandardMCVariables():
    t.column(name+'_'+var, vars[var])

def addMomentumVariables(t, name, momentum):
  vars = getMomentumVariables(momentum)
  for var in getMomentumVariables():
    t.column(name+'_'+var, vars[var])


class GenericArrayVariables():
  def __init__(self, maxentries, generator):
    self.maxentries = maxentries
    self.generator = generator
    self.variables = {}
    for variable in self.generator():
      self.variables[variable] = std.vector('double')()

  def add(self, p, **options):
    if self.variables.itervalues().next().size() >= self.maxentries:
      return
    vars = self.generator(p, **options)
    for var in vars:
      self.variables[var].push_back(vars[var])

  def addToTuple(self, tup, name, Nname=None):
    for variable in self.generator():
      if not Nname:
        self.variables[variable].resize(self.maxentries, -1000.)
        tup.array(name+'_'+variable, self.variables[variable])
      else:
        tup.farray(name+'_'+variable, self.variables[variable], Nname, self.maxentries)

class StandardParticleArrayVaraiables(GenericArrayVariables):
  def __init__(self, maxentries):
    GenericArrayVariables.__init__(self, maxentries, getStandardVariables)

class StandardMCParticleArrayVaraiables(GenericArrayVariables):
  def __init__(self, maxentries):
    GenericArrayVariables.__init__(self, maxentries, getStandardMCVariables)

class StandardTrackArrayVaraiables(GenericArrayVariables):
  def __init__(self, maxentries):
    GenericArrayVariables.__init__(self, maxentries, getTrackVariables)

class StandardVertexArrayVaraiables(GenericArrayVariables):
  def __init__(self, maxentries):
    GenericArrayVariables.__init__(self, maxentries, getVertexVariables)

class StandardResidualArrayVariables(GenericArrayVariables):
  def __init__(self, maxentries):
    GenericArrayVariables.__init__(self, maxentries, getClusterResidualVariables)


