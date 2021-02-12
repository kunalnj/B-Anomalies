###############################################################################
#                                                                             #
# This file contains the code for an algorithm to run over Bs->phi3 tau tau   #
# MC. It first finds the true MC particles. It then finds the reconstructed   #
# particles the MC particles have produced. Then it combines the two kaons    #
# into a reconstructed phi3, with a vertex fit. It doesn't do anything        #
# complicated with the tau, it just 'reconstructs' them as the muon.          #
#                                                                             #
###############################################################################


from Bender.Main import *

from HistoUtils import book
from TupleUtils import nTuple, releaseTuples
from GaudiPython.Bindings import gbl

from RecoUtilities import *
from TupleUtilities import *
import numpy as np

class Bs2Phi3TauTau(Algo):
  """Analyse the tau tau
  """
  def initialize ( self ):
    # This just sets up the output objects, and all the tools needed to run the reconstruction
    sc = Algo.initialize ( self )
    if sc.isFailure(): return sc

    self.nsignal = 0
    self.ev = 0
    self.histograms = {}
    # how you would book a histogram
    #self.histograms['nPVs'] = book('nPVs', ';nPVs', 100, 0, 100)

    self.tuple = nTuple('DecayTuple','DecayTuple')

    gaudi = appMgr() 
    self.trigtool = gaudi.toolsvc().create('TriggerTisTos', interface='ITriggerTisTos')
    self.disttool = gaudi.toolsvc().create('LoKi::DistanceCalculator',interface='IDistanceCalculator')
    self.vfittool = gaudi.toolsvc().create('LoKi::VertexFitter', interface='IVertexFit')
    self.combtool = gaudi.toolsvc().create('OfflineVertexFitter', interface='IParticleCombiner')
    self.pvtool   = gaudi.toolsvc().create('P2PVWithIPChi2',interface='IRelatedPVFinder')

    return sc

  def finalize(self):
    # Cleans up
    releaseTuples()
    return SUCCESS

  def analyse(self):
    """The main 'analysis' method
    """
    # This is the important function. It is called once per event

    self.ev += 1
    print self.ev
        
    particles = get('/Event/MC/Particles')        #all MC particles in the event
    protos    = get('/Event/Rec/ProtoP/Charged')  #all recontructed particles in the event
    pvs       = get('/Event/Rec/Vertex/Primary')  #all Primary vertices in an event
    print '# particles:', len(particles)
    print '# protos:', len(protos)
    print '# pvs:', len(pvs)

    #Locate all the interesting particles in the decay chain to confirm this is
    #the right decay. In almost all cases the first MC particle is the one we
    #care about, but loop over all particles until we fine it, just in case.
    Bs       = None
    phi3     = None
    tauplus  = None
    tauminus = None
    for particle in particles:
      if particle.particleID().abspid() == 531: #The Bs
        for vertex in particle.endVertices():
          if vertex.isDecay():
            for d in vertex.products():
              if d.particleID().abspid() in [333, 337]: #The phi or phi3
                phi3 = d
              elif d.particleID().pid() == 15:   #the tau-
                tauminus = d
              elif d.particleID().pid() == -15: #the tau+
                tauplus = d
            Bs = particle
            if Bs and phi3 and tauplus and tauminus:
              break
            else:
              Bs       = None
              phi3     = None
              tauplus  = None
              tauminus = None


    if not (Bs and phi3 and tauplus and tauminus):
      print 'No Bs?? That is strange'
      return SUCCESS      ## IMPORTANT!!!

    #Now we know we have a Bs -> phi(3) tau+ tau- decay


    # build the reconstructed phi(3) out of the two kaons.
    Kplus  = None
    Kminus = None
    for vertex in phi3.endVertices():
      if vertex.isDecay():
        for d in vertex.products():
          if d.particleID().pid() == 321:
            Kplus = d
          elif d.particleID().pid() == -321:
            Kminus = d

    phi3Cand = []
    Kplusarrayadder  = StandardParticleArrayVaraiables(20)
    Kminusarrayadder = StandardParticleArrayVaraiables(20)
    phi3arrayadder   = StandardParticleArrayVaraiables(20)

    #Find the reconstructcted protoparticles that correspond to the true MC particles
    recoKp =  [protoToParticle(p, Kplus.particleID().pid())  for p in findProtos(Kplus, protos)]
    recoKm =  [protoToParticle(p, Kminus.particleID().pid()) for p in findProtos(Kminus, protos)]

    if not len(recoKp):
      print 'No reconstructed K+ found, skipping'
      return SUCCESS

    if not len(recoKm):
      print 'No reconstructed K- found, skipping'
      return SUCCESS

    for recoK in recoKp:
      Kplusarrayadder.add(recoK, pvs=pvs, disttool=self.disttool, triggertool=self.trigtool, pvtool=self.pvtool)

    for recoK in recoKm:
      Kminusarrayadder.add(recoK, pvs=pvs, disttool=self.disttool, triggertool=self.trigtool, pvtool=self.pvtool)

    #Do a particle fit for all combinations of K+ and K- that were reconstructed from the true K+ and K-
    for recoP1 in recoKp:
      for recoP2 in recoKm:
        phi3vx   = gbl.LHCb.Vertex()
        recophi3 = gbl.LHCb.Particle(phi3.particleID())
        children = gbl.LHCb.Particle.ConstVector()
        children.push_back(recoP1)
        children.push_back(recoP2)
        try:
          self.combtool.combine(children, recophi3, phi3vx)
          recophi3.setEndVertex(phi3vx)
        except:
          print 'Failed to make particle'
          continue
        phi3Cand.append(recophi3)
        phi3arrayadder.add(recophi3, pvs=pvs, disttool=self.disttool, vfittool = self.vfittool, triggertool=self.trigtool, pvtool=self.pvtool)

    #get the tau+ info
    tauplusCand       = []
    tauplusarrayadder = StandardParticleArrayVaraiables(20)
    muplusarrayadder = StandardParticleArrayVaraiables(20)
    recomuplus = []
    muplus = None

    for vertex in tauplus.endVertices():
      if vertex.isDecay():
        for d in vertex.products():
          if d.particleID().abspid() == 13:
            muplus = d
            break


    recomuplus = [protoToParticle(p, muplus.particleID().pid()) for p in findProtos(muplus, protos)]
    for p in recomuplus:
      muplusarrayadder.add(p, pvs=pvs, disttool=self.disttool, vfittool = self.vfittool, triggertool=self.trigtool, pvtool=self.pvtool)

    #This is where you want to do interesting things to the tau....
    persist = []
    for rtc in recomuplus:
      recotauplus  = reconstructTau(rtc, tauplus.particleID(), persist)
      tauplusCand.append(recotauplus)
      tauplusarrayadder.add(recotauplus, pvs=pvs, disttool=self.disttool, triggertool=self.trigtool, pvtool=self.pvtool)

    #get the tau- info
    tauminusCand       = []
    tauminusarrayadder = StandardParticleArrayVaraiables(20)
    muminusarrayadder = StandardParticleArrayVaraiables(20)
    recomuminus = []
    muminus = None

    for vertex in tauminus.endVertices():
      if vertex.isDecay():
        for d in vertex.products():
          if d.particleID().abspid() == 13:
            muminus = d
            break


    recomuminus = [protoToParticle(p, muminus.particleID().pid()) for p in findProtos(muminus, protos)]
    for p in recomuminus:
      muminusarrayadder.add(p, pvs=pvs, disttool=self.disttool, vfittool = self.vfittool, triggertool=self.trigtool, pvtool=self.pvtool)

    for rtc in recomuminus:
      recotauminus  = reconstructTau(rtc, tauminus.particleID(), persist)
      tauminusCand.append(recotauminus)
      tauminusarrayadder.add(recotauminus, pvs=pvs, disttool=self.disttool, triggertool=self.trigtool, pvtool=self.pvtool)


    if not len(recomuplus):
      print 'No reconstructed mu+ found, skipping'
      return SUCCESS

    if not len(recomuminus):
      print 'No reconstructed mu- found, skipping'
      return SUCCESS


    #Now make the Bs candidate from the phi(3) and the 'taus' (which are really the muons)
    BsCand       = []
    Bsarrayadder = StandardParticleArrayVaraiables(20)
    for recophi3 in phi3Cand:
      for recotauplus in tauplusCand:
        for recotauminus in tauminusCand:
          bsvx     = gbl.LHCb.Vertex()
          recoBs   = gbl.LHCb.Particle(Bs.particleID())
          children = gbl.LHCb.Particle.ConstVector()
          children.push_back(recophi3)
          children.push_back(recotauplus)
          children.push_back(recotauminus)
          try:
            self.combtool.combine(children, recoBs, bsvx)
            recoBs.setEndVertex(bsvx)
          except:
            print 'Failed to make particle'
            continue
          BsCand.append(recoBs)
          Bsarrayadder.add(recoBs, pvs=pvs, disttool=self.disttool, triggertool=self.trigtool, pvtool=self.pvtool)


    #DOCA, IP, FD, VFIT:

    _phi = recophi3[0]

    mu_plus = recomuplus[0]
    mu_minus = recomuminus[0]

    doca_mu_plus = getdoca(self.disttool, mu_plus, _phi)[0]
    doca_mu_minus = getdoca(self.disttool, mu_minus, _phi)[0]
    doca_chi2_mu_plus = getdoca(self.disttool, mu_plus, _phi)[1]
    doca_chi2_mu_minus = getdoca(self.disttool, mu_minus, _phi)[1]

    ip_mu_plus = getIp(mu_plus, _phi.endVertex(), self.disttool)[0]
    ip_mu_minus = getIp(mu_minus, _phi.endVertex(), self.disttool)[0]
    ip_chi2_mu_plus = getIp(mu_plus, _phi.endVertex(), self.disttool)[1]
    ip_chi2_mu_minus = getIp(mu_minus, _phi.endVertex(), self.disttool)[1]
                           
    vertexfit_plus = getvertex(self.vfittool, mu_plus, _phi)[0]
    vertexfit_minus = getvertex(self.vfittool, mu_minus, _phi)[0]

    vfit_plus_x = vertexfit_plus.position().x()
    vfit_plus_y = vertexfit_plus.position().y()
    vfit_plus_z = vertexfit_plus.position().z()

    vfit_minus_x = vertexfit_minus.position().x()
    vfit_minus_y = vertexfit_minus.position().y()
    vfit_minus_z = vertexfit_minus.position().z()

    tau_plus = recotauplus[0]
    tau_minus = recotauminus[0]

    fd_tau_plus = getFd(vertexfit_plus, _phi.endVertex(), self.disttool)[0]
    fd_tau_minus = getFd(vertexfit_minus, _phi.endVertex(), self.disttool)[0]
    fd_chi2_tau_plus = getFd(vertexfit_plus, _phi.endVertex(), self.disttool)[1]
    fd_chi2_tau_minus = getFd(vertexfit_minus, _phi.endVertex(), self.disttool)[1]

    print 'Reconstruction over. Found:'
    print len(recoKp), 'K+'
    print len(recoKm), 'K-'
    print len(phi3Cand), 'phi(1080) or phi3(1850)0'
    print len(tauplusCand), 'tau+'
    print len(tauminusCand), 'tau-'
    print len(BsCand), 'B_s'

    #DOCA:
    
    print doca_mu_plus, "doca mu +"
    print doca_chi2_mu_plus, "doca chi2 mu +"
    print doca_mu_minus, "doca mu -"
    print doca_chi2_mu_minus, "doca chi2 mu -"
    
    #IP:

    print ip_mu_plus, "ip mu +"
    print ip_chi2_mu_plus, "ip chi2 mu +"
    print ip_mu_minus, "ip mu -"
    print ip_chi2_mu_minus, "ip chi2 mu -"

    # #FD:

    print fd_tau_plus, "fd mu +"
    print fd_chi2_tau_plus, "fd chi2 mu +"
    print fd_tau_minus, "fd mu -"
    print fd_chi2_tau_minus, "fd chi2 mu -"

    #Now to add everything to the output tuple

    addStandardMCVariables(self.tuple, 'Bs_TRUE', Bs)
    Bsarrayadder.addToTuple(self.tuple, 'Bs', 'nBs')

    addStandardMCVariables(self.tuple, 'phi3_TRUE', phi3)
    phi3arrayadder.addToTuple(self.tuple, 'phi3', 'nphi3')
    addStandardMCVariables(self.tuple, 'Kplus_TRUE', Kplus)
    Kplusarrayadder.addToTuple(self.tuple, 'Kplus', 'nKplus')
    addStandardMCVariables(self.tuple, 'Kminus_TRUE', Kminus)
    Kminusarrayadder.addToTuple(self.tuple, 'Kminus', 'nKminus')

    addStandardMCVariables(self.tuple, 'tauplus_TRUE', tauplus)
    tauplusarrayadder.addToTuple(self.tuple, 'tauplus', 'ntauplus')
    muplusarrayadder.addToTuple(self.tuple, 'muplus', 'nmuplus')

    addStandardMCVariables(self.tuple, 'tauminus_TRUE', tauminus)
    tauminusarrayadder.addToTuple(self.tuple, 'tauminus', 'ntauminus')
    muminusarrayadder.addToTuple(self.tuple, 'muminus', 'nmuminus')

    pvarrayadder = StandardVertexArrayVaraiables(10)
    for pv in pvs:
      pvarrayadder.add(pv)
    pvarrayadder.addToTuple(self.tuple, 'PV', 'nPV')

    #append DOCA:

    self.tuple.column("DOCA_mu_plus", doca_mu_plus)
    self.tuple.column("DOCA_chi2_mu_plus", doca_chi2_mu_plus)
    self.tuple.column("DOCA_mu_minus", doca_mu_minus)
    self.tuple.column("DOCA_chi2_mu_minus", doca_chi2_mu_minus)

    #append IP:

    self.tuple.column("IP_mu_plus", ip_mu_plus)
    self.tuple.column("IP_chi2_mu_plus", ip_chi2_mu_plus)
    self.tuple.column("IP_mu_minus", ip_mu_minus)
    self.tuple.column("IP_chi2_mu_minus", ip_chi2_mu_minus)

    #append FD:

    self.tuple.column("FD_tau_plus", fd_tau_plus)
    self.tuple.column("FD_chi2_tau_plus", fd_chi2_tau_plus)
    self.tuple.column("FD_tau_minus", fd_tau_minus)
    self.tuple.column("FD_chi2_tau_minus", fd_chi2_tau_minus)
    
    # append vfits:
    
    self.tuple.column("VFIT_plus_x", vfit_plus_x)
    self.tuple.column("VFIT_plus_y", vfit_plus_y)
    self.tuple.column("VFIT_plus_z", vfit_plus_z)
    
    self.tuple.column("VFIT_minus_x", vfit_minus_x)
    self.tuple.column("VFIT_minus_y", vfit_minus_y)
    self.tuple.column("VFIT_minus_z", vfit_minus_z)
    
    
    ##################################################################################

    self.tuple.column('eventNumber', int(get('/Event/DAQ/ODIN').eventNumber()))
    self.tuple.column('runNumber', int(get('/Event/DAQ/ODIN').runNumber()))

    # write out all the reconstruction summary information, so number of tracks of each type, etc
    recsummary = get('/Event/Rec/Summary')
    self.tuple.column('nPVs',recsummary.info(recsummary.nPVs, 0))
    self.tuple.column('nLongTracks',recsummary.info(recsummary.nLongTracks, 0))
    self.tuple.column('nDownstreamTracks',recsummary.info(recsummary.nDownstreamTracks, 0))
    self.tuple.column('nUpstreamTracks',recsummary.info(recsummary.nUpstreamTracks, 0))
    self.tuple.column('nVeloTracks',recsummary.info(recsummary.nVeloTracks, 0))
    self.tuple.column('nTTracks',recsummary.info(recsummary.nTTracks, 0))
    self.tuple.column('nBackTracks',recsummary.info(recsummary.nBackTracks, 0))
    self.tuple.column('nTracks',recsummary.info(recsummary.nTracks, 0))
    self.tuple.column('nGhosts',recsummary.info(recsummary.nGhosts, 0))
    self.tuple.column('nRich1Hits',recsummary.info(recsummary.nRich1Hits, 0))
    self.tuple.column('nRich2Hits',recsummary.info(recsummary.nRich2Hits, 0))
    self.tuple.column('nVeloClusters',recsummary.info(recsummary.nVeloClusters, 0))
    self.tuple.column('nITClusters',recsummary.info(recsummary.nITClusters, 0))
    self.tuple.column('nTTClusters',recsummary.info(recsummary.nTTClusters, 0))
    self.tuple.column('nUTClusters',recsummary.info(recsummary.nUTClusters, 0))
    self.tuple.column('nOTClusters',recsummary.info(recsummary.nOTClusters, 0))
    self.tuple.column('nFTClusters',recsummary.info(recsummary.nFTClusters, 0))
    self.tuple.column('nSPDhits',recsummary.info(recsummary.nSPDhits, 0))
    self.tuple.column('nMuonCoordsS0',recsummary.info(recsummary.nMuonCoordsS0, 0))
    self.tuple.column('nMuonCoordsS1',recsummary.info(recsummary.nMuonCoordsS1, 0))
    self.tuple.column('nMuonCoordsS2',recsummary.info(recsummary.nMuonCoordsS2, 0))
    self.tuple.column('nMuonCoordsS3',recsummary.info(recsummary.nMuonCoordsS3, 0))
    self.tuple.column('nMuonCoordsS4',recsummary.info(recsummary.nMuonCoordsS4, 0))
    self.tuple.column('nMuonTracks',recsummary.info(recsummary.nMuonTracks, 0))

    # write out the number of proto particles created
    self.tuple.column('nProtoCharged',len(get('/Event/Rec/ProtoP/Charged')))
    self.tuple.column('nProtoNeutral',len(get('/Event/Rec/ProtoP/Neutrals')))
    self.tuple.column('nCaloElectrons',len(get('/Event/Rec/Calo/Electrons')))
    self.tuple.column('nCaloPhotons',len(get('/Event/Rec/Calo/Photons')))

    self.tuple.write()
    return SUCCESS      ## IMPORTANT!!!




  analyse.counter = 0
# =============================================================================

