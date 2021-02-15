###############################################################################
#                                                                             #
# This file contains a number of utility functions used to help with the      #
# reconstruction of the particles.                                            #
#                                                                             #
###############################################################################

from Bender.Main import *

from ROOT import *
import math
from LoKiCore.basic import LHCb
import LoKiPhysMC.Track2MC_Configuration
#import Relations.Rels
import LoKiMC.MC
from GaudiPython.Bindings import gbl

idtoname = {
               11:'e-',
              -11:'e+',
               12:'nu_e',
              -12:'nubar_e',
               13:'mu-',
              -13:'mu+',
               14:'nu_mu',
              -14:'nubar_mu',
               15:'tau-',
              -15:'tau+',
               16:'nu_tau',
              -16:'nubar_tau',
               22:'gamma',
              111:'pi0',
              211:'pi+',
             -211:'pi-',
              113:'rho0',
              213:'rho+',
             -213:'rho-',
              130:'K0-Long',
              310:'K0-Short',
              311:'K0',
             -311:'K~0',
              321:'K+',
             -321:'K-',
              411:'D+',
             -411:'D-',
              413:'D*(2010)+',
             -413:'D*(2010)-',
              421:'D0',
             -421:'D~0',
              423:'D*(2007)0',
             -423:'D*~(2007)0',
              431:'Ds+',
             -431:'Ds-',
              511:'B0',
             -511:'B~0',
              521:'B+',
             -521:'B-',
           }



# just returns the name of a particle if it's in the list above
def getName(pid):
  try:
    return idtoname[pid]
  except:
    return pid

#Produces a Particle from a track fit
def trackToParticle(track, pid):
  proto = LHCb.ProtoParticle()
  proto.setTrack(track)
  particle = LHCb.Particle(LHCb.ParticleID(pid))

  mass = {13:105.6583715, 211:139.57018, 321:493.679, 2212:938.272081}[abs(pid)]
  particle.setProto(proto)
  pvec = gbl.Gaudi.XYZPoint()
  pcov = gbl.Gaudi.SymMatrix3x3()
  track.position(pvec, pcov)
  particle.setReferencePoint(pvec)
  particle.setPosCovMatrix(pcov)
  p = track.momentum()
  momentum = gbl.Gaudi.LorentzVector()
  momentum.SetPxPyPzE(p.x(), p.y(), p.z(), math.sqrt(p.Mag2() + mass*mass))
  particle.setMomentum(momentum)
  particle.setMeasuredMass(mass)
  return particle

#Produces a Particle from a protoparticle
def protoToParticle(proto, pid):
  particle = LHCb.Particle(LHCb.ParticleID(pid))
  track = proto.track()

  mass = {13:105.6583715, 211:139.57018, 321:493.679, 2212:938.272081}[abs(pid)]
  particle.setProto(proto)
  pvec = gbl.Gaudi.XYZPoint()
  pcov = gbl.Gaudi.SymMatrix3x3()
  track.position(pvec, pcov)
  particle.setReferencePoint(pvec)
  particle.setPosCovMatrix(pcov)
  p = track.momentum()
  momentum = gbl.Gaudi.LorentzVector()
  momentum.SetPxPyPzE(p.x(), p.y(), p.z(), math.sqrt(p.Mag2() + mass*mass))
  particle.setMomentum(momentum)
  particle.setMeasuredMass(mass)
  return particle

#Searches for all protoparticles that we generated from a give MC particle
def findProtos(mcparticle, protos):
  retlist = []
  table = get('Relations/Rec/Track/Default') #tracks to particles
  for p in protos:
    mcps = table.relations(p.track())
    for mcp in mcps:
      if mcparticle == mcp._to():
        retlist.append(p)
  return retlist

#Reconstruct the tau, this doesn't do anything fun yet
def reconstructTau(rtc, pid, persist):
  tauvx    = gbl.LHCb.Vertex()
  if rtc.referencePoint():
    tauvx.setPosition(rtc.referencePoint())
  recotau  = rtc.clone()
  recotau.setParticleID(pid)
  recotau.addToDaughters(rtc)
  recotau.setEndVertex(tauvx)
  persist.append(recotau)
  persist.append(tauvx)
  return recotau

  

