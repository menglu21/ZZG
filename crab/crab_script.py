import os
import sys
import optparse
import ROOT
import re

from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from PhysicsTools.NanoAODTools.postprocessing.modules.common.countHistogramsModule import *
from PhysicsTools.NanoAODTools.postprocessing.analysis.modules.eleRECOSFProducer import *
from PhysicsTools.NanoAODTools.postprocessing.analysis.modules.eleIDSFProducer import *
from PhysicsTools.NanoAODTools.postprocessing.analysis.modules.photonRECOSFProducer import *
from PhysicsTools.NanoAODTools.postprocessing.analysis.modules.photonIDSFProducer import *
from PhysicsTools.NanoAODTools.postprocessing.analysis.modules.muonScaleResProducer import *
from PhysicsTools.NanoAODTools.postprocessing.analysis.modules.muonIDISOSFProducer import *
from PhysicsTools.NanoAODTools.postprocessing.analysis.modules.ZZGProducer import *
#from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jetmetUncertainties import *
from PhysicsTools.NanoAODTools.postprocessing.framework.crabhelper import inputFiles, runsAndLumis
### main python file to run ###

def main():

  usage = 'usage: %prog [options]'
  parser = optparse.OptionParser(usage)
  parser.add_option('--year', dest='year', help='which year sample', default='2018', type='string')
  parser.add_option('-m', dest='ismc', help='to apply sf correction or not', default=True, action='store_true')
  parser.add_option('-d', dest='ismc', help='to apply sf correction or not', action='store_false')
  (opt, args) = parser.parse_args()

  if opt.ismc:
    if opt.year == "2016a":
      p = PostProcessor(".", inputFiles(), modules=[countHistogramsModule(),muonIDISOSF2016(),muonScaleRes2016a(),eleRECOSF2016(),eleIDSF2016(),photonRECOSF2016(),photonIDSF2016(), ZZG2016()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis())
    if opt.year == "2016b":
      p = PostProcessor(".", inputFiles(), modules=[countHistogramsModule(),muonIDISOSF2016(),muonScaleRes2016b(),eleRECOSF2016(),eleIDSF2016(),photonRECOSF2016(),photonIDSF2016(), ZZG2016()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis())
    if opt.year == "2017":
      p = PostProcessor(".", inputFiles(), modules=[countHistogramsModule(),muonIDISOSF2017(),muonScaleRes2017(),eleRECOSF2017(),eleIDSF2017(),photonRECOSF2017(),photonIDSF2017(), ZZG2017()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis())
    if opt.year == "2018":
      p = PostProcessor(".", inputFiles(), modules=[countHistogramsModule(),muonIDISOSF2018(),muonScaleRes2018(),eleRECOSF2018(),eleIDSF2018(),photonRECOSF2018(),photonIDSF2018(), ZZG2018()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis())


# Sequence for data
  if not (opt.ismc):
    if opt.year == "2016b" or opt.year == "2016c" or opt.year == "2016d":
      p = PostProcessor(".", inputFiles(), modules=[muonScaleRes2016a(),ZZG2016()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis())
    if opt.year == "2016e" or opt.year == "2016f":
      p = PostProcessor(".", inputFiles(), modules=[muonScaleRes2016a(),ZZG2016()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis())
    if opt.year == "2016g" or opt.year == "2016h":
      p = PostProcessor(".", inputFiles(), modules=[muonScaleRes2016b(),ZZG2016()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis())
    if opt.year == "2017b":
      p = PostProcessor(".", inputFiles(), modules=[muonScaleRes2017(),ZZG2017()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis())
    if opt.year == "2017c":
      p = PostProcessor(".", inputFiles(), modules=[muonScaleRes2017(),ZZG2017()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis())
    if opt.year == "2017d":
      p = PostProcessor(".", inputFiles(), modules=[muonScaleRes2017(),ZZG2017()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis())
    if opt.year == "2017e":
      p = PostProcessor(".", inputFiles(), modules=[muonScaleRes2017(),ZZG2017()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis())
    if opt.year == "2017f":
      p = PostProcessor(".", inputFiles(), modules=[muonScaleRes2017(),ZZG2017()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis())
    if opt.year == "2018a":
      p = PostProcessor(".", inputFiles(), modules=[muonScaleRes2018(),ZZG2018()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis())
    if opt.year == "2018b":
      p = PostProcessor(".", inputFiles(), modules=[muonScaleRes2018(),ZZG2018()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis())
    if opt.year == "2018c":
      p = PostProcessor(".", inputFiles(), modules=[muonScaleRes2018(),ZZG2018()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis())
    if opt.year == "2018d":
      p = PostProcessor(".", inputFiles(), modules=[muonScaleRes2018(),ZZG2018()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis())
  p.run()

if __name__ == "__main__":
    sys.exit(main())
