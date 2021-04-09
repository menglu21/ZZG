import os
import sys
import optparse
import ROOT
import re

from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from PhysicsTools.NanoAODTools.postprocessing.analysis.modules.eleRECOSFProducer import *
from PhysicsTools.NanoAODTools.postprocessing.analysis.modules.eleIDSFProducer import *
from PhysicsTools.NanoAODTools.postprocessing.analysis.modules.photonRECOSFProducer import *
from PhysicsTools.NanoAODTools.postprocessing.analysis.modules.photonIDSFProducer import *
from PhysicsTools.NanoAODTools.postprocessing.analysis.modules.muonScaleResProducer import *
from PhysicsTools.NanoAODTools.postprocessing.analysis.modules.muonIDISOSFProducer import *
from PhysicsTools.NanoAODTools.postprocessing.analysis.modules.ZZGProducer import *
#from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jetmetUncertainties import *

### main python file to run ###

def main():

  usage = 'usage: %prog [options]'
  parser = optparse.OptionParser(usage)
  parser.add_option('--year', dest='year', help='which year sample', default='2018', type='string')
  parser.add_option('-m', dest='ismc', help='to apply sf correction or not', default=True, action='store_true')
  parser.add_option('-d', dest='ismc', help='to apply sf correction or not', action='store_false')
  parser.add_option('-i', '--in', dest='inputs', help='input directory with files', default=None, type='string')
  parser.add_option('-o', '--out', dest='output', help='output directory with files', default=None, type='string')
  (opt, args) = parser.parse_args()
  print 'year:', opt.year

  if opt.ismc:
    if opt.year == "2016a":
      p = PostProcessor(opt.output, [opt.inputs], branchsel=None, modules=[muonIDISOSF2016(),muonScaleRes2016a(),eleRECOSF2016(),eleIDSF2016(),photonRECOSF2016(),photonIDSF2016()], provenance=True)
    if opt.year == "2016b":
      p = PostProcessor(opt.output, [opt.inputs], branchsel=None, modules=[muonIDISOSF2016(),muonScaleRes2016b(),eleRECOSF2016(),eleIDSF2016(),photonRECOSF2016(),photonIDSF2016()], provenance=True)
    if opt.year == "2017":
      p = PostProcessor(opt.output, [opt.inputs], branchsel=None, modules=[muonIDISOSF2017(),muonScaleRes2017(),eleRECOSF2017(),eleIDSF2017(),photonRECOSF2017(),photonIDSF2017(), ZZG2017()], provenance=True)
    if opt.year == "2018":
      p = PostProcessor(opt.output, [opt.inputs], branchsel=None, modules=[muonIDISOSF2018(),muonScaleRes2018(),eleRECOSF2018(),eleIDSF2018(),photonRECOSF2018(),photonIDSF2018()], provenance=True)


# Sequence for data
  if not (opt.ismc):
    if opt.year == "2016b" or opt.year == "2016c" or opt.year == "2016d":
      p = PostProcessor(opt.output, [opt.inputs], branchsel=None, modules=[muonScaleRes2016a()], provenance=True)
    if opt.year == "2016e" or opt.year == "2016f":
      p = PostProcessor(opt.output, [opt.inputs], branchsel=None, modules=[muonScaleRes2016a()], provenance=True)
    if opt.year == "2016g" or opt.year == "2016h":
      p = PostProcessor(opt.output, [opt.inputs], branchsel=None, modules=[muonScaleRes2016b()], provenance=True)
    if opt.year == "2017b":
      p = PostProcessor(opt.output, [opt.inputs], branchsel=None, modules=[muonScaleRes2017()], provenance=True)
    if opt.year == "2017c":
      p = PostProcessor(opt.output, [opt.inputs], branchsel=None, modules=[muonScaleRes2017()], provenance=True)
    if opt.year == "2017d":
      p = PostProcessor(opt.output, [opt.inputs], branchsel=None, modules=[muonScaleRes2017()], provenance=True)
    if opt.year == "2017e":
      p = PostProcessor(opt.output, [opt.inputs], branchsel=None, modules=[muonScaleRes2017()], provenance=True)
    if opt.year == "2017f":
      p = PostProcessor(opt.output, [opt.inputs], branchsel=None, modules=[muonScaleRes2017()], provenance=True)
    if opt.year == "2018a":
      p = PostProcessor(opt.output, [opt.inputs], branchsel=None, modules=[muonScaleRes2018()], provenance=True)
    if opt.year == "2018b":
      p = PostProcessor(opt.output, [opt.inputs], branchsel=None, modules=[muonScaleRes2018()], provenance=True)
    if opt.year == "2018c":
      p = PostProcessor(opt.output, [opt.inputs], branchsel=None, modules=[muonScaleRes2018()], provenance=True)
    if opt.year == "2018d":
      p = PostProcessor(opt.output, [opt.inputs], branchsel=None, modules=[muonScaleRes2018()], provenance=True)
  p.run()

if __name__ == "__main__":
    sys.exit(main())
