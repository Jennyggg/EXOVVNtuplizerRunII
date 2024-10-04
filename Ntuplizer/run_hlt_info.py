# import ROOT in batch mode
import sys
oldargv = sys.argv[:]
sys.argv = [ '-b-' ]
import ROOT
ROOT.gROOT.SetBatch(True)
sys.argv = oldargv

# load FWLite C++ libraries
ROOT.gSystem.Load("libFWCoreFWLite.so");
ROOT.gSystem.Load("libDataFormatsFWLite.so");
ROOT.FWLiteEnabler.enable()

# load FWlite python libraries
from DataFormats.FWLite import Handle, Events

triggerBits, triggerBitLabel = Handle("edm::TriggerResults"), ("TriggerResults","","HLT")
triggerObjects, triggerObjectLabel  = Handle("std::vector<pat::TriggerObjectStandAlone>"), "slimmedPatTrigger"
triggerPrescales, triggerPrescaleLabel  = Handle("pat::PackedTriggerPrescales"), "patTrigger"

events = Events(sys.argv[1])
passedtriggers = {}

for iev,event in enumerate(events):
    event.getByLabel(triggerBitLabel, triggerBits)
    event.getByLabel(triggerObjectLabel, triggerObjects)
    event.getByLabel(triggerPrescaleLabel, triggerPrescales)
    #print "\nEvent %d: run %6d, lumi %4d, event %12d" % (iev,event.eventAuxiliary().run(), event.eventAuxiliary().luminosityBlock(),event.eventAuxiliary().event())
    #print "\n === TRIGGER PATHS ==="
    names = event.object().triggerNames(triggerBits.product())
    for i in xrange(triggerBits.product().size()):
      if triggerBits.product().accept(i):
         #print "Fired trigger: ",names.triggerName(i), ", prescale ", triggerPrescales.product().getPrescaleForIndex(i)
         if not (names.triggerName(i) in passedtriggers):
           passedtriggers[names.triggerName(i)]=1
         else:
           passedtriggers[names.triggerName(i)]+=1
#        print "Trigger ", names.triggerName(i), ", prescale ", triggerPrescales.product().getPrescaleForIndex(i), ": ", ("PASS" if triggerBits.product().accept(i) else "fail (or not run)")

    if iev > 1000: break

for name in passedtriggers.keys():
  print name, ': ',passedtriggers[name]
