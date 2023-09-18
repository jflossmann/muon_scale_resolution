import ROOT

def rochester():
  ROOT.gROOT.Reset()
  ROOT.gROOT.SetBatch()
  ROOT.gStyle.SetOptStat(0)

  Anno_2018=True
  lumi=86.6855
  if Anno_2018: lumi=373

  reweight=1

  crossSektion=13730 #pb (from DAS)

  
