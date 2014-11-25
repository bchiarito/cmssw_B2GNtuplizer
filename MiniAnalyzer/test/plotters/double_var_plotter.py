import os
import glob
import math
import ROOT
from ROOT import *
import sys

from optparse import OptionParser

parser = OptionParser()

parser.add_option('--cut', metavar='F', type='string', action='store',
                  dest='cut',
                  help='')

parser.add_option('--var1', metavar='F', type='string', action='store',
                  dest='var1',
                  help='')

parser.add_option('--var2', metavar='F', type='string', action='store',
                  dest='var2',
                  help='')

parser.add_option('--Min', metavar='F', type='float', action='store',
                  dest='Min',
                  help='')

parser.add_option('--Max', metavar='F', type='float', action='store',
                  dest='Max',
                  help='')

parser.add_option('--name', metavar='F', type='string', action='store',
	    	  default = "blahblahblah",
                  dest='name',
                  help='')

parser.add_option('--log', metavar='F', type='string', action='store',
                  default='no',
                  dest='log',
                  help='')

parser.add_option('--scale', metavar='F', type='float', action='store',
                  default='1.0',
                  dest='scale',
                  help='')

parser.add_option('--bin', metavar='F', type='int', action='store',
                  default=100,
                  dest='bin',
                  help='')

parser.add_option('--file', metavar='F', type='string', action='store',
                  default='no',
                  dest='fi',
                  help='')

(options, args) = parser.parse_args()

scale = options.scale
cut = options.cut
var1 = options.var1
var2 = options.var2
x = options.Min
y = options.Max
log = options.log
bin = options.bin
fi = options.fi
name = options.name

#f = ROOT.TFile( options.name + ".root", "recreate" )
#f.cd()

chain = ROOT.TChain("tree")
chain.Add(fi)
newhist = ROOT.TH1F(name, name, bin, x, y)	
chain.Draw(var1+">>"+name,""+ cut, "goff")
#newhist.Scale(1/newhist.Integral())
newhist.SetLineColor(ROOT.kBlue)
newhist.SetFillColor(0)
newhist.SetLineWidth(2)
newhist.SetLineStyle(2)	
newhist.SetStats(0)
#f.Write()

name2 = name + "_two"
newhist2 = ROOT.TH1F(name2, name2, bin, x, y)	
chain.Draw(var2+">>"+name2,""+ cut, "goff")
#newhist2.Scale(1/newhist2.Integral())
newhist2.SetLineColor(ROOT.kRed)
newhist2.SetFillColor(0)
newhist2.SetLineWidth(2)
newhist2.SetLineStyle(2)	
newhist2.SetStats(0)


c = TCanvas()
c.cd()
newhist.SetTitle(name)
newhist.GetXaxis().SetTitle(var1 + "  w/  " + cut)
if cut == "":
	newhist.GetXaxis().SetTitle(var1)
newhist.GetYaxis().SetTitle("events")

if log == "yes":
	c.SetLogy()

newhist.SetMaximum(max(newhist.GetMaximum(), newhist2.GetMaximum()) * 1.05)

newhist.Draw()
newhist2.Draw("same")
c.SaveAs(name + ".png")

print str(newhist.GetEntries())

raw_input()

