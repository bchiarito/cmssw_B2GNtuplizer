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

parser.add_option('--Min1', metavar='F', type='float', action='store',
                  dest='Min1',
                  help='')

parser.add_option('--Max1', metavar='F', type='float', action='store',
                  dest='Max1',
                  help='')

parser.add_option('--Min2', metavar='F', type='float', action='store',
                  dest='Min2',
                  help='')

parser.add_option('--Max2', metavar='F', type='float', action='store',
                  dest='Max2',
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

parser.add_option('--bin1', metavar='F', type='int', action='store',
                  default=100,
                  dest='bin1',
                  help='')

parser.add_option('--bin2', metavar='F', type='int', action='store',
                  default=100,
                  dest='bin2',
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
x = options.Min1
y = options.Max1
x2 = options.Min2
y2 = options.Max2
log = options.log
bin = options.bin1
bin2 = options.bin2
fi = options.fi
name = options.name

#f = ROOT.TFile( options.name + ".root", "recreate" )
#f.cd()

chain = ROOT.TChain("tree")
chain.Add(fi)
newhist = ROOT.TH2F(name, name, bin, x, y, bin2, x2, y2)	
chain.Draw(var2+":"+var1+">>"+name,""+ cut, "Colz")
#newhist.Scale(1/newhist.Integral())
newhist.SetMarkerColor(ROOT.kBlue)
newhist.SetFillColor(0)
newhist.SetLineWidth(2)
newhist.SetLineStyle(2)	
newhist.SetStats(0)
#f.Write()

c = TCanvas()
c.cd()
newhist.SetTitle(name)
newhist.GetXaxis().SetTitle(var1)
newhist.GetYaxis().SetTitle(var2)

if log == "yes":
	c.SetLogy()

newhist.Draw()
c.SaveAs(name + ".png")

print str(newhist.GetEntries())

raw_input()

