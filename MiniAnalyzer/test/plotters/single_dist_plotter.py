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

parser.add_option('--var', metavar='F', type='string', action='store',
                  dest='var',
                  help='')

parser.add_option('--Min', metavar='F', type='float', action='store',
                  dest='Min',
                  help='')

parser.add_option('--Max', metavar='F', type='float', action='store',
                  dest='Max',
                  help='')

parser.add_option('--name', metavar='F', type='string', action='store',
	    	  default = "blank",
                  dest='name',
                  help='')

parser.add_option('--log', action='store_true', default=False,
                  dest='log',
                  help='log sacle on y')

parser.add_option('--scale', action='store_true', default=False,
                  dest='scale',
                  help='scale to integral = 1')

parser.add_option('--bin', metavar='F', type='int', action='store',
                  default=100,
                  dest='bin',
                  help='')

parser.add_option('--file', metavar='F', type='string', action='store',
                  default='no',
                  dest='fi',
                  help='')

parser.add_option('--save', action='store_true', default=False,
                  dest='save',
                  help='save plot')

parser.add_option('--title', metavar='F', type='string', action='store',
	    	  default = "blank",
                  dest='title',
                  help='')

parser.add_option('--noplot', action='store_true', default=False,
                  dest='noplot',
                  help='Do no plot anything')



(options, args) = parser.parse_args()

scale = options.scale
cut = options.cut
var = options.var
x = options.Min
y = options.Max
log = options.log
bin = options.bin
fi = options.fi
name = options.name
title = options.title

#f = ROOT.TFile( options.name + ".root", "recreate" )
#f.cd()

chain = ROOT.TChain("tree")
chain.Add(fi)
newhist = ROOT.TH1F(name, name, bin, x, y)	
chain.Draw(var+">>"+name,""+ cut, "goff")
if scale:
  newhist.Scale(1/newhist.Integral())
newhist.SetLineColor(ROOT.kBlue)
newhist.SetFillColor(0)
newhist.SetLineWidth(2)
newhist.SetLineStyle(2)	
newhist.SetStats(0)
#f.Write()

if not options.noplot:
  c = TCanvas()
  c.cd()
  newhist.SetTitle(title)
  newhist.GetXaxis().SetTitle(var + "  w/  " + cut)
  if cut == "":
    newhist.GetXaxis().SetTitle(var)
  newhist.GetYaxis().SetTitle("Events")
  if log:
    c.SetLogy()
  newhist.Draw()
  leg = ROOT.TLegend(0.55, 0.85, 0.9, 0.9)
  leg.AddEntry(newhist, "T1tttt_mGo1300_mStop300_mCh285_mChi280", "l")
  leg.Draw("same")
  if options.save == True:
    c.SaveAs(name + ".png")

print "file: " + str(newhist.GetEntries())

if options.save == False and options.noplot == False:
  raw_input()

