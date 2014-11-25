#! /usr/bin/env python
import os
import glob
import sys
from DataFormats.FWLite import Events, Handle
import ROOT
from optparse import OptionParser
from ttreemaker import treemaker

parser = OptionParser()

parser.add_option('-d', '--data', action='store_true', default=False,
                  dest='isData',
                  help='Set to data')
parser.add_option('-m', '--mc', action='store_false', default=False,
                  dest='isData',
                  help='Set to MC')
parser.add_option('-n', '--name', default='my_tree',
                  dest='name',
                  help='Name of TTree')
parser.add_option('-f', '--file',
                  dest='files',
                  help='File or group of files using a wildcard (remember to use \\ to input a wildcard)')
parser.add_option('-M', '--max', type='int', default=-1,
		  dest='max',
		  help='Maximum number of events to process')
parser.add_option('-w', '--weight', default=1.0,
		  dest='weight',
		  help='Weight of each event as it is added to the TTree')

(options, args) = parser.parse_args()

print options.name  
files = glob.glob( options.files )
print len(files)
events = Events (files)
ntotal = events.size()
filename = options.name
analyzer = treemaker(filename, options.isData, options.weight)
count = 0
print "Start event loop"
for event in events:
	count = count + 1
	if count % 10000 == 0 or count == 1:
      		percentDone = float(count) / float(ntotal) * 100.0
       		print 'Processing {0:10.0f}/{1:10.0f} : {2:5.2f} %'.format(count, ntotal, percentDone )
	num = analyzer.analyze(event)
	analyzer.reset(num)
	if count == options.max:
      		percentDone = float(count) / float(ntotal) * 100.0
       		print 'Processed  {0:10.0f}/{1:10.0f} : {2:5.2f} %'.format(count, ntotal, percentDone )
		print 'Hit Max, Done'
		break
del analyzer
