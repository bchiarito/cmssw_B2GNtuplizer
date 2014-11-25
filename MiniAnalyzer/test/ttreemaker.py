# Treemaker

# Error codes
# 0 - no error
# 1 - invalid handle

import ROOT
from DataFormats.FWLite import Events, Handle
from array import array
from math import *
import sys
sys.path.insert(0, '/uscms_data/d2/bchiari1/CMSSW_5_3_13/test/ttree')
from JetTools import *
from lepWmaker import *
from operator import itemgetter

#Returns the index of the jet (from a collection "jets") closest to the given four-vector
def ClosestJet(jets, fourvec):
        DR = 9999.
        index = -1
        for j in range(0,len(jets)):
            if jets[j].Pt() > 0 :
                dR = fourvec.DeltaR(jets[j])
                if dR < DR :
                        DR = fourvec.DeltaR(jets[j])
                        index = j
        return index

#Match by dR
def MatchCol(Col, jet, maxDR=0.4):
        j = -1
        dr = maxDR
        for i in range(len(Col)):
                C = ROOT.TLorentzVector()
                C.SetPtEtaPhiM( Col[i].Pt(), Col[i].eta(), Col[i].phi(), Col[i].mass() )
                dr = abs(jet.DeltaR(C))
                if dr < maxDR :
                        j = i
                        break
        if dr > maxDR:
                return -1
        return j

#Make sure our jets are strictly sorted in pt, note that we start with a vector, but end up with a list.
def ReorderByPt(UnsortedCol):
        ColByPt = []
        ColPts = [(temp_col_index, temp_col.Pt()) for temp_col_index, temp_col in enumerate(UnsortedCol)]
        ColPts.sort(key=itemgetter(1), reverse=True)
        ColIndex = [sorted_col_index[0] for sorted_col_index in ColPts]
        for k in range(len(ColIndex)):
                OrderedJet = ROOT.TLorentzVector()
                OrderedJet.SetPtEtaPhiM( UnsortedCol[ColIndex[k]].Pt(), UnsortedCol[ColIndex[k]].eta(), UnsortedCol[ColIndex[k]].phi(), UnsortedCol[ColIndex[k]].mass())
                ColByPt.append(OrderedJet)
        return ColByPt

#fits for the MET's missing parameters to make a W with a given lepton.
def make_lepW(met, lep):
        newmet = ROOT.TLorentzVector()
        newmet_m = ROOT.TLorentzVector()
        newmet_p = ROOT.TLorentzVector()
        newmet.SetPtEtaPhiM(met.Pt(),0,met.Phi(),0)
        newmet_m.SetPtEtaPhiM(met.Pt(),0,met.Phi(),0)
        newmet_p.SetPtEtaPhiM(met.Pt(),0,met.Phi(),0)
        phivec = [math.cos(met.Phi()), math.sin(met.Phi())]
        P_lep = math.sqrt((lep.Px()*lep.Px())+(lep.Py()*lep.Py())+(lep.Pz()*lep.Pz()))
        P_phi = (lep.Px()*phivec[0])+(lep.Py()*phivec[1])
        b = (80.4*80.4) + (P_lep*P_lep) - (lep.E()*lep.E()) + (2*met.Pt()*P_phi)
        arg = (lep.E()*lep.E()) * ((4*met.Pt()*met.Pt()*((lep.Pz()*lep.Pz())-(lep.E()*lep.E())))+(b*b))
        if arg <= 0:
                Pz_met = lep.Pz()*b/(2*((lep.E()*lep.E()) -(lep.Pz()*lep.Pz())))
                newmet.SetPz(Pz_met)
                newmet.SetE(math.sqrt(newmet.Px()*newmet.Px()+newmet.Py()*newmet.Py()+newmet.Pz()*newmet.Pz()))
                return [newmet, newmet]
        else:
                Pz_met_p = ((lep.Pz()*b)+math.sqrt(arg))/(2*((lep.E()*lep.E()) -(lep.Pz()*lep.Pz())))
                Pz_met_m = ((lep.Pz()*b)-math.sqrt(arg))/(2*((lep.E()*lep.E()) -(lep.Pz()*lep.Pz())))
                newmet_p.SetPz(Pz_met_p)
                newmet_p.SetE(math.sqrt(newmet_p.Px()*newmet_p.Px()+newmet_p.Py()*newmet_p.Py()+newmet_p.Pz()*newmet_p.Pz()))
                newmet_m.SetPz(Pz_met_m)
                newmet_m.SetE(math.sqrt(newmet_m.Px()*newmet_m.Px()+newmet_m.Py()*newmet_m.Py()+newmet_m.Pz()*newmet_m.Pz()))
                return [newmet_p, newmet_m]

class treemaker:
	def __init__(self, name, isdata, weight):
		# variables we'll need:
		self.w = weight
		self.name = name
		if isdata == False:
			print "this is an MC file"
			self.isMC = 1
		else:
			print "this is a DATA file"
			self.isMC = 0
		# handles/labels
		self.p4vec = "vector<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > > "
		# gen:
		self.genHan = Handle("vector<reco::GenParticle>")
		self.genLab = ("prunedGenParticles")
		# event:
		self.metHan = Handle(self.p4vec)
		self.metLab = ("b2g", "met")
		# ca8 jets:
		self.jetHan = Handle(self.p4vec)
		self.jetLab = ("b2g", "CA8prunedjets")
		self.unjetHan = Handle(self.p4vec)
		self.unjetLab = ("b2g", "CA8jets")
		self.t1Han = Handle( "std::vector<float>" )
		self.t1Lab = ("b2g", "CA8tau1")
		self.t2Han = Handle( "std::vector<float>" )
		self.t2Lab = ("b2g", "CA8tau2")		
		self.t3Han = Handle( "std::vector<float>" )
		self.t3Lab = ("b2g", "CA8tau3")
		self.t4Han = Handle( "std::vector<float>" )
		self.t4Lab = ("b2g", "CA8tau4")
		# ak4 jets
		self.ak4jetHan = Handle(self.p4vec)
		self.ak4jetLab = ("b2g", "AK4jets")
		# top tag
		self.toptagjetHan = Handle(self.p4vec)
		self.toptagjetLab = ("b2g", "TopTagjets")
		self.toptagminmassHan = Handle("std::vector<float>")
		self.toptagminmassLab = ("b2g", "TopTagminmass")
		self.toptagtopmassHan = Handle("std::vector<float>")
		self.toptagtopmassLab = ("b2g", "TopTagtopmass")
		self.toptagnsubHan = Handle("std::vector<unsigned int>")
		self.toptagnsubLab = ("b2g", "TopTagnsub")
		self.toptagWmassHan = Handle("std::vector<float>")
		self.toptagWmassLab = ("b2g", "TopTagWmass")
		#leptons:
		self.muonHan = Handle(self.p4vec)
		self.muonLab = ("b2g","muons")
		self.muonRelIsoHan = Handle("std::vector<float>")
		self.muonRelIsoLab = ("b2g","muonsRelIso")
		self.muonIsTightHan = Handle("std::vector<bool>")
		self.muonIsTightLab = ("b2g","muonsIsTight")
		self.elecHan = Handle(self.p4vec)
		self.elecLab = ("b2g","electrons")
		self.elecRelIsoHan = Handle("std::vector<float>")
		self.elecRelIsoLab = ("b2g","electronsRelIso")
		self.__book__()
	def __book__(self):
		self.f = ROOT.TFile( self.name + ".root", "recreate" )
        	self.f.cd()
		self.tree = ROOT.TTree("tree", "tree")
		self.ErrHist = ROOT.TH1F('Errors', 'Error Codes', 10, 0, 5)
		self.jetlist = []
		# General
		self.weight = array('f', [self.w])
		self.addBranch('weight', self.weight)
		self.numjets = array('f', [0.0])
		self.addBranch('numjets', self.numjets)
		self.numak4jets = array('f', [0.0])
		self.addBranch('numak4jets', self.numak4jets)
		self.numtoptagjets = array('f', [0.0])
		self.addBranch('numtoptagjets', self.numtoptagjets)
		self.Ht = array('f', [0.0])
		self.addBranch('Ht', self.Ht)
		self.eventType = array('f', [-1.0])
		self.addBranch('eventType', self.eventType)
		self.metpt = array('f', [-1.0])
		self.addBranch('metpt', self.metpt)
		self.metphi = array('f', [100.0])
		self.addBranch('metphi', self.metphi)
		self.metetaP = array('f', [100.0])
		self.addBranch('metetaP', self.metetaP)
		self.metetaN = array('f', [100.0])
		self.addBranch('metetaN', self.metetaN)
		# Jets
		self.jet1mass = array('f', [-1.0])
		self.addBranch('jet1mass', self.jet1mass)
		self.jet1pt = array('f', [-1.0])
		self.addBranch('jet1pt', self.jet1pt)
		self.jet1phi = array('f', [100.0])
		self.addBranch('jet1phi', self.jet1phi)
		self.jet1eta = array('f', [100.0])
		self.addBranch('jet1eta', self.jet1eta)
		self.jet1csv = array('f', [0.0])
		self.addBranch('jet1csv', self.jet1csv)
		self.jet1tau1 = array('f', [1.0])
		self.addBranch('jet1tau1', self.jet1tau1)
		self.jet1tau2 = array('f', [1.0])
		self.addBranch('jet1tau2', self.jet1tau2)
		self.jet1tau3 = array('f', [1.0])
		self.addBranch('jet1tau3', self.jet1tau3)
		self.jet1tau4 = array('f', [1.0])
		self.addBranch('jet1tau4', self.jet1tau4)
		self.jet2mass = array('f', [-1.0])
		self.addBranch('jet2mass', self.jet2mass)
		self.jet2pt = array('f', [-1.0])
		self.addBranch('jet2pt', self.jet2pt)
		self.jet2phi = array('f', [100.0])
		self.addBranch('jet2phi', self.jet2phi)
		self.jet2eta = array('f', [100.0])
		self.addBranch('jet2eta', self.jet2eta)
		self.jet2csv = array('f', [0.0])
		self.addBranch('jet2csv', self.jet2csv)
		self.jet2tau1 = array('f', [1.0])
		self.addBranch('jet2tau1', self.jet2tau1)
		self.jet2tau2 = array('f', [1.0])
		self.addBranch('jet2tau2', self.jet2tau2)
		self.jet2tau3 = array('f', [1.0])
		self.addBranch('jet2tau3', self.jet2tau3)
		self.jet2tau4 = array('f', [1.0])
		self.addBranch('jet2tau4', self.jet2tau4)
		self.jet3mass = array('f', [-1.0])
		self.addBranch('jet3mass', self.jet3mass)
		self.jet3pt = array('f', [-1.0])
		self.addBranch('jet3pt', self.jet3pt)
		self.jet3phi = array('f', [100.0])
		self.addBranch('jet3phi', self.jet3phi)
		self.jet3eta = array('f', [100.0])
		self.addBranch('jet3eta', self.jet3eta)
		self.jet3csv = array('f', [0.0])
		self.addBranch('jet3csv', self.jet3csv)
		self.jet3tau1 = array('f', [1.0])
		self.addBranch('jet3tau1', self.jet3tau1)
		self.jet3tau2 = array('f', [1.0])
		self.addBranch('jet3tau2', self.jet3tau2)
		self.jet3tau3 = array('f', [1.0])
		self.addBranch('jet3tau3', self.jet3tau3)
		self.jet3tau4 = array('f', [1.0])
		self.addBranch('jet3tau4', self.jet3tau4)
		self.jet4mass = array('f', [-1.0])
		self.addBranch('jet4mass', self.jet4mass)
		self.jet4pt = array('f', [-1.0])
		self.addBranch('jet4pt', self.jet4pt)
		self.jet4phi = array('f', [100.0])
		self.addBranch('jet4phi', self.jet4phi)
		self.jet4eta = array('f', [100.0])
		self.addBranch('jet4eta', self.jet4eta)
		self.jet4csv = array('f', [0.0])
		self.addBranch('jet4csv', self.jet4csv)
		self.jet4tau1 = array('f', [1.0])
		self.addBranch('jet4tau1', self.jet4tau1)
		self.jet4tau2 = array('f', [1.0])
		self.addBranch('jet4tau2', self.jet4tau2)
		self.jet4tau3 = array('f', [1.0])
		self.addBranch('jet4tau3', self.jet4tau3)
		self.jet4tau4 = array('f', [1.0])
		self.addBranch('jet4tau4', self.jet4tau4)
		self.jet5mass = array('f', [-1.0])
		self.addBranch('jet5mass', self.jet5mass)
		self.jet5pt = array('f', [-1.0])
		self.addBranch('jet5pt', self.jet5pt)
		self.jet5phi = array('f', [100.0])
		self.addBranch('jet5phi', self.jet5phi)
		self.jet5eta = array('f', [100.0])
		self.addBranch('jet5eta', self.jet5eta)
		self.jet5csv = array('f', [0.0])
		self.addBranch('jet5csv', self.jet5csv)
		self.jet5tau1 = array('f', [1.0])
		self.addBranch('jet5tau1', self.jet5tau1)
		self.jet5tau2 = array('f', [1.0])
		self.addBranch('jet5tau2', self.jet5tau2)
		self.jet5tau3 = array('f', [1.0])
		self.addBranch('jet5tau3', self.jet5tau3)
		self.jet5tau4 = array('f', [1.0])
		self.addBranch('jet5tau4', self.jet5tau4)
		self.jet6mass = array('f', [-1.0])
		self.addBranch('jet6mass', self.jet6mass)
		self.jet6pt = array('f', [-1.0])
		self.addBranch('jet6pt', self.jet6pt)
		self.jet6phi = array('f', [100.0])
		self.addBranch('jet6phi', self.jet6phi)
		self.jet6eta = array('f', [100.0])
		self.addBranch('jet6eta', self.jet6eta)
		self.jet6csv = array('f', [0.0])
		self.addBranch('jet6csv', self.jet6csv)
		self.jet6tau1 = array('f', [1.0])
		self.addBranch('jet6tau1', self.jet6tau1)
		self.jet6tau2 = array('f', [1.0])
		self.addBranch('jet6tau2', self.jet6tau2)
		self.jet6tau3 = array('f', [1.0])
		self.addBranch('jet6tau3', self.jet6tau3)
		self.jet6tau4 = array('f', [1.0])
		self.addBranch('jet6tau4', self.jet6tau4)
		self.jet7mass = array('f', [-1.0])
		self.addBranch('jet7mass', self.jet7mass)
		self.jet7pt = array('f', [-1.0])
		self.addBranch('jet7pt', self.jet7pt)
		self.jet7phi = array('f', [100.0])
		self.addBranch('jet7phi', self.jet7phi)
		self.jet7eta = array('f', [100.0])
		self.addBranch('jet7eta', self.jet7eta)
		self.jet7csv = array('f', [0.0])
		self.addBranch('jet7csv', self.jet7csv)
		self.jet7tau1 = array('f', [1.0])
		self.addBranch('jet7tau1', self.jet7tau1)
		self.jet7tau2 = array('f', [1.0])
		self.addBranch('jet7tau2', self.jet7tau2)
		self.jet7tau3 = array('f', [1.0])
		self.addBranch('jet7tau3', self.jet7tau3)
		self.jet7tau4 = array('f', [1.0])
		self.addBranch('jet7tau4', self.jet7tau4)
		self.jet8mass = array('f', [-1.0])
		self.addBranch('jet8mass', self.jet8mass)
		self.jet8pt = array('f', [-1.0])
		self.addBranch('jet8pt', self.jet8pt)
		self.jet8phi = array('f', [100.0])
		self.addBranch('jet8phi', self.jet8phi)
		self.jet8eta = array('f', [100.0])
		self.addBranch('jet8eta', self.jet8eta)
		self.jet8csv = array('f', [0.0])
		self.addBranch('jet8csv', self.jet8csv)
		self.jet8tau1 = array('f', [1.0])
		self.addBranch('jet8tau1', self.jet8tau1)
		self.jet8tau2 = array('f', [1.0])
		self.addBranch('jet8tau2', self.jet8tau2)
		self.jet8tau3 = array('f', [1.0])
		self.addBranch('jet8tau3', self.jet8tau3)
		self.jet8tau4 = array('f', [1.0])
		self.addBranch('jet8tau4', self.jet8tau4)
		self.jet9mass = array('f', [-1.0])
		self.addBranch('jet9mass', self.jet9mass)
		self.jet9pt = array('f', [-1.0])
		self.addBranch('jet9pt', self.jet9pt)
		self.jet9phi = array('f', [100.0])
		self.addBranch('jet9phi', self.jet9phi)
		self.jet9eta = array('f', [100.0])
		self.addBranch('jet9eta', self.jet9eta)
		self.jet9csv = array('f', [0.0])
		self.addBranch('jet9csv', self.jet9csv)
		self.jet9tau1 = array('f', [1.0])
		self.addBranch('jet9tau1', self.jet9tau1)
		self.jet9tau2 = array('f', [1.0])
		self.addBranch('jet9tau2', self.jet9tau2)
		self.jet9tau3 = array('f', [1.0])
		self.addBranch('jet9tau3', self.jet9tau3)
		self.jet9tau4 = array('f', [1.0])
		self.addBranch('jet9tau4', self.jet9tau4)
		self.jet10mass = array('f', [-1.0])
		self.addBranch('jet10mass', self.jet10mass)
		self.jet10pt = array('f', [-1.0])
		self.addBranch('jet10pt', self.jet10pt)
		self.jet10phi = array('f', [100.0])
		self.addBranch('jet10phi', self.jet10phi)
		self.jet10eta = array('f', [100.0])
		self.addBranch('jet10eta', self.jet10eta)
		self.jet10csv = array('f', [0.0])
		self.addBranch('jet10csv', self.jet10csv)
		self.jet10tau1 = array('f', [1.0])
		self.addBranch('jet10tau1', self.jet10tau1)
		self.jet10tau2 = array('f', [1.0])
		self.addBranch('jet10tau2', self.jet10tau2)
		self.jet10tau3 = array('f', [1.0])
		self.addBranch('jet10tau3', self.jet10tau3)
		self.jet10tau4 = array('f', [1.0])
		self.addBranch('jet10tau4', self.jet10tau4)
		# ak4 Jets
		self.ak4jet1mass = array('f', [-1.0])
		self.addBranch('ak4jet1mass', self.ak4jet1mass)
		self.ak4jet1pt = array('f', [-1.0])
		self.addBranch('ak4jet1pt', self.ak4jet1pt)
		self.ak4jet1phi = array('f', [100.0])
		self.addBranch('ak4jet1phi', self.ak4jet1phi)
		self.ak4jet1eta = array('f', [100.0])
		self.addBranch('ak4jet1eta', self.ak4jet1eta)
		self.ak4jet1csv = array('f', [0.0])
		self.addBranch('ak4jet1csv', self.ak4jet1csv)
		self.ak4jet2mass = array('f', [-1.0])
		self.addBranch('ak4jet2mass', self.ak4jet2mass)
		self.ak4jet2pt = array('f', [-1.0])
		self.addBranch('ak4jet2pt', self.ak4jet2pt)
		self.ak4jet2phi = array('f', [100.0])
		self.addBranch('ak4jet2phi', self.ak4jet2phi)
		self.ak4jet2eta = array('f', [100.0])
		self.addBranch('ak4jet2eta', self.ak4jet2eta)
		self.ak4jet2csv = array('f', [0.0])
		self.addBranch('ak4jet2csv', self.ak4jet2csv)
		self.ak4jet3mass = array('f', [-1.0])
		self.addBranch('ak4jet3mass', self.ak4jet3mass)
		self.ak4jet3pt = array('f', [-1.0])
		self.addBranch('ak4jet3pt', self.ak4jet3pt)
		self.ak4jet3phi = array('f', [100.0])
		self.addBranch('ak4jet3phi', self.ak4jet3phi)
		self.ak4jet3eta = array('f', [100.0])
		self.addBranch('ak4jet3eta', self.ak4jet3eta)
		self.ak4jet3csv = array('f', [0.0])
		self.addBranch('ak4jet3csv', self.ak4jet3csv)
		self.ak4jet4mass = array('f', [-1.0])
		self.addBranch('ak4jet4mass', self.ak4jet4mass)
		self.ak4jet4pt = array('f', [-1.0])
		self.addBranch('ak4jet4pt', self.ak4jet4pt)
		self.ak4jet4phi = array('f', [100.0])
		self.addBranch('ak4jet4phi', self.ak4jet4phi)
		self.ak4jet4eta = array('f', [100.0])
		self.addBranch('ak4jet4eta', self.ak4jet4eta)
		self.ak4jet4csv = array('f', [0.0])
		self.addBranch('ak4jet4csv', self.ak4jet4csv)
		self.ak4jet5mass = array('f', [-1.0])
		self.addBranch('ak4jet5mass', self.ak4jet5mass)
		self.ak4jet5pt = array('f', [-1.0])
		self.addBranch('ak4jet5pt', self.ak4jet5pt)
		self.ak4jet5phi = array('f', [100.0])
		self.addBranch('ak4jet5phi', self.ak4jet5phi)
		self.ak4jet5eta = array('f', [100.0])
		self.addBranch('ak4jet5eta', self.ak4jet5eta)
		self.ak4jet5csv = array('f', [0.0])
		self.addBranch('ak4jet5csv', self.ak4jet5csv)
		self.ak4jet6mass = array('f', [-1.0])
		self.addBranch('ak4jet6mass', self.ak4jet6mass)
		self.ak4jet6pt = array('f', [-1.0])
		self.addBranch('ak4jet6pt', self.ak4jet6pt)
		self.ak4jet6phi = array('f', [100.0])
		self.addBranch('ak4jet6phi', self.ak4jet6phi)
		self.ak4jet6eta = array('f', [100.0])
		self.addBranch('ak4jet6eta', self.ak4jet6eta)
		self.ak4jet6csv = array('f', [0.0])
		self.addBranch('ak4jet6csv', self.ak4jet6csv)
		self.ak4jet7mass = array('f', [-1.0])
		self.addBranch('ak4jet7mass', self.ak4jet7mass)
		self.ak4jet7pt = array('f', [-1.0])
		self.addBranch('ak4jet7pt', self.ak4jet7pt)
		self.ak4jet7phi = array('f', [100.0])
		self.addBranch('ak4jet7phi', self.ak4jet7phi)
		self.ak4jet7eta = array('f', [100.0])
		self.addBranch('ak4jet7eta', self.ak4jet7eta)
		self.ak4jet7csv = array('f', [0.0])
		self.addBranch('ak4jet7csv', self.ak4jet7csv)
		self.ak4jet8mass = array('f', [-1.0])
		self.addBranch('ak4jet8mass', self.ak4jet8mass)
		self.ak4jet8pt = array('f', [-1.0])
		self.addBranch('ak4jet8pt', self.ak4jet8pt)
		self.ak4jet8phi = array('f', [100.0])
		self.addBranch('ak4jet8phi', self.ak4jet8phi)
		self.ak4jet8eta = array('f', [100.0])
		self.addBranch('ak4jet8eta', self.ak4jet8eta)
		self.ak4jet8csv = array('f', [0.0])
		self.addBranch('ak4jet8csv', self.ak4jet8csv)
		self.ak4jet9mass = array('f', [-1.0])
		self.addBranch('ak4jet9mass', self.ak4jet9mass)
		self.ak4jet9pt = array('f', [-1.0])
		self.addBranch('ak4jet9pt', self.ak4jet9pt)
		self.ak4jet9phi = array('f', [100.0])
		self.addBranch('ak4jet9phi', self.ak4jet9phi)
		self.ak4jet9eta = array('f', [100.0])
		self.addBranch('ak4jet9eta', self.ak4jet9eta)
		self.ak4jet9csv = array('f', [0.0])
		self.addBranch('ak4jet9csv', self.ak4jet9csv)
		self.ak4jet10mass = array('f', [-1.0])
		self.addBranch('ak4jet10mass', self.ak4jet10mass)
		self.ak4jet10pt = array('f', [-1.0])
		self.addBranch('ak4jet10pt', self.ak4jet10pt)
		self.ak4jet10phi = array('f', [100.0])
		self.addBranch('ak4jet10phi', self.ak4jet10phi)
		self.ak4jet10eta = array('f', [100.0])
		self.addBranch('ak4jet10eta', self.ak4jet10eta)
		self.ak4jet10csv = array('f', [0.0])
		self.addBranch('ak4jet10csv', self.ak4jet10csv)
		# Top Tag
		self.toptagjet1mass = array('f', [-1.0])
		self.addBranch('toptagjet1mass', self.toptagjet1mass)
		self.toptagjet1pt = array('f', [-1.0])
		self.addBranch('toptagjet1pt', self.toptagjet1pt)
		self.toptagjet1phi = array('f', [100.0])
		self.addBranch('toptagjet1phi', self.toptagjet1phi)
		self.toptagjet1eta = array('f', [100.0])
		self.addBranch('toptagjet1eta', self.toptagjet1eta)
		self.toptagjet1minmass = array('f', [-1.0])
		self.addBranch('toptagjet1minmass', self.toptagjet1minmass)
		self.toptagjet1topmass = array('f', [-1.0])
		self.addBranch('toptagjet1topmass', self.toptagjet1topmass)
		self.toptagjet1Wmass = array('f', [-1.0])
		self.addBranch('toptagjet1Wmass', self.toptagjet1Wmass)
		self.toptagjet1nsub = array('f', [-1.0])
		self.addBranch('toptagjet1nsub', self.toptagjet1nsub)
		self.toptagjet2mass = array('f', [-1.0])
		self.addBranch('toptagjet2mass', self.toptagjet2mass)
		self.toptagjet2pt = array('f', [-1.0])
		self.addBranch('toptagjet2pt', self.toptagjet2pt)
		self.toptagjet2phi = array('f', [100.0])
		self.addBranch('toptagjet2phi', self.toptagjet2phi)
		self.toptagjet2eta = array('f', [100.0])
		self.addBranch('toptagjet2eta', self.toptagjet2eta)
		self.toptagjet2minmass = array('f', [-1.0])
		self.addBranch('toptagjet2minmass', self.toptagjet2minmass)
		self.toptagjet2topmass = array('f', [-1.0])
		self.addBranch('toptagjet2topmass', self.toptagjet2topmass)
		self.toptagjet2Wmass = array('f', [-1.0])
		self.addBranch('toptagjet2Wmass', self.toptagjet2Wmass)
		self.toptagjet2nsub = array('f', [-1.0])
		self.addBranch('toptagjet2nsub', self.toptagjet2nsub)
		self.toptagjet3mass = array('f', [-1.0])
		self.addBranch('toptagjet3mass', self.toptagjet3mass)
		self.toptagjet3pt = array('f', [-1.0])
		self.addBranch('toptagjet3pt', self.toptagjet3pt)
		self.toptagjet3phi = array('f', [100.0])
		self.addBranch('toptagjet3phi', self.toptagjet3phi)
		self.toptagjet3eta = array('f', [100.0])
		self.addBranch('toptagjet3eta', self.toptagjet3eta)
		self.toptagjet3minmass = array('f', [-1.0])
		self.addBranch('toptagjet3minmass', self.toptagjet3minmass)
		self.toptagjet3topmass = array('f', [-1.0])
		self.addBranch('toptagjet3topmass', self.toptagjet3topmass)
		self.toptagjet3Wmass = array('f', [-1.0])
		self.addBranch('toptagjet3Wmass', self.toptagjet3Wmass)
		self.toptagjet3nsub = array('f', [-1.0])
		self.addBranch('toptagjet3nsub', self.toptagjet3nsub)
		self.toptagjet4mass = array('f', [-1.0])
		self.addBranch('toptagjet4mass', self.toptagjet4mass)
		self.toptagjet4pt = array('f', [-1.0])
		self.addBranch('toptagjet4pt', self.toptagjet4pt)
		self.toptagjet4phi = array('f', [100.0])
		self.addBranch('toptagjet4phi', self.toptagjet4phi)
		self.toptagjet4eta = array('f', [100.0])
		self.addBranch('toptagjet4eta', self.toptagjet4eta)
		self.toptagjet4minmass = array('f', [-1.0])
		self.addBranch('toptagjet4minmass', self.toptagjet4minmass)
		self.toptagjet4topmass = array('f', [-1.0])
		self.addBranch('toptagjet4topmass', self.toptagjet4topmass)
		self.toptagjet4Wmass = array('f', [-1.0])
		self.addBranch('toptagjet4Wmass', self.toptagjet4Wmass)
		self.toptagjet4nsub = array('f', [-1.0])
		self.addBranch('toptagjet4nsub', self.toptagjet4nsub)
		self.toptagjet5mass = array('f', [-1.0])
		self.addBranch('toptagjet5mass', self.toptagjet5mass)
		self.toptagjet5pt = array('f', [-1.0])
		self.addBranch('toptagjet5pt', self.toptagjet5pt)
		self.toptagjet5phi = array('f', [100.0])
		self.addBranch('toptagjet5phi', self.toptagjet5phi)
		self.toptagjet5eta = array('f', [100.0])
		self.addBranch('toptagjet5eta', self.toptagjet5eta)
		self.toptagjet5minmass = array('f', [-1.0])
		self.addBranch('toptagjet5minmass', self.toptagjet5minmass)
		self.toptagjet5topmass = array('f', [-1.0])
		self.addBranch('toptagjet5topmass', self.toptagjet5topmass)
		self.toptagjet5Wmass = array('f', [-1.0])
		self.addBranch('toptagjet5Wmass', self.toptagjet5Wmass)
		self.toptagjet5nsub = array('f', [-1.0])
		self.addBranch('toptagjet5nsub', self.toptagjet5nsub)
		self.toptagjet6mass = array('f', [-1.0])
		self.addBranch('toptagjet6mass', self.toptagjet6mass)
		self.toptagjet6pt = array('f', [-1.0])
		self.addBranch('toptagjet6pt', self.toptagjet6pt)
		self.toptagjet6phi = array('f', [100.0])
		self.addBranch('toptagjet6phi', self.toptagjet6phi)
		self.toptagjet6eta = array('f', [100.0])
		self.addBranch('toptagjet6eta', self.toptagjet6eta)
		self.toptagjet6minmass = array('f', [-1.0])
		self.addBranch('toptagjet6minmass', self.toptagjet6minmass)
		self.toptagjet6topmass = array('f', [-1.0])
		self.addBranch('toptagjet6topmass', self.toptagjet6topmass)
		self.toptagjet6Wmass = array('f', [-1.0])
		self.addBranch('toptagjet6Wmass', self.toptagjet6Wmass)
		self.toptagjet6nsub = array('f', [-1.0])
		self.addBranch('toptagjet6nsub', self.toptagjet6nsub)
		self.toptagjet7mass = array('f', [-1.0])
		self.addBranch('toptagjet7mass', self.toptagjet7mass)
		self.toptagjet7pt = array('f', [-1.0])
		self.addBranch('toptagjet7pt', self.toptagjet7pt)
		self.toptagjet7phi = array('f', [100.0])
		self.addBranch('toptagjet7phi', self.toptagjet7phi)
		self.toptagjet7eta = array('f', [100.0])
		self.addBranch('toptagjet7eta', self.toptagjet7eta)
		self.toptagjet7minmass = array('f', [-1.0])
		self.addBranch('toptagjet7minmass', self.toptagjet7minmass)
		self.toptagjet7topmass = array('f', [-1.0])
		self.addBranch('toptagjet7topmass', self.toptagjet7topmass)
		self.toptagjet7Wmass = array('f', [-1.0])
		self.addBranch('toptagjet7Wmass', self.toptagjet7Wmass)
		self.toptagjet7nsub = array('f', [-1.0])
		self.addBranch('toptagjet7nsub', self.toptagjet7nsub)
		self.toptagjet8mass = array('f', [-1.0])
		self.addBranch('toptagjet8mass', self.toptagjet8mass)
		self.toptagjet8pt = array('f', [-1.0])
		self.addBranch('toptagjet8pt', self.toptagjet8pt)
		self.toptagjet8phi = array('f', [100.0])
		self.addBranch('toptagjet8phi', self.toptagjet8phi)
		self.toptagjet8eta = array('f', [100.0])
		self.addBranch('toptagjet8eta', self.toptagjet8eta)
		self.toptagjet8minmass = array('f', [-1.0])
		self.addBranch('toptagjet8minmass', self.toptagjet8minmass)
		self.toptagjet8topmass = array('f', [-1.0])
		self.addBranch('toptagjet8topmass', self.toptagjet8topmass)
		self.toptagjet8Wmass = array('f', [-1.0])
		self.addBranch('toptagjet8Wmass', self.toptagjet8Wmass)
		self.toptagjet8nsub = array('f', [-1.0])
		self.addBranch('toptagjet8nsub', self.toptagjet8nsub)
		self.toptagjet9mass = array('f', [-1.0])
		self.addBranch('toptagjet9mass', self.toptagjet9mass)
		self.toptagjet9pt = array('f', [-1.0])
		self.addBranch('toptagjet9pt', self.toptagjet9pt)
		self.toptagjet9phi = array('f', [100.0])
		self.addBranch('toptagjet9phi', self.toptagjet9phi)
		self.toptagjet9eta = array('f', [100.0])
		self.addBranch('toptagjet9eta', self.toptagjet9eta)
		self.toptagjet9minmass = array('f', [-1.0])
		self.addBranch('toptagjet9minmass', self.toptagjet9minmass)
		self.toptagjet9topmass = array('f', [-1.0])
		self.addBranch('toptagjet9topmass', self.toptagjet9topmass)
		self.toptagjet9Wmass = array('f', [-1.0])
		self.addBranch('toptagjet9Wmass', self.toptagjet9Wmass)
		self.toptagjet9nsub = array('f', [-1.0])
		self.addBranch('toptagjet9nsub', self.toptagjet9nsub)
		self.toptagjet10mass = array('f', [-1.0])
		self.addBranch('toptagjet10mass', self.toptagjet10mass)
		self.toptagjet10pt = array('f', [-1.0])
		self.addBranch('toptagjet10pt', self.toptagjet10pt)
		self.toptagjet10phi = array('f', [100.0])
		self.addBranch('toptagjet10phi', self.toptagjet10phi)
		self.toptagjet10eta = array('f', [100.0])
		self.addBranch('toptagjet10eta', self.toptagjet10eta)
		self.toptagjet10minmass = array('f', [-1.0])
		self.addBranch('toptagjet10minmass', self.toptagjet10minmass)
		self.toptagjet10topmass = array('f', [-1.0])
		self.addBranch('toptagjet10topmass', self.toptagjet10topmass)
		self.toptagjet10Wmass = array('f', [-1.0])
		self.addBranch('toptagjet10Wmass', self.toptagjet10Wmass)
		self.toptagjet10nsub = array('f', [-1.0])
		self.addBranch('toptagjet10nsub', self.toptagjet10nsub)
		# Leps
		self.Htlep = array('f', [-1.0])
		self.addBranch('Htlep', self.Htlep)
		self.numleptons = array('f', [-1.0])
		self.addBranch('numleptons', self.numleptons)
		self.numgoodleptons = array('f', [-1.0])
		self.addBranch('numgoodleptons', self.numgoodleptons)
		self.lep1id = array('f', [-1.0])
		self.addBranch('lep1id', self.lep1id)
		self.lep1pt = array('f', [-1.0])
		self.addBranch('lep1pt', self.lep1pt)
		self.lep1phi = array('f', [100.0])
		self.addBranch('lep1phi', self.lep1phi)
		self.lep1eta = array('f', [100.0])
		self.addBranch('lep1eta', self.lep1eta)
		self.lep1mass = array('f', [0.0])
		self.addBranch('lep1mass', self.lep1mass)
		self.lep1reliso = array('f', [0.0])
		self.addBranch('lep1reliso', self.lep1reliso)	
		self.lep2id = array('f', [-1.0])
		self.addBranch('lep2id', self.lep2id)
		self.lep2pt = array('f', [-1.0])
		self.addBranch('lep2pt', self.lep2pt)
		self.lep2phi = array('f', [100.0])
		self.addBranch('lep2phi', self.lep2phi)
		self.lep2eta = array('f', [100.0])
		self.addBranch('lep2eta', self.lep2eta)
		self.lep2mass = array('f', [0.0])
		self.addBranch('lep2mass', self.lep2mass)
		self.lep2reliso = array('f', [0.0])
		self.addBranch('lep2reliso', self.lep2reliso)	
		self.lep3id = array('f', [-1.0])
		self.addBranch('lep3id', self.lep3id)
		self.lep3pt = array('f', [-1.0])
		self.addBranch('lep3pt', self.lep3pt)
		self.lep3phi = array('f', [100.0])
		self.addBranch('lep3phi', self.lep3phi)
		self.lep3eta = array('f', [100.0])
		self.addBranch('lep3eta', self.lep3eta)
		self.lep3mass = array('f', [0.0])
		self.addBranch('lep3mass', self.lep3mass)
		self.lep3reliso = array('f', [0.0])
		self.addBranch('lep3reliso', self.lep3reliso)	
		self.lep4id = array('f', [-1.0])
		self.addBranch('lep4id', self.lep4id)
		self.lep4pt = array('f', [-1.0])
		self.addBranch('lep4pt', self.lep4pt)
		self.lep4phi = array('f', [100.0])
		self.addBranch('lep4phi', self.lep4phi)
		self.lep4eta = array('f', [100.0])
		self.addBranch('lep4eta', self.lep4eta)
		self.lep4mass = array('f', [0.0])
		self.addBranch('lep4mass', self.lep4mass)
		self.lep4reliso = array('f', [0.0])
		self.addBranch('lep4reliso', self.lep4reliso)	
		self.lep2D_rel = array('f', [-1.0])
		self.addBranch('lep2D_rel', self.lep2D_rel)
		self.lep2D_dr = array('f', [-1.0])
		self.addBranch('lep2D_dr', self.lep2D_dr)
		self.lep2D_mass = array('f', [-1.0])
		self.addBranch('lep2D_mass', self.lep2D_mass)
		self.leptri_lep = array('f', [-1.0])
		self.addBranch('leptri_lep', self.leptri_lep)
		self.leptri_jet = array('f', [-1.0])
		self.addBranch('leptri_jet', self.leptri_jet)
		self.WcandPpt = array('f', [-1.0])
                self.addBranch('WcandPpt', self.WcandPpt)
                self.WcandNpt = array('f', [-1.0])
                self.addBranch('WcandNpt', self.WcandNpt)
                self.WcandPphi = array('f', [100.0])
                self.addBranch('WcandPphi', self.WcandPphi)
                self.WcandNphi = array('f', [100.0])
                self.addBranch('WcandNphi', self.WcandNphi)
                self.WcandPeta = array('f', [100.0])
                self.addBranch('WcandPeta', self.WcandPeta)
                self.WcandNeta = array('f', [100.0])
                self.addBranch('WcandNeta', self.WcandNeta)
                self.WcandPm = array('f', [0.0])
                self.addBranch('WcandPm', self.WcandPm)
                self.WcandNm = array('f', [0.0])
                self.addBranch('WcandNm', self.WcandNm)
	def analyze(self, event):
		self.weight[0] = self.w
		# met
		event.getByLabel(self.metLab, self.metHan)
		if not (self.metHan.isValid()):
			return 1
		self.metpt[0] = self.metHan.product()[0].Pt()
		self.metphi[0] = self.metHan.product()[0].Phi()
		self.unfitmet = ROOT.TLorentzVector()
		self.unfitmet.SetPtEtaPhiM(self.metpt[0], 0, self.metphi[0], 0)
		# jets setup
		event.getByLabel (self.jetLab, self.jetHan)
	    	event.getByLabel (self.unjetLab, self.unjetHan)
		event.getByLabel (self.ak4jetLab, self.ak4jetHan)
		event.getByLabel (self.toptagjetLab, self.toptagjetHan)
		event.getByLabel (self.toptagnsubLab, self.toptagnsubHan)
		event.getByLabel (self.toptagminmassLab, self.toptagminmassHan)
		event.getByLabel (self.toptagWmassLab, self.toptagWmassHan)
		event.getByLabel (self.toptagtopmassLab, self.toptagtopmassHan)
	        event.getByLabel (self.t3Lab, self.t3Han)
        	self.Tau3  =  self.t3Han.product() 
	        event.getByLabel (self.t2Lab, self.t2Han)
        	self.Tau2  =  self.t2Han.product()
	        event.getByLabel (self.t1Lab, self.t1Han)
        	self.Tau1  =  self.t1Han.product()
		if not (self.jetHan.isValid() and self.unjetHan.isValid()) :
			return 1
		self.unjets = self.unjetHan.product()
		self.jets = self.jetHan.product()
		self.ak4jets = self.ak4jetHan.product()
		self.toptagjets = self.toptagjetHan.product()	
		self.toptagnsubs = self.toptagnsubHan.product()	
		self.toptagminmasses = self.toptagminmassHan.product()	
		self.toptagtopmasses = self.toptagtopmassHan.product()	
		self.toptagWmasses = self.toptagWmassHan.product()	
		self.numjets[0] = float(len(self.jets))
		self.numak4jets[0] = float(len(self.ak4jets))
		self.numtoptagjets[0] = float(len(self.toptagjets))	
		# ca8 jets
		if self.numjets[0] > 0:
			self.jet1 = ROOT.TLorentzVector()
			self.jet1.SetPtEtaPhiM(self.jets[0].Pt(), self.jets[0].Eta(), self.jets[0].Phi(), self.jets[0].M())
			self.jetlist.append(self.jet1)
			self.jet1pt[0] = self.jet1.Pt()
			self.jet1mass[0] = self.jet1.M()
			self.jet1phi[0] = self.jet1.Phi()
			self.jet1eta[0] = self.jet1.Eta()
			unjet_index_jet1 = MatchCol(self.unjets, self.jet1)
			if unjet_index_jet1 >= 0:
				self.jet1tau1[0] = self.Tau1[unjet_index_jet1]
				self.jet1tau2[0] = self.Tau2[unjet_index_jet1]
				self.jet1tau3[0] = self.Tau3[unjet_index_jet1]		
		if self.numjets[0] > 1:
			self.jet2 = ROOT.TLorentzVector()
			self.jet2.SetPtEtaPhiM(self.jets[1].Pt(), self.jets[1].Eta(), self.jets[1].Phi(), self.jets[1].M())
			self.jetlist.append(self.jet2)
			self.jet2pt[0] = self.jet2.Pt()
			self.jet2mass[0] = self.jet2.M()
			self.jet2phi[0] = self.jet2.Phi()
			self.jet2eta[0] = self.jet2.Eta()
			unjet_index_jet2 = MatchCol(self.unjets, self.jet2)
			if unjet_index_jet2 >= 0:
				self.jet2tau1[0] = self.Tau1[unjet_index_jet2]
				self.jet2tau2[0] = self.Tau2[unjet_index_jet2]
				self.jet2tau3[0] = self.Tau3[unjet_index_jet2]		
		if self.numjets[0] > 2:
			self.jet3 = ROOT.TLorentzVector()
			self.jet3.SetPtEtaPhiM(self.jets[2].Pt(), self.jets[2].Eta(), self.jets[2].Phi(), self.jets[2].M())
			self.jetlist.append(self.jet3)
			self.jet3pt[0] = self.jet3.Pt()
			self.jet3mass[0] = self.jet3.M()
			self.jet3phi[0] = self.jet3.Phi()
			self.jet3eta[0] = self.jet3.Eta()
			unjet_index_jet3 = MatchCol(self.unjets, self.jet3)
			if unjet_index_jet3 >= 0:
				self.jet3tau1[0] = self.Tau1[unjet_index_jet3]
				self.jet3tau2[0] = self.Tau2[unjet_index_jet3]
				self.jet3tau3[0] = self.Tau3[unjet_index_jet3]		
		if self.numjets[0] > 3:
			self.jet4 = ROOT.TLorentzVector()
			self.jet4.SetPtEtaPhiM(self.jets[3].Pt(), self.jets[3].Eta(), self.jets[3].Phi(), self.jets[3].M())
			self.jetlist.append(self.jet4)
			self.jet4pt[0] = self.jet4.Pt()
			self.jet4mass[0] = self.jet4.M()
			self.jet4phi[0] = self.jet4.Phi()
			self.jet4eta[0] = self.jet4.Eta()
			unjet_index_jet4 = MatchCol(self.unjets, self.jet4)
			if unjet_index_jet4 >= 0:
				self.jet4tau1[0] = self.Tau1[unjet_index_jet4]
				self.jet4tau2[0] = self.Tau2[unjet_index_jet4]
				self.jet4tau3[0] = self.Tau3[unjet_index_jet4]		
		if self.numjets[0] > 4:
			self.jet5 = ROOT.TLorentzVector()
			self.jet5.SetPtEtaPhiM(self.jets[4].Pt(), self.jets[4].Eta(), self.jets[4].Phi(), self.jets[4].M())
			self.jetlist.append(self.jet5)
			self.jet5pt[0] = self.jet5.Pt()
			self.jet5mass[0] = self.jet5.M()
			self.jet5phi[0] = self.jet5.Phi()
			self.jet5eta[0] = self.jet5.Eta()
			unjet_index_jet5 = MatchCol(self.unjets, self.jet5)
			if unjet_index_jet5 >= 0:
				self.jet5tau1[0] = self.Tau1[unjet_index_jet5]
				self.jet5tau2[0] = self.Tau2[unjet_index_jet5]
				self.jet5tau3[0] = self.Tau3[unjet_index_jet5]		
		if self.numjets[0] > 5:
			self.jet6 = ROOT.TLorentzVector()
			self.jet6.SetPtEtaPhiM(self.jets[5].Pt(), self.jets[5].Eta(), self.jets[5].Phi(), self.jets[5].M())
			self.jetlist.append(self.jet6)
			self.jet6pt[0] = self.jet6.Pt()
			self.jet6mass[0] = self.jet6.M()
			self.jet6phi[0] = self.jet6.Phi()
			self.jet6eta[0] = self.jet6.Eta()
			unjet_index_jet6 = MatchCol(self.unjets, self.jet6)
			if unjet_index_jet6 >= 0:
				self.jet6tau1[0] = self.Tau1[unjet_index_jet6]
				self.jet6tau2[0] = self.Tau2[unjet_index_jet6]
				self.jet6tau3[0] = self.Tau3[unjet_index_jet6]		
		if self.numjets[0] > 6:
			self.jet7 = ROOT.TLorentzVector()
			self.jet7.SetPtEtaPhiM(self.jets[6].Pt(), self.jets[6].Eta(), self.jets[6].Phi(), self.jets[6].M())
			self.jetlist.append(self.jet7)
			self.jet7pt[0] = self.jet7.Pt()
			self.jet7mass[0] = self.jet7.M()
			self.jet7phi[0] = self.jet7.Phi()
			self.jet7eta[0] = self.jet7.Eta()
			unjet_index_jet7 = MatchCol(self.unjets, self.jet7)
			if unjet_index_jet7 >= 0:
				self.jet7tau1[0] = self.Tau1[unjet_index_jet7]
				self.jet7tau2[0] = self.Tau2[unjet_index_jet7]
				self.jet7tau3[0] = self.Tau3[unjet_index_jet7]
		if self.numjets[0] > 7:
			self.jet8 = ROOT.TLorentzVector()
			self.jet8.SetPtEtaPhiM(self.jets[7].Pt(), self.jets[7].Eta(), self.jets[7].Phi(), self.jets[7].M())
			self.jetlist.append(self.jet8)
			self.jet8pt[0] = self.jet8.Pt()
			self.jet8mass[0] = self.jet8.M()
			self.jet8phi[0] = self.jet8.Phi()
			self.jet8eta[0] = self.jet8.Eta()
			unjet_index_jet8 = MatchCol(self.unjets, self.jet8)
			if unjet_index_jet8 >= 0:
				self.jet8tau1[0] = self.Tau1[unjet_index_jet8]
				self.jet8tau2[0] = self.Tau2[unjet_index_jet8]
				self.jet8tau3[0] = self.Tau3[unjet_index_jet8]
		if self.numjets[0] > 8:
			self.jet9 = ROOT.TLorentzVector()
			self.jet9.SetPtEtaPhiM(self.jets[8].Pt(), self.jets[8].Eta(), self.jets[8].Phi(), self.jets[8].M())
			self.jetlist.append(self.jet9)
			self.jet9pt[0] = self.jet9.Pt()
			self.jet9mass[0] = self.jet9.M()
			self.jet9phi[0] = self.jet9.Phi()
			self.jet9eta[0] = self.jet9.Eta()
			unjet_index_jet9 = MatchCol(self.unjets, self.jet9)
			if unjet_index_jet9 >= 0:
				self.jet9tau1[0] = self.Tau1[unjet_index_jet9]
				self.jet9tau2[0] = self.Tau2[unjet_index_jet9]
				self.jet9tau3[0] = self.Tau3[unjet_index_jet9]
		if self.numjets[0] > 9:
			self.jet10 = ROOT.TLorentzVector()
			self.jet10.SetPtEtaPhiM(self.jets[9].Pt(), self.jets[9].Eta(), self.jets[9].Phi(), self.jets[9].M())
			self.jetlist.append(self.jet10)
			self.jet10pt[0] = self.jet10.Pt()
			self.jet10mass[0] = self.jet10.M()
			self.jet10phi[0] = self.jet10.Phi()
			self.jet10eta[0] = self.jet10.Eta()
			unjet_index_jet10 = MatchCol(self.unjets, self.jet10)
			if unjet_index_jet10 >= 0:
				self.jet10tau1[0] = self.Tau1[unjet_index_jet10]
				self.jet10tau2[0] = self.Tau2[unjet_index_jet10]
				self.jet10tau3[0] = self.Tau3[unjet_index_jet10]
		# top tag
		if self.numtoptagjets[0] > 0:
			self.toptagjet1 = ROOT.TLorentzVector()
			self.toptagjet1.SetPtEtaPhiM(self.toptagjets[0].Pt(), self.toptagjets[0].Eta(), self.toptagjets[0].Phi(), self.toptagjets[0].M())
			self.toptagjet1pt[0] = self.toptagjet1.Pt()
			self.toptagjet1mass[0] = self.toptagjet1.M()
			self.toptagjet1phi[0] = self.toptagjet1.Phi()
			self.toptagjet1eta[0] = self.toptagjet1.Eta()
			self.toptagjet1minmass[0] = self.toptagminmasses[0]
			self.toptagjet1topmass[0] = self.toptagtopmasses[0]
			self.toptagjet1Wmass[0] = self.toptagWmasses[0]
			self.toptagjet1nsub[0] = self.toptagnsubs[0]
		if self.numtoptagjets[0] > 1:
			self.toptagjet2 = ROOT.TLorentzVector()
			self.toptagjet2.SetPtEtaPhiM(self.toptagjets[1].Pt(), self.toptagjets[1].Eta(), self.toptagjets[1].Phi(), self.toptagjets[1].M())
			self.toptagjet2pt[0] = self.toptagjet2.Pt()
			self.toptagjet2mass[0] = self.toptagjet2.M()
			self.toptagjet2phi[0] = self.toptagjet2.Phi()
			self.toptagjet2eta[0] = self.toptagjet2.Eta()
			self.toptagjet2minmass[0] = self.toptagminmasses[1]
			self.toptagjet2topmass[0] = self.toptagtopmasses[1]
			self.toptagjet2Wmass[0] = self.toptagWmasses[1]
			self.toptagjet2nsub[0] = self.toptagnsubs[1]
		if self.numtoptagjets[0] > 2:
			self.toptagjet3 = ROOT.TLorentzVector()
			self.toptagjet3.SetPtEtaPhiM(self.toptagjets[2].Pt(), self.toptagjets[2].Eta(), self.toptagjets[2].Phi(), self.toptagjets[2].M())
			self.toptagjet3pt[0] = self.toptagjet3.Pt()
			self.toptagjet3mass[0] = self.toptagjet3.M()
			self.toptagjet3phi[0] = self.toptagjet3.Phi()
			self.toptagjet3eta[0] = self.toptagjet3.Eta()
			self.toptagjet3minmass[0] = self.toptagminmasses[2]
			self.toptagjet3topmass[0] = self.toptagtopmasses[2]
			self.toptagjet3Wmass[0] = self.toptagWmasses[2]
			self.toptagjet3nsub[0] = self.toptagnsubs[2]
		if self.numtoptagjets[0] > 3:
			self.toptagjet4 = ROOT.TLorentzVector()
			self.toptagjet4.SetPtEtaPhiM(self.toptagjets[3].Pt(), self.toptagjets[3].Eta(), self.toptagjets[3].Phi(), self.toptagjets[3].M())
			self.toptagjet4pt[0] = self.toptagjet4.Pt()
			self.toptagjet4mass[0] = self.toptagjet4.M()
			self.toptagjet4phi[0] = self.toptagjet4.Phi()
			self.toptagjet4eta[0] = self.toptagjet4.Eta()
			self.toptagjet4minmass[0] = self.toptagminmasses[3]
			self.toptagjet4topmass[0] = self.toptagtopmasses[3]
			self.toptagjet4Wmass[0] = self.toptagWmasses[3]
			self.toptagjet4nsub[0] = self.toptagnsubs[3]
		if self.numtoptagjets[0] > 4:
			self.toptagjet5 = ROOT.TLorentzVector()
			self.toptagjet5.SetPtEtaPhiM(self.toptagjets[4].Pt(), self.toptagjets[4].Eta(), self.toptagjets[4].Phi(), self.toptagjets[4].M())
			self.toptagjet5pt[0] = self.toptagjet5.Pt()
			self.toptagjet5mass[0] = self.toptagjet5.M()
			self.toptagjet5phi[0] = self.toptagjet5.Phi()
			self.toptagjet5eta[0] = self.toptagjet5.Eta()
			self.toptagjet5minmass[0] = self.toptagminmasses[4]
			self.toptagjet5topmass[0] = self.toptagtopmasses[4]
			self.toptagjet5Wmass[0] = self.toptagWmasses[4]
			self.toptagjet5nsub[0] = self.toptagnsubs[4]
		if self.numtoptagjets[0] > 5:
			self.toptagjet6 = ROOT.TLorentzVector()
			self.toptagjet6.SetPtEtaPhiM(self.toptagjets[5].Pt(), self.toptagjets[5].Eta(), self.toptagjets[5].Phi(), self.toptagjets[5].M())
			self.toptagjet6pt[0] = self.toptagjet6.Pt()
			self.toptagjet6mass[0] = self.toptagjet6.M()
			self.toptagjet6phi[0] = self.toptagjet6.Phi()
			self.toptagjet6eta[0] = self.toptagjet6.Eta()
			self.toptagjet6minmass[0] = self.toptagminmasses[5]
			self.toptagjet6topmass[0] = self.toptagtopmasses[5]
			self.toptagjet6Wmass[0] = self.toptagWmasses[5]
			self.toptagjet6nsub[0] = self.toptagnsubs[5]
		if self.numtoptagjets[0] > 6:
			self.toptagjet7 = ROOT.TLorentzVector()
			self.toptagjet7.SetPtEtaPhiM(self.toptagjets[6].Pt(), self.toptagjets[6].Eta(), self.toptagjets[6].Phi(), self.toptagjets[6].M())
			self.toptagjet7pt[0] = self.toptagjet7.Pt()
			self.toptagjet7mass[0] = self.toptagjet7.M()
			self.toptagjet7phi[0] = self.toptagjet7.Phi()
			self.toptagjet7eta[0] = self.toptagjet7.Eta()
			self.toptagjet7minmass[0] = self.toptagminmasses[6]
			self.toptagjet7topmass[0] = self.toptagtopmasses[6]
			self.toptagjet7Wmass[0] = self.toptagWmasses[6]
			self.toptagjet7nsub[0] = self.toptagnsubs[6]
		if self.numtoptagjets[0] > 7:
			self.toptagjet8 = ROOT.TLorentzVector()
			self.toptagjet8.SetPtEtaPhiM(self.toptagjets[7].Pt(), self.toptagjets[7].Eta(), self.toptagjets[7].Phi(), self.toptagjets[7].M())
			self.toptagjet8pt[0] = self.toptagjet8.Pt()
			self.toptagjet8mass[0] = self.toptagjet8.M()
			self.toptagjet8phi[0] = self.toptagjet8.Phi()
			self.toptagjet8eta[0] = self.toptagjet8.Eta()
			self.toptagjet8minmass[0] = self.toptagminmasses[7]
			self.toptagjet8topmass[0] = self.toptagtopmasses[7]
			self.toptagjet8Wmass[0] = self.toptagWmasses[7]
			self.toptagjet8nsub[0] = self.toptagnsubs[7]
		if self.numtoptagjets[0] > 8:
			self.toptagjet9 = ROOT.TLorentzVector()
			self.toptagjet9.SetPtEtaPhiM(self.toptagjets[8].Pt(), self.toptagjets[8].Eta(), self.toptagjets[8].Phi(), self.toptagjets[8].M())
			self.toptagjet9pt[0] = self.toptagjet9.Pt()
			self.toptagjet9mass[0] = self.toptagjet9.M()
			self.toptagjet9phi[0] = self.toptagjet9.Phi()
			self.toptagjet9eta[0] = self.toptagjet9.Eta()
			self.toptagjet9minmass[0] = self.toptagminmasses[8]
			self.toptagjet9topmass[0] = self.toptagtopmasses[8]
			self.toptagjet9Wmass[0] = self.toptagWmasses[8]
			self.toptagjet9nsub[0] = self.toptagnsubs[8]
		if self.numtoptagjets[0] > 9:
			self.toptagjet10 = ROOT.TLorentzVector()
			self.toptagjet10.SetPtEtaPhiM(self.toptagjets[9].Pt(), self.toptagjets[9].Eta(), self.toptagjets[9].Phi(), self.toptagjets[9].M())
			self.toptagjet10pt[0] = self.toptagjet10.Pt()
			self.toptagjet10mass[0] = self.toptagjet10.M()
			self.toptagjet10phi[0] = self.toptagjet10.Phi()
			self.toptagjet10eta[0] = self.toptagjet10.Eta()
			self.toptagjet10minmass[0] = self.toptagminmasses[9]
			self.toptagjet10topmass[0] = self.toptagtopmasses[9]
			self.toptagjet10Wmass[0] = self.toptagWmasses[9]
			self.toptagjet10nsub[0] = self.toptagnsubs[9]
		
		# leptons
		event.getByLabel (self.muonLab, self.muonHan)
		event.getByLabel (self.muonRelIsoLab, self.muonRelIsoHan)
		event.getByLabel (self.muonIsTightLab, self.muonIsTightHan)
		event.getByLabel (self.elecLab, self.elecHan)
		event.getByLabel (self.elecRelIsoLab, self.elecRelIsoHan)
		if not (self.muonHan.isValid() and self.elecHan.isValid()):
			return 1
		self.muons = self.muonHan.product()
		self.muonsRelIso = self.muonRelIsoHan.product()
		self.muonsIsTight = self.muonIsTightHan.product()
		self.electrons = self.elecHan.product()
		self.electronsRelIso = self.elecRelIsoHan.product()
		# Classify event according to number of good leptons
		n_good_electrons = 0
		n_good_muons = 0
		goodleps = []
		goodleptons = []
		for i in range(0, len(self.electrons)):
			if self.electrons[i].Pt() > 35 and fabs(self.electrons[i].Eta()) < 2.5:
				n_good_electrons += 1
				goodlepton = ROOT.TLorentzVector()
				goodlepton.SetPtEtaPhiM(self.electrons[i].Pt(), self.electrons[i].Eta(), self.electrons[i].Phi(), self.electrons[i].M())
				goodleps.append( (goodlepton, 1, self.electronsRelIso[i]) )		
				goodleptons.append(goodlepton)		
		for i in range(0, len(self.muons)):
			if self.muons[i].Pt() > 45 and fabs(self.muons[i].Eta()) < 2.1 and self.muonsIsTight[i]:
				n_good_muons += 1
				goodlepton = ROOT.TLorentzVector()
				goodlepton.SetPtEtaPhiM(self.muons[i].Pt(), self.muons[i].Eta(), self.muons[i].Phi(), self.muons[i].M())
				goodleps.append( (goodlepton, 2, self.muonsRelIso[i]) )
				goodleptons.append(goodlepton)		
		self.numleptons[0] = self.muons.size() + self.electrons.size()
		self.numgoodleptons[0] = len(goodleps)
		Ht_lep = 0
		for lep in goodleps:
			Ht_lep += lep[0].Pt()
		self.Htlep[0] = Ht_lep
		# Keep four highest pt leptons
		pt1 = 0
		i1 = -1
		for i in range(len(goodleps)):
			pt = goodleps[i][0].Pt()
			if pt > pt1:
				pt1 = pt
				i1 = i
		if i1 != -1:
			self.lep1id[0] = goodleps[i1][1]
			self.lep1pt[0] = goodleps[i1][0].Pt()
			self.lep1eta[0] = goodleps[i1][0].Eta()
			self.lep1phi[0] = goodleps[i1][0].Phi()
			self.lep1mass[0] = goodleps[i1][0].M()
			self.lep1reliso[0] = goodleps[i1][2]	
		pt2 = 0
		i2 = -1
		for i in range(len(goodleps)):
			if i == i1:
				continue
			pt = goodleps[i][0].Pt()
			if pt > pt2:
				pt2 = pt
				i2 = i
		if i2 != -1:
			self.lep2id[0] = goodleps[i2][1]
			self.lep2pt[0] = goodleps[i2][0].Pt()
			self.lep2eta[0] = goodleps[i2][0].Eta()
			self.lep2phi[0] = goodleps[i2][0].Phi()
			self.lep2mass[0] = goodleps[i2][0].M()
			self.lep2reliso[0] = goodleps[i2][2]	
		pt3 = 0
		i3 = -1
		for i in range(len(goodleps)):
			if i == i1 or i == i2:
				continue
			pt = goodleps[i][0].Pt()
			if pt > pt3:
				pt3 = pt
				i3 = i
		if i3 != -1:
			self.lep3id[0] = goodleps[i3][1]
			self.lep3pt[0] = goodleps[i3][0].Pt()
			self.lep3eta[0] = goodleps[i3][0].Eta()
			self.lep3phi[0] = goodleps[i3][0].Phi()
			self.lep3mass[0] = goodleps[i3][0].M()
			self.lep3reliso[0] = goodleps[i3][2]	
		pt4 = 0
		i4 = -1
		for i in range(len(goodleps)):
			if i == i1 or i == i2 or i == i3:
				continue
			pt = goodleps[i][0].Pt()
			if pt > pt4:
				pt4 = pt
				i4 = i
		if i4 != -1:
			self.lep4id[0] = goodleps[i4][1]
			self.lep4pt[0] = goodleps[i4][0].Pt()
			self.lep4eta[0] = goodleps[i4][0].Eta()
			self.lep4phi[0] = goodleps[i4][0].Phi()
			self.lep4mass[0] = goodleps[i4][0].M()
			self.lep4reliso[0] = goodleps[i4][2]	
		self.leadlep = ROOT.TLorentzVector()
		if (n_good_electrons == 0 and n_good_muons == 0) : 
			self.eventType[0] = 0.0 # full hadronic
		if (n_good_electrons == 1 and n_good_muons == 0) : 
			self.eventType[0] = 1.0 # semileptonic with electron
			self.leadlep = goodleptons[0]
			self.doLepStuff()
		if (n_good_electrons == 0 and n_good_muons == 1) :
			self.eventType[0] = 2.0 # semileptonic with muon
			self.leadlep = goodleptons[0]
			self.doLepStuff()
		if (n_good_electrons + n_good_muons >= 2) :
			self.eventType[0] = 3.0 # full leptonic or multi-leptons
			leading_lep = goodleptons[0]
			for test_lep in goodleptons:
				if test_lep.Pt() > leading_lep.Pt():
					leading_lep = test_lep
			self.leadlep = leading_lep
			self.doLepStuff()
		# Reco
		Ht_event = 0.0
		for i in range(len(self.jets)):
			Ht_event += self.jets[i].Pt()
		self.Ht[0] = Ht_event
		# Fill tree
		self.tree.Fill()
		return 0

	def doLepStuff(self):
		# MET
		self.newmetvals= make_lepW(self.unfitmet, self.leadlep)
		self.metetaP[0] = self.newmetvals[0].Eta()
		self.metetaN[0] = self.newmetvals[1].Eta()
		# W reco
		self.WcandP = self.newmetvals[0] + self.leadlep
		self.WcandN = self.newmetvals[1] + self.leadlep
		self.WcandPpt[0] = self.WcandP.Pt()
		self.WcandNpt[0] = self.WcandN.Pt()
		self.WcandPphi[0] = self.WcandP.Phi()
		self.WcandNphi[0] = self.WcandN.Phi()
		self.WcandPeta[0] = self.WcandP.Eta()
		self.WcandNeta[0] = self.WcandN.Eta()
		self.WcandPm[0] = self.WcandP.M()
		self.WcandNm[0] = self.WcandN.M()
		# Selection variables on lepton
		jet25list = []
		for jet in self.jetlist:
			if jet.Pt() > 25:
				jet25list.append(jet)
		if len(jet25list) > 0:
			index = ClosestJet(jet25list, self.leadlep)
	                self.lep2D_rel[0] = self.leadlep.Perp(jet25list[index].Vect())
	                self.lep2D_dr[0] = self.leadlep.DeltaR(jet25list[index])
			self.lep2D_mass[0] = (self.leadlep + jet25list[index]).M()
	                self.leptri_lep[0] = abs(self.leadlep.Phi() - self.metphi[0])
	                self.leptri_jet[0] = abs(self.jets[0].Phi() - self.metphi[0])
		
	def reset(self, err):
		if err != 0:
			self.ErrHist.Fill(err)
		self.jetlist = []
		self.weight[0] = self.w
		# General
		self.numjets[0] = 0.0	
		self.numak4jets[0] = 0.0	
		self.numtoptagjets[0] = 0.0	
		self.eventType[0] = -1.0		
		self.Ht[0] = 0.0
		self.metpt[0] = -1.0
		self.metphi[0] = 100.0
		self.metetaP[0] = 100.0
		self.metetaN[0] = 100.0
		# lep
		self.Htlep[0] = 0
		self.numleptons[0] = 0
		self.numgoodleptons[0] = 0
		self.lep1id[0] = -1.0
		self.lep1pt[0] = -1.0
		self.lep1phi[0] = 100.0
		self.lep1eta[0] = 100.0
		self.lep1mass[0] = 100.0
		self.lep1reliso[0] = 100.0
		self.lep2id[0] = -1.0
		self.lep2pt[0] = -1.0
		self.lep2phi[0] = 100.0
		self.lep2eta[0] = 100.0
		self.lep2mass[0] = 100.0
		self.lep2reliso[0] = 100.0
		self.lep3id[0] = -1.0
		self.lep3pt[0] = -1.0
		self.lep3phi[0] = 100.0
		self.lep3eta[0] = 100.0
		self.lep3mass[0] = 100.0
		self.lep3reliso[0] = 100.0
		self.lep4id[0] = -1.0
		self.lep4pt[0] = -1.0
		self.lep4phi[0] = 100.0
		self.lep4eta[0] = 100.0
		self.lep4mass[0] = 100.0
		self.lep1reliso[0] = 100.0
		self.lep2D_rel[0] = -1.0
		self.lep2D_dr[0] = -1.0
		self.lep2D_mass[0] = -1.0
		self.leptri_lep[0] = -1.0
		self.leptri_jet[0] = -1.0
		# W
		self.WcandPpt[0] = -1.0
		self.WcandNpt[0] = -1.0
		self.WcandPphi[0] = 100.0
		self.WcandNphi[0] = 100.0
		self.WcandPeta[0] = 100.0
		self.WcandNeta[0] = 100.0
		self.WcandPm[0] = 0.0
		self.WcandNm[0] = 0.0
		# jet
		self.jet1mass[0] = -1.0
		self.jet1pt[0] = -1.0
		self.jet1eta[0] = 100.0
		self.jet1phi[0] = 100.0
		self.jet1csv[0] = 0.0
		self.jet1tau1[0] = 1.0
		self.jet1tau2[0] = 1.0
		self.jet1tau3[0] = 1.0
		self.jet1tau4[0] = 1.0
		self.jet2mass[0] = -1.0
		self.jet2pt[0] = -1.0
		self.jet2eta[0] = 100.0
		self.jet2phi[0] = 100.0
		self.jet2csv[0] = 0.0
		self.jet2tau1[0] = 1.0
		self.jet2tau2[0] = 1.0
		self.jet2tau3[0] = 1.0
		self.jet2tau4[0] = 1.0
		self.jet3mass[0] = -1.0
		self.jet3pt[0] = -1.0
		self.jet3eta[0] = 100.0
		self.jet3phi[0] = 100.0
		self.jet3csv[0] = 0.0
		self.jet3tau1[0] = 1.0
		self.jet3tau2[0] = 1.0
		self.jet3tau3[0] = 1.0
		self.jet3tau4[0] = 1.0
		self.jet4mass[0] = -1.0
		self.jet4pt[0] = -1.0
		self.jet4eta[0] = 100.0
		self.jet4phi[0] = 100.0
		self.jet4csv[0] = 0.0
		self.jet4tau1[0] = 1.0
		self.jet4tau2[0] = 1.0
		self.jet4tau3[0] = 1.0
		self.jet4tau4[0] = 1.0
		self.jet5mass[0] = -1.0
		self.jet5pt[0] = -1.0
		self.jet5eta[0] = 100.0
		self.jet5phi[0] = 100.0
		self.jet5csv[0] = 0.0
		self.jet5tau1[0] = 1.0
		self.jet5tau2[0] = 1.0
		self.jet5tau3[0] = 1.0
		self.jet5tau4[0] = 1.0
		self.jet6mass[0] = -1.0
		self.jet6pt[0] = -1.0
		self.jet6eta[0] = 100.0
		self.jet6phi[0] = 100.0
		self.jet6csv[0] = 0.0
		self.jet6tau1[0] = 1.0
		self.jet6tau2[0] = 1.0
		self.jet6tau3[0] = 1.0
		self.jet6tau4[0] = 1.0
		self.jet7mass[0] = -1.0
		self.jet7pt[0] = -1.0
		self.jet7eta[0] = 100.0
		self.jet7phi[0] = 100.0
		self.jet7csv[0] = 0.0
		self.jet7tau1[0] = 1.0
		self.jet7tau2[0] = 1.0
		self.jet7tau3[0] = 1.0
		self.jet7tau4[0] = 1.0
		self.jet8mass[0] = -1.0
		self.jet8pt[0] = -1.0
		self.jet8eta[0] = 100.0
		self.jet8phi[0] = 100.0
		self.jet8csv[0] = 0.0
		self.jet8tau1[0] = 1.0
		self.jet8tau2[0] = 1.0
		self.jet8tau3[0] = 1.0
		self.jet8tau4[0] = 1.0
		self.jet9mass[0] = -1.0
		self.jet9pt[0] = -1.0
		self.jet9eta[0] = 100.0
		self.jet9phi[0] = 100.0
		self.jet9csv[0] = 0.0
		self.jet9tau1[0] = 1.0
		self.jet9tau2[0] = 1.0
		self.jet9tau3[0] = 1.0
		self.jet9tau4[0] = 1.0
		self.jet10mass[0] = -1.0
		self.jet10pt[0] = -1.0
		self.jet10eta[0] = 100.0
		self.jet10phi[0] = 100.0
		self.jet10csv[0] = 0.0
		self.jet10tau1[0] = 1.0
		self.jet10tau2[0] = 1.0
		self.jet10tau3[0] = 1.0
		self.jet10tau4[0] = 1.0
		# toptag
		self.toptagjet1mass[0] = -1.0
		self.toptagjet1pt[0] = -1.0
		self.toptagjet1phi[0] = 100.0
		self.toptagjet1eta[0] = 100.0
		self.toptagjet1minmass[0] = -1.0
		self.toptagjet1topmass[0] = -1.0
		self.toptagjet1Wmass[0] = -1.0
		self.toptagjet1nsub[0] = -1.0
		self.toptagjet2mass[0] = -1.0
		self.toptagjet2pt[0] = -1.0
		self.toptagjet2phi[0] = 100.0
		self.toptagjet2eta[0] = 100.0
		self.toptagjet2minmass[0] = -1.0
		self.toptagjet2topmass[0] = -1.0
		self.toptagjet2Wmass[0] = -1.0
		self.toptagjet2nsub[0] = -1.0
		self.toptagjet3mass[0] = -1.0
		self.toptagjet3pt[0] = -1.0
		self.toptagjet3phi[0] = 100.0
		self.toptagjet3eta[0] = 100.0
		self.toptagjet3minmass[0] = -1.0
		self.toptagjet3topmass[0] = -1.0
		self.toptagjet3Wmass[0] = -1.0
		self.toptagjet3nsub[0] = -1.0
		self.toptagjet4mass[0] = -1.0
		self.toptagjet4pt[0] = -1.0
		self.toptagjet4phi[0] = 100.0
		self.toptagjet4eta[0] = 100.0
		self.toptagjet4minmass[0] = -1.0
		self.toptagjet4topmass[0] = -1.0
		self.toptagjet4Wmass[0] = -1.0
		self.toptagjet4nsub[0] = -1.0
		self.toptagjet5mass[0] = -1.0
		self.toptagjet5pt[0] = -1.0
		self.toptagjet5phi[0] = 100.0
		self.toptagjet5eta[0] = 100.0
		self.toptagjet5minmass[0] = -1.0
		self.toptagjet5topmass[0] = -1.0
		self.toptagjet5Wmass[0] = -1.0
		self.toptagjet5nsub[0] = -1.0
		self.toptagjet6mass[0] = -1.0
		self.toptagjet6pt[0] = -1.0
		self.toptagjet6phi[0] = 100.0
		self.toptagjet6eta[0] = 100.0
		self.toptagjet6minmass[0] = -1.0
		self.toptagjet6topmass[0] = -1.0
		self.toptagjet6Wmass[0] = -1.0
		self.toptagjet6nsub[0] = -1.0
		self.toptagjet7mass[0] = -1.0
		self.toptagjet7pt[0] = -1.0
		self.toptagjet7phi[0] = 100.0
		self.toptagjet7eta[0] = 100.0
		self.toptagjet7minmass[0] = -1.0
		self.toptagjet7topmass[0] = -1.0
		self.toptagjet7Wmass[0] = -1.0
		self.toptagjet7nsub[0] = -1.0
		self.toptagjet8mass[0] = -1.0
		self.toptagjet8pt[0] = -1.0
		self.toptagjet8phi[0] = 100.0
		self.toptagjet8eta[0] = 100.0
		self.toptagjet8minmass[0] = -1.0
		self.toptagjet8topmass[0] = -1.0
		self.toptagjet8Wmass[0] = -1.0
		self.toptagjet8nsub[0] = -1.0
		self.toptagjet9mass[0] = -1.0
		self.toptagjet9pt[0] = -1.0
		self.toptagjet9phi[0] = 100.0
		self.toptagjet9eta[0] = 100.0
		self.toptagjet9minmass[0] = -1.0
		self.toptagjet9topmass[0] = -1.0
		self.toptagjet9Wmass[0] = -1.0
		self.toptagjet9nsub[0] = -1.0
		self.toptagjet10mass[0] = -1.0
		self.toptagjet10pt[0] = -1.0
		self.toptagjet10phi[0] = 100.0
		self.toptagjet10eta[0] = 100.0
		self.toptagjet10minmass[0] = -1.0
		self.toptagjet10topmass[0] = -1.0
		self.toptagjet10Wmass[0] = -1.0
		self.toptagjet10nsub[0] = -1.0

	def addBranch(self, name, var):
		self.tree.Branch(name, var, name+'/F')

	def __del__(self):
	        self.f.cd()
	        self.f.Write()
	        self.f.Close()




