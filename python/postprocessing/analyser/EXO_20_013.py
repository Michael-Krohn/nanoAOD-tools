from PhysicsTools.NanoAODTools.postprocessing.tools import *
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
import os
import sys
import ROOT
import argparse
import math
ROOT.PyConfig.IgnoreCommandLineOptions = True


class EXO_20_013(Module):
    def __init__(self):
        self.writeHistFile = True

    def beginJob(self, histFile=None, histDirName=None):
        Module.beginJob(self, histFile, histDirName)
	# define histograms
	self.addObject(ROOT.TH1F('Cutflow','Cutflow',6,0,6))
	self.addObject(ROOT.TH1F('DeltaR_ll','DeltaR_ll',100,0,6))
	self.addObject(ROOT.TH1F('m_ll','m_ll',100,12,1012))
	self.addObject(ROOT.TH1F('mT','mT',100,0,1000))
	self.addObject(ROOT.TH1F('pT_l1','pT_l1',100,0,150))
	self.addObject(ROOT.TH1F('pT_l2','pT_l2',100,0,150))
	self.addObject(ROOT.TH1F('eta_l1','eta_l1',50,-3,3))
	self.addObject(ROOT.TH1F('eta_l2','eta_l2',50,-3,3))
	self.addObject(ROOT.TH1F('pT_ll','pT_ll',150,0,300))
	self.addObject(ROOT.TH1F('met','met',150,0,300))

#    def endJob(self):
#        pass

#    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
#        self.out = wrappedOutputTree
#        self.out.branch("EventMass", "F")

#    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
#        pass

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        rawElectrons = Collection(event, "Electron")
        rawMuons = Collection(event, "Muon")
        rawJets = Collection(event, "Jet")
	met = Object(event, "MET") 

	electrons = []; muons = []; jets = []
	leptons = []

	self.Cutflow.Fill(0)

        eventSum = ROOT.TLorentzVector()
        for lep in rawMuons:
	    if lep.tightId and lep.pfIsoId > 3 and abs(lep.dxy) < 0.02 and abs(lep.dz) < 0.1:
		muons.append(lep)
		leptons.append(lep)
        for lep in rawElectrons:
	    if lep.mvaFall17V2Iso_WP90:
		electrons.append(lep)
		leptons.append(lep)

	#pT sort all objects
	electrons.sort(key=lambda x : x.p4().Pt()) ; muons.sort(key=lambda x : x.p4().Pt()) ; leptons.sort(key=lambda x : x.p4().Pt())
	if len(electrons) + len(muons) < 2: return False

        for j in rawJets:
	    if j.jetId is 4:
		if deltaR(leptons[0],j) > 0.4 and deltaR(leptons[1],j) > 0.4:
		    jets.append(j)

	jets.sort(key=lambda x : x.p4().Pt())
	# EVENT SELECTION

	# Preselection:
	# - At least 2 leptons, opposite flavor, opposite sign
	# - pT(l1)/pT(l2) > 25/20 GeV
	# - Veto additional loose leptons w/ pT(l3) > 10 GeV
	# - pT(ll) > 30 GeV
	# - m(ll) > 12 GeV
	# - MET > 20 GeV

	if len(electrons) + len(muons) < 2: return False
	if len(electrons) < 1 and len(muons) < 1: return False
	if leptons[0].charge == leptons[1].charge: return False
	if leptons[0].pt < 25: return False
	if leptons[1].pt < 20: return False
	self.Cutflow.Fill(1)
	if len(leptons) > 2:
	    if leptons[2].pt > 10: return False
	self.Cutflow.Fill(2)

	if (leptons[0].p4() + leptons[1].p4()).Pt() < 30: return False
	self.Cutflow.Fill(3)
	if (leptons[0].p4() + leptons[1].p4()).M() < 12: return False
	self.Cutflow.Fill(4)
	
	if met.pt < 20: return False
	self.Cutflow.Fill(5)

	mT = math.sqrt(2*leptons[1].pt*met.pt*(1-math.cos(deltaPhi(leptons[1].phi,met.phi))))
	self.DeltaR_ll.Fill(deltaR(leptons[0],leptons[1]))
	self.m_ll.Fill((leptons[0].p4() + leptons[1].p4()).M())
	self.mT.Fill(mT)
	self.pT_l1.Fill(leptons[0].pt)
	self.pT_l2.Fill(leptons[1].pt)
	self.eta_l1.Fill(leptons[0].eta)
	self.eta_l2.Fill(leptons[1].eta)
	self.met.Fill(met.pt)
	self.pT_ll.Fill((leptons[0].p4() + leptons[1].p4()).Pt())

        return True


# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed

exampleModuleConstr = lambda: exampleProducer(jetSelection=lambda j: j.pt > 30)

parser = argparse.ArgumentParser(description='FastSim EXO-20-013 Analysis Validation Tool')
parser.add_argument('-d', dest='directory',default="")
parser.add_argument('-o', dest='output',default="histOut.root")
args = parser.parse_args()

directory = "" ; outputname = "" ; preselection = ""

if args.directory == "" :
	print("Directory should be given")
	quit()

if args.directory != "" :
	directory = args.directory

outputname = args.output 
files = [] 

for filename in os.listdir(directory) :
	if filename.endswith(".root") :
		files.append(os.path.join(directory, filename))
		print(os.path.join(directory,filename))
	else : continue


p = PostProcessor(".", files, cut=preselection, branchsel=None, modules=[
                  EXO_20_013()], noOut=True, histFileName=outputname, histDirName="plots")
p.run()
