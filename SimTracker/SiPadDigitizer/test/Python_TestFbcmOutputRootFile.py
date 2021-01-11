#from ROOT import * 
import ROOT
import os
import sys
import array
import math

from DataFormats.FWLite import Events, Handle

#ROOT.gROOT.SetBatch(True)

ROOT.gSystem.Load("libFWCoreFWLite")
ROOT.gSystem.Load("libDataFormatsFbcmDigi")
ROOT.FWLiteEnabler.enable()

#edm::DetSetVector<SiPadDigiData>      "simFbcmDigis"              "SiPad"  

events = Events ("../../../p1/p2/GEN_SIM_DIGI.root")
handle  = Handle ("edm::DetSetVector<SiPadDigiData>")
label = ("simFbcmDigis","SiPad")
outputfileName = 'CheckEventNo_'
ii=1
for event in events:
	event.getByLabel (label, handle)
	fbcmDets = handle.product()
	fileName = outputfileName + str(ii) + '.txt'
	f= open(fileName,"w")
	for Det in fbcmDets:
		for Digidata in Det:
			Side=Digidata.SideIndex()
			Station=Digidata.StationIndex()
			Die=Digidata.SiliconDieIndex() 
			SiPadNo=Digidata.SiPadIndex()
			Radius=Digidata.Radius() 
			phi=Digidata.Phi_Degrees() 
			Area_=Digidata.Area() 
			ampl=Digidata.ChargeSum()
			BVect = Digidata.BxSlotHitAnalysisVector()
			Ch_P = Digidata.CahrgePsimVector()
			str1 = '\nSide:' + str(Side) + ', Station:' + str(Station) +', Die:' + str(Die) +', SiPadNo:' + str(SiPadNo) +', Radius:' +  '%.2f' % (Radius) +', phi:' +  '%.2f' % (phi) +', Area:' +  '%.4f' % (Area_) +', ChargeSum:' +  '%.1f' % (ampl) + ', BVect_Size:'+str(BVect.size())+ ', Ch_P_Size:'+str(Ch_P.size()) +'\n'
			f.write(str1)
			for chp in Ch_P:
				charge = chp.first
				PSim = chp.second
				evntID = PSim.eventId()
				str2 = '\tCharge:'+ '%8.1f' % (charge) +', evntID:' + str(evntID.event()) + ', EventBx:' + str(evntID.bunchCrossing()) +', EventRawId:' + str(evntID.rawId()) +', trackId:' + str(PSim.trackId()) +', hitIndex:' + str(PSim.hitIndex()) +', tofBin:' + str(PSim.tofBin()) +', time:' + '%.3f' % (PSim.time()) + ', Tof:'+ '%.3f' % (PSim.Tof()) + ', Pabs:'+ '%.3e' % (PSim.Pabs()) +', EnergyLoss:' + '%.3e' % (PSim.EnergyLoss()) + ', ParticleType:'+str(PSim.ParticleType()) + ', DetUnitId:'+str(PSim.DetUnitId()) + ', ProcessType:'+str(PSim.ProcessType()) + ', Bx:'+str(PSim.BunchCrossing()) + '\n'
				f.write(str2)

			for hitA in BVect:
				ToTToAVect=hitA.TotToaVectort()
				str3 = '\t\tBxSlotNo:'+ str(hitA.BxSlotNo()) +', AlignerDelay:' +  '%.1f' % (hitA.lpGBTAlignerDelay()) + ', Bx_HitStatus:' + str(hitA.Bx_HitStatusInt()) +', nbrOfRecognizedHitsInBx:' + str(hitA.nbrOfRecognizedHitsInBx()) + ', ToTToAVect_size:' + str(ToTToAVect.size()) +'\n'
				f.write(str3)
				for toa in ToTToAVect: 
					str4 = '\t\t\tToA:'+  '%.1f' % (toa.ToA()) +', ToT:' +  '%.1f' % (toa.ToT()) + ', ToAStatus:' + str(toa.ToAToAStatusInt()) +', SubBxBinNumber:' + str(toa.SubBxBinNumber()) + ', IsToTValid:' + str(toa.IsToTValid()) + '\n'
					f.write(str4)

			str5 = '=======================================================================================================\n'
			f.write(str5)
				
	ii=ii+1
	f.close()
ii=ii-1
print 'Done: '+str(ii) + ' files were written for '+ str(ii) +' events.'
