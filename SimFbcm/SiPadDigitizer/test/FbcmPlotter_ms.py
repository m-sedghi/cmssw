#! /usr/bin/env python

import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(True)

import os 
import sys
import argparse 
import math

#Run sammple :
#./FbcmPlotter_ms -s 0 -p 100


def puCase(x):
    return {
        '0p5': '0p5/210216_122628',
        '1': '1/210216_122729',
		'1p5': '1p5/210216_122817',
		'10': '10/210216_122528',
		'50': '50/210216_122429',
		'100': '100/210216_122326',
		'140': '140/210216_122237',
		'200': '200/210216_122146',
    }.get(x, '100/210216_122326')

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument( '-i' , '--infile' , dest='infile' , default='out_pu' , help='the name of the input file' , type=str )
    parser.add_argument( '-s' , '--SensorGroupIndex' , dest='SensorGroupIndex' , default=0 , help='the sensor group index' , type=int , choices=[i for i in range(8)])
    parser.add_argument( '-o' , '--outfile' , dest='outfile' , default='auto' , help='the name of the output file' , type=str )
    parser.add_argument( '-p' , '--pileup' , dest='pileup' , default=None , help='enter the pileup' , type=str )
    parser.add_argument( '-n' , '--number' , dest='num' , default='1' , help='enter the file number' , type=str )

    opt = parser.parse_args()
    print(puCase('1p5'))

#    if ((not opt.infile) and (not opt.pileup)):
        #print('please specify the input file name using -i option or the pileup using -p')
		
    if 	not opt.pileup :
		print('please specify the input file pileup name using -p option')
		return 1
    
    inFileName = opt.infile + opt.pileup + "_" + opt.num + ".root"
    opt.infile=inFileName
    #print(opt.infile)
    if opt.outfile == 'auto':
        opt.outfile = "{0}_{1}.root".format(opt.infile.split('.')[0] , opt.SensorGroupIndex)
		

    dir_infile = "/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/FBCM/NuGun/crab_FBCMnTuplePU{0}/0000/".format(puCase(opt.pileup))
    print(dir_infile+opt.infile)
	
    fIn = ROOT.TFile.Open( dir_infile + opt.infile)

    nEvents = fIn.Get("FbcmNtuple/hNEvents").GetBinContent( 1 )
    print(nEvents)
    #print(type(opt.SensorGroupIndex))
    geometryInfo = fIn.Get('FbcmNtuple/GeometryInfo')

    sensorsPerRho = ROOT.TH1I("hSensorsPerRho" , "Sensors per row" , 30 , 0 , 30 )
    geometryInfo.Draw( "SensorRho >> hSensorsPerRho" , "SensorGroupIndex=={0}".format( opt.SensorGroupIndex ) )

    totalAreaPerRho = ROOT.TH1D("hTotalAreaPerRho" , "" , 30 , 0 , 30 )
    geometryInfo.Draw("SensorRho >> hTotalAreaPerRho" , "SensorX*SensorY*(SensorGroupIndex=={0})".format( opt.SensorGroupIndex ) )

    totalEffectiveAreaPerRho = totalAreaPerRho.Clone("totalEffectiveAreaPerRho")
    totalEffectiveAreaPerRho.Scale( 3.0*nEvents)

    allHits = fIn.Get('FbcmNtuple/PU{0}'.format( opt.pileup ) )
    simhitsPerRho = ROOT.TH1D("hSimHitsPerRho" , "Sim Hits per rho" , 30 , 0 , 30 )
    digiHitsPerRho = ROOT.TH1D("hDigiHitsPerRho" , "Digi Hits per rho" , 30 , 0 , 30 )
	
    SensorArea = 0.0
    nNonZeros = ROOT.TH1D("hNonZeroPerRho" , "number of nonZeros per row" , 30 , 0 , 30 ) 
		
    for hit in allHits: 
		if hit.SensorGroupIndex== int(opt.SensorGroupIndex):
			simhitsPerRho.Fill( hit.SensorRho, hit.nSimParticles )
			digiHitsPerRho.Fill(hit.SensorRho , hit.nValidDigiToAs )
			SensorArea = hit.SensorArea
			for bx in range(3):
				if  hit.DigiHitStatus[bx]==-1:
					totalEffectiveAreaPerRho.Fill( hit.SensorRho , -hit.SensorArea )
				elif hit.DigiHitStatus[bx]==1:
					nNonZeros.Fill(hit.SensorRho)
    
    simhitsPerRho.Divide( totalAreaPerRho )
    simhitsPerRho.Scale( 1.0/7.0/nEvents )
    Effective_M = totalEffectiveAreaPerRho.Clone("Effective_M")
    Effective_M.Scale(1.0/SensorArea)
    nZeros=Effective_M.Clone("nZeros")
    nZeros.Add(nNonZeros,-1)
    nZeros.Divide( Effective_M )
    for b in range( nZeros.GetNbinsX() ) :
		bc = nZeros.GetBinContent(b+1)
		if bc > 0:
			nZeros.SetBinContent( b + 1 , -math.log( bc ) )
			
    ZeroCountingPerCm=nZeros.Clone("ZeroCont")
    ZeroCountingPerCm.Scale(1.0/SensorArea)
    totalEffectiveAreaPerRho.Scale( 1.0/nEvents/3.0 )
    digiHitsPerRho2 = digiHitsPerRho.Clone("digiHitsPerRho2")
    digiHitsPerRho.Divide( totalEffectiveAreaPerRho )
    digiHitsPerRho.Scale( 1.0/3.0/nEvents )
    digiHitsPerRho2.Divide(totalAreaPerRho)
    digiHitsPerRho2.Scale( 1.0/3.0/nEvents )
    fout = ROOT.TFile.Open( opt.outfile , "recreate")
    simhitsPerRho.Write()
    digiHitsPerRho.Write()
    digiHitsPerRho2.Write()
    Effective_M.Write()
    totalEffectiveAreaPerRho.Write()
    totalAreaPerRho.Write()
    nZeros.Write()
    ZeroCountingPerCm.Write()
    fout.Close()
    return 0

if __name__ == "__main__":
    sys.exit( main() )
