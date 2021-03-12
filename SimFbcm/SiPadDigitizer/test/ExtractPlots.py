#! /usr/bin/env python

import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(True)

import os 
import sys
import argparse 
import math

class SensorGroupInformation:
    def CopyH1IToH1F(self , h ):
        hret = ROOT.TH1F()
        h.Copy( hret )
        return hret

    def __init__(self , sensorGroup , pu , fIn):
        self.SensorGroup = sensorGroup
        self.PU = pu

        self.InDir = fIn.GetDirectory( 'SensorGroup{0}'.format( self.SensorGroup ) )
        self.nSimHits = self.CopyH1IToH1F( self.InDir.Get( 'h{0}_SG{1}'.format( 'nSimHits' , self.SensorGroup ) ) )
        self.nDigiHits = self.CopyH1IToH1F( self.InDir.Get( 'h{0}_SG{1}'.format( 'nDigiHits' , self.SensorGroup ) ) )
        self.nOnes = self.CopyH1IToH1F( self.InDir.Get( 'h{0}_SG{1}'.format( 'nOnes' , self.SensorGroup ) ) )
        self.nUnknowns = self.CopyH1IToH1F( self.InDir.Get( 'h{0}_SG{1}'.format( 'nUnknowns' , self.SensorGroup ) ) )
        self.nTotals = self.CopyH1IToH1F( self.InDir.Get( 'h{0}_SG{1}'.format( 'TotalNumbers' , self.SensorGroup ) ) )

        self.nSimHits.Divide( self.nTotals )
        self.nSimHits.Scale( 1.0/7.0 )
        self.nDigiHits.Divide( self.nTotals )
        self.nDigiHits.Scale( 1.0 / 3.0 )
        
        self.nTotals.Scale( 3.0 )
        self.nZeros = self.nTotals.Clone("nZeros_SG{0}".format(self.SensorGroup))
        self.nZeros.Add( self.nOnes , -1.0 )
        self.nZeros.Add( self.nUnknowns , -1.0 )
        self.nTotals.Add( self.nUnknowns , -1.0 )
        self.nZeros.Divide( self.nTotals )
        for b in range( self.nZeros.GetNbinsX() ) :
            bc = self.nZeros.GetBinContent(b+1)
            if bc > 0:
                self.nZeros.SetBinContent( b + 1 , -math.log( bc ) )

        

    def Write(self ):
        self.c1 = ROOT.TCanvas('cPU{0}SG{1}'.format( self.PU , self.SensorGroup ), 'PU{0}SG{1}'.format( self.PU , self.SensorGroup ) )

        self.nZeros.SetLineColor( 3 )
        self.nZeros.SetTitle( '<Digi>-ZeroCounting' )
        self.nZeros.SetStats( False )
        self.nZeros.Draw(  )

        self.nDigiHits.SetLineColor( 2 )
        self.nDigiHits.SetTitle( '<Digi>-counting' )
        self.nDigiHits.SetStats(False)
        self.nDigiHits.Draw('same')

        self.nSimHits.SetLineColor( 4 )
        self.nSimHits.SetTitle('<sim-hits>')
        self.nSimHits.SetStats(False)
        self.nSimHits.Draw("SAME")

        self.c1.BuildLegend()
        self.c1.SaveAs( '{0}.png'.format( self.c1.GetName() ) )
        

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument( '-i' , '--infile' , dest='infile' , help='the name of the input file' , type=str  )
    parser.add_argument( '-p' , '--pu' , dest='PU' , default=None , help='the pu of the file' , type=str , choices=['0p5', '1' , '1p5' , '10' , '50' , '100' , '140' , '200'])

    opt = parser.parse_args()

    if not opt.infile:
        print('please specify the input file name using -i option')
        return 1


    ROOT.gStyle.SetOptTitle(0)

    fIn = ROOT.TFile.Open( opt.infile ) 
    for sg in range( 8 ):
        s = SensorGroupInformation(sg , opt.PU , fIn )
        s.Write()
        

    return 0

if __name__ == "__main__":
    sys.exit( main() )
