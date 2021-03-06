#! /usr/bin/env python

import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(True)

from ROOT import TCanvas, TGraph

import os 
import sys
import argparse 
import math
import copy

# import matplotlib
# import matplotlib.pyplot as plt
# from matplotlib.colors import BoundaryNorm
# from matplotlib.ticker import MaxNLocator
import numpy as np


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
        self.nDigiHitsV2 = self.CopyH1IToH1F( self.InDir.Get( 'h{0}_SG{1}'.format( 'nDigiHitsV2' , self.SensorGroup ) ) )
        self.nOnes = self.CopyH1IToH1F( self.InDir.Get( 'h{0}_SG{1}'.format( 'nOnes' , self.SensorGroup ) ) )
        self.nUnknowns = self.CopyH1IToH1F( self.InDir.Get( 'h{0}_SG{1}'.format( 'nUnknowns' , self.SensorGroup ) ) )
        self.nTotals = self.CopyH1IToH1F( self.InDir.Get( 'h{0}_SG{1}'.format( 'TotalNumbers' , self.SensorGroup ) ) ) 
        self.TotRho = self.InDir.Get( 'h{0}_SG{1}'.format( 'TotRho' , self.SensorGroup ) )
        self.ToaRho = self.InDir.Get( 'h{0}_SG{1}'.format( 'ToaRho' , self.SensorGroup ) )
        self.TofRho = self.InDir.Get( 'h{0}_SG{1}'.format( 'TofRho' , self.SensorGroup ) )
        self.BxTofRho = self.InDir.Get( 'h{0}_SG{1}'.format( 'BxTofRho' , self.SensorGroup ) ) 
        
        self.numOfBins = self.nSimHits.GetSize() # this includes Underflow and Overflow, so it is binNum+2
        RhoAxisVect = self.nSimHits.GetXaxis()
        RhoXBins = RhoAxisVect.GetXbins()
        #.GetArray()
        #self.RhoBinVect = self.nSimHits.GetXaxis.GetXbins()
        
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

        self.nUnknownsRatio = self.nUnknowns.Clone("nUnknownsRatio_SG{0}".format(self.SensorGroup))
        self.nUnknownsRatio.Divide( self.nTotals )

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
        self.c1.SaveAs( './results/Figures/{0}.png'.format( self.c1.GetName() ) )
        
    def hZeroCounting(self ):
        return self.nZeros
        
    def hCounting(self ):
        #return self.nSimHits
        return self.nDigiHits
    
    def hSimHits(self ):
        return self.nSimHits
    
    def hUnknownRatio(self ):
        return self.nUnknownsRatio
        
    def WriteToA(self, Rho):
        self.c2 = ROOT.TCanvas('cPU{0}SG{1}'.format( self.PU , self.SensorGroup ), 'PU{0}SG{1}'.format( self.PU , self.SensorGroup ) )
        #Tpadd=ROOT.TPad("H1","H2",)
        b = self.nSimHits.FindBin(Rho)
        y_proj = self.ToaRho.ProjectionY("py",b)
        y_proj.SetName(self.ToaRho.GetName()+"_projY{0}".format(self.SensorGroup))
        y_proj.Draw("same, hist")
        ToA_1FH = ROOT.TH1F()
        y_proj.Copy(ToA_1FH)
        ToA_1FH.SetTitle( 'BIB_ToA' )
        ToA_1FH.SetStats( False )
        #ToA_1FH.GetXaxis().SetRangeUser(-15.,15.)
        ToA_1FH.Rebin(2)
        ToA_1FH.Draw()
        self.c2.BuildLegend()
        self.c2.SaveAs( './{0}.png'.format( self.c2.GetName() ) )
        
        print("Rho is bin No: {0} ".format(b) )
        #return 0
        
        # self.nZeros.SetLineColor( 3 )
        # self.nZeros.SetTitle( '<Digi>-ZeroCounting' )
        # self.nZeros.SetStats( False )
        # self.nZeros.Draw(  )

        # self.nDigiHits.SetLineColor( 2 )
        # self.nDigiHits.SetTitle( '<Digi>-counting' )
        # self.nDigiHits.SetStats(False)
        # self.nDigiHits.Draw('same')

        # self.nSimHits.SetLineColor( 4 )
        # self.nSimHits.SetTitle('<sim-hits>')
        # self.nSimHits.SetStats(False)
        # self.nSimHits.Draw("SAME")

        # self.c1.BuildLegend()
        # self.c1.SaveAs( './results/Figures/{0}.png'.format( self.c1.GetName() ) )
        
        
        
            
            
# def puCase(x):
    # return {
        # '0p5': './results/0p5/outPU',
        # '1': './results/1/outPU',
        # '1p5': './results/1p5/outPU',
        # '10': './results/10/outPU',
        # '50': './results/50/outPU',
        # '100': './results/100/outPU',
        # '140': './results/140/outPU',
        # '200': './results/200/outPU',
    # }.get(x, './results/100/outPU')
    
            
def puCase(x):
    return {
        '0' : '/histBibSelfMixOut_pu',
        '0p5': '----------',
        '1': '----------',
        '1p5': '----------',
        '10': '----------',
        '50': '----------',
        '100': '----------',
        '140': '----------',
        '200': '----------',
    }.get(x, '/---------')    

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    #parser.add_argument( '-i' , '--infile' , dest='infile' , help='the name of the input file' , type=str  )
    parser.add_argument( '-s' , '--src' , dest='srcDir' ,default='.', help='the name of the source base directory' , type=str  )
    parser.add_argument( '-p' , '--pu' , dest='PU' , default=None , help='the pu of the file' , type=str , choices=['0', '0p5', '1' , '1p5' , '10' , '50' , '100' , '140' , '200', 'all'])
    
    opt = parser.parse_args()
    destDir=opt.srcDir +'/'
    # if not opt.infile:
        # print('please specify the input file name using -i option')
        # return 1
    
    if not opt.PU:
        print('please specify the pu using -p option')
        return 1
    
    if opt.PU=='all':
        opt.PU = ['0p5', '1' , '1p5' , '10' , '50' , '100' , '140' , '200']
    else :
        opt.PU =[opt.PU]
    
    
    
    #print(opt.PU)
    
    ROOT.gStyle.SetOptTitle(0)
    puDict = {}
    #puDict = dict.fromkeys(opt.PU, list())
    for pu in opt.PU:
        opt_infile = opt.srcDir + puCase(pu) + pu +".root"
        fIn = ROOT.TFile.Open( opt_infile ) 
        print(opt_infile)
        sg_list=list()
        for sg in range( 8 ):
            s = SensorGroupInformation(sg , pu , fIn )
            sg_list.append(copy.deepcopy(s))
            puDict[pu] = sg_list

    print("------------")
    


    for pu in opt.PU:
        print(destDir + "out2D_ZeroCounting_pu{0}.txt".format(pu))
        fz = open(destDir + "out2D_ZeroCounting_pu{0}.txt".format(pu), "w")
        fc = open(destDir + "out2D_DigiCounting_pu{0}.txt".format(pu), "w")
        fs = open(destDir + "out2D_Psim_pu{0}.txt".format(pu), "w")
        fu = open(destDir + "out2D_UnKnownRatio_pu{0}.txt".format(pu), "w")
        print(pu)
        for sg in [4]: #range( 8):
            #puDict[pu][sg].Write()
            #sg=4
            puDict[pu][sg].WriteToA(20.)
            hZc = puDict[pu][sg].hZeroCounting()
            hC = puDict[pu][sg].hCounting()
            hS = puDict[pu][sg].hSimHits()
            hU = puDict[pu][sg].hUnknownRatio()
            for b in range( hZc.GetNbinsX() ) :
                bc = hZc.GetBinContent(b+1)
                fz.write("{0}".format(bc)+", ")
            fz.write("\n")
            for b in range( hC.GetNbinsX() ) :
                bc = hC.GetBinContent(b+1)
                fc.write("{0}".format(bc)+", ")
            fc.write("\n")
            
            for b in range( hS.GetNbinsX() ) :
                bc = hS.GetBinContent(b+1)
                fs.write("{0}".format(bc)+", ")
            fs.write("\n")
            for b in range( hU.GetNbinsX() ) :
                bc = hU.GetBinContent(b+1)
                fu.write("{0}".format(bc)+", ")
            fu.write("\n")
            
            
    fc.close()
    fz.close()
    fs.close()
    fs.close()
    
    
    # np.random.seed(19680801)
    # Z = np.random.rand(6, 10)
    # x = np.arange(-0.5, 10, 1)  # len = 11
    # y = np.arange(4.5, 11, 1)  # len = 7

    # fig, ax = plt.subplots()
    # ax.pcolormesh(x, y, Z)
    # plt.show()
    return 0

if __name__ == "__main__":
    sys.exit( main() )
