import ROOT
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from array import array
import rootnotes
import time
import os
import commands as cmd
#ROOT.gROOT.LoadMacro('tdrStyle.C')

def DPHI(x):
    if x>=np.pi: return x-2*np.pi
    elif x<-np.pi: return x+2*np.pi
    else: return x

#Style functions
def SetAxis(Histo,Axis,TOffset,TSize,LOffset,LSize,Ndiv):
    """Sets offset and size of an axis in the histogram. Axis must be 'X' or 'Y', and Histo should be a valid root histogram"""
    if Axis=='X':
        Histo.GetXaxis().SetTitleSize(TSize)
        Histo.GetXaxis().SetTitleOffset(TOffset)
        Histo.GetXaxis().SetLabelSize(LSize)
        Histo.GetXaxis().SetLabelOffset(LOffset)
        Histo.GetXaxis().SetNdivisions(Ndiv)
    elif Axis=='Y':
        Histo.GetYaxis().SetTitleSize(TSize)
        Histo.GetYaxis().SetTitleOffset(TOffset)
        Histo.GetYaxis().SetLabelSize(LSize)
        Histo.GetYaxis().SetLabelOffset(LOffset)
        Histo.GetYaxis().SetNdivisions(Ndiv)
    else: print "Please correct axis selection: Valid values are 'X' or 'Y'"
def SetCos(Hist,FillColor,FillStyle,LineColor,LineWidth,LineStyle,MarkerStyle):
    """Hist, FillColor, FillStyle, LineColor, LineWidth, LineStyle, MarkerStyle"""
    Hist.SetLineStyle(LineStyle); Hist.SetLineWidth(LineWidth); Hist.SetLineColor(LineColor)
    Hist.SetFillStyle(FillStyle); Hist.SetFillColor(FillColor)
    Hist.SetMarkerStyle(MarkerStyle)

#Propagation of error functions
def DivE(x,y,dx,dy):
    if x!=0 and y!=0: return (x/y)*((dx/x)+(dy/y))
    else: return 0.0

def MulE(x,y,dx,dy):
    if x!=0 and y!=0: return (x*y)*((dx/x)+(dy/y))
    else: return 0.0

def SqrtE(x,dx):
    if x!=0: return np.sqrt(x)*0.5*(dx/x)
    else: return 0.0

def EffE(eff,N):
    if N!=0: return np.sqrt((eff*(1-eff))/N)
    else: return 0.0

def EffV(a,b):
    if b!=0: return a/b
    else: return 0.0

def SysErr(n,N,dn,dN):
    if N!=0: return (abs(N-n)/N)*np.sqrt((((dn+dN)/abs(N-n))**2+(dN/N)**2)) #DivE(abs(N-n),N,(dn+dN),dN)
    else: return 0.0

#Getting Info strings
def GetMR(Histo):
    return "Mean={0:.2f}".format(Histo.GetMean())+"#pm{0:.2f}".format(Histo.GetMeanError())+" RMS={0:.2f}".format(Histo.GetRMS())+"#pm{0:.2f}".format(Histo.GetRMSError())

def GetEWI(Histo):
    INT=Histo.Integral(0,Histo.GetNbinsX()+1)
    ENT=Histo.GetEntries()
    if ENT!=0: W=INT/ENT
    else: W=0
    #return "Entries={0:.2f}".format(ENT)+" W={0:.2f}".format(W)+" Int={0:.2f}".format(INT)+"+-{0:.2f}".format(W*np.sqrt(ENT))
    return "Entries={0:.2f}".format(ENT)+" Int={0:.2f}".format(INT)+"#pm{0:.2f}".format(W*np.sqrt(ENT))

#Function to calculate efficiencies from cuts on a varialbe
def GetEffHisto(CutType,CenValue,Tree,SampleName,Var,BinsLim,CutApp):
    """CutType: 
                'w'-> window --------> CenValue must be also set
                'g'-> great than
                'l'-> less than
    """
    HistName=Var.replace("(","").replace(")","").replace("[","").replace("]","").replace(":","").replace(".","").replace("-","_")+"BaseHist"
    EffHistName=Var.replace("(","").replace(")","").replace("[","").replace("]","").replace(":","").replace(".","").replace("-","_")+"EffHist"+SampleName
    print Var+" >> "+HistName+BinsLim
    if CutType=="w":
        WVar="TMath::Abs("+Var+"-"+CenValue+")"
        Tree.Draw(WVar+" >> "+HistName+BinsLim,CutApp)
        TemHisto=ROOT.gDirectory.Get(HistName)
        Eff=TemHisto.Clone(EffHistName)
        for j in xrange(1,TemHisto.GetXaxis().GetNbins()+2):
            EffBin=1.0; EffErr=1.0
            if TemHisto.Integral()!=0:
                EffBin=TemHisto.Integral(1,j)/TemHisto.Integral()
                if EffBin>=0.99: EffErr=0.0
                else: EffErr=EffE(EffBin,TemHisto.Integral())
            Eff.SetBinContent(j,EffBin); Eff.SetBinError(j,EffErr)
    else:
        Tree.Draw(Var+" >> "+HistName+BinsLim,CutApp)
        TemHisto=ROOT.gDirectory.Get(HistName)
        Eff=TemHisto.Clone(EffHistName)
        for j in xrange(1,TemHisto.GetXaxis().GetNbins()+2):
            EffBin=1.0; EffErr=1.0
            if TemHisto.Integral()!=0:
                if CutType=="g": EffBin=TemHisto.Integral(j,TemHisto.GetXaxis().GetNbins()+2)/TemHisto.Integral(0,TemHisto.GetXaxis().GetNbins()+2)
                elif CutType=="l": EffBin=TemHisto.Integral(1,j)/TemHisto.Integral(0,TemHisto.GetXaxis().GetNbins()+2)
                if EffBin>=0.99: EffErr=0.0
                else: EffErr=EffE(EffBin,TemHisto.Integral())
            Eff.SetBinContent(j,EffBin); Eff.SetBinError(j,EffErr)
    return Eff

#Functions to calculate the correct error for very low statistics (specially bins with zero entries)
def CorrectErrorLowStats(H1):
    H1.SetBinErrorOption(ROOT.TH1.kPoisson)
    for j in xrange(1,H1.GetXaxis().GetNbins()+1):
        H1.GetBinErrorLow(j)
        H1.GetBinErrorUp(j)

def AddCorrectError(Histo1,Histo2,Weight1,Weight2):
    "Adds Histo2 to Histo1"
    for j in xrange(1,Histo1.GetXaxis().GetNbins()+1):
        BinC1=Histo1.GetBinContent(j)
        BinC2=Histo2.GetBinContent(j)
        if BinC1!=0: 
            BinE1=Histo1.GetBinError(j)
        else:
            BinE1=1.841*Weight1
        if BinC2!=0: 
            BinE2=Histo2.GetBinError(j)
        else:
            BinE2=1.841*Weight2
        BinCT=BinC1+BinC2
        BinET=np.sqrt(BinE1**2+BinE2**2)
        print "histo 1:", BinC1, BinE1
        print "histo 2:", BinC2, BinE2
        print "sum histo:", BinCT, BinET
        Histo1.SetBinContent(j,BinCT)
        Histo1.SetBinError(j,BinET)
