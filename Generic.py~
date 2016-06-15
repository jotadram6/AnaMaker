import ROOT
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from array import array
import rootnotes
ROOT.gROOT.LoadMacro('tdrStyle.C')
ROOT.gROOT.LoadMacro('PUR.C')
def DPHI(x):
    if x>=np.pi: return x-2*np.pi
    elif x<-np.pi: return x+2*np.pi
    else: return x
VersionProc="V121"
Base="/home/jruizalv/"
InclusiveVersions=["Version121/","Version122/","Version123/"]

DATA_F="Full_Data.root"

datafiles = [Base+InclusiveVersions[0]+DATA_F,Base+InclusiveVersions[1]+DATA_F,Base+InclusiveVersions[2]+DATA_F]

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

#Corrections
PU_weight=[0.0, 0.283904, 0.223747, 0.155297, 0.26765, 0.298877, 0.723858, 0.544046, 0.536757, 0.733619, 1.0669, 1.48245, 1.80291, 1.80134, 1.58708, 1.34094, 1.15851, 1.0611, 1.02807, 1.04437, 1.08871, 1.12881, 1.15127, 1.15997, 1.15499, 1.12965, 1.07905, 1.00491, 0.914114, 0.813124, 0.707008, 0.599841, 0.495265, 0.39636, 0.306662, 0.22925, 0.166144, 0.118378, 0.0850749, 0.063948, 0.0523125, 0.0476093, 0.0478425, 0.0516782, 0.0583517, 0.0676538, 0.0795915, 0.0944514, 0.112777, 0.135238, 0.162647, 0.196085, 0.236688, 0.286083, 0.34616, 0.418813, 0.506749, 0.612955, 0.741067, 0.894819, 2.26873]
def PUR_function(TI):
  if (TI>=61) or (TI<0): return 1
  else: return PU_weight[TI]

#Declaring samples, tree structure, and ntuples
NsignalFiles=9

files = ["TpJetM600_5318_Full_analyzed.root",
         "TpJetM650_5318_Full_analyzed.root",
         "TpJetM700_5318_Full_analyzed.root",
         "TpJetM750_5318_Full_analyzed.root",
         "TpJetM800_5318_Full_analyzed.root",
         "TpJetM850_5318_Full_analyzed.root",
         "TpJetM900_5318_Full_analyzed.root",
         "TpJetM950_5318_Full_analyzed.root",
         "TpJetM1000_5318_Full_analyzed.root",
         "TTJetsF_5318_Full_analyzed.root",
         "QCD_HT_500_1000_5318_Full_analyzed.root",
         "QCD_HT_1000_Inf_5318_Full_analyzed.root"]

TreeStructureA="M5J:DRHJ:DRWH:RelHT:M2HP:DRTp6thJ:HM:chi2:MTHAsym:NTI:weight"
TreeStructure="M5J:DRHJ:DRWH:RelHT:M2HP:DRTp6thJ:HM:chi2:MTHAsym"

#Signal Sample
TntupleA = ROOT.TNtuple("TntupleA","TntupleA",TreeStructureA)
Q500ntupleA = ROOT.TNtuple("Q500ntupleA","Q500ntupleA",TreeStructureA)
Q1000ntupleA = ROOT.TNtuple("Q1000ntupleA","Q1000ntupleA",TreeStructureA)
QntuplesA=[Q500ntupleA,Q1000ntupleA]
S600ntupleA = ROOT.TNtuple("S600ntupleA","S600ntupleA",TreeStructureA)
S650ntupleA = ROOT.TNtuple("S650ntupleA","S650ntupleA",TreeStructureA)
S700ntupleA = ROOT.TNtuple("S700ntupleA","S700ntupleA",TreeStructureA)
S750ntupleA = ROOT.TNtuple("S750ntupleA","S750ntupleA",TreeStructureA)
S800ntupleA = ROOT.TNtuple("S800ntupleA","S800ntupleA",TreeStructureA)
S850ntupleA = ROOT.TNtuple("S850ntupleA","S850ntupleA",TreeStructureA)
S900ntupleA = ROOT.TNtuple("S900ntupleA","S900ntupleA",TreeStructureA)
S950ntupleA = ROOT.TNtuple("S950ntupleA","S950ntupleA",TreeStructureA)
S1000ntupleA = ROOT.TNtuple("S1000ntupleA","S1000ntupleA",TreeStructureA)
NMassPointsA=[S600ntupleA,S650ntupleA,S700ntupleA,S750ntupleA,S800ntupleA,S850ntupleA,S900ntupleA,S950ntupleA,S1000ntupleA]
ListTreesNamesA=["S600ntupleA","S650ntupleA","S700ntupleA","S750ntupleA","S800ntupleA","S850ntupleA","S900ntupleA","S950ntupleA","S1000ntupleA","TntupleA","Q500ntupleA","Q1000ntupleA"]

#Control Sample
TntupleE = ROOT.TNtuple("TntupleE","TntupleE",TreeStructure)
Q500ntupleE = ROOT.TNtuple("Q500ntupleE","Q500ntupleE",TreeStructure)
Q1000ntupleE = ROOT.TNtuple("Q1000ntupleE","Q1000ntupleE",TreeStructure)
QntuplesE=[Q500ntupleE,Q1000ntupleE]
S600ntupleE = ROOT.TNtuple("S600ntupleE","S600ntupleE",TreeStructure)
S650ntupleE = ROOT.TNtuple("S650ntupleE","S650ntupleE",TreeStructure)
S700ntupleE = ROOT.TNtuple("S700ntupleE","S700ntupleE",TreeStructure)
S750ntupleE = ROOT.TNtuple("S750ntupleE","S750ntupleE",TreeStructure)
S800ntupleE = ROOT.TNtuple("S800ntupleE","S800ntupleE",TreeStructure)
S850ntupleE = ROOT.TNtuple("S850ntupleE","S850ntupleE",TreeStructure)
S900ntupleE = ROOT.TNtuple("S900ntupleE","S900ntupleE",TreeStructure)
S950ntupleE = ROOT.TNtuple("S950ntupleE","S950ntupleE",TreeStructure)
S1000ntupleE = ROOT.TNtuple("S1000ntupleE","S1000ntupleE",TreeStructure)
NMassPointsE=[S600ntupleE,S650ntupleE,S700ntupleE,S750ntupleE,S800ntupleE,S850ntupleE,S900ntupleE,S950ntupleE,S1000ntupleE]
ListTreesNamesE=["S600ntupleE","S650ntupleE","S700ntupleE","S750ntupleE","S800ntupleE","S850ntupleE","S900ntupleE","S950ntupleE","S1000ntupleE","TntupleE","Q500ntupleE","Q1000ntupleE"]

#Sntuple.Print()
#"M5J:DRHJ:DRWH:RelHT:Eta6thJ:M2HP:DRTp6thJ:HM:chi2:MTHAsym"
CutHT = ROOT.TCut("THT>550")
Cut1 = ROOT.TCut("chi2<8")
CutDRbb = ROOT.TCut("DRHJ<=1.2")
#CutDRWH = ROOT.TCut("DRWH>=1.6 && DRWH<=4.0")
CutDRWH = ROOT.TCut("DRWH>=-1.0")
CutHM = ROOT.TCut("HM>=105 && HM<=145")
CutM2HP = ROOT.TCut("M2HP>6.8")
CutRelHT = ROOT.TCut("RelHT>=0.67")
CutDRTp6thJ = ROOT.TCut("DRTp6thJ>4.8")
MTHAsymWin="10.50"; MTHAsymMean="0.17"
CutMTHAsym=ROOT.TCut("MTHAsym>=("+MTHAsymMean+"-"+MTHAsymWin+") && MTHAsym<=("+MTHAsymMean+"+"+MTHAsymWin+")")
#CutMTHAsym = ROOT.TCut("MTHAsym>=0.1 && MTHAsym<=0.3")

Lumi=19694.513

SWeight=Lumi*(0.577/0.94)*0.1437/99464
TWeight=Lumi*247.74/61941200
Q500Weight=Lumi*8426.0/29019521
Q1000Weight=Lumi*204.0/13438527
QWeights=[Q500Weight,Q1000Weight]
NevtsMassPoints=[95167.,98817.,99464.,99375.,95801.,99257.,99174.,95960.,99078.]
XSNevtsMassPoints=[0.2154,0.1778,0.1437,0.1186,0.1,0.0843,0.0726,0.0626,0.0539]
NWeights=Lumi*(0.577/0.94)*np.array(XSNevtsMassPoints)/np.array(NevtsMassPoints)

NevtsDiboson=[9989803.,9989440.,9789108.]
XSDiboson=[33.6, 56.0, 7.6]
DibosonWeights=Lumi*np.array(XSDiboson)/np.array(NevtsDiboson)

NevtsSingleT=[259575.,3752921.,496681.,139803.,1932775.,492545.]
XSSingleT=[3.79, 54.87, 11.1, 1.76, 29.74, 11.1]
SingleTWeights=Lumi*np.array(XSSingleT)/np.array(NevtsSingleT)

NevtsDY=[2001966.,1997817.]
XSDY=[3840.86, 3060.099]
DYWeights=Lumi*np.array(XSDY)/np.array(NevtsDY)

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

def GetEffHisto(CutType,CenValue,Tree,SampleName,Var,BinsLim,CutApp):
    """CutType: 
                'w'-> window --------> CenValue must be also set
                'g'-> great than
                'l'-> less than
    """
    HistName=Var+"BaseHist"
    EffHistName=Var+"EffHist"+SampleName
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

def NormFunc(Cut,InvCut,Var,BinsLim,Tree1,Tree2,CutVar,CutVarBinsLim,CutNm1,Chi2TestBool,RFileName):
    """Tree1: Signal Sample
    Tree2: Control Sample
    Output: N_cs_in,N_cs_out,N_ss_in,N_ss_out,R_cs,R_ss,Np_ss_in
    """
    NHistoSSIn=Var+"SSIn"
    NHistoCSIn=Var+"CSIn"
    NHistoSSOut=Var+"SSOut"
    NHistoCSOut=Var+"CSOut"
    Tree1.Draw(Var+" >> "+NHistoSSIn+BinsLim,Cut)
    Tree1.Draw(Var+" >> "+NHistoSSOut+BinsLim,InvCut)
    Tree2.Draw(Var+" >> "+NHistoCSIn+BinsLim,Cut)
    Tree2.Draw(Var+" >> "+NHistoCSOut+BinsLim,InvCut)
    HistoSSIn=ROOT.gDirectory.Get(NHistoSSIn)
    HistoSSOut=ROOT.gDirectory.Get(NHistoSSOut)
    HistoCSIn=ROOT.gDirectory.Get(NHistoCSIn)
    HistoCSOut=ROOT.gDirectory.Get(NHistoCSOut)
    N_cs_in=HistoCSIn.GetEntries()
    N_cs_out=HistoCSOut.GetEntries()
    N_ss_in=HistoSSIn.GetEntries()
    N_ss_out=HistoSSOut.GetEntries()
    R_cs=N_cs_in/N_cs_out
    R_ss=N_ss_in/N_ss_out
    Np_ss_in=R_cs*N_ss_out
    if Chi2TestBool:
        RFile = ROOT.TFile("HiggsMass"+RFileName+".root", "recreate")
        Tree1.Draw(CutVar+" >> SS_CutVar"+CutVarBinsLim,CutNm1)
        Tree2.Draw(CutVar+" >> CS_CutVar"+CutVarBinsLim,CutNm1)
        CutVarHistoSS=ROOT.gDirectory.Get("SS_CutVar")
        CutVarHistoCS=ROOT.gDirectory.Get("CS_CutVar")
        CutVarHistoSS.Write()
        CutVarHistoCS.Write()
        Chi2TestRes = CutVarHistoCS.Chi2Test(CutVarHistoSS,"UU CHI2/NDF")
        RFile.Close()
        return N_cs_in,N_cs_out,N_ss_in,N_ss_out,R_cs,R_ss,Np_ss_in,Chi2TestRes
    else: return N_cs_in,N_cs_out,N_ss_in,N_ss_out,R_cs,R_ss,Np_ss_in

def AddingTriggerSys(Cut,Var,BinsLim,Tree1,HistoName):
    TriggerSysError=[0.05, 0.03, 0.01, 0.01, 0.01] #as function of jet6 pt: [20,40), [40,60), [60,80), [80,100), [100,120)
    MinPT=[30,50,70,90,110]; MaxPT=[50,70,90,110,130]
    M5JHistos=[]

    RestOfjet6Cut= ROOT.TCut("jet6pt>=120")
    Tree1.Draw(Var+" >> "+HistoName+BinsLim,Cut*RestOfjet6Cut)
    FinalHisto=ROOT.gDirectory.Get(HistoName); FinalHisto.Sumw2()
    print "Entries for M5J with jet6pt>120:", FinalHisto.GetEntries()

    for i in xrange(5):
        CutJet6Pt = ROOT.TCut("(jet6pt>="+str(MinPT[i])+" && jet6pt<"+str(MaxPT[i])+")")
        Tree1.Draw(Var+" >> "+HistoName+str(i)+BinsLim,Cut*CutJet6Pt)
        M5JHistos.append(ROOT.gDirectory.Get(HistoName+str(i)))
        M5JHistos[-1].Sumw2()
        print "Entries for M5J with jet6pt>:{0:.0f}".format(MinPT[i]), M5JHistos[-1].GetEntries()
        if M5JHistos[-1].GetEntries()!=0:
            for j in xrange(1,M5JHistos[-1].GetXaxis().GetNbins()+1):
                BinCon=M5JHistos[-1].GetBinContent(j); SysErr=BinCon*TriggerSysError[i]
                BinErr=M5JHistos[-1].GetBinError(j); TotErr=np.sqrt(BinErr**2+SysErr**2)
                print "For bin=", j, "Base error=", BinErr, "Sys error=", SysErr, "Total Error=", TotErr
                M5JHistos[-1].SetBinError(j,TotErr)
        FinalHisto.Add(M5JHistos[-1])

    return FinalHisto
            

def AddingReWTriggerSys(Cut,Var,BinsLim,Tree1,HistoName):
    TriggerSysError=[1.04, 0.82, 1.04, 1.0, 1.0] #as function of jet6 pt: [20,40), [40,60), [60,80), [80,100), [100,120)
    MinPT=[30,50,70,90,110]; MaxPT=[50,70,90,110,130]
    M5JHistos=[]

    RestOfjet6Cut= ROOT.TCut("jet6pt>=130")
    Tree1.Draw(Var+" >> "+HistoName+BinsLim,Cut*RestOfjet6Cut)
    FinalHisto=ROOT.gDirectory.Get(HistoName); FinalHisto.Sumw2()
    print "Entries for M5J with jet6pt>130:", FinalHisto.GetEntries()

    for i in xrange(5):
        CutJet6Pt = ROOT.TCut("(jet6pt>="+str(MinPT[i])+" && jet6pt<"+str(MaxPT[i])+")")
        Tree1.Draw(Var+" >> "+HistoName+str(i)+BinsLim,Cut*CutJet6Pt)
        M5JHistos.append(ROOT.gDirectory.Get(HistoName+str(i)))
        M5JHistos[-1].Sumw2(); M5JHistos[-1].Scale(TriggerSysError[i])
        print "Entries and weight for M5J with jet6pt>{0:.0f}: {1:.0f} {2:.2f}".format(MinPT[i], M5JHistos[-1].GetEntries(), TriggerSysError[i])
        #if M5JHistos[-1].GetEntries()!=0:
        #    for j in xrange(1,M5JHistos[-1].GetXaxis().GetNbins()+1):
        #        BinCon=M5JHistos[-1].GetBinContent(j); SysErr=BinCon*TriggerSysError[i]
        #        BinErr=M5JHistos[-1].GetBinError(j); TotErr=np.sqrt(BinErr**2+SysErr**2)
        #        print "For bin=", j, "Base error=", BinErr, "Sys error=", SysErr, "Total Error=", TotErr
        #        M5JHistos[-1].SetBinError(j,TotErr)
        FinalHisto.Add(M5JHistos[-1])

    return FinalHisto

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
