import ROOT as rt
from array import *
from sms import *
from color import *

class smsPlotABS(object):
    # modelname is the sms name (see sms.py)
    # histo is the 2D xsec map
    # obsLimits is a list of opbserved limits [NOMINAL, +1SIGMA, -1SIGMA]
    # expLimits is a list of expected limits [NOMINAL, +1SIGMA, -1SIGMA]
    # label is a label referring to the analysis (e.g. RA1, RA2, RA2b, etc)

    def __init__(self, modelname, histo, obsLimits, expLimits, energy, lumi, preliminary, label):
        self.standardDef(modelname, histo, obsLimits, expLimits, energy, lumi, preliminary)
        self.LABEL = label
        self.c = rt.TCanvas("cABS_%s" %label,"cABS_%s" %label,300,300)
        self.histo = histo

    def standardDef(self, modelname, histo, obsLimits, expLimits, energy, lumi, preliminary):
        # which SMS?
        #print modelname
        self.model = sms(modelname)
        self.OBS = obsLimits
        self.EXP = expLimits
        self.lumi = lumi
        self.energy = energy
        self.preliminary = preliminary
        # create the reference empty histo
        self.emptyhisto = self.emptyHistogramFromModel()

    def emptyHistogramFromModel(self):
        self.emptyHisto = rt.TH2D("emptyHisto", "", 1, self.model.Xmin, self.model.Xmax, 1, self.model.Ymin, self.model.Ymax)
        
    # define the plot canvas
    def setStyle(self):
        # canvas style
        rt.gStyle.SetOptStat(0)
        rt.gStyle.SetOptTitle(0)        

        self.c.SetLogz()
        self.c.SetTickx(1)
        self.c.SetTicky(1)

        self.c.SetRightMargin(0.19)
        self.c.SetTopMargin(0.08)
        self.c.SetLeftMargin(0.14)
        self.c.SetBottomMargin(0.14)

        # set x axis
        self.emptyHisto.GetXaxis().SetLabelFont(42)
        self.emptyHisto.GetXaxis().SetLabelOffset(0.015)
        self.emptyHisto.GetXaxis().SetLabelSize(0.04)
        self.emptyHisto.GetXaxis().SetTitleFont(42)
        self.emptyHisto.GetXaxis().SetTitleSize(0.05)
        self.emptyHisto.GetXaxis().SetTitleOffset(1.2)
        self.emptyHisto.GetXaxis().SetTitle(self.model.sParticle)
        #self.emptyHisto.GetXaxis().CenterTitle(True)

        # set y axis
        self.emptyHisto.GetYaxis().SetLabelFont(42)
        self.emptyHisto.GetYaxis().SetLabelOffset(0.01)
        self.emptyHisto.GetYaxis().SetLabelSize(0.04)
        self.emptyHisto.GetYaxis().SetTitleFont(42)
        self.emptyHisto.GetYaxis().SetTitleSize(0.05)
        self.emptyHisto.GetYaxis().SetTitleOffset(1.3)
        self.emptyHisto.GetYaxis().SetTitle(self.model.LSP)
        #self.emptyHisto.GetYaxis().CenterTitle(True)
                
    def DrawText(self):
        #redraw axes
        self.c.RedrawAxis()
        # white background
        graphWhite = rt.TGraph(5)
        graphWhite.SetName("white")
        graphWhite.SetTitle("white")
        graphWhite.SetFillColor(rt.kWhite)
        graphWhite.SetFillStyle(1001)
        graphWhite.SetLineColor(rt.kBlack)
        graphWhite.SetLineStyle(1)
        graphWhite.SetLineWidth(3)
        graphWhite.SetPoint(0,self.model.Xmin, self.model.Ymax)
        graphWhite.SetPoint(1,self.model.Xmax, self.model.Ymax)
        graphWhite.SetPoint(2,self.model.Xmax, self.model.Ymax*0.69)
        graphWhite.SetPoint(3,self.model.Xmin, self.model.Ymax*0.69)
        graphWhite.SetPoint(4,self.model.Xmin, self.model.Ymax)
        graphWhite.Draw("FSAME")
        graphWhite.Draw("LSAME")
        self.c.graphWhite = graphWhite
        
        # CMS LABEL
##        textCMS = rt.TLatex(0.22,0.98,"CMS %s, %s fb^{-1}, #sqrt{s} = %s TeV" %(self.preliminary, self.lumi, self.energy))
##        textCMS = rt.TLatex(0.22,0.98,"CMS %s, 18.1 fb^{-1}, #sqrt{s} = %s TeV" %(self.preliminary, self.energy))
        textCMS = rt.TLatex(0.15,0.98,  "CMS                                     18.1 fb^{-1} (%s TeV)" %(self.energy))
##        textCMS = rt.TLatex(0.17,0.98,"CMS %s, 18.1-19.6 fb^{-1}, #sqrt{s} = %s TeV" %(self.preliminary, self.energy))
##        textCMS = rt.TLatex(0.15,0.98,  "CMS                             18.1-19.6 fb^{-1} (%s TeV)" %(self.energy))
##        textCMS = rt.TLatex(0.15,0.98,"CMS             18.1-19.6 fb^{-1}           #sqrt{s} = %s TeV" %(self.energy))
        textCMS.SetNDC()
        textCMS.SetTextAlign(13)
        textCMS.SetTextFont(42)
        textCMS.SetTextSize(0.038)
        textCMS.Draw()
        self.c.textCMS = textCMS

        # Masses LABEL
        textMass = rt.TLatex(0.05,0.07,"m_{#tilde{#tau}} = m_{#tilde{#nu}_{#tau}} = 0.5 (m_{#tilde{#chi}^{#pm}_{1}} + m_{#tilde{#chi}^{0}_{1}})")
        textMass.SetNDC()
        textMass.SetTextAlign(13)
        textMass.SetTextFont(42)
        textMass.SetTextSize(0.038)
        #textMass.Draw()
        self.c.textMass = textMass

        # MODEL LABEL
        textModelLabel= rt.TLatex(0.16,0.90,"%s NLO+NLL exclusion" %self.model.label)
        textModelLabel.SetNDC()
        textModelLabel.SetTextAlign(13)
        textModelLabel.SetTextFont(42)
        textModelLabel.SetTextSize(0.040)
        #textModelLabel.Draw()
        self.c.textModelLabel = textModelLabel

        # NLO NLL XSEC
        textNLONLL= rt.TLatex(0.16,0.32,"NLO-NLL exclusion")
        textNLONLL.SetNDC()
        textNLONLL.SetTextAlign(13)
        textNLONLL.SetTextFont(42)
        textNLONLL.SetTextSize(0.040)
        textNLONLL.Draw()
        #self.c.textNLONLL = textNLONLL
 
        # "LEP excluded
        textLEPExcl= rt.TLatex(0.16,0.27,"LEP excluded")
        textLEPExcl.SetNDC()
        textLEPExcl.SetTextAlign(13)
        textLEPExcl.SetTextFont(42)
        textLEPExcl.SetTextSize(0.0275)
        textLEPExcl.SetTextAngle(-45.);
        textLEPExcl.Draw()
        self.c.textLEPExcl = textLEPExcl

    def Save(self,label):
        # save the output
        self.c.SaveAs("%s.pdf" %label)
        #self.c.SaveAs("%s.C" %label)
        
    def DrawLegend(self):
        xRange = self.model.Xmax-self.model.Xmin
        yRange = self.model.Ymax-self.model.Ymin

##        textModelLabelMinus1 = rt.TLatex(self.model.Xmin+3*xRange/100, self.model.Ymax-0.40*yRange/100*10, "e#tau_{h}, #mu#tau_{h} and #tau_{h}#tau_{h} combined")
        textModelLabelMinus1 = rt.TLatex(self.model.Xmin+3*xRange/100, self.model.Ymax-0.40*yRange/100*10, "#tau_{h}#tau_{h} channel")
        textModelLabelMinus1.SetTextFont(42)
        textModelLabelMinus1.SetTextSize(0.04)
        textModelLabelMinus1.Draw()
        self.c.textModelLabelMinus1 = textModelLabelMinus1

        textModelLabel0 = rt.TLatex(self.model.Xmin+43*xRange/100, self.model.Ymax-1.0*yRange/100*10, "NLO+NLL exclusion")
        textModelLabel0.SetTextFont(42)
        textModelLabel0.SetTextSize(0.030)
        textModelLabel0.Draw()
        self.c.textModelLabel0 = textModelLabel0
        textModelLabel1 = rt.TLatex(self.model.Xmin+3*xRange/100, self.model.Ymax-1.0*yRange/100*10, "pp #rightarrow #tilde{#chi}^{+}_{1} #tilde{#chi}^{-}_{1}")
        textModelLabel1.SetTextFont(42)
        textModelLabel1.SetTextSize(0.030)
        textModelLabel1.Draw()
        self.c.textModelLabel1 = textModelLabel1
        textModelLabel2 = rt.TLatex(self.model.Xmin+3*xRange/100, self.model.Ymax-1.60*yRange/100*10, "#tilde{#chi}^{#pm}_{1} #rightarrow #tilde{#tau}^{#pm} #nu_{#tau},#tilde{#nu}_{#tau} #tau^{#pm}")
        textModelLabel2.SetTextFont(42)
        textModelLabel2.SetTextSize(0.030)
        textModelLabel2.Draw()
        self.c.textModelLabel2 = textModelLabel2
        textModelLabel3 = rt.TLatex(self.model.Xmin+3*xRange/100, self.model.Ymax-2.2*yRange/100*10, "#tilde{#tau}^{#pm} #rightarrow #tau^{#pm} #tilde{#chi}_{1}^{0}, #tilde{#nu} #rightarrow #nu #tilde{#chi}_{1}^{0}")
        textModelLabel3.SetTextFont(42)
        textModelLabel3.SetTextSize(0.030)
        textModelLabel3.Draw()
        self.c.textModelLabel3 = textModelLabel3
 
        textMassLabel = rt.TLatex(self.model.Xmin+3*xRange/100, self.model.Ymax-2.7*yRange/100*10,"m_{#tilde{#tau}} = m_{#tilde{#nu}_{#tau}} = 0.5 (m_{#tilde{#chi}^{#pm}_{1}}+ m_{#tilde{#chi}^{0}_{1}})")
        textMassLabel.SetTextFont(42)
        textMassLabel.SetTextSize(0.030)
        textMassLabel.Draw()
        self.c.textMassLabel = textMassLabel


        
        LObs = rt.TGraph(2)
        LObs.SetName("LObs")
        LObs.SetTitle("LObs")
        LObs.SetLineColor(color(self.OBS['colorLine']))
        LObs.SetLineStyle(1)
        LObs.SetLineWidth(4)
        LObs.SetMarkerStyle(20)
        LObs.SetPoint(0,self.model.Xmin+43*xRange/100, self.model.Ymax-1.35*yRange/100*10)
        LObs.SetPoint(1,self.model.Xmin+50*xRange/100, self.model.Ymax-1.35*yRange/100*10)

        LObsP = rt.TGraph(2)
        LObsP.SetName("LObsP")
        LObsP.SetTitle("LObsP")
        LObsP.SetLineColor(color(self.OBS['colorLine']))
        LObsP.SetLineStyle(1)
        LObsP.SetLineWidth(2)
        LObsP.SetMarkerStyle(20)
        LObsP.SetPoint(0,self.model.Xmin+43*xRange/100, self.model.Ymax-1.20*yRange/100*10)
        LObsP.SetPoint(1,self.model.Xmin+50*xRange/100, self.model.Ymax-1.20*yRange/100*10)

        LObsM = rt.TGraph(2)
        LObsM.SetName("LObsM")
        LObsM.SetTitle("LObsM")
        LObsM.SetLineColor(color(self.OBS['colorLine']))
        LObsM.SetLineStyle(1)
        LObsM.SetLineWidth(2)
        LObsM.SetMarkerStyle(20)
        LObsM.SetPoint(0,self.model.Xmin+43*xRange/100, self.model.Ymax-1.60*yRange/100*10)
        LObsM.SetPoint(1,self.model.Xmin+50*xRange/100, self.model.Ymax-1.60*yRange/100*10)

        textObs = rt.TLatex(self.model.Xmin+51*xRange/100, self.model.Ymax-1.60*yRange/100*10, "Observed #pm1 #sigma_{theory}")
        textObs.SetTextFont(42)
        textObs.SetTextSize(0.030)
        textObs.Draw()
        self.c.textObs = textObs

        LExpP = rt.TGraph(2)
        LExpP.SetName("LExpP")
        LExpP.SetTitle("LExpP")
        LExpP.SetLineColor(color(self.EXP['colorLine']))
        LExpP.SetLineStyle(7)
        LExpP.SetLineWidth(2)  
        LExpP.SetPoint(0,self.model.Xmin+43*xRange/100, self.model.Ymax-1.85*yRange/100*10)
        LExpP.SetPoint(1,self.model.Xmin+50*xRange/100, self.model.Ymax-1.85*yRange/100*10)

        LExp = rt.TGraph(2)
        LExp.SetName("LExp")
        LExp.SetTitle("LExp")
        LExp.SetLineColor(color(self.EXP['colorLine']))
        LExp.SetLineStyle(7)
        LExp.SetLineWidth(4)
        LExp.SetPoint(0,self.model.Xmin+43*xRange/100, self.model.Ymax-2.00*yRange/100*10)
        LExp.SetPoint(1,self.model.Xmin+50*xRange/100, self.model.Ymax-2.00*yRange/100*10)
        
        LExpM = rt.TGraph(2)
        LExpM.SetName("LExpM")
        LExpM.SetTitle("LExpM")
        LExpM.SetLineColor(color(self.EXP['colorLine']))
        LExpM.SetLineStyle(7)
        LExpM.SetLineWidth(2)  
        LExpM.SetPoint(0,self.model.Xmin+43*xRange/100, self.model.Ymax-2.2*yRange/100*10)
        LExpM.SetPoint(1,self.model.Xmin+50*xRange/100, self.model.Ymax-2.2*yRange/100*10)

        textExp = rt.TLatex(self.model.Xmin+51*xRange/100, self.model.Ymax-2.2*yRange/100*10, "Expected #pm1 #sigma_{experiment}")
        textExp.SetTextFont(42)
        textExp.SetTextSize(0.030)
        textExp.Draw()
        self.c.textExp = textExp

        LObs.Draw("LSAME")
        LObsM.Draw("LSAME")
        LObsP.Draw("LSAME")
        LExp.Draw("LSAME")
        LExpM.Draw("LSAME")
        LExpP.Draw("LSAME")
        
        self.c.LObs = LObs
        self.c.LObsM = LObsM
        self.c.LObsP = LObsP
        self.c.LExp = LExp
        self.c.LExpM = LExpM
        self.c.LExpP = LExpP

    def DrawDiagonal(self):
        diagonal = rt.TGraph(3, self.model.diagX, self.model.diagY)
        diagonal.SetName("diagonal")
        diagonal.SetFillColor(rt.kWhite)
        diagonal.SetLineColor(rt.kGray)
        diagonal.SetLineStyle(2)
        #diagonal.Draw("FSAME")
        diagonal.Draw("LSAME")
        self.c.diagonal = diagonal
        
    def DrawLines(self):
        # observed
        self.OBS['nominal'].SetLineColor(color(self.OBS['colorLine']))
        self.OBS['nominal'].SetLineStyle(1)
        self.OBS['nominal'].SetLineWidth(4)
        # observed + 1sigma
        self.OBS['plus'].SetLineColor(color(self.OBS['colorLine']))
        self.OBS['plus'].SetLineStyle(1)
        self.OBS['plus'].SetLineWidth(2)        
        # observed - 1sigma
        self.OBS['minus'].SetLineColor(color(self.OBS['colorLine']))
        self.OBS['minus'].SetLineStyle(1)
        self.OBS['minus'].SetLineWidth(2)        
        # expected + 1sigma
        self.EXP['plus'].SetLineColor(color(self.EXP['colorLine']))
        self.EXP['plus'].SetLineStyle(7)
        self.EXP['plus'].SetLineWidth(2)                
        # expected
        self.EXP['nominal'].SetLineColor(color(self.EXP['colorLine']))
        self.EXP['nominal'].SetLineStyle(7)
        self.EXP['nominal'].SetLineWidth(4)        
        # expected - 1sigma
        self.EXP['minus'].SetLineColor(color(self.EXP['colorLine']))
        self.EXP['minus'].SetLineStyle(7)
        self.EXP['minus'].SetLineWidth(2)                        
        # DRAW LINES
        self.EXP['nominal'].Draw("LSAME")
        self.EXP['plus'].Draw("LSAME")
        self.EXP['minus'].Draw("LSAME")
        self.OBS['nominal'].Draw("LSAME")
        self.OBS['plus'].Draw("LSAME")
        self.OBS['minus'].Draw("LSAME")        

        
