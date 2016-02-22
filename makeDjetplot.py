#!/usr/bin/env python

import array
import collections
import loadlib
import os
import ROOT
import style
from math import *

#plotloc='"~/www/VBF/Djet/'
plotloc=''

class Plot(object):
    maindir = "root://lxcms03://data3/Higgs/160203/"
    basename = "ZZ4lAnalysis.root"
    min = 0.
    max = 1.
    bins = 100
    units = ""

    def __init__(self, title, color, style=1001, *CJLSTdirs, **kwargs):
        for kwarg in kwargs:
            if kwarg == "maindir":
                self.maindir = kwargs[kwarg]
            elif kwarg == "basename":
                self.basename = kwargs[kwarg]
            else:
                raise TypeError("Unknown kwarg: %s=%s" % (kwarg, kwargs[kwarg]))
        self.filenames = [os.path.join(self.maindir, dir, self.basename) for dir in CJLSTdirs]
        self.title = title
        self.color = color
        self.style = style
        self._h = None

    def __hash__(self):
        return hash((tuple(self.filenames), self.title, self.color))
    def __str__(self):
        return self.title

    def addtolegend(self, tlegend, option = "l"):
        tlegend.AddEntry(self.h(), self.title, option)

class TreePlot(Plot):
    def __init__(self, title, color, style, *CJLSTdirs, **kwargs):
        Plot.__init__(self, title, color, style,  *CJLSTdirs, **kwargs)

    def h(self, bins = None):
        if self._h is not None and bins is None:
            return self._h
        if bins is None:
            bins = [Bin(-1, float("inf"))]
        t = ROOT.TChain("ZZTree/candTree")
        for filename in self.filenames:
            t.Add(filename)

        h = {}
        sumofweights = {}
        for bin in bins:
            h[bin] = ROOT.TH1F("h"+self.title+str(bin), "D_{jet}", self.bins, self.min, self.max)
            h[bin].SetLineColor(self.color)
            h[bin].SetMarkerColor(self.color)
            h[bin].SetLineWidth(3)
            h[bin].Sumw2()
            sumofweights[bin] = 0

        length = t.GetEntries()

        for i, entry in enumerate(t):
            t.GetEntry(i)

            wt = entry.genHEPMCweight

            try:
                Djet = (entry.pvbf_VAJHU_old / (entry.pvbf_VAJHU_old + entry.phjj_VAJHU_old))
            except ZeroDivisionError:
                pass

            for bin in bins:
                if bin.min < entry.ZZMass < bin.max:
                    if entry.pvbf_VAJHU_old >= 0 and entry.phjj_VAJHU_old >= 0 and not (entry.pvbf_VAJHU_old == entry.phjj_VAJHU_old == 0):
                        h[bin].Fill(Djet, wt)

                    sumofweights[bin] += wt

            if (i+1) % 10000 == 0 or i+1 == length:
                print (i+1), "/", length

        for bin in bins:
            try:
                h[bin].Scale(1.0/sumofweights[bin])
            except ZeroDivisionError:
                pass

        if len(h) == 1:
            self._h = h.values()[0]
        else:
            self._h = h
        return self._h

    def vbf_eff_jec(self, bins=None):
        t = ROOT.TChain("ZZTree/candTree")
        for filename in self.filenames:
            t.Add(filename)

        n_pass_ct = {}
        n_fail_ct = {}
        n_pass_ct_sumw2 = {}
        n_fail_ct_sumw2 = {}
        n_pass_up = {}
        n_pass_dn = {}
        n_total = {}
        eff_jec = {}
        for bin in bins:
            n_pass_ct[bin] = 0.0
            n_fail_ct[bin] = 0.0
            n_pass_ct_sumw2[bin] = 0.0
            n_fail_ct_sumw2[bin] = 0.0
            n_pass_up[bin] = 0.0
            n_pass_dn[bin] = 0.0
            n_total[bin] = 0.0
            eff_jec[bin] = [0.0, 0.0, 0.0, 0.0] # 4 values: jec center, jec up, jec down, statistical error

        print 'Processing ',self.title
        length = t.GetEntries()

        for i, entry in enumerate(t):
            #t.GetEntry(i)

            wt = entry.genHEPMCweight

            try:
                pass_ct = entry.nCleanedJetsPt30>=2 and 1.0/(1.0+entry.phjj_VAJHU_old/entry.pvbf_VAJHU_old)>0.5
            except ZeroDivisionError:
                pass_ct = False
                pass
            try:
                pass_up = entry.nCleanedJetsPt30_jecUp>=2 and 1.0/(1.0+entry.phjj_VAJHU_old_up/entry.pvbf_VAJHU_old_up)>0.5
            except ZeroDivisionError:
                pass_up = False
                pass
            try:
                pass_dn = entry.nCleanedJetsPt30_jecDn>=2 and 1.0/(1.0+entry.phjj_VAJHU_old_dn/entry.pvbf_VAJHU_old_dn)>0.5
            except ZeroDivisionError:
                pass_dn = False
                pass

            if (i+1) % 10000 == 0 or i+1 == length:
                print (i+1), "/", length

            for bin in bins:
                if bin.min < entry.ZZMass < bin.max:
                    n_total[bin] += wt
                    if pass_ct:
                        n_pass_ct[bin] += wt
                        n_pass_ct_sumw2[bin] +=wt*wt
                    else:
                        n_fail_ct[bin] += wt
                        n_fail_ct_sumw2[bin] +=wt*wt
                    if pass_up: n_pass_up[bin] += wt
                    if pass_dn: n_pass_dn[bin] += wt

        for bin in bins:
            try:
                eff_jec[bin][0] = n_pass_ct[bin]/n_total[bin] # jec center
                eff_jec[bin][1] = n_pass_up[bin]/n_total[bin] # jec up
                eff_jec[bin][2] = n_pass_dn[bin]/n_total[bin] # jec down
                eff_jec[bin][3] = sqrt((n_fail_ct[bin]**2*n_pass_ct_sumw2[bin]+n_pass_ct[bin]**2*n_fail_ct_sumw2[bin])/(n_fail_ct[bin]**4)) # statistical uncertainty deff^2=(nf^2*dnp^2+np^2*dnf^2)/nf^4
            except ZeroDivisionError:
                eff_jec[bin][0] = 0.0
                eff_jec[bin][1] = 0.0
                eff_jec[bin][2] = 0.0
                eff_jec[bin][3] = 0.0
                pass

        return eff_jec


class ZXPlot(Plot):
    def __init__(self, color, style):
        Plot.__init__(self, "Z+X", color, style, "AllData")

    def h(self, bins = None):
        if self._h is not None and bins is None:
            return self._h
        if bins is None:
            bins = [Bin(-1, float("inf"))]
        t = ROOT.TChain("CRZLLTree/candTree")
        for filename in self.filenames:
            t.Add(filename)

        h = {}
        sumofweights = {}
        for bin in bins:
            h[bin] = ROOT.TH1F("h"+self.title+str(bin), "D_{jet}", self.bins, self.min, self.max)
            h[bin].SetLineColor(self.color)
            h[bin].SetMarkerColor(self.color)
            h[bin].SetLineWidth(3)
            h[bin].Sumw2()
            sumofweights[bin] = 0

        length = t.GetEntries()

        for i, entry in enumerate(t):
            wt = ROOT.fakeRate13TeV(entry.LepPt.at(2),entry.LepEta.at(2),entry.LepLepId.at(2)) * ROOT.fakeRate13TeV(entry.LepPt.at(3),entry.LepEta.at(3),entry.LepLepId.at(3))

            try:
                Djet = (entry.pvbf_VAJHU_old / (entry.pvbf_VAJHU_old + entry.phjj_VAJHU_old))
            except ZeroDivisionError:
                pass

            for bin in bins:
                if bin.min < entry.ZZMass < bin.max:
                    if entry.pvbf_VAJHU_old >= 0 and entry.phjj_VAJHU_old >= 0 and not (entry.pvbf_VAJHU_old == entry.phjj_VAJHU_old == 0):
                        h[bin].Fill(Djet, wt)

                    sumofweights[bin] += wt

            if (i+1) % 10000 == 0 or i+1 == length:
                print (i+1), "/", length

        for bin in bins:
            try:
                h[bin].Scale(1.0/sumofweights[bin])
            except ZeroDivisionError:
                pass

        if len(h) == 1:
            self._h = h.values()[0]
        else:
            self._h = h
        return self._h

    def vbf_eff_jec(self, bins=None):
        t = ROOT.TChain("CRZLLTree/candTree")
        for filename in self.filenames:
            t.Add(filename)

        n_pass_ct = {}
        n_fail_ct = {}
        n_pass_ct_sumw2 = {}
        n_fail_ct_sumw2 = {}
        n_pass_up = {}
        n_pass_dn = {}
        n_total = {}
        eff_jec = {}
        for bin in bins:
            n_pass_ct[bin] = 0.0
            n_fail_ct[bin] = 0.0
            n_pass_ct_sumw2[bin] = 0.0
            n_fail_ct_sumw2[bin] = 0.0
            n_pass_up[bin] = 0.0
            n_pass_dn[bin] = 0.0
            n_total[bin] = 0.0
            eff_jec[bin] = [0.0, 0.0, 0.0, 0.0] # 4 values: jec center, jec up, jec down, statistical error

        print 'Processing ',self.title
        length = t.GetEntries()

        for i, entry in enumerate(t):
            #t.GetEntry(i)

            wt = ROOT.fakeRate13TeV(entry.LepPt.at(2),entry.LepEta.at(2),entry.LepLepId.at(2)) * ROOT.fakeRate13TeV(entry.LepPt.at(3),entry.LepEta.at(3),entry.LepLepId.at(3))

            try:
                pass_ct = entry.nCleanedJetsPt30>=2 and 1.0/(1.0+entry.phjj_VAJHU_old/entry.pvbf_VAJHU_old)>0.5
            except ZeroDivisionError:
                pass_ct = False
                pass
            try:
                pass_up = entry.nCleanedJetsPt30_jecUp>=2 and 1.0/(1.0+entry.phjj_VAJHU_old_up/entry.pvbf_VAJHU_old_up)>0.5
            except ZeroDivisionError:
                pass_up = False
                pass
            try:
                pass_dn = entry.nCleanedJetsPt30_jecDn>=2 and 1.0/(1.0+entry.phjj_VAJHU_old_dn/entry.pvbf_VAJHU_old_dn)>0.5
            except ZeroDivisionError:
                pass_dn = False
                pass

            if (i+1) % 10000 == 0 or i+1 == length:
                print (i+1), "/", length

            for bin in bins:
                if bin.min < entry.ZZMass < bin.max:
                    n_total[bin] += wt
                    if pass_ct: 
                        n_pass_ct[bin] += wt
                        n_pass_ct_sumw2[bin] +=wt*wt
                    else:
                        n_fail_ct[bin] += wt
                        n_fail_ct_sumw2[bin] +=wt*wt
                    if pass_up: n_pass_up[bin] += wt
                    if pass_dn: n_pass_dn[bin] += wt

        for bin in bins:
            try:
                eff_jec[bin][0] = n_pass_ct[bin]/n_total[bin] # jec center
                eff_jec[bin][1] = n_pass_up[bin]/n_total[bin] # jec up
                eff_jec[bin][2] = n_pass_dn[bin]/n_total[bin] # jec down
                eff_jec[bin][3] = sqrt((n_fail_ct[bin]**2*n_pass_ct_sumw2[bin]+n_pass_ct[bin]**2*n_fail_ct_sumw2[bin])/(n_fail_ct[bin]**4)) # statistical uncertainty deff^2=(nf^2*dnp^2+np^2*dnf^2)/nf^4
            except ZeroDivisionError:
                pass

        return eff_jec

def makeDjetplots(*plots):
    c1 = ROOT.TCanvas()
    legend = ROOT.TLegend(0.6, 0.5, 0.9, 0.9)
    legend.SetLineStyle(0)
    legend.SetLineColor(0)
    legend.SetFillStyle(0)
    hstack = ROOT.THStack("hstack", "D_{jet}")
    max, min, bins, units = None, None, None, None
    for plot in plots:
        hstack.Add(plot.h())
        plot.addtolegend(legend)
        if max is None:
            max, min, bins, units = plot.max, plot.min, plot.bins, plot.units
        assert (max, min, bins, units) == (plot.max, plot.min, plot.bins, plot.units)
    hstack.Draw("nostackhist")
    hstack.GetXaxis().SetTitle("D_{jet}")
    hstack.GetYaxis().SetTitle("fraction of events / %s%s" % ((max-min)/bins, " "+units if units else ""))
    legend.Draw()
    c1.SaveAs(plotloc+"Djet.png")
    c1.SaveAs(plotloc+"Djet.eps")
    c1.SaveAs(plotloc+"Djet.root")
    c1.SaveAs(plotloc+"Djet.pdf")

class Bin(object):
    def __init__(self, min, max):
        self.min = min
        self.max = max
        self.center = (max+min)*.5
        self.error = (max-min)*.5
    def __str__(self):
        return "%s-%s GeV" % (self.min, self.max)

def makeDjettable(massbins, *plots):
    print massbins
    bins = [Bin(massbins[i], massbins[i+1]) for i in range(len(massbins)-1)]
    for a in bins: print a

    fraction = collections.OrderedDict()
    x = {}
    y = {}
    ex = {}
    ey = {}
    nbins = {}
    g = {}
    mg = ROOT.TMultiGraph()
    legend = ROOT.TLegend(0.6, 0.4, 0.9, 0.8)
    legend.SetLineStyle(0)
    legend.SetLineColor(0)
    legend.SetFillStyle(0)

    for plot in plots:
        x[plot] = array.array("d")
        y[plot] = array.array("d")
        ex[plot] = array.array("d")
        ey[plot] = array.array("d")
        nbins[plot] = 0
        fraction[plot] = collections.OrderedDict()

        h = plot.h(bins)
        for bin in bins:
            integralerror = array.array("d", [0])
            fraction[plot][bin] = h[bin].IntegralAndError(51, 100, integralerror)
            if plot.title == "ttH" and bin.min >= 500: continue
            nbins[plot] += 1
            x[plot].append(bin.center)
            y[plot].append()
            ex[plot].append(bin.error)
            ey[plot].append(integralerror[0])
        g[plot] = ROOT.TGraphErrors(nbins[plot], x[plot], y[plot], ex[plot], ey[plot])
        mg.Add(g[plot])
        g[plot].SetLineColor(plot.color)
        g[plot].SetMarkerColor(plot.color)
        legend.AddEntry(g[plot], plot.title, "lp")

    c1 = ROOT.TCanvas()
    mg.Draw("AP")
    mg.GetXaxis().SetTitle("m_{4l}")
    mg.GetYaxis().SetTitle("fraction of events with D_{jet}>0.5")
    legend.Draw()
    c1.SaveAs(plotloc+"fraction.png")
    c1.SaveAs(plotloc+"fraction.eps")
    c1.SaveAs(plotloc+"fraction.root")
    c1.SaveAs(plotloc+"fraction.pdf")

    print r"\begin{center}"
    print r"\begin{tabular}{ |%s| }" % ("|".join("c" * (len(plots)+1)))

    #http://stackoverflow.com/a/9536084
    header_format = " & ".join(["{:>15}"] * (len(plots) + 1)) + r" \\"
    row_format = " & ".join(["{:>15}"] + [r"{:14.2f}\%"] * (len(plots))) + r"\\"
    print r"\hline"
    print header_format.format("", *plots)
    for bin in bins:
        print r"\hline"
        print row_format.format(bin, *(fraction[plot][bin]*100 for plot in plots))
    print r"\hline"
    print r"\end{tabular}"
    print r"\end{center}"

def makeJECTable(massbins, *plots, output='makeJEC.root'):

    outfile = ROOT.TFile(output, "RECREATE") 
    print massbins
    bins = [Bin(massbins[i], massbins[i+1]) for i in range(len(massbins)-1)]
    for a in bins: print a

    eff = collections.OrderedDict()
    x = {}
    y = {}
    ex = {}
    ey = {}
    nbins = {}
    g = {}
    mg = ROOT.TMultiGraph()
    mg.SetName('mg')

    eff_unc_up = collections.OrderedDict()
    eff_unc_dn = collections.OrderedDict()
    jec_ex_up = {}
    jec_ex_dn = {}
    jec_ey_up = {}
    jec_ey_dn = {}
    jec_g = {}
    jec_mg = ROOT.TMultiGraph()
    jec_mg.SetName('mg_jec')

    legend = ROOT.TLegend(0.6, 0.4, 0.9, 0.8)
    legend.SetName("legend")
    legend.SetLineStyle(0)
    legend.SetLineColor(0)
    legend.SetFillStyle(0)

    for plot in plots:
        x[plot] = array.array("d")
        y[plot] = array.array("d")
        ex[plot] = array.array("d")
        ey[plot] = array.array("d")
        nbins[plot] = 0
        eff[plot] = collections.OrderedDict()
        eff_unc_up[plot] = collections.OrderedDict()
        eff_unc_dn[plot] = collections.OrderedDict()
        jec_ex_up[plot] = array.array("d")
        jec_ex_dn[plot] = array.array("d")
        jec_ey_up[plot] = array.array("d")
        jec_ey_dn[plot] = array.array("d")

        eff_jec = plot.vbf_eff_jec(bins)
        for bin in bins:
            eff[plot][bin] = eff_jec[bin][0] # eff at jec center
            try:
                eff_unc_up[plot][bin] = abs(eff_jec[bin][1]-eff_jec[bin][0])/eff_jec[bin][0] # eff jec unc up
                eff_unc_dn[plot][bin] = abs(eff_jec[bin][2]-eff_jec[bin][0])/eff_jec[bin][0] # eff jec unc dn
            except ZeroDivisionError:
                eff_unc_up[plot][bin] = 0.0
                eff_unc_dn[plot][bin] = 0.0
                pass
            if plot.title == "ttH" and bin.min >= 500: continue
            nbins[plot] += 1
            x[plot].append(bin.center)
            y[plot].append(eff[plot][bin])
            ex[plot].append(bin.error)
            ey[plot].append(eff_jec[bin][3]) # eff statistical error at jec center
            jec_ex_up[plot].append(bin.error)
            jec_ex_dn[plot].append(bin.error)
            jec_ey_up[plot].append(abs(eff_jec[bin][1]-eff_jec[bin][0]))
            jec_ey_dn[plot].append(abs(eff_jec[bin][2]-eff_jec[bin][0]))
 
        g[plot] = ROOT.TGraphErrors(nbins[plot], x[plot], y[plot], ex[plot], ey[plot])
        g[plot].SetName('gr_'+plot.title)
        mg.Add(g[plot])
        g[plot].SetLineColor(plot.color)
        g[plot].SetMarkerColor(plot.color)
        legend.AddEntry(g[plot], plot.title, "lp")
        if plot.title != "Z+X" : 
            jec_g[plot] = ROOT.TGraphAsymmErrors(nbins[plot], x[plot], y[plot], jec_ex_dn[plot], jec_ex_up[plot], jec_ey_dn[plot], jec_ey_up[plot])
            jec_g[plot].SetName('gr_jec_'+plot.title)
            jec_g[plot].SetMarkerColorAlpha(plot.color, 0)
            jec_g[plot].SetLineColorAlpha(plot.color, 0)
            jec_g[plot].SetFillColor(plot.color)
            jec_g[plot].SetFillStyle(plot.style)
            jec_mg.Add(jec_g[plot])

    legend.AddEntry(jec_g[plots[0]],"JEC uncertainties", "f")

    c1 = ROOT.TCanvas('c1','c1')
    pt =ROOT.TPaveText(0.1577181,0.9562937,0.9580537,0.9947552,"brNDC")
    pt.SetName('pt_cms')
    pt.SetBorderSize(0)
    pt.SetTextAlign(12)
    pt.SetFillStyle(0)
    pt.SetTextFont(42)
    pt.SetTextSize(0.03)
    text = pt.AddText(0.15,0.3,"CMS Simulation       #sqrt{s} = 13 TeV")
    #text = pt.AddText(0.15,0.3,"CMS Preliminary")
    #text = pt.AddText(0.55,0.3,"#sqrt{s} = 13 TeV, L = 2.26 fb^{-1}")

    mg.Draw("AP")
    jec_mg.Draw("E2")
    legend.Draw()
    mg.GetXaxis().SetTitle("m_{4l}")
    mg.GetYaxis().SetTitle("VBF-tag Efficiency")
    mg.GetYaxis().SetRangeUser(0,0.45)
    pt.Draw()
    c1.Update()
    c1.SaveAs(plotloc+"VBFeffJec.png")
    c1.SaveAs(plotloc+"VBFeffJec.eps")
    c1.SaveAs(plotloc+"VBFeffJec.root")
    c1.SaveAs(plotloc+"VBFeffJec.pdf")

    # write to file
    outfile.cd()
    c1.Write()
    mg.Write()
    mg_jec.Write()
    legend.Write()
    pt.Write()
    for plot in plots:
        mg[plot].Write()
        mg_jec[plot].Write()

    c2 = ROOT.TCanvas('c2','c2')
    mg.Draw("AP")
    jec_mg.Draw("E2")
    legend.Draw()
    mg.GetXaxis().SetTitle("m_{4l}")
    mg.GetYaxis().SetTitle("VBF-tag Efficiency")
    mg.GetYaxis().SetRangeUser(0,0.45)
    mg.GetXaxis().SetRangeUser(100.0, 200.0)
    pt.Draw()
    c2.Update()
    c2.SaveAs(plotloc+"VBFeffJec100.png")
    c2.SaveAs(plotloc+"VBFeffJec100.eps")
    c2.SaveAs(plotloc+"VBFeffJec100.root")
    c2.SaveAs(plotloc+"VBFeffJec100.pdf")


    print r"%%%%%%%%%%%%%% VBF-tag Efficiency %%%%%%%%%%%%%%%%%%%%%%%"
    print r"\begin{center}"
    print r"\begin{tabular}{ |%s| }" % ("|".join("c" * (len(plots)+1)))

    #http://stackoverflow.com/a/9536084
    header_format = " & ".join(["{:>15}"] * (len(plots) + 1)) + r" \\"
    row_format = " & ".join(["{:>15}"] + [r"{:14.2f}\%"] * (len(plots))) + r"\\"
    print r"\hline"
    print header_format.format("", *plots)
    for bin in bins:
        print r"\hline"
        print row_format.format(bin, *(eff[plot][bin]*100 for plot in plots))
    print r"\hline"
    print r"\end{tabular}"
    print r"\end{center}"
    print r"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"


    print r"%%%%%%%%%%%%%% VBF-tag Efficiency with JEC %%%%%%%%%%%%%%%"
    print r"\begin{center}"
    print r"\begin{tabular}{ |%s| }" % ("|".join("c" * (len(plots)+1)))

    #http://stackoverflow.com/a/9536084
    header_format = " & ".join(["{:>15}"] * (len(plots) + 1)) + r" \\"
    row_format = " & ".join(["{:>15}"] + [r"{:14.2f}\%"] * (len(plots))) + r"\\"
    print r"\hline"
    print header_format.format("", *plots)
    for bin in bins:
        print r"\hline"
        print row_format.format(bin, *(eff[plot][bin]*100 for plot in plots))
        print row_format.format('$epsilon(JEC-up)$', *((eff[plot][bin]+eff_unc_up[plot][bin]*eff[plot][bin])*100 for plot in plots))
        print row_format.format('$epsilon(JEC-dn)$', *((eff[plot][bin]-eff_unc_dn[plot][bin]*eff[plot][bin])*100 for plot in plots))
        print row_format.format('$\delta\epsilon(JEC)/\epsilon$', *((eff_unc_up[plot][bin]+eff_unc_dn[plot][bin])/2*100 for plot in plots))
    print r"\hline"
    print r"\end{tabular}"
    print r"\end{center}"
    print r"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"

if __name__ == "__main__":
    forplot = False
    fortable = False
    doJEC = True
    if forplot:
        plots = (
                 TreePlot("VBF",  1,              "VBFH125"),
                 TreePlot("H+jj", 2,              "ggH125"),
                 TreePlot("ZH",   ROOT.kGreen-6,  "ZH125"),
                 TreePlot("WH",   3,              "WplusH125"),
                 TreePlot("ttH",  4,              "ttH125",     maindir = "root://lxcms03://data3/Higgs/160111_ggZZincomplete/"),
                 TreePlot("qqZZ", 6,              "ZZTo4l"),
                 TreePlot("ggZZ", ROOT.kViolet-1, "ggZZ2e2mu", "ggZZ2e2tau", "ggZZ2mu2tau", 
                                              #"ggZZ4e", "ggZZ4mu", "ggZZ4tau"
                                            ),
                 ZXPlot(7),
                )
        makeDjetplots(*plots)
    elif fortable or doJEC:
        plots = (

                 TreePlot("VBF",  1, 3003, 
                    "VBFH125",
#                    "VBFH125","VBFH124", "VBFH125", "VBFH126", "VBFH130", "VBFH135", "VBFH140", "VBFH155", "VBFH160", "VBFH165", "VBFH170", "VBFH175", "VBFH200", "VBFH210", "VBFH230", "VBFH250", "VBFH270", "VBFH300", "VBFH350", "VBFH400", "VBFH450", "VBFH500", "VBFH550", "VBFH600", "VBFH700", "VBFH750", "VBFH800", "VBFH900", "VBFH1000", 
                         ),
#                 TreePlot("ggH", 2, 3004, 
#                    "ggH125",
#                    "ggH115", "ggH120", "ggH124", "ggH125", "ggH126", "ggH130", "ggH135", "ggH140", "ggH145", "ggH150", "ggH155", "ggH160", "ggH165", "ggH170", "ggH175", "ggH180", "ggH190", "ggH210", "ggH230", "ggH250", "ggH270", "ggH300", "ggH350", "ggH400", "ggH450", "ggH500", "ggH550", "ggH600", "ggH700", "ggH800", "ggH900", "ggH1000",
#                          ),
#                 TreePlot("ZH",   ROOT.kGreen-6, 3005, 
#                   "ZH125", 
#                    "ZH120", "ZH124", "ZH125", "ZH145", "ZH150", "ZH165", "ZH180", "ZH200", "ZH300", "ZH400", 
#                           ),
#                 TreePlot("WH",   3,  3011,  
#                  "WplusH125", 
#                    "WplusH115", "WplusH120","WplusH125", "WplusH130", "WplusH135", "WplusH140", "WplusH145", "WplusH150", "WplusH155", "WplusH160", "WplusH165", "WplusH175", "WplusH180", "WplusH190", "WplusH210", "WplusH230", "WplusH250", "WplusH270", "WplusH300", "WplusH350", "WplusH400", "WminusH115", "WminusH120", "WminusH124", "WminusH125", "WminusH126", "WminusH130", "WminusH135", "WminusH140", "WminusH145", "WminusH150", "WminusH155", "WminusH160", "WminusH165", "WminusH170", "WminusH175", "WminusH180", "WminusH190", "WminusH210", "WminusH230", "WminusH250", "WminusH270", "WminusH300", "WminusH350", "WminusH400",
#                          ),
#                 TreePlot("ttH",  4, 3012,  "ttH125"),
#                 TreePlot("qqZZ", 6, 3013,  "ZZTo4l"),
#                 TreePlot("ggZZ", ROOT.kViolet-1, 3016, #"ggZZ2e2mu", 
#                    "ggZZ2e2mu", "ggZZ2e2tau", "ggZZ2mu2tau", "ggZZ4e", "ggZZ4mu", "ggZZ4tau", 
#                           ),
#                 ZXPlot(7, 3017),
                )
        if fortable:
            makeDjettable([100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 300, 400, 500, 600, 700, 800, 900, 1000], *plots)
        elif doJEC:
#            makeJECTable([100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 300, 400, 500, 600, 700, 800, 900, 1000], *plots)
            makeJECTable([100, 150], *plots, output='outjec.root')

