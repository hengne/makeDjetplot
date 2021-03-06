import array
import collections
import os
import ROOT
import style
import makeDjetplot

#rootfile = ROOT.TFile.Open("fromUlascan/HZZ4l-DjetCutShapes.root")

class AnalyticFunction(object):
    def __init__(self, name, evalstr, low, hi):
        self.evalstr = evalstr
        self.low = low
        self.hi = hi
        self.name = name
        self.Cstr  = "if (channel == %(name)s)\n"
        self.Cstr += "{\n"
        if low is not None:
            self.Cstr += "if (m4l < %(low)f) m4l = %(low)f;\n" 
        if hi is not None:
            self.Cstr += "if (m4l > %(hi)f) m4l = %(hi)f;\n"
        self.Cstr += "return %(evalstr)s;\n"
        self.Cstr += "}\n"
        self.Cstr %= self.__dict__

def getfunction(name):

    fname = {
             "ggZZ": "MINLO_Djetcutshape",
             "H+jj": "MINLO_Djetcutshape",
             "qqZZ": "qqZZ_Djetcutshape",
             "VBF": "PhantomSig_Djetcutshape",
             "Z+X": "ZX_Djetcutshape",
            }
    fother = {
              "ZH": ROOT.TF1("ZH_Djetcutshape", "x < [2] ? [0]+[1]*[2] : (x < [3] ? [0]+[1]*x : [0]+[1]*[3])", 0, 5000),
              "WH": ROOT.TF1("WH_Djetcutshape", "x < [2] ? [0]+[1]*[2] : (x < [3] ? [0]+[1]*x : [0]+[1]*[3])", 0, 5000),
              "ttH": ROOT.TF1("ttH_Djetcutshape", "x < [2] ? [0]+[1]*[2]+[4]*[2]*[2] : (x < [3] ? [0]+[1]*x+[4]*x*x : [0]+[1]*[3]+[4]*[3]*[3])", 0, 5000),
              "VBF": ROOT.TF1("VBF_Djetcutshape", "x < [2] ? [0]+[1]*[2]+[4]*[2]*[2] : (x < [3] ? [0]+[1]*x+[4]*x*x : [0]+[1]*[3]+[4]*[3]*[3])", 0, 5000),
              "H+jj": ROOT.TF1("ggZZ_Djetcutshape", "[0]", 0, 5000),
              "qqZZ": ROOT.TF1("qqZZ_Djetcutshape", "([0]-[1]*x*TMath::Gaus((x-[2])/[3]))", 0, 5000),
              "Z+X": ROOT.TF1("ZX_Djetcutshape", "[0]", 0, 5000),
             }

    #the numbers are not really this precise, they are just copy and pasted
    fother["ZH"].SetParameters(3.149234e-02, -9.108965e-05, 100, 200)
    fother["WH"].SetParameters(3.363341e-02, -9.065518e-05, 100, 200)
    fother["ttH"].SetParameters(1.067331e-01, -2.617962e-04, 100, 500, 2.580946e-07)
    fother["VBF"].SetParameters(3.850116e-01, 8.321654e-05, 100, 1000, -1.062607e-07)
    fother["H+jj"].SetParameters(2.368989e-02, 0)
    fother["qqZZ"].SetParameters(6.54811139624252893e-03, 5.86652284998493653e-06, 2.43263229325644204e+02, 2.27247741344343623e+01)
    fother["Z+X"].SetParameters(1.736278e-02, 0)
    print fother["Z+X"].GetParameter(0)

    try:
        return fother[name]
    except KeyError:
        #f = rootfile.Get(fname[name])
        assert f
        return f

def getplotsfromcanvas(canvas):
    legend = canvas.GetListOfPrimitives().At(2)
    plots = {}
    for entry in legend.GetListOfPrimitives():
        graph = entry.GetObject()
        name = entry.GetLabel()
        plots[name] = graph
    return plots

def draw(filename):
    f = ROOT.TFile(filename)
    if not f:
        raise IOError("No file %s!" % filename)
    c = f.GetListOfKeys().At(0).ReadObj()
    if not c or type(c) != ROOT.TCanvas:
        raise IOError("no canvas in file " + filename + "!")
    multigraph = c.GetListOfPrimitives().At(1)
    legend = c.GetListOfPrimitives().At(2)

    plots = getplotsfromcanvas(c)
    functions = {}
    for name in plots:
        if name == "ggZZ":
            multigraph.RecursiveRemove(plots[name])
            for entry in legend.GetListOfPrimitives():
                if entry.GetObject() == plots[name]:
                    legend.GetListOfPrimitives().RecursiveRemove(entry)
            legend.RecursiveRemove(plots[name])
            continue
        plot = plots[name]
        f = getfunction(name)
        if f is not None:
            functions[name] = f
            f.SetLineColor(plot.GetLineColor())
            f.SetLineWidth(3)
            f.Draw("same")

    legend.Draw()
    #functions["Z+X"] = getfunction("Z+X")
    #functions["Z+X"].SetLineColor(7)
    #functions["Z+X"].SetLineWidth(3)
    #functions["Z+X"].Draw("same")


    c.SaveAs("/afs/cern.ch/user/h/hroskes/www/VBF/Djet/fits.png")
    c.SaveAs("/afs/cern.ch/user/h/hroskes/www/VBF/Djet/fits.eps")
    c.SaveAs("/afs/cern.ch/user/h/hroskes/www/VBF/Djet/fits.root")
    c.SaveAs("/afs/cern.ch/user/h/hroskes/www/VBF/Djet/fits.pdf")

if __name__ == '__main__':
    draw("/afs/cern.ch/user/h/hroskes/www/VBF/Djet/fraction.root")
