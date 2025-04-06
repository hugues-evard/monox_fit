import ROOT as r
import sys
import array
import re
from HiggsAnalysis.CombinedLimit.ModelTools import *

MAXBINS = 100


def naming_convention(id, catid, convention="BU"):
    if convention == "BU":
        return "model_mu_cat_%s_bin_%d" % (catid, id)
    elif convention == "IC":
        m = re.match(".*(201\d).*", catid)
        if not m or len(m.groups()) > 1:
            raise RuntimeError("Cannot derive year from category ID: " + catid)
        year = m.groups()[0]
        return "MTR_%s_QCDZ_SR_bin%d" % (year, id + 1)
    else:
        raise RuntimeError("Unknown naming convention: " + convention)


def getNormalizedHist(hist):
    thret = hist.Clone()
    nb = hist.GetNbinsX()
    for b in range(1, nb + 1):
        sfactor = 1.0 / hist.GetBinWidth(b)
        thret.SetBinContent(b, hist.GetBinContent(b) * sfactor)
        thret.SetBinError(b, hist.GetBinError(b) * sfactor)
        # thret.GetYaxis().SetTitle("Events")
        thret.GetYaxis().SetTitle("Events/GeV")
    return thret
