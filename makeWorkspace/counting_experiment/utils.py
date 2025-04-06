import re
from HiggsAnalysis.CombinedLimit.ModelTools import *  # type: ignore


def naming_convention(id, catid, convention="BU"):
    if convention == "BU":
        return f"model_mu_cat_{catid}_bin_{id}"
    elif convention == "IC":
        m = re.match(".*(201\d).*", catid)
        if not m or len(m.groups()) > 1:
            raise RuntimeError("Cannot derive year from category ID: " + catid)
        year = m.groups()[0]
        return f"MTR_{year}_QCDZ_SR_bin{id+1}"
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
