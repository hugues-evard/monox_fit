import ROOT  # type: ignore
import sys
from HiggsAnalysis.CombinedLimit.ModelTools import SafeWorkspaceImporter  # type: ignore


class Channel:
    # This class holds a "channel" which is as dumb as saying it holds a dataset and scale factors
    def __init__(self, cname, wspace, wspace_out, catid, scalefactors, convention="BU"):
        self.catid = catid
        self.chid = cname
        self.scalefactors = scalefactors
        self.chname = f"ControlRegion_{self.chid}"
        self.backgroundname = ""
        self.wspace_out = wspace_out
        self.wspace_out._safe_import = SafeWorkspaceImporter(self.wspace_out)
        # self.wspace_out = getattr(self.wspace_out, "import")
        self.set_wspace(wspace)
        self.nuisances = []
        self.bkg_nuisances = []
        self.systematics = {}
        self.crname = cname
        self.nbins = scalefactors.GetNbinsX()
        self.convention = convention

    def ret_title(self):
        return self.crname

    def add_systematic_shape(self, sys, file):
        # Unused
        sys.exit("Nothing Will Happen with add_systematic, use add_nuisance")
        sfup = f"{self.scalefactors.GetName()}_{sys}_Up"
        sfdn = f"{self.scalefactors.GetName()}_{sys}_Down"
        print(f"Looking for systematic shapes ... {sfup}, {sfdn}")
        try:
            print(file.Get(sfup).GetName())
            print(file.Get(sfdn).GetName())
        except AttributeError:
            print("Missing one of ", sfup, sfdn, " in ", file.GetName())
            print("Following is in directory ")
            file.Print()
            sys.exit()
        self.systematics[sys] = [file.Get(sfup), file.Get(sfdn)]

    def add_systematic_yield(self, syst, kappa):
        # Unused
        sys.exit("Nothing Will Happen with add_systematic, use add_nuisance")
        sfup = f"{self.scalefactors.GetName()}_{syst}_Up"
        sfdn = f"{self.scalefactors.GetName()}_{syst}_Down"
        sfup = self.scalefactors.Clone()
        sfup.SetName(sfup)
        sfdn = self.scalefactors.Clone()
        sfdn.SetName(sfdn)
        # log-normal scalefactors
        sfup.Scale(1 + kappa)
        sfdn.Scale(1.0 / (1 + kappa))
        self.systematics[syst] = [sfup, sfdn]

    def add_nuisance(self, name, size, bkg=False):
        # print "Error, Nuisance parameter model not supported fully for shape variations, dont use it!"
        if not (self.wspace_out.var(name)):
            # nuis = ROOT.RooRealVar("nuis_%s"%name,"Nuisance - %s"%name,0,-3,3);
            nuis = ROOT.RooRealVar(name, f"Nuisance - {name}", 0, -3, 3)
            nuis.setAttribute("NuisanceParameter_EXTERNAL", True)
            if bkg:
                nuis.setAttribute("BACKGROUND_NUISANCE", True)
            self.wspace_out._import(nuis)
            cont = ROOT.RooGaussian(
                f"const_{name}", f"Constraint - {name}", self.wspace_out.var(nuis.GetName()), ROOT.RooFit.RooConst(0), ROOT.RooFit.RooConst(1)
            )
            self.wspace_out._import(cont)

        # run through all of the bins in the control regions and create a function to interpolate
        for b in range(self.nbins):
            if self.convention == "BU":
                fname = f"sys_function_{name}_cat_{self.catid}_ch_{self.chid}_bin_{b}"
            else:
                fname = f"sys_function_{name}_cat_{self.catid}_ch_{self.chid}_bin{b+1}"

            func = ROOT.RooFormulaVar(fname, "Systematic Varation", f"@0*{size}", ROOT.RooArgList(self.wspace_out.var(name)))
            if not self.wspace_out.function(func.GetName()):
                self.wspace_out._import(func)
        # else
        #  nuis = self.wspace_out.var("nuis_%s"%name)
        if bkg:
            self.bkg_nuisances.append(name)
        else:
            self.nuisances.append(name)

    def add_nuisance_shape(self, name, file, setv="", functype="quadratic"):
        if not (self.wspace_out.var(name)):
            nuis = ROOT.RooRealVar(name, f"Nuisance - {name}", 0, -3, 3)
            nuis.setAttribute("NuisanceParameter_EXTERNAL", True)
            self.wspace_out._import(nuis)
            nuis_IN = ROOT.RooRealVar(f"nuis_IN_{name}", f"Constraint Mean - {name}", 0, -10, 10)
            nuis_IN.setConstant()
            self.wspace_out._import(nuis_IN)

            cont = ROOT.RooGaussian(
                f"const_{name}",
                f"Constraint - {name}",
                self.wspace_out.var(nuis.GetName()),
                self.wspace_out.var(nuis_IN.GetName()),
                ROOT.RooFit.RooConst(1),
            )
            self.wspace_out._import(cont)

        sfup = f"{self.scalefactors.GetName()}_{name}_Up"
        sfdn = f"{self.scalefactors.GetName()}_{name}_Down"
        print(f"Looking for systematic shapes ... {sfup}, {sfdn}")
        sysup, sysdn = file.Get(sfup), file.Get(sfdn)
        try:
            sysup.GetName()
            sysdn.GetName()
        except ReferenceError:
            print("Missing one of ", sfup, sfdn, " in ", file.GetName())
            print("Following is in directory ")
            file.ls()
            sys.exit()
        # Now we loop through each bin and construct a polynomial function per bin
        for b in range(self.nbins):
            # Name of the function depends on naming scheme
            if self.convention == "BU":
                fname = f"sys_function_{name}_cat_{self.catid}_ch_{self.chid}_bin_{b}"
            else:
                fname = f"sys_function_{name}_cat_{self.catid}_ch_{self.chid}_bin{b+1}"
            if functype == "quadratic":
                if self.scalefactors.GetBinContent(b + 1) == 0:
                    nsf = 0
                    vu = 0
                    vd = 0
                else:
                    nsf = 1.0 / (self.scalefactors.GetBinContent(b + 1))
                    vu = 1.0 / (sysup.GetBinContent(b + 1)) - nsf

                    if sysdn.GetBinContent(b + 1) == 0:
                        vd = 0
                    else:
                        # Note this should be <ve if down is lower, its not a bug
                        vd = 1.0 / (sysdn.GetBinContent(b + 1)) - nsf
                coeff_a = 0.5 * (vu + vd)
                coeff_b = 0.5 * (vu - vd)

                # this is now relative deviation, SF-SF_0 = func => SF = SF_0*(1+func/SF_0)
                func = ROOT.RooFormulaVar(fname, "Systematic Varation", f"({coeff_a}*@0*@0+{coeff_b}*@0)/{nsf}", ROOT.RooArgList(self.wspace_out.var(name)))

                if coeff_a == 0 and coeff_b == 0:
                    func.setAttribute("temp", True)
            elif functype == "lognorm":
                n0 = self.scalefactors.GetBinContent(b + 1)
                if n0 == 0:
                    sigma = 0
                    direction = 1
                else:
                    sfmax = max(0, sysup.GetBinContent(b + 1))
                    sfmin = max(0, sysdn.GetBinContent(b + 1))
                    sigma = 0.5 * abs(sfmax - sfmin)

                    direction = 1 if sfmax > sfmin else -1

                func = ROOT.RooFormulaVar(
                    fname,
                    "Systematic Variation",
                    "({N} * (1+{SIGMA}/{N})**({DIRECTION}*@0) - {N}) / {N}".format(N=n0, SIGMA=sigma, DIRECTION=direction),
                    ROOT.RooArgList(self.wspace_out.var(name)),
                )
                if sigma == 0:
                    func.setAttribute("temp", True)
            self.wspace_out.var(name).setVal(0)
            if not self.wspace_out.function(func.GetName()):
                self.wspace_out._import(func)
        if setv != "":
            if "SetTo" in setv:
                vv = float(setv.split("=")[1])
                self.wspace_out.var(f"nuis_IN_{name}").setVal(vv)
                self.wspace_out.var(name).setVal(vv)
            else:
                print(f"DIRECTIVE {setv} IN SYSTEMATIC {name}, NOT UNDERSTOOD!")
                sys.exit()
        self.nuisances.append(name)

    def set_wspace(self, w):
        # Only used once, in __init__
        self.wspace = w
        self.wspace._safe_import = SafeWorkspaceImporter(self.wspace)
        # self.wspace._import = getattr(self.wspace,"import") # workaround: import is a python keyword

    def ret_bkg_nuisances(self):
        return self.bkg_nuisances

    def ret_nuisances(self):
        return self.nuisances

    def ret_name(self):
        return self.chname

    def ret_chid(self):
        # Unused
        return self.chid

    def ret_sfactor(self, i, syst="", direction=1):
        if self.scalefactors.GetBinContent(i + 1) == 0:
            return 0
        if syst and syst in self.systematics.keys():
            if direction > 0:
                index = 0
            else:
                index = 1
            return 1.0 / (self.systematics[syst][index].GetBinContent(i + 1))
        else:
            return 1.0 / (self.scalefactors.GetBinContent(i + 1))

    def ret_background(self):
        return self.backgroundname

    def has_background(self):
        return len(self.backgroundname)
