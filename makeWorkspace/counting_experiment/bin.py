import ROOT as r
import sys
import array
import re
from HiggsAnalysis.CombinedLimit.ModelTools import *
from .utils import *

MAXBINS = 100


class Bin:
    def __init__(self, category, catid, chid, id, var, datasetname, wspace, wspace_out, xmin, xmax, convention):
        self.category = category
        self.chid = chid  # This is the thing that links two bins from different controls togeher
        self.id = id
        self.catid = catid
        # self.type_id   = 10*MAXBINS*catid+MAXBINS*chid+id

        self.convention = convention

        if self.convention == "BU":
            self.binid = "cat_%s_ch_%s_bin_%d" % (catid, chid, id)
        elif self.convention == "IC":
            self.binid = "cat_%s_ch_%s_bin%d" % (catid, chid, id + 1)

        self.wspace_out = wspace_out
        self.wspace_out._safe_import = SafeWorkspaceImporter(self.wspace_out)

        self.set_wspace(wspace)

        self.var = self.wspace_out.var(var.GetName())
        # self.dataset   = self.wspace.data(datasetname)

        self.rngename = "rnge_%s" % self.binid
        self.var.setRange(self.rngename, xmin, xmax)
        self.xmin = xmin
        self.xmax = xmax
        self.cen = (xmax + xmin) / 2

        self.initY = 0
        self.initE = 0
        self.initE_precorr = 0
        self.initB = 0
        self.binerror = 0
        self.binerror_m = 0

        # self.dataset.sumEntries("%s>=%g && %s<%g "%(var.GetName(),xmin,var.GetName(),xmax))
        self.o = 1
        # r.RooRealVar("observed","Observed Events bin",1)
        self.obs = self.wspace_out.var("observed")

        # <-------------------------- Check this is cool
        self.argset = r.RooArgSet(wspace.var(self.var.GetName()))
        self.obsargset = r.RooArgSet(self.wspace_out.var("observed"), self.wspace_out.cat("bin_number"))

        self.b = 0
        # self.constBkg = True

    def add_background(self, bkg):
        if "Purity" in bkg:
            tmp_pfunc = r.TF1("tmp_bkg_%s" % self.id, bkg.split(":")[-1])  # ?
            b = self.o * (1 - tmp_pfunc.Eval(self.cen))
            # self.constBkg = False
        else:
            bkg_set = self.wspace.data(bkg)
            # if not self.wspace_out.data(bkg): self.wspace_out._import(bkg)
            b = bkg_set.sumEntries("%s>=%g && %s<%g " % (self.var.GetName(), self.xmin, self.var.GetName(), self.xmax))

        # Now model nuisances for background
        nuisances = self.cr.ret_bkg_nuisances()
        if len(nuisances) > 0:
            prod = 0
            print("Is this really true? How many nuisance:", len(nuisances))

            if len(nuisances) > 1:
                nuis_args = r.RooArgList()
                for nuis in nuisances:
                    print("Adding Background Nuisance ", nuis)
                    # Nuisance*Scale is the model
                    # form_args = r.RooArgList(self.wspace_out.var("nuis_%s"%nuis),self.wspace_out.function("sys_function_%s_%s"%(nuis,self.binid)))
                    print("Trying to continue", self.wspace_out.function("sys_function_%s_%s" % (nuis, self.binid)).GetName())
                    print("Does it have an attribute:", self.wspace_out.function("sys_function_%s_%s" % (nuis, self.binid)).getAttribute("temp"))
                    if self.wspace_out.function("sys_function_%s_%s" % (nuis, self.binid)).getAttribute("temp"):
                        continue
                    form_args = r.RooArgList(self.wspace_out.function("sys_function_%s_%s" % (nuis, self.binid)))
                    delta_nuis = r.RooFormulaVar("delta_bkg_%s_%s" % (self.binid, nuis), "Delta Change from %s" % nuis, "1+@0", form_args)
                    self.wspace_out._import(delta_nuis, r.RooFit.RecycleConflictNodes())
                    nuis_args.add(self.wspace_out.function(delta_nuis.GetName()))
                prod = r.RooProduct("prod_background_%s" % self.binid, "Nuisance Modifier", nuis_args)
            else:
                print("Adding Background Nuisance ", nuisances[0])
                # if (self.wspace_out.function.getAttribute("temp")):
                ##  prod = r.RooFormulaVar("prod_background_%s"%self.binid,"Delta Change in Background from %s"%nuisances[0],"1",r.RooArgList())
                # else:
                prod = r.RooFormulaVar(
                    "prod_background_%s" % self.binid,
                    "Delta Change in Background from %s" % nuisances[0],
                    "1+@0",
                    r.RooArgList(self.wspace_out.function("sys_function_%s_%s" % (nuisances[0], self.binid))),
                )

            self.b = r.RooFormulaVar("background_%s" % self.binid, "Number of expected background events in %s" % self.binid, "@0*%f" % b, r.RooArgList(prod))
        else:
            self.b = r.RooFormulaVar(
                "background_%s" % self.binid, "Number of expected background events in %s" % self.binid, "@0", r.RooArgList(r.RooFit.RooConst(b))
            )
        self.wspace_out._import(self.b)
        self.b = self.wspace_out.function(self.b.GetName())

    def ret_initY(self):
        return self.initY

    def set_initY(self, mcdataset):
        print(
            "INIT Y: ",
            "%s>=%g && %s<%g" % (self.var.GetName(), self.xmin, self.var.GetName(), self.xmax),
            self.rngename,
            self.wspace,
            self.wspace.data(mcdataset),
            mcdataset,
        )
        self.initY = self.wspace.data(mcdataset).sumEntries("%s>=%g && %s<%g" % (self.var.GetName(), self.xmin, self.var.GetName(), self.xmax), self.rngename)

    def set_initE_precorr(self):
        return 0
        self.initE_precorr = (
            self.wspace_out.var(naming_convention(self.id, self.catid, self.convention)).getVal() * self.wspace_out.var(self.sfactor.GetName()).getVal()
        )

    def set_initE(self):
        return 0
        self.initE = self.ret_expected()
        self.initB = self.ret_background()
        self.set_initE_precorr()

    def set_label(self, cat):
        self.categoryname = cat.GetName()
        # self.wspace._import(cat,r.RooFit.RecycleConflictNodes())

    def set_wspace(self, w):
        self.wspace = w
        # self.wspace._import = getattr(self.wspace,"import") # workaround: import is a python keyword
        self.wspace._safe_import = SafeWorkspaceImporter(self.wspace)

    def set_sfactor(self, val):
        # print "Scale Factor for " ,self.binid,val
        if self.wspace_out.var("sfactor_%s" % self.binid):
            self.sfactor.setVal(val)
            self.wspace_out.var(self.sfactor.GetName()).setVal(val)
        else:
            self.sfactor = r.RooRealVar("sfactor_%s" % self.binid, "Scale factor for bin %s" % self.binid, val, 0.00001, 10000)
            self.sfactor.removeRange()
            self.sfactor.setConstant()
            self.wspace_out._import(self.sfactor, r.RooFit.RecycleConflictNodes())

    def setup_expect_var(self, functionalForm=""):
        print(functionalForm)
        if not len(functionalForm):
            if not self.wspace_out.var(naming_convention(self.id, self.catid, self.convention)):
                self.model_mu = r.RooRealVar(
                    naming_convention(self.id, self.catid, self.convention), "Model of N expected events in %d" % self.id, self.initY, 0, 3 * self.initY
                )
                # self.model_mu.removeMax() TODO
            else:
                self.model_mu = self.wspace_out.var(naming_convention(self.id, self.catid, self.convention))
        else:
            print("Setting up dependence!!")
            if self.convention == "BU":
                DEPENDANT = "%s_bin_%d" % (functionalForm, self.id)
            else:
                DEPENDANT = "%s_bin%d" % (functionalForm, self.id + 1)

            self.model_mu = self.wspace_out.function("pmu_%s" % (DEPENDANT))

        arglist = r.RooArgList((self.model_mu), self.wspace_out.var(self.sfactor.GetName()))

        # Multiply by each of the uncertainties in the control region, dont alter the Poisson pdf, we will add the constraint at the end. Actually we won't use this right now.
        nuisances = self.cr.ret_nuisances()
        if len(nuisances) > 0:
            prod = 0
            if len(nuisances) > 1:
                nuis_args = r.RooArgList()
                for nuis in nuisances:
                    if self.wspace_out.function("sys_function_%s_%s" % (nuis, self.binid)).getAttribute("temp"):
                        continue

                    print("Adding Nuisance ", nuis)
                    # Nuisance*Scale is the model
                    # form_args = r.RooArgList(self.wspace_out.var("nuis_%s"%nuis),self.wspace_out.function("sys_function_%s_%s"%(nuis,self.binid)))
                    form_args = r.RooArgList(self.wspace_out.function("sys_function_%s_%s" % (nuis, self.binid)))
                    delta_nuis = r.RooFormulaVar("delta_%s_%s" % (self.binid, nuis), "Delta Change from %s" % nuis, "1+@0", form_args)
                    self.wspace_out._import(delta_nuis, r.RooFit.RecycleConflictNodes())
                    nuis_args.add(self.wspace_out.function(delta_nuis.GetName()))

                prod = r.RooProduct("prod_%s" % self.binid, "Nuisance Modifier", nuis_args)
            else:
                print("Adding Nuisance ", nuisances[0])
                prod = r.RooFormulaVar(
                    "prod_%s" % self.binid,
                    "Delta Change from %s" % nuisances[0],
                    "1+@0",
                    r.RooArgList(self.wspace_out.function("sys_function_%s_%s" % (nuisances[0], self.binid))),
                )
            arglist.add(prod)
            self.pure_mu = r.RooFormulaVar("pmu_%s" % self.binid, "Number of expected (signal) events in %s" % self.binid, "(@0*@1)*@2", arglist)
        else:
            self.pure_mu = r.RooFormulaVar("pmu_%s" % self.binid, "Number of expected (signal) events in %s" % self.binid, "(@0*@1)", arglist)
        # Finally we add in the background
        bkgArgList = r.RooArgList(self.pure_mu)
        # if self.constBkg: self.mu = r.RooFormulaVar("mu_%s"%self.binid,"Number of expected events in %s"%self.binid,"%f+@0"%self.b,bkgArgList)
        # else : self.mu = r.RooFormulaVar("mu_%s"%self.binid,"Number of expected events in %s"%self.binid,"@0/%f"%self.b,bkgArgList)
        self.mu = r.RooFormulaVar("mu_%s" % self.binid, "Number of expected events in %s" % self.binid, "@0", bkgArgList)

        # self.mu = r.RooFormulaVar("mu_%s"%self.binid,"Number of expected events in %s"%self.binid,"@0/(@1*@2)",r.RooArgList(self.integral,self.sfactor,self.pdfFullInt))
        self.wspace_out._import(self.mu, r.RooFit.RecycleConflictNodes())
        self.wspace_out._import(self.obs, r.RooFit.RecycleConflictNodes())
        self.wspace_out.factory("Poisson::pdf_%s(observed,mu_%s)" % (self.binid, self.binid))

    def add_to_dataset(self):
        return
        # create a dataset called observed
        # self.wspace_out.var("observed").setVal(self.o)
        # self.wspace_out.cat(self.categoryname).setIndex(self.type_id)
        lv = self.wspace_out.var("observed")
        lc = self.wspace_out.cat("bin_number")
        local_obsargset = r.RooArgSet(lv, lc)
        if not self.wspace_out.data("combinedData"):
            obsdata = r.RooDataSet("combinedData", "Data in all Bins", local_obsargset)
            self.wspace_out._import(obsdata)
        obsdata = self.wspace_out.data("combinedData")
        obsdata.addFast(local_obsargset)

    def set_control_region(self, control):
        self.cr = control

    def ret_binid(self):
        return self.binid

    def ret_observed_dset(self):
        return self.wspace_out.data(dsname)

    def ret_observed(self):
        return self.o

    def ret_err(self):
        return self.binerror

    def add_err(self, e):
        self.binerror = (self.binerror**2 + e**2) ** 0.5

    def add_model_err(self, e):
        self.binerror_m = (self.binerror_m**2 + e**2) ** 0.5

    def ret_expected(self):
        return self.wspace_out.function(self.mu.GetName()).getVal()

    def ret_expected_err(self):
        return self.wspace_out.function(self.mu.GetName()).getError()

    def ret_model_err(self):
        return self.binerror_m

    def ret_background(self):
        # if self.constBkg: return self.b
        # else: return (1-self.b)*(self.ret_expected())
        return 0  # self.wspace_out.function(self.b.GetName()).getVal()

    def ret_correction(self):
        return (self.wspace_out.var(self.model_mu.GetName()).getVal()) / self.initY

    def ret_correction_err(self):
        return self.ret_model_err() / self.initY

    def ret_model(self):
        return self.wspace_out.var(self.model_mu.GetName()).getVal()

    def Print(self):
        print(
            "Channel/Bin -> ",
            self.chid,
            self.binid,
            ", Var -> ",
            self.var.GetName(),
            ", Range -> ",
            self.xmin,
            self.xmax,
            "MODEL MU (prefit/current state)= ",
            self.initY,
            "/",
            self.ret_model(),
        )
        print(
            " .... observed = ",
            self.o,
            ", expected = ",
            self.wspace_out.function(self.mu.GetName()).getVal(),
            " (of which %f is background)" % self.ret_background(),
            ", scale factor = ",
            self.wspace_out.function(self.sfactor.GetName()).getVal(),
        )
        print(", Pre-corrections (nuisance at 0) expected (-bkg) ", self.initE_precorr)
