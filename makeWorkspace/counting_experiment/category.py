import ROOT  # type: ignore
import sys
import array
from HiggsAnalysis.CombinedLimit.ModelTools import SafeWorkspaceImporter  # type: ignore
from .bin import Bin
from .utils import *

# MAXBINS = 100


class Category:
    # This class holds a "category" which contains a bunch of channels
    # It needs to hold a combined_pdf object, a combined_dataset object and
    # the target dataset for this channel
    def __init__(
        self,
        corrname,
        catid,
        cname,  # name for the parametric variation templates
        _fin,  # TDirectory
        _fout,  # and output file
        _wspace,  # RooWorkspace (in)
        _wspace_out,  # RooWorkspace (out)
        _bins,  # just get the bins
        _varname,  # name of the variale
        _target_datasetname,  # only for initial fit values
        _control_regions,  # CRs constructed
        diag,  # a diagonalizer object
        convention="BU",  # Naming convention to use: either BU or IC
    ):
        self.GNAME = corrname
        self.cname = cname
        self.category = catid
        self.catid = catid + "_" + corrname
        # A crappy way to store canvases to be saved in the end
        self.canvases = {}
        self.histograms = []
        self.model_hist = 0
        self._fin = _fin
        self._fout = _fout

        self._wspace = _wspace
        self._wspace_out = _wspace_out

        self._wspace_out._safe_import = SafeWorkspaceImporter(self._wspace_out)
        self._wspace._safe_import = SafeWorkspaceImporter(self._wspace)

        # self.diag = diag
        self.additional_vars = {}
        self.additional_targets = []

        self.channels = []
        self.all_hists = []
        self.cr_prefit_hists = []
        # Setup a bunch of the attributes for this category
        self._var = _wspace.var(_varname)
        self._varname = _varname
        self._bins = _bins[:]
        self._control_regions = _control_regions
        # self._data_mc  = _wspace.data(_target_datasetname)
        self._target_datasetname = _target_datasetname
        self.sample = self._wspace_out.cat("bin_number")
        self._obsvar = self._wspace_out.var("observed")
        # self._obsdata = self._wspace_out.data("combinedData")
        if self._wspace_out.var(self._var.GetName()):
            a = 1

        else:
            self._wspace_out._import(self._var, ROOT.RooFit.RecycleConflictNodes())
        self._var = self._wspace_out.var(self._var.GetName())
        self.isSecondDependant = False
        # for j,cr in enumerate(self._control_regions):
        # for i,bl in enumerate(self._bins):
        #  if i >= len(self._bins)-1 : continue
        #  self.sample.defineType("cat_%s_ch_%s_bin_%d"%(self.catid,j,i),10*MAXBINS*catid+MAXBINS*j+i)
        #  self.sample.setIndex(10*MAXBINS*catid+MAXBINS*j+i)

        self.convention = convention

    def setDependant(self, BASE, CONTROL):
        self.isSecondDependant = True
        self.BASE = BASE
        self.CONTROL = CONTROL

    def addTarget(self, vn, CR, correct=True):
        self.additional_targets.append([vn, CR, correct])

    def addVar(self, vnam, n, xmin, xmax):
        self.additional_vars[vnam] = [n, xmin, xmax]

    def fillExpectedHist(self, cr, expected_hist):
        bc = 0
        for i, ch in enumerate(self.channels):
            if ch.chid == cr.chid:
                bc += 1
                expected_hist.SetBinContent(bc, ch.ret_expected())
                expected_hist.SetBinError(bc, ch.ret_err())

    def fillExpectedCorr(self, cr, expected_hist, regen=False):
        bc = 0
        for i, ch in enumerate(self.channels):
            if ch.chid == cr.chid:
                bc += 1
                prefitValue = ch.initE_precorr if regen else ch.initE - ch.initB
                expected_hist.SetBinContent(bc, (ch.ret_expected() - ch.ret_background()) / (prefitValue))
                expected_hist.SetBinError(bc, ch.ret_err() / (ch.initE - ch.initB))

    def fillObservedHist(self, cr, observed_hist):
        bc = 0
        for i, ch in enumerate(self.channels):
            if ch.chid == cr.chid:
                bc += 1
                observed_hist.SetBinContent(bc, ch.ret_observed())
                observed_hist.SetBinError(bc, (ch.ret_observed()) ** 0.5)

    def fillBackgroundHist(self, cr, background_hist):
        bc = 0
        for i, ch in enumerate(self.channels):
            if ch.chid == cr.chid:
                bc += 1
                background_hist.SetBinContent(bc, ch.ret_background())

    def fillModelHist(self, model_hist):
        for i, ch in enumerate(self.channels):
            if i >= len(self._bins) - 1:
                break
            model_hist.SetBinContent(i + 1, ch.ret_model())

    def makeWeightHists(self, cr_i=-1, regen=False):
        hist = ROOT.TH1F("control_Region_weights", "Expected Post-fit/Pre-fit", len(self._bins) - 1, array.array("d", self._bins))
        if cr_i == -1:
            for i, ch in enumerate(self.channels):
                if i >= len(self._bins) - 1:
                    break
                hist.SetBinContent(i + 1, ch.ret_correction())
                hist.SetBinError(i + 1, ch.ret_correction_err())
        elif cr_i == -2:  # no correction
            for i, ch in enumerate(self.channels):
                if i >= len(self._bins) - 1:
                    break
                hist.SetBinContent(i + 1, 1)
                hist.SetBinError(i + 1, 0)
        else:
            self.fillExpectedCorr(self._control_regions[cr_i], hist, regen)

        return hist.Clone()

    def init_channels(self):
        # print "self._wspace_out.Print(V)", self._wspace_out.Print("V")
        # ROOT.RooCategory("bin_number","bin_number")
        sample = self._wspace_out.cat("bin_number")
        # print "zeynep sample", sample, self._wspace_out.cat("bin_number")

        # for j,cr in enumerate(self._control_regions):
        for j, cr in enumerate(self._control_regions):
            for i, bl in enumerate(self._bins):
                if i >= len(self._bins) - 1:
                    continue
                xmin, xmax = bl, self._bins[i + 1]
                if i == len(self._bins) - 2:
                    xmax = 999999.0
                ch = Bin(self.category, self.catid, cr.chid, i, self._var, "", self._wspace, self._wspace_out, xmin, xmax, convention=self.convention)
                ch.set_control_region(cr)
                if cr.has_background():
                    ch.add_background(cr.ret_background())
                ch.set_label(sample)  # should import the sample category label
                ch.set_initY(self._target_datasetname)
                ch.set_sfactor(cr.ret_sfactor(i))
                # This has to the the last thing
                # Note, we can have an expected value which is itself a RooFormulaVar

                # model_mu added to ws somewhere setup_expect_var (in else block)
                if self.isSecondDependant:
                    ch.setup_expect_var(f"cat_{self.category}_{self.BASE}_ch_{self.CONTROL}")
                else:
                    ch.setup_expect_var()

                # initialise expected  (but this will be somewhat a "post" state), i.e after fiddling with the nuisance parameters.
                ch.set_initE()
                ch.add_to_dataset()
                self.channels.append(ch)
        # fit is buggered so need to scale by 1.1

        for j, cr in enumerate(self._control_regions):
            # save the prefit histos
            cr_pre_hist = ROOT.TH1F(
                f"control_region_{cr.ret_name()}", f"Expected {cr.ret_name()} control region", len(self._bins) - 1, array.array("d", self._bins)
            )
            self.fillExpectedHist(cr, cr_pre_hist)
            cr_pre_hist.SetLineWidth(2)
            cr_pre_hist.SetLineColor(ROOT.kRed)
            self.all_hists.append(cr_pre_hist.Clone())
            self.cr_prefit_hists.append(cr_pre_hist.Clone())

        for i, bl in enumerate(self.channels):
            if i >= len(self._bins) - 1:
                break
            model_mu = self._wspace_out.var(naming_convention(bl.id, bl.catid, self.convention))
            # self._wspace_out.var(model_mu.GetName()).setVal(1.2*model_mu.getVal())

    def ret_control_regions(self):
        return self._control_regions

    def ret_channels(self):
        return self.channels

    def generate_systematic_templates(self, diag, npars):
        if self.model_hist == 0:
            sys.exit(
                "Error in generate_systematic_templates: cannot generate template variations before nominal model is created, first run Category.save_model() !!!! "
            )

        # First store nominal values in control regions (to make error bands)
        nominals = []
        for j, cr in enumerate(self._control_regions):
            nominal_values = []
            for i, ch in enumerate(self.channels):
                if ch.chid != cr.chid:
                    continue
                nominal_values.append(ch.ret_expected())
            nominals.append(nominal_values)

        # The parameters have changed so re-generate the templates
        # We also re-calculate the expectations in each CR to update the errors for the plotting
        leg_var = ROOT.TLegend(0.56, 0.1, 0.89, 0.91)
        leg_var.SetFillColor(0)
        leg_var.SetTextFont(42)

        # We will make a plot along the way
        canvr = ROOT.TCanvas("canv_variations_ratio")
        canv = ROOT.TCanvas("canv_variations")
        model_hist_spectrum = getNormalizedHist(self.model_hist)
        model_hist_spectrum.SetTitle("")
        model_hist_spectrum.GetXaxis().SetTitle("E_{T}^{miss} (GeV)")
        model_hist_spectrum.Draw()
        self.all_hists.append(model_hist_spectrum)

        sys_c = 0
        systrats = []

        for par in range(npars):
            hist_up = ROOT.TH1F(
                f"{self.GNAME}_combined_model_par_{par}_Up" f"combined_model par {par} Up 1 sigma - {self.cname} ",
                len(self._bins) - 1,
                array.array("d", self._bins),
            )
            hist_dn = ROOT.TH1F(
                f"{self.GNAME}_combined_model_par_{par}_Down" f"combined_model par {par} Down 1 sigma - {self.cname} ",
                len(self._bins) - 1,
                array.array("d", self._bins),
            )

            diag.setEigenset(par, 1)  # up variation
            # fillModelHist(hist_up,channels)
            histW = self.makeWeightHists()
            diag.generateWeightedTemplate(hist_up, histW, self._varname, self._varname, self._wspace.data(self._target_datasetname))

            # Also want to calculate for each control region an error per bin associated, its very easy to do, but only do it for "Up" variation and the error will symmetrize itself
            for j, cr in enumerate(self._control_regions):
                chi = 0
                for i, ch in enumerate(self.channels):
                    if ch.chid != cr.chid:
                        continue
                    derr = abs(ch.ret_expected() - nominals[j][chi])
                    ch.add_err(derr)
                    chi += 1

            # also add in signalregion the errors
            for i, ch in enumerate(self.channels):
                derr = abs(hist_up.GetBinContent(i + 1) - self.model_hist.GetBinContent(i + 1))
                ch.add_model_err(derr)
                if i > len(self._bins) - 1:
                    break

            diag.setEigenset(par, -1)  # up variation
            # fillModelHist(hist_dn,channels)
            histW = self.makeWeightHists()
            diag.generateWeightedTemplate(hist_dn, histW, self._varname, self._varname, self._wspace.data(self._target_datasetname))

            # Reset parameter values
            diag.resetPars()

            # make the plots
            canv.cd()
            hist_up.SetLineWidth(2)
            hist_dn.SetLineWidth(2)
            if sys_c + 2 == 10:
                sys_c += 1
            hist_up.SetLineColor(sys_c + 2)
            hist_dn.SetLineColor(sys_c + 2)
            hist_dn.SetLineStyle(2)

            # _fout.WriteTObject(hist_up)
            # _fout.WriteTObject(hist_dn)
            self.histograms.append(hist_up.Clone())
            self.histograms.append(hist_dn.Clone())

            hist_up = getNormalizedHist(hist_up)
            hist_dn = getNormalizedHist(hist_dn)
            self.all_hists.append(hist_up)
            self.all_hists.append(hist_dn)

            hist_up.Draw("samehist")
            hist_dn.Draw("samehist")

            flat = self.model_hist.Clone()
            hist_up_cl = hist_up.Clone()
            hist_up_cl.SetName(hist_up_cl.GetName() + "_ratio")
            hist_dn_cl = hist_dn.Clone()
            hist_dn_cl.SetName(hist_dn_cl.GetName() + "_ratio")
            hist_up_cl.Divide(model_hist_spectrum)
            hist_dn_cl.Divide(model_hist_spectrum)
            flat.Divide(self.model_hist)

            # ratio plot
            canvr.cd()
            flat.SetTitle("")
            flat.GetXaxis().SetTitle("E_{T}^{miss} (GeV)")
            # flat.GetYaxis().SetRangeUser(0.85,1.2)
            if par == 0:
                flat.Draw("hist")
            self.all_hists.append(flat)
            self.all_hists.append(hist_up_cl)
            self.all_hists.append(hist_dn_cl)
            systrats.append(hist_up_cl.Clone())
            systrats.append(hist_dn_cl.Clone())
            # hist_up_cl.Draw('histsame')
            # hist_dn_cl.Draw('histsame')
            leg_var.AddEntry(hist_up_cl, f"Parameter {par}", "L")
            sys_c += 1

        # find maximum
        maxdiff = 0
        for syst in systrats:
            max_local = max([syst.GetBinContent(b + 1) for b in range(syst.GetNbinsX())])
            # print max_local, maxdiff, syst.GetName()
            if max_local > maxdiff:
                maxdiff = max_local
        print("MaxDiff = ", maxdiff)
        maxdiff -= 1
        canvr.cd()
        dHist = ROOT.TH1F("dummy", ";E_{T}^{miss};Variation/Nominal", 1, self._bins[0], self._bins[-1])
        dHist.SetBinContent(1, 1)
        dHist.SetMaximum(1 + 1.1 * maxdiff)
        dHist.SetMinimum(1 - 1.1 * maxdiff)
        dHist.Draw("AXIS")
        for isy, syst in enumerate(systrats):
            syst.Draw("histsame")

        canv.cd()
        leg_var.Draw()
        canvr.cd()
        leg_var.Draw()
        self._fout.WriteTObject(canv)
        self._fout.WriteTObject(canvr)

    def save_model(self, diag):
        # Need to make ratio
        self.model_hist = ROOT.TH1F(f"{self.cname}_combined_model", f"combined_model - {self.cname}", len(self._bins) - 1, array.array("d", self._bins))
        # fillModelHist(model_hist,channels)

        histW = self.makeWeightHists()
        diag.generateWeightedTemplate(self.model_hist, histW, self._varname, self._varname, self._wspace.data(self._target_datasetname))
        self.model_hist.SetLineWidth(2)
        self.model_hist.SetLineColor(1)

        # _fout = ROOT.TFile("combined_model.root","RECREATE")
        # _fout.WriteTObject(self.model_hist)
        self.model_hist.SetName(f"{self.GNAME}_combined_model")
        histW.SetName(f"{self.GNAME}_correction_weights_{self.cname}")
        histW.SetLineWidth(2)
        histW.SetLineColor(4)
        self.histograms.append(histW)

    def save_all_models_internal(self, diag):
        # First we make errors for the nominal model histogram
        error_hist_F = ROOT.TH1F(f"{self.cname}_combined_model_ERRORS", "combined_model - {self.cname}", len(self._bins) - 1, array.array("d", self._bins))
        histW = self.makeWeightHists()
        histW_U = self.makeWeightHists()
        for b in range(histW_U.GetNbinsX()):
            # now its ~the default correction +1 sigma
            histW_U.SetBinContent(b + 1, histW_U.GetBinContent(b + 1) + histW_U.GetBinError(b + 1))
        diag.generateWeightedTemplate(error_hist_F, histW_U, self._varname, self._varname, self._wspace.data(self._target_datasetname))
        for b in range(error_hist_F.GetNbinsX()):
            sterr = error_hist_F.GetBinError(b + 1)
            self.model_hist.SetBinError(b + 1, (sterr**2 + (abs(error_hist_F.GetBinContent(b + 1) - self.model_hist.GetBinContent(b + 1))) ** 2) ** 0.5)

        # First, I think we want to pull in all of the "pre-fit" targets and make
        # them as the denominator, not necessary for the signal region

        for tg_v in self.additional_targets:
            tg = tg_v[0]
            cr_i = tg_v[1]
            histW = self.makeWeightHists(cr_i, True)
            histW_U = self.makeWeightHists(cr_i, True)
            histW.SetName(f"{self.GNAME}_{tg}_combined_model_WEIGHTS_CR_FORTARGET")
            self.histograms.append(histW.Clone())
            for b in range(histW_U.GetNbinsX()):
                # now its ~the default correction +1 sigma
                histW_U.SetBinContent(b + 1, histW_U.GetBinContent(b + 1) + histW_U.GetBinError(b + 1))
            model_tg = ROOT.TH1F(f"{self.GNAME}_{tg}_combined_model", f"combined_model - {self.cname}", len(self._bins) - 1, array.array("d", self._bins))
            diag.generateWeightedTemplate(model_tg, histW, self._varname, self._varname, self._wspace.data(tg))
            model_tg_errs = ROOT.TH1F(
                f"{self.GNAME}_{tg}_combined_model_ERRORS", "combined_model - {self.cname}", len(self._bins) - 1, array.array("d", self._bins)
            )
            diag.generateWeightedTemplate(model_tg_errs, histW_U, self._varname, self._varname, self._wspace.data(tg))
            # Errors are set as
            for b in range(model_tg_errs.GetNbinsX()):
                # add statisticsl part
                sterr = model_tg.GetBinError(b + 1)
                model_tg.SetBinError(b + 1, (sterr**2 + (abs(model_tg_errs.GetBinContent(b + 1) - model_tg.GetBinContent(b + 1))) ** 2) ** 0.5)
            self.histograms.append(model_tg.Clone())

        # Also make a weighted version of each other variable
        for varx in self.additional_vars.keys():
            histW = self.makeWeightHists()
            histW_U = self.makeWeightHists()
            for b in range(histW_U.GetNbinsX()):
                # now its ~the default correction +1 sigma
                histW_U.SetBinContent(b + 1, histW_U.GetBinContent(b + 1) + histW_U.GetBinError(b + 1))
            nb = self.additional_vars[varx][0]
            min = self.additional_vars[varx][1]
            max = self.additional_vars[varx][2]
            model_hist_vx = ROOT.TH1F(f"{self.GNAME}_combined_model{varx}", f"combined_model - {self.cname}", nb, min, max)
            model_hist_vx_errs = ROOT.TH1F(f"{self.GNAME}_combined_model{varx}_ERRORS", f"combined_model - {self.cname}", nb, min, max)
            diag.generateWeightedTemplate(model_hist_vx, histW, self._varname, varx, self._wspace.data(self._target_datasetname))
            diag.generateWeightedTemplate(model_hist_vx_errs, histW_U, self._varname, varx, self._wspace.data(self._target_datasetname))
            for b in range(model_hist_vx_errs.GetNbinsX()):
                sterr = model_hist_vx.GetBinError(b + 1)
                model_hist_vx.SetBinError(b + 1, (sterr**2 + (abs(model_hist_vx_errs.GetBinContent(b + 1) - model_hist_vx.GetBinContent(b + 1))) ** 2) ** 0.5)
            self.histograms.append(model_hist_vx.Clone())

            for tg_v in self.additional_targets:
                tg = tg_v[0]
                cr_i = tg_v[1]
                histW = self.makeWeightHists(cr_i, True)
                histW_U = self.makeWeightHists(cr_i, True)
                for b in range(histW_U.GetNbinsX()):
                    # now its ~the default correction +1 sigma
                    histW_U.SetBinContent(b + 1, histW_U.GetBinContent(b + 1) + histW_U.GetBinError(b + 1))
                model_hist_vx_tg = ROOT.TH1F(f"{self.GNAME}_{tg}_combined_model{varx}", "combined_model - {self.cname}", nb, min, max)
                model_hist_vx_tg_errs = ROOT.TH1F(f"{self.GNAME}_{tg}_combined_model_ERRORS{varx}", f"combined_model - {self.cname}", nb, min, max)
                diag.generateWeightedTemplate(model_hist_vx_tg, histW, self._varname, varx, self._wspace.data(tg))
                diag.generateWeightedTemplate(model_hist_vx_tg_errs, histW_U, self._varname, varx, self._wspace.data(tg))
                for b in range(model_hist_vx_tg_errs.GetNbinsX()):
                    sterr = model_hist_vx_tg.GetBinError(b + 1)
                    model_hist_vx_tg.SetBinError(
                        b + 1, (sterr**2 + (abs(model_hist_vx_tg_errs.GetBinContent(b + 1) - model_hist_vx_tg.GetBinContent(b + 1))) ** 2) ** 0.5
                    )
                self.histograms.append(model_hist_vx_tg.Clone())

    def make_post_fit_plots(self):
        c = ROOT.TCanvas(f"{self._target_datasetname}region_mc_fit_before_after")
        hist_original = ROOT.TH1F(f"{self.cname}_OriginalZvv", "", len(self._bins) - 1, array.array("d", self._bins))
        hist_post = ROOT.TH1F(f"{self.cname}_NewZvv", "", len(self._bins) - 1, array.array("d", self._bins))
        for i, ch in enumerate(self.channels):
            if i >= len(self._bins) - 1:
                break
            hist_original.SetBinContent(i + 1, ch.ret_initY())
            hist_post.SetBinContent(i + 1, ch.ret_model())
        hist_original.GetXaxis().SetTitle("fake MET (GeV)")
        hist_original.GetXaxis().SetTitle("Events/GeV")
        hist_original.SetLineWidth(2)
        hist_post.SetLineWidth(2)
        hist_original.SetLineColor(2)
        hist_post.SetLineColor(4)
        hist_original = getNormalizedHist(hist_original)
        hist_post = getNormalizedHist(hist_post)
        hist_original.Draw("hist")
        hist_post.Draw("histsame")
        self._fout.WriteTObject(c)

        lat = ROOT.TLatex()
        lat.SetNDC()
        lat.SetTextSize(0.04)
        lat.SetTextFont(42)

        # now build post fit plots in each control region with some indication of systematic variations from fit?
        for j, cr in enumerate(self._control_regions):
            c = ROOT.TCanvas(f"c_{cr.ret_name()}", "", 800, 800)
            cr_hist = ROOT.TH1F(
                f"{self.cname}_control_region_{cr.ret_name()}",
                f"Expected {cr.ret_name()} control region",
                len(self._bins) - 1,
                array.array("d", self._bins),
            )
            da_hist = ROOT.TH1F(
                f"{self.cname}_data_control_region_{cr.ret_name()}",
                f"data {cr.ret_name()} control region",
                len(self._bins) - 1,
                array.array("d", self._bins),
            )
            mc_hist = ROOT.TH1F(
                f"{self.cname}_mc_control_region_{cr.ret_name()}",
                f"Background {cr.ret_name()} control region",
                len(self._bins) - 1,
                array.array("d", self._bins),
            )
            self.fillObservedHist(cr, da_hist)
            self.fillBackgroundHist(cr, mc_hist)
            self.fillExpectedHist(cr, cr_hist)
            da_hist.SetTitle("")
            cr_hist.SetFillColor(ROOT.kBlue - 10)
            mc_hist.SetFillColor(ROOT.kGray)

            cr_hist = getNormalizedHist(cr_hist)
            da_hist = getNormalizedHist(da_hist)
            mc_hist = getNormalizedHist(mc_hist)
            pre_hist = getNormalizedHist(self.cr_prefit_hists[j])

            cr_hist.SetLineColor(ROOT.kBlue)
            cr_hist.SetLineWidth(2)
            mc_hist.SetLineColor(1)
            mc_hist.SetLineWidth(2)
            da_hist.SetMarkerColor(1)
            da_hist.SetLineColor(1)
            da_hist.SetLineWidth(2)
            da_hist.SetMarkerStyle(20)
            self.histograms.append(da_hist)
            self.histograms.append(cr_hist)
            self.histograms.append(mc_hist)
            self.histograms.append(pre_hist)

            c.cd()
            pad1 = ROOT.TPad("p1", "p1", 0, 0.28, 1, 1)
            pad1.SetBottomMargin(0.01)
            pad1.SetCanvas(c)
            pad1.Draw()
            pad1.cd()
            tlg = ROOT.TLegend(0.54, 0.53, 0.89, 0.89)
            tlg.SetFillColor(0)
            tlg.SetTextFont(42)
            tlg.AddEntry(da_hist, f"Data - {cr.ret_title()}", "PEL")
            tlg.AddEntry(cr_hist, "Expected (post-fit)", "FL")
            tlg.AddEntry(mc_hist, "Backgrounds Component", "F")
            tlg.AddEntry(pre_hist, "Expected (pre-fit)", "L")
            da_hist.GetYaxis().SetTitle("Events/GeV")
            da_hist.GetXaxis().SetTitle("fake MET (GeV)")
            da_hist.Draw("Pe")
            cr_hist.Draw("sameE2")
            cr_line = cr_hist.Clone()
            cr_line.SetFillColor(0)
            self.all_hists.append(cr_line)
            pre_hist.Draw("samehist")
            cr_line.Draw("histsame")
            mc_hist.Draw("samehist")
            da_hist.Draw("Pesame")
            tlg.Draw()
            lat.DrawLatex(0.1, 0.92, "#bf{CMS} #it{Preliminary}")
            pad1.SetLogy()
            pad1.RedrawAxis()

            # Ratio plot
            c.cd()
            pad2 = ROOT.TPad("p2", "p2", 0, 0.068, 1, 0.285)
            pad2.SetTopMargin(0.04)
            pad2.SetCanvas(c)
            pad2.Draw()
            pad2.cd()
            # Need to make sure cr hist has no errors for when we divide
            cr_hist_noerr = cr_hist.Clone()
            cr_hist_noerr.SetName(cr_hist.GetName() + "noerr")
            for b in range(cr_hist_noerr.GetNbinsX()):
                cr_hist_noerr.SetBinError(b + 1, 0)
            pre_hist_noerr = pre_hist.Clone()
            pre_hist_noerr.SetName(pre_hist.GetName() + "noerr")
            for b in range(pre_hist_noerr.GetNbinsX()):
                pre_hist_noerr.SetBinError(b + 1, 0)

            ratio = da_hist.Clone()
            ratio_pre = da_hist.Clone()
            ratio.GetYaxis().SetRangeUser(0.21, 1.79)
            ratio.Divide(cr_hist_noerr)
            ratio_pre.Divide(pre_hist_noerr)
            ratio.GetYaxis().SetTitle("Data/Bkg")
            ratio.GetYaxis().SetNdivisions(5)
            ratio.GetYaxis().SetLabelSize(0.1)
            ratio.GetYaxis().SetTitleSize(0.12)
            ratio.GetXaxis().SetTitleSize(0.085)
            ratio.GetXaxis().SetLabelSize(0.12)
            self.all_hists.append(ratio)
            self.all_hists.append(ratio_pre)
            ratio.GetXaxis().SetTitle("")
            ratio.SetLineColor(cr_hist.GetLineColor())
            ratio.SetMarkerColor(cr_hist.GetLineColor())
            ratio.Draw()
            ratio_pre.SetLineColor(pre_hist.GetLineColor())
            ratio_pre.SetMarkerColor(pre_hist.GetLineColor())
            ratio_pre.SetLineWidth(2)
            eline = ratio.Clone()
            eline.SetName(f"OneWithError_{ratio.GetName()}")
            self.all_hists.append(eline)
            for b in range(ratio.GetNbinsX()):
                eline.SetBinContent(b + 1, 1)
                if cr_hist.GetBinContent(b + 1) > 0:
                    eline.SetBinError(b + 1, cr_hist.GetBinError(b + 1) / cr_hist.GetBinContent(b + 1))
            eline.SetFillColor(ROOT.kBlue - 10)
            eline.SetLineColor(ROOT.kBlue - 10)
            eline.SetMarkerSize(0)
            eline.Draw("sameE2")
            line = ROOT.TLine(da_hist.GetXaxis().GetXmin(), 1, da_hist.GetXaxis().GetXmax(), 1)
            line.SetLineColor(1)
            line.SetLineStyle(2)
            line.SetLineWidth(2)
            line.Draw()
            ratio.Draw("same")
            ratio_pre.Draw("pelsame")
            ratio.Draw("samepel")
            self.all_hists.append(line)
            pad2.RedrawAxis()
            self._fout.WriteTObject(c)

    def save(self):
        # for canv in self.canvases.keys():
        #  self._fout.WriteTObject(self.canvases[canv])
        # finally THE model
        self._fout.WriteTObject(self.model_hist)
        print("Saving hitograms")
        for hist in self.histograms:
            print("Saving - ", hist.GetName())
            self._fout.WriteTObject(hist)
