#!/usr/bin/env python3

import os
import re
import pdb  # TODO
import argparse
from math import sqrt
from typing import Optional
from collections import defaultdict
from collections.abc import Callable
import ROOT  # type: ignore
from HiggsAnalysis.CombinedLimit.ModelTools import *  # type: ignore

from utils.generic.general import is_MC_bkg
from utils.generic.logger import initialize_colorized_logger
from utils.generic.colors import green

logger = initialize_colorized_logger(log_level="INFO")

# Load the Combine library (required for RooWorkspace manipulation)
ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit")


def parse_args():
    """Parse and validate command-line arguments."""
    parser = argparse.ArgumentParser(description="Convert input histograms from a ROOT file into a RooWorkspace.")
    parser.add_argument("--input_filename", type=str, required=True, help="Path to the input ROOT file (must end with .root).")
    parser.add_argument("--output_filename", type=str, required=True, help="Path to save the output RooWorkspace ROOT file.")
    parser.add_argument("--category", type=str, required=True, help="Analysis category (e.g., 'vbf_2017').")
    parser.add_argument("--variable", type=str, default=None, help="Variable name to extract (default: 'mjj' for VBF, otherwise 'met').")
    parser.add_argument("--root_folder", type=str, default=None, help="Optional folder path inside the input ROOT file.")

    args = parser.parse_args()

    # Validation
    if not os.path.isfile(args.input_filename):
        logger.critical(f"Input file not found: {args.input_filename}", exception_cls=IOError)

    if not args.input_filename.endswith(".root"):
        logger.critical(f"Input file must be a ROOT file: {args.input_filename}", exception_cls=IOError)

    return args


def ensure_nonzero_integral(hist: ROOT.TH1) -> None:
    """Ensure the histogram has a non-zero integral.

    This is required by the Combine framework to avoid invalid likelihoods.
    If the histogram is empty, it sets the first bin to a small positive value.
    """
    if hist.Integral() <= 0:
        hist.SetBinContent(1, 1e-4)


def merge_overflow_into_last_bin(hist: ROOT.TH1) -> None:
    """Move the overflow content into the last visible bin of the histogram."""
    n_bins = hist.GetNbinsX()
    overflow_content = hist.GetBinContent(n_bins + 1)
    overflow_error = hist.GetBinError(n_bins + 1)
    if overflow_content == 0 and overflow_error == 0:
        return

    pdb.set_trace()
    last_bin_content = hist.GetBinContent(n_bins)
    last_bin_error = hist.GetBinError(n_bins)
    new_content = last_bin_content + overflow_content
    new_error = sqrt(last_bin_error**2 + overflow_error**2)
    hist.SetBinContent(n_bins, new_content)
    hist.SetBinError(n_bins, new_error)
    hist.SetBinContent(n_bins + 1, 0)
    hist.SetBinError(n_bins + 1, 0)


def multiply_histogram_by_function(histogram: ROOT.TH1, function: Callable[[float], float]) -> None:
    """Apply a scaling function to the content and error of each bin in a histogram."""
    for bin_idx in range(1, histogram.GetNbinsX() + 1):
        bin_center = histogram.GetBinCenter(bin_idx)
        scale = function(bin_center)
        content = histogram.GetBinContent(bin_idx)
        error = histogram.GetBinError(bin_idx)
        histogram.SetBinContent(bin_idx, content * scale)
        histogram.SetBinError(bin_idx, error * scale)


def get_jes_file(category: str) -> ROOT.TFile:
    """Get the JES shape uncertainty ROOT file."""
    jer_suffix = "jer_smeared" if "vbf" in category else "not_jer_smeared"
    jes_files = {
        r"(monoj|monov).*": f"inputs/sys/{category}/monoj_monov_shape_jes_uncs_smooth_{jer_suffix}.root",
        r"vbf.*": f"inputs/sys/{category}/vbf_shape_jes_uncs_smooth_{jer_suffix}.root",
    }
    for pattern, filepath in jes_files.items():
        if re.match(pattern, category):
            jes_file = ROOT.TFile(filepath, "READ")
            logger.debug(f"Using JES/JER uncertainty file: {jes_file.GetName()}")
            return jes_file
    logger.critical(f"No JES file found for category: {category}", exception_cls=RuntimeError)


def get_jes_variations(hist: ROOT.TH1, category: str) -> dict[str, ROOT.TH1]:
    """Build varied histograms by applying JES variations from an external file."""
    f_jes = get_jes_file(category)
    keys = [key.GetName() for key in f_jes.GetListOfKeys()]
    f_jes.Close()

    # Use QCD Z(vv) shapes from the source file
    # Save varied histograms for all JES variations and the histogram names in this dictionary
    variations: dict[str, ROOT.TH1] = {}
    if "vbf" in category:
        tag = "ZJetsToNuNu"
        is_valid = lambda k: tag in k and "jesTotal" not in k
        prefix = f"{tag}20\\d\\d_"
    else:
        channel = "monov" if "monov" in category else "monojet"
        is_valid = lambda k: k.startswith(channel)
        prefix = f"{channel}_20\\d\\d_"

    for key in filter(is_valid, keys):
        variation = re.sub(prefix, "", key)
        variation_name = f"{hist.GetName()}_{variation}"

        # Skip invalid combinations
        if ("up" in variation_name and "Down" in variation) or ("down" in variation_name and "Up" in variation):
            continue

        # Clean up duplicated suffix
        variation_name = variation_name.replace("_jec_Total_up", "").replace("_jec_Total_down", "")
        varied_hist = hist.Clone(variation_name)
        varied_hist.SetDirectory(0)
        # TODO: old method. multiply by JES factor to get the varied yields
        # varied_hist.Multiply(f_jes.Get(key))
        variations[variation_name] = varied_hist

    return variations


def get_photon_id_variations(hist: ROOT.TH1, category: str) -> dict[str, ROOT.TH1]:
    """Get photon ID variations from file, returns all the varied histograms stored in a dictionary."""
    match = re.match(r".*(201[6-8]).*", category)
    if not match:
        logger.critical(f"Could not extract year from category: {category}", exception_cls=ValueError)
    year = match.group(1)

    name_map = {
        f"CMS_eff{year}_phoUp": f"monojet_{year}_photon_id_up",
        f"CMS_eff{year}_phoDown": f"monojet_{year}_photon_id_dn",
        f"CMS_eff{year}_pho_extrapUp": f"monojet_{year}_photon_id_extrap_up",
        f"CMS_eff{year}_pho_extrapDown": f"monojet_{year}_photon_id_extrap_dn",
    }

    variations = {}
    f_pho = ROOT.TFile(f"inputs/sys/{category}/photon_id_unc.root", "READ")
    for variation, histo_name in name_map.items():
        variation_name = f"{hist.GetName()}_{variation}"
        varied_hist = hist.Clone(variation_name)
        varied_hist.Multiply(f_pho.Get(histo_name))
        varied_hist.SetDirectory(0)
        variations[variation_name] = varied_hist
    f_pho.Close()
    return variations


def get_photon_qcd_variations(hist: ROOT.TH1, category: str) -> dict[str, ROOT.TH1]:
    """Create photon QCD purity fit variations based on category/year.

    Args:
        hist (ROOT.TH1): The central histogram to vary.
        category (str): Category string containing the year (2016/17/18).

    Returns:
        dict[str, ROOT.TH1]: Dictionary with Up and Down varied histograms.
    """
    match = re.match(r".*(201[6-8]).*", category)
    if not match:
        logger.critical(f"Could not extract year from category: {category}", exception_cls=ValueError)
    year = match.group(1)

    unc = 1.05 if year == "2018" else 1.10
    tag = f"purity_fit_{year}"
    func_up = lambda x: 1 + (unc - 1) / 550 * (x - 250)
    func_dn = lambda x: 1 - (unc - 1) / 550 * (x - 250)

    variations: dict[str, ROOT.TH1] = {}
    for direction, func in [("Up", func_up), ("Down", func_dn)]:
        variation_name = f"{hist.GetName()}_{tag}{direction}"
        varied_hist = hist.Clone(variation_name)
        varied_hist.SetDirectory(0)
        multiply_histogram_by_function(histogram=varied_hist, function=func)
        variations[variation_name] = varied_hist

    return variations


def get_diboson_variations(hist: ROOT.TH1, category: str, process: str) -> dict[str, ROOT.TH1]:
    """Apply shape variations from a diboson systematic file."""
    channel = re.sub(r"(loose|tight|_201\d)", "", category)
    variations: dict[str, ROOT.TH1] = {}
    shape_file = ROOT.TFile(f"inputs/sys/{category}/shape_diboson_unc.root", "READ")
    for key in shape_file.GetListOfKeys():
        key_name = key.GetName()
        if process not in key_name or channel not in key_name:
            continue
        variation_label = key_name.replace(f"{channel}_{process}_", "")
        variation_name = f"{hist.GetName()}_{variation_label}"
        varied_hist = hist.Clone(variation_name)
        varied_hist.SetDirectory(0)
        varied_hist.Multiply(shape_file.Get(key_name))
        variations[variation_name] = varied_hist
    shape_file.Close()
    return variations


def get_signal_theory_variations(hist: ROOT.TH1, category: str) -> dict[str, ROOT.TH1]:
    """Retrieve signal theory systematic variations from external file."""
    central_name = hist.GetName()
    if not central_name.startswith("signal_"):
        logger.debug(f"Skipping non-signal histogram: {central_name}")
        return {}

    channel = re.sub(r"(loose|tight|_201\d)", "", category)
    # Skip vbf for now (TODO)
    if channel == "vbf":
        return {}

    real_process = central_name.replace("signal_", "")
    # Match process type
    pattern_map = [
        (r"(vbf|ggh|ggzh|zh|wh)(\d+)?", lambda m: m.group(1)),
        (r"(vector|axial|pseudoscalar|scalar)_monow_.*", lambda _: "wh"),
        (r"(vector|axial|pseudoscalar|scalar)_(monoz|monov)_.*", lambda _: "zh"),
        (r"(vector|axial|pseudoscalar|scalar)_monojet_.*", lambda _: "ggh"),
        (r"add_md\\d+_d\\d", lambda _: "ggh"),
        (r"lq_m\\d+_d[\\d,p]+", lambda _: "ggh"),
        (r".*S3D.*", lambda _: "ggh"),
        (r".*svj.*", lambda _: "ggh"),
    ]
    process_for_unc = None
    for pattern, extractor in pattern_map:
        match = re.match(pattern, real_process)
        if match:
            process_for_unc = extractor(match)
            break

    if not process_for_unc:
        logger.warning(f"No signal theory mapping for process: {real_process}")
        return {}

    theory_file = ROOT.TFile(f"inputs/sys/{category}/signal_theory_unc.root", "READ")
    variations: dict[str, ROOT.TH1] = {}
    unc_types = ["pdf", "scale"]
    directions = ["Up", "Down"]
    for unc_type in unc_types:
        for direction in directions:
            if unc_type == "scale":
                hist_name = f"signal_{real_process}_QCDscale_{real_process}_ACCEPT{direction}"
            elif unc_type == "pdf":
                hist_name = f"signal_{real_process}_pdf_{real_process}_ACCEPT{direction}"

            variation_key = f"{channel}_{process_for_unc}_{unc_type}{direction}"
            variation_hist = theory_file.Get(variation_key)
            logger.info(variation_key, variation_hist)
            if not variation_hist:
                logger.warning(f"Missing variation histogram: {variation_key}")
                continue

            varied_hist = hist.Clone(hist_name)
            varied_hist.SetDirectory(0)
            varied_hist.Multiply(variation_hist)
            variations[hist_name] = varied_hist

    theory_file.Close()
    return variations


def add_histograms(histograms: list[ROOT.TH1], new_name: str) -> ROOT.TH1:
    """Add a list of histograms into a new histogram with a given name."""
    if not histograms:
        logger.critical("Empty list provided to add_histograms", exception_cls=ValueError)

    summed_hist = histograms[0].Clone(new_name)
    summed_hist.SetDirectory(0)

    for hist in histograms[1:]:
        summed_hist.Add(hist)

    return summed_hist


def get_mergedMC_stat_variations(to_merge_mc_bkgs: dict[str, list[ROOT.TH1]], category: str) -> dict[str, ROOT.TH1]:
    """Create autoMCstats-like per-bin statistical variation histograms for merged MC backgrounds.

    TODO: method not understood yet.
    Args:
        to_merge_mc_bkgs (dict): A mapping of region name to a list of MC background histograms.
        category (str): Analysis category name (used for naming).

    Returns:
        dict: A dictionary mapping histogram names to ROOT.TH1 up/down variation histograms.
    """
    variations: dict[str, ROOT.TH1] = {}

    for region, hists in to_merge_mc_bkgs.items():
        merged_name = f"{region}_mergedMCBkg"
        merged_hist = add_histograms(histograms=hists, new_name=merged_name)

        for bin_idx in range(1, merged_hist.GetNbinsX() + 1):
            bin_content = merged_hist.GetBinContent(bin_idx)
            bin_error = merged_hist.GetBinError(bin_idx)
            if bin_error <= 0:
                logger.critical(f"Please check why there is no error in the background prediction in bin {bin_idx} for {region}.", exception_cls=ValueError)

            if bin_content <= 0:
                continue  # Skip variation if central value is zero or negative

            ratio_up = 1.0 + bin_error / bin_content
            ratio_dn = max(0.0, 1.0 - bin_error / bin_content)

            variation_base = f"{merged_name}_{category}_stat_bin{bin_idx}"

            for hist in hists:
                base_name = hist.GetName()
                name_up = f"{base_name}_{variation_base}Up"
                name_dn = f"{base_name}_{variation_base}Down"

                h_up = hist.Clone(name_up)
                h_up.SetDirectory(0)
                content = hist.GetBinContent(bin_idx)
                h_up.SetBinContent(bin_idx, max(0.0, content * ratio_up))
                variations[name_up] = h_up
                h_dn = hist.Clone(name_dn)
                h_dn.SetDirectory(0)
                h_dn.SetBinContent(bin_idx, max(0.0, content * ratio_dn))
                variations[name_dn] = h_dn

    return variations


def write_histogram_to_workspace(
    hist: ROOT.TH1,
    name: str,
    category: str,
    workspace: ROOT.RooWorkspace,
    output_dir: ROOT.TDirectory,
    observable: ROOT.RooRealVar,
) -> None:
    """Convert histogram to RooDataHist and import into the workspace and ROOT output file."""
    logger.debug(f"Creating RooDataHist for {name}")
    roo_hist = ROOT.RooDataHist(name, f"DataSet - {category}, {name}", ROOT.RooArgList(observable), hist)
    msg_service = ROOT.RooMsgService.instance()
    prev_level = msg_service.globalKillBelow()
    msg_service.setGlobalKillBelow(ROOT.RooFit.WARNING)
    workspace._import(roo_hist)
    msg_service.setGlobalKillBelow(prev_level)

    # Write the individual histograms for easy transfer factor calculation later on
    hist.SetDirectory(0)
    hist.Sumw2()
    output_dir.cd()
    output_dir.WriteTObject(hist)


def write_variations_to_workspace(
    variations: dict,
    category: str,
    workspace: ROOT.RooWorkspace,
    output_dir: ROOT.TDirectory,
    observable: ROOT.RooRealVar,
) -> None:
    """Write multiple histograms from a dictionary of variations into the workspace."""
    for name, hist in variations.items():
        if not hist:
            logger.critical(f"Null histogram for {name}", exception_cls=RuntimeError)
        write_histogram_to_workspace(hist=hist, name=name, category=category, workspace=workspace, output_dir=output_dir, observable=observable)


def process_histogram(
    hist: ROOT.TH1,
    category: str,
    workspace: ROOT.RooWorkspace,
    output_dir: ROOT.TDirectory,
    observable: ROOT.RooRealVar,
    to_merge_mc_bkgs: dict[str, list[ROOT.TH1]],
) -> None:
    """Process a single histogram by applying systematic variations and importing it into the workspace."""
    name = hist.GetName()
    logger.debug(f"Processing histogram: {name}")
    common_kwargs = {"category": category, "workspace": workspace, "output_dir": output_dir, "observable": observable}

    ensure_nonzero_integral(hist=hist)
    merge_overflow_into_last_bin(hist=hist)
    write_histogram_to_workspace(hist=hist, name=name, **common_kwargs)

    if "data" in name:
        return

    # MC stat
    if is_MC_bkg(name) and not any(x in name.lower() for x in ["up", "down"]):
        # for MC-based background, merge the stat unc into single nuisance
        region = name.split("_")[0]
        to_merge_mc_bkgs[region].append(hist)

    # JES variations: Get them from the source file and save them to workspace
    if "jec_Total" in name:
        jes_variations = get_jes_variations(hist, category)
        write_variations_to_workspace(variations=jes_variations, **common_kwargs)

    return
    # TODO: import shapes for photon id
    if "gjets" in name:
        pdb.set_trace()
        photon_id_vars = get_photon_id_variations(hist, category)
        write_variations_to_workspace(variations=photon_id_vars, **common_kwargs)

    return
    if name == "gjets_qcd":
        photon_qcd_vars = get_photon_qcd_variations(hist, category)
        write_variations_to_workspace(variations=photon_qcd_vars, **common_kwargs)

    # Diboson variations
    diboson_processes = ["wz", "ww", "zz", "zgamma", "wgamma"]
    process = "_".join(name.split("_")[1:])
    if process in diboson_processes:
        pdb.set_trace()
        diboson_vars = get_diboson_variations(hist, category, process)
        write_variations_to_workspace(variations=diboson_vars, **common_kwargs)

    return
    # Signal theory variations
    signal_theory_vars = get_signal_theory_variations(hist, category)
    write_variations_to_workspace(variations=signal_theory_vars, **common_kwargs)


def create_workspace(
    input_filename: str,
    output_filename: str,
    category: str,
    variable: str,
    root_folder: Optional[str] = None,
) -> None:
    """Create a RooWorkspace and fill it with histograms from the input ROOT file.

    Args:
        input_filename (str): Path to the input ROOT file.
        output_filename (str): Path to save the RooWorkspace.
        category (str): Analysis category (e.g., "vbf_2017").
        variable (str): Observable variable (e.g., mjj, met).
        root_folder (Optional[str]): Subdirectory inside the input ROOT file.
    """
    os.makedirs(os.path.dirname(output_filename), exist_ok=True)

    input_file = ROOT.TFile(input_filename, "READ")
    output_file = ROOT.TFile(output_filename, "RECREATE")

    input_dir = input_file if root_folder is None else input_file.Get(root_folder)
    output_dir = output_file.mkdir(f"category_{category}")

    workspace = ROOT.RooWorkspace(f"wspace_{category}", f"wspace_{category}")
    workspace._safe_import = SafeWorkspaceImporter(workspace)  # type: ignore
    logger.info(green(f"Creating main observable: {variable}"))
    observable = ROOT.RooRealVar(variable, variable, 0, 10000)

    # Loop through all histograms in the input file and add them to the work space.
    logger.info(green("Adding histograms to workspace..."))
    to_merge_mc_bkgs: dict[str, list[ROOT.TH1]] = defaultdict(list)
    for key in input_dir.GetListOfKeys():
        obj = key.ReadObj()
        if not isinstance(obj, (ROOT.TH1D, ROOT.TH1F)):
            continue
        process_histogram(hist=obj, category=category, workspace=workspace, output_dir=output_dir, observable=observable, to_merge_mc_bkgs=to_merge_mc_bkgs)

    # now do the merging of MC-based bkg
    stat_variations = get_mergedMC_stat_variations(to_merge_mc_bkgs, category)
    write_variations_to_workspace(
        variations=stat_variations,
        category=category,
        workspace=workspace,
        output_dir=output_dir,
        observable=observable,
    )

    # Finalize workspace and close files
    output_dir.cd()
    output_dir.WriteTObject(workspace)
    output_dir.Write()
    output_file.Write()
    input_file.Close()
    output_file.Close()


def main() -> None:
    """Main entry point for the script."""
    args = parse_args()
    create_workspace(
        input_filename=args.input_filename,
        output_filename=args.output_filename,
        category=args.category,
        root_folder=args.root_folder,
        variable=args.variable,
    )


if __name__ == "__main__":
    main()
