#!/usr/bin/env python3

import os
import re
import argparse
import ROOT  # type: ignore
from HiggsAnalysis.CombinedLimit.ModelTools import *  # type: ignore

from utils.generic.logger import initialize_colorized_logger
from utils.workspace.convert_to_combine_workspace import convert_to_combine_workspace

logger = initialize_colorized_logger(log_level="INFO")

# Setup ROOT and external C++ dependencies
ROOT.gSystem.AddIncludePath("-I$CMSSW_BASE/src/")
ROOT.gSystem.AddIncludePath("-I$ROOFITSYS/include")
ROOT.gSystem.Load("libRooFit.so")
ROOT.gSystem.Load("libRooFitCore.so")
ROOT.gROOT.SetBatch(True)
# ROOT.gSystem.SetBuildDir("/tmp/ROOT_build", True)
ROOT.gROOT.ProcessLine(".L makeWorkspace/diagonalizer.cc+")
from ROOT import diagonalizer  # type: ignore


def parse_args() -> argparse.Namespace:
    """Parse and validate command-line arguments."""
    parser = argparse.ArgumentParser(description="Construct fit model from RooWorkspace.")
    parser.add_argument("--input_filename", type=str, required=True, help="Input file containing histograms and workspace.")
    parser.add_argument("--output_filename", type=str, default="combined_model.root", help="Path to the output ROOT file.")
    parser.add_argument("--category", type=str, required=True, help="Analysis category, e.g., 'vbf_2017'.")
    parser.add_argument("--rename", type=str, default="", help="Optional new name for the observable variable.")

    args = parser.parse_args()

    args.input_filename = os.path.abspath(args.input_filename)
    args.output_filename = os.path.abspath(args.output_filename)

    if not os.path.isfile(args.input_filename):
        logger.critical(f"Input file not found: {args.input_filename}", exception_cls=IOError)
    if not args.input_filename.endswith(".root"):
        logger.critical(f"Input file must be a ROOT file: {args.input_filename}", exception_cls=IOError)

    return args


def generate_combine_model(
    input_filename: str,
    output_filename: str,
    category: str,
    rename: str = "",
) -> None:
    """Generate a Combine RooWorkspace with control region models."""
    # Determine CR configurations based on category
    if "mono" in category:
        controlregions_def = ["Z_constraints", "W_constraints"]
    elif "vbf" in category:
        controlregions_def = [
            "Z_constraints_qcd_withphoton",
            "W_constraints_qcd",
            "Z_constraints_ewk_withphoton",
            "W_constraints_ewk",
        ]
    else:
        logger.critical("Could not infer control region definitions from category.", exception_cls=RuntimeError)

    os.makedirs(os.path.dirname(output_filename), exist_ok=True)
    input_file = ROOT.TFile.Open(input_filename)
    output_file = ROOT.TFile(output_filename, "RECREATE")

    workspace = ROOT.RooWorkspace("combinedws")
    workspace._safe_import = SafeWorkspaceImporter(workspace)  # type: ignore

    logger.info("Creating global observables")
    sample_type = ROOT.RooCategory("bin_number", "Bin Number")
    observed = ROOT.RooRealVar("observed", "Observed Events bin", 1)
    workspace._import(sample_type)  # Global variables for dataset
    workspace._import(observed)

    # Loop over control region definitions, and load their model definitions
    diag = diagonalizer(workspace)
    cmb_categories = []
    for cr_name in controlregions_def:
        module = __import__(cr_name)
        match = re.match(r".*201(6|7|8).*", category)
        if not match or len(match.groups()) != 1:
            logger.critical(f"Cannot determine year from category: {category}", exception_cls=RuntimeError)
        year = int("201" + match.group(1))

        cr_dir = output_file.mkdir(f"{cr_name}_category_{category}")
        # This `cmodel` functions is where models are created.
        # This function does two things:
        # 1) Compute transfer factors
        #   Processes are express different processes as a ratio with of QCD Znunu in SR
        #   For the qcd_zjets model, this is directly the ratio of the two processes
        #   For other models, processes are first taken as a ratio with either EWK Znunu, QCD Wjets and EWK Wjets,
        #   but these are later expressed as transfer factors to make QCD Znunu appear at the `init_channels` step
        # 2) Add nuisances and add them to the workspace
        #   for veto, JES/JER, theory and statistical uncertainties
        #   for each transfer factor in the model
        if "MTR" in rename:
            model = module.cmodel(category, cr_name, input_file, cr_dir, workspace, diag, year, convention="IC")
        else:
            model = module.cmodel(category, cr_name, input_file, cr_dir, workspace, diag, year)

        cmb_categories.append(model)

    # Initialize model channels: model_mu_cat_vbf_2017_qcd_zjets_bin_0 TODO
    for model in cmb_categories:
        logger.info(f"Initializing model channels for: {model}")
        # This is where the actual model distributions as a function of QCD Znunu in SR are made for all processes.
        # Processes modeled with `qcd_zjets` are expressed as:
        #   (process yield) * [transfer factor = (QCD Znunu in SR) / (process yield)] * Product of all nuisances
        # Processes using the other models will first fetch the above model, before multiplying the transfer factor and nuisances
        # For instance, EWK Zll in the diMuon region, modeled with `ewk_zjets` will fetch
        #   (EWK Znunu in SR) * [transfer factor = (QCD Znunu in SR) / (EWK Znunu in SR)] * Product of all nuisances (for EWK Znunu in SR)
        # and multiply by
        #   [transfer factor = (EWK Znunu in SR) / (EWK Zll in diMuon)] * Product of all nuisances (for EWK Zll in diMuon)
        # These models are made for each bin of the mjj distribution and saved to the workspace
        model.init_channels()
        # _ = model.ret_channels() TODO not used

    # Save pre-fit snapshot
    workspace.saveSnapshot("PRE_EXT_FIT_Clean", workspace.allVars())

    # Convert workspace to Combine workspace format
    # This actually builds the histograms used in the fit
    # It fetches the modeled distribution for every bin of every process computed at the above `init_channels` step
    # and stores them as the `RooParametricHist` that will be used in the datacard
    convert_to_combine_workspace(
        wsin_combine=workspace,
        f_simple_hists=input_file,
        category=category,
        cmb_categories=cmb_categories,
        controlregions_def=controlregions_def,
        rename_variable=rename,
    )

    output_file.WriteTObject(workspace)
    logger.info(f"--> Produced constraints model in Combine workspace: {output_file.GetName()}")


def main() -> None:
    """Main entry point for the script."""
    args = parse_args()
    generate_combine_model(
        input_filename=args.input_filename,
        output_filename=args.output_filename,
        category=args.category,
        rename=args.rename,
    )


if __name__ == "__main__":
    main()
