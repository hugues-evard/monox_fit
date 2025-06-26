import ROOT  # type:ignore
from typing import Any
from counting_experiment import Category, Channel
from utils.generic.logger import initialize_colorized_logger
from utils.workspace.jes_utils import get_jes_variations_names, get_jes_file
from utils.workspace.flat_uncertainties import get_veto_unc

logger = initialize_colorized_logger(log_level="INFO")


def define_model(
    category_id: str,
    category_name: str,
    input_file: ROOT.TFile,
    output_file: ROOT.TFile,
    output_workspace: ROOT.RooWorkspace,
    diagonalizer: Any,
    year: int,
    variable: str,
    convention: str,
    model_name: str,
    target_name: str,
    samples_map: dict[str, str],
    channel_names: dict[str, str],
    veto_channel_list: list[str],
    jes_jer_channel_list: list[str],
    jes_jer_process: str,
    theory_channel_list: list[str],
    region_names: dict[str, str],
):
    """
    Defines a statistical model for a given category using transfer factors.

    This function constructs a model by:
    - Fetching input histograms for the target and control samples.
    - Computing transfer factors (ratios of target to control samples).
    - Creating `Channel` objects for each transfer factor.
    - Applying systematic uncertainties such as veto nuisances, JES/JER variations, and theory uncertainties.
    - Adding bin-by-bin statistical uncertainties.
    - Storing results in an output ROOT file and workspace.
    - Returning a `Category` object encapsulating the model details.

    Args:
        category_id (str): Unique identifier for the category.
        category_name (str): Name of the category.
        input_file (ROOT.TFile): Input ROOT file containing histograms.
        output_file (ROOT.TFile): Output ROOT file to store results.
        output_workspace (ROOT.RooWorkspace): Output workspace for the statistical model.
        diagonalizer (Any): Object for diagonalizing correlations.
        year (int): Data-taking year.
        convention (str): Naming convention for systematic uncertainties.
        model_name (str): Name of the model.
        target_name (str): Name of the target MC sample in the input ROOT file.
        samples_map (dict[str, str]): Mapping of control MC sample names to their ROOT file entries.
        channel_names (dict[str, str]): Mapping of transfer factor labels to channel names.
        veto_channel_list (list[str]): Channels where veto uncertainties are applied.
        veto_dict (dict[str, float]): Dictionary of veto nuisance values.
        jes_jer_channel_list (list[str]): Channels where JES/JER uncertainties are applied.
        jes_jer_process (str): Process label for JES/JER uncertainties.
        theory_channel_list (list[str]): Channels where theory uncertainties are applied.
        region_names (dict[str, str]): Mapping of transfer factor labels to region names.

    Returns:
        Category: A `Category` object encapsulating the defined model.
    """

    # Some setup
    input_tdir = input_file.Get(f"category_{category_id}")
    input_wspace = input_tdir.Get(f"wspace_{category_id}")

    # Defining the nominal transfer factors
    # Nominal MC process to model
    target = input_tdir.Get(target_name)
    # Control MC samples
    control_samples = fetch_control_samples(input_tdir=input_tdir, samples_map=samples_map)

    # Compute and save a copy of the transfer factors (target divided by control)
    transfer_factors = define_transfer_factors(
        control_samples=control_samples,
        category_id=category_id,
        target_sample=target,
        output_file=output_file,
    )

    # Create a `Channel` object for each transfer factor
    CRs = define_channels(
        transfer_factors=transfer_factors,
        category_id=category_id,
        input_wspace=input_wspace,
        output_workspace=output_workspace,
        convention=convention,
        channel_names=channel_names,
        model_name=model_name,
    )

    add_veto_nuisances(
        channel_objects=CRs,
        channel_list=veto_channel_list,
        model_name=model_name,
        year=year,
    )
    add_jes_jer_uncertainties(
        transfer_factors=transfer_factors,
        channel_objects=CRs,
        channel_list=jes_jer_channel_list,
        year=year,
        category_id=category_id,
        output_file=output_file,
        process=jes_jer_process,
        production_mode=model_name.split("_")[0],
        syst_folder=f"inputs/sys/{variable}",
    )
    add_theory_uncertainties(
        control_samples=control_samples,
        transfer_factors=transfer_factors,
        target_sample=target,
        channel_objects=CRs,
        channel_list=theory_channel_list,
        year=year,
        category_id=category_id,
        output_file=output_file,
        production_mode=model_name.split("_")[0],
        syst_folder=f"inputs/sys/{variable}",
    )

    # Add Bin by bin nuisances to cover statistical uncertainties
    for sample, transfer_factor in transfer_factors.items():
        do_stat_unc(transfer_factor, proc=sample, region=region_names[sample], CR=CRs[sample], cid=category_id, outfile=output_file)

    # Extract the bin edges of the distribution
    bin_edges = [target.GetBinLowEdge(b + 1) for b in range(target.GetNbinsX() + 1)]

    # Create and return `Category` object
    return Category(
        corrname=model_name,
        catid=category_id,
        cname=category_name,
        _fin=input_tdir,
        _fout=output_file,
        _wspace=input_wspace,
        _wspace_out=output_workspace,
        _bins=bin_edges,
        _varname=variable,
        _target_datasetname=target.GetName(),
        _control_regions=list(CRs.values()),
        diag=diagonalizer,
        convention=convention,
    )


def fetch_control_samples(
    input_tdir: ROOT.TDirectoryFile,
    samples_map: dict[str, str],
) -> dict[str : ROOT.TH1]:
    """
    Fetch all needed control MC samples used by the model.

    Args:
        input_tdir (ROOT.TDirectoryFile): Directory where samples are stored
        samples_map (dict[str, str]): Dictionnary mapping the name of the region / transfer factor to the name of the MC sample in the directory.
    """

    return {region_name: input_tdir.Get(sample_name) for region_name, sample_name in samples_map.items()}


def define_transfer_factors(
    control_samples: dict[str : ROOT.TH1],
    category_id: str,
    target_sample: ROOT.TH1,
    output_file: ROOT.TFile,
) -> dict[str : ROOT.TH1]:
    """
    Compute the transfer factors for each MC sample.

    Args:
        control_samples (dict[str, ROOT.TH1]): Dictionary mapping transfer factors labels to their control MC sample.
        category_id (str): Unique identifier for the category.
        target_sample (Any): Histogram of the target process.
        output_file (ROOT.TFile): Output ROOT file for a copy of the transfer factors.
    """
    transfer_factors = {}
    for label, sample in control_samples.items():
        factor = target_sample.Clone(f"{label}_weights_{category_id}")
        factor.Divide(sample)
        # transfer_factors[label] = factor
        transfer_factors.update({label: factor})
        output_file.WriteTObject(factor)
    return transfer_factors


def define_channels(
    transfer_factors: dict[str, ROOT.TH1],
    category_id: str,
    model_name: str,
    input_wspace: ROOT.RooWorkspace,
    output_workspace: ROOT.RooWorkspace,
    convention: str,
    channel_names: dict[str, str],
) -> dict[str, Channel]:
    """
    Stores transfer factors in Channel objects.

    Args:
        transfer_factors (dict[str, ROOT.TH1]): Dictionary mapping transfer factors labels to their distributions.
        category_id (str): Unique identifier for the category.
        model_name (str): Unique identifier for the model.
        input_wspace (ROOT.RooWorkspace): The input ROOT workspace containing the model parameters.
        output_workspace (ROOT.RooWorkspace): The output ROOT workspace used for importing variables and functions.
        convention (str): A string defining the naming convention used for systematic uncertainties ("IC" or "BU").
        channel_names (dict[str, str]): Dictionary mapping transfer factors labels to the name of their channels.
    """
    return {
        sample: Channel(
            cname=channel_names[sample],
            wspace=input_wspace,
            wspace_out=output_workspace,
            catid=category_id + "_" + model_name,
            scalefactors=transfer_factor,
            convention=convention,
        )
        for sample, transfer_factor in transfer_factors.items()
    }


def add_veto_nuisances(channel_objects: dict[str, Channel], channel_list: list[str], model_name: str, year: int) -> None:
    """
    Adds veto systematic uncertainties to the specified control regions.

    Args:
        channel_objects (dict[str, Channel]): Dictionary mapping control region names to `Channel` objects.
        channel_list (list[str]): List of control regions to apply veto uncertainties.
        veto_dict (dict[str, float]): Dictionnary mapping the name of the nuissance to add and its value.
    """
    veto_dict = {f"CMS_veto_{key}_{year}": value for key, value in get_veto_unc(model=model_name).items()}
    for channel in channel_list:
        for veto_name, veto_value in veto_dict.items():
            channel_objects[channel].add_nuisance(veto_name, veto_value)


def add_jes_jer_uncertainties(
    transfer_factors: dict[str, Any],
    channel_objects: dict[str, Channel],
    channel_list: list[str],
    year: str,
    category_id: str,
    output_file: ROOT.TFile,
    process: str,
    production_mode: str,
    syst_folder: str,
) -> None:
    """
    Adds JES and JER uncertainties to transfer factors.

    This function:
    - Retrieves JES/JER uncertainty variations from an external file.
    - Applies up/down variations to transfer factors.
    - Stores the modified transfer factors in the output file.
    - Adds nuisance parameters for JES/JER uncertainties to the corresponding channels.

    Args:
        transfer_factors (dict[str, Any]): Dictionary of transfer factors.
        channel_objects (dict[str, Channel]): Dictionary of `Channel` objects.
        channel_list (list[str]): List of control regions to apply JES/JER uncertainties.
        year (str): Data-taking year.
        category_id (str): Unique identifier for the category.
        output_file (ROOT.TFile): Output ROOT file for storing variations.
        process (str): Label indicating the process being modeled, either "znunu" or "wlnu".
        production_mode (str): Label indicating the if the production_mode is strong or electroweak, either "qcd" or "ewk".
    """

    jes_region_labels = {
        "qcd_w": "wlnu",
        "qcd_zmm": "zmumu",
        "qcd_zee": "zee",
        "qcd_photon": "gjets",
        "qcd_wmn": "wmunu",
        "qcd_wen": "wenu",
        "ewk_w": "wlnu",
        "ewk_zmm": "zmumu",
        "ewk_zee": "zee",
        "ewk_photon": "gjets",
        "ewk_wmn": "wmunu",
        "ewk_wen": "wenu",
    }
    # Get the JES/JER uncertainty file for transfer factors
    # Read the split uncertainties from there
    jet_variations = get_jes_variations_names(year=year)

    for sample in channel_list:
        for var in jet_variations:
            fjes = get_jes_file(category=category_id, source=var, syst_folder=syst_folder)
            for var_direction in ["Up", "Down"]:
                # Scale transfer factor by relative variation and write to output file
                add_variation(
                    nominal=transfer_factors[sample],
                    unc_file=fjes,
                    unc_name=f"{process}_over_{jes_region_labels[sample]}_{production_mode}_{var}{var_direction}",
                    new_name=f"{sample}_weights_{category_id}_{var}_{var_direction}",
                    outfile=output_file,
                )
            # Add function (quadratic) to model the nuisance
            channel_objects[sample].add_nuisance_shape(var, output_file)


def add_theory_uncertainties(
    control_samples: dict[str, Any],
    transfer_factors: dict[str, ROOT.TH1],
    target_sample: Any,
    channel_objects: dict[str, Channel],
    channel_list: list[str],
    year: str,
    category_id: str,
    output_file: ROOT.TFile,
    production_mode: str,
    syst_folder: str,
) -> None:
    """
    Adds theoretical uncertainties (scale, PDF, and EWK corrections) to transfer factors.

    This function:
    - Saves copies of control samples used to derive theory variations.
    - Retrieves theoretical uncertainty histograms from an external file.
    - Computes and stores up/down variations for QCD scale, PDF, and EWK uncertainties.
    - Applies bin-by-bin decorrelated EWK uncertainties.
    - Adds nuisance parameters for theoretical uncertainties to the corresponding channels.

    Args:
        control_samples (dict[str, Any]): Dictionary of control region histograms.
        target_sample (Any): Histogram of the target process.
        channel_objects (dict[str, Channel]): Dictionary of `Channel` objects.
        channel_list (list[str]): List of control regions to apply theoretical uncertainties.
        year (str): Data-taking year.
        category_id (str): Unique identifier for the category.
        output_file (ROOT.TFile): Output ROOT file for storing variations.
    """

    # Save a (renamed) copy of samples used to derive theory variations
    # Done to perfectly mirrors what is done in Z_constraints_qcd_withphoton
    spectrum_label = {
        "qcd_w": "qcd_w",
        "qcd_photon": "qcd_gjets",
        "ewk_w": "ewk_w",
        "ewk_photon": "ewk_photon",
    }
    # TODO: also write Z->nunu spectrum
    spectrums = {region: control_samples[region].Clone() for region in channel_list}
    for region, sample in spectrums.items():
        sample.SetName(f"{spectrum_label[region]}_spectrum_{category_id}_")
        output_file.WriteTObject(sample)

    nbins = target_sample.GetNbinsX()

    # different labels to convert naming scheme between the different histogram and nuisances to read and write
    label_dict = {
        "qcd_w": ("zoverw", "ZnunuWJets_QCD", "qcd_ewk"),
        "qcd_photon": ("goverz", "Photon_QCD", "qcd_photon_ewk"),
    }
    if "vbf" in category_id:
        label_dict["ewk_w"] = ("zoverw", "ZnunuWJets_EWK", "ewk_ewk")
        label_dict["ewk_photon"] = ("goverz", "Photon_EWK", "ewkphoton_ewk")

    for region in channel_list:
        ratio, qcd_label, ewk_label = label_dict[region]

        # Add QCD and PDF uncertainties
        # TODO follow https://cms-analysis.docs.cern.ch/guidelines/systematics/systematics/#pdf-uncertainties
        # QCD_ren_scale_<process> QCD_fac_scale_<process>
        for var in [("mur", "renscale"), ("muf", "facscale"), ("pdf", "pdf")]:
            vbf_sys = ROOT.TFile.Open(f"{syst_folder}/{category_id}/systematics_{var[0]}.root", "READ")
            for var_direction in ["Up", "Down"]:
                add_variation(
                    nominal=transfer_factors[region],
                    unc_file=vbf_sys,
                    # unc_name=f"uncertainty_ratio_{denom_label}_mjj_unc_{ratio}_nlo_{var[0]}_{dir}_{year}",
                    unc_name=f"{ratio}_over_{region}_{production_mode}_{var[0]}{var_direction}",
                    new_name=f"{region}_weights_{category_id}_{qcd_label}_{var[1]}_vbf_{var_direction}",
                    outfile=output_file,
                )
            vbf_sys.Close()

            # Add function (quadratic) to model the nuisance
            channel_objects[region].add_nuisance_shape(f"{qcd_label}_{var[1]}_vbf", output_file)

        if "vbf" in category_id:
            # EWK uncertainty (decorrelated among bins)
            ewk_sys = ROOT.TFile.Open(f"{syst_folder}/{category_id}/systematics_pdf.root", "READ")  # TODO
            for dir in ["Up", "Down"]:
                ratio_ewk = transfer_factors[region].Clone(f"{region}_weights_{category_id}_ewk_{dir}")
                # ratio_ewk.Multiply(vbf_sys.Get(f"uncertainty_ratio_{denom_label}_mjj_unc_w_ewkcorr_overz_common_{dir}_{year}"))
                ratio_ewk.Multiply(ewk_sys.Get(f"signal_ewkzjets_over_signal_ewkwjets_pdf{dir}"))

                for b in range(nbins):
                    new_name = f"{region}_weights_{category_id}_{ewk_label}_{category_id.replace(f'_{year}', '')}_bin{b}_{dir}"
                    ewk_w = transfer_factors[region].Clone(new_name)
                    ewk_w.SetBinContent(b + 1, ratio_ewk.GetBinContent(b + 1))
                    output_file.WriteTObject(ewk_w)
            ewk_sys.Close()

            for b in range(nbins):
                # Add function (quadratic) to model the nuisance
                channel_objects[region].add_nuisance_shape(f"{ewk_label}_{category_id.replace(f'_{year}', '')}_bin{b}", output_file)


# Ported from W_constraints, WIP
def do_stat_unc(
    histogram,
    proc,
    cid,
    region,
    CR,
    outfile,
):
    """Add stat. unc. variations to the workspace"""

    # Add one variation per bin
    for b in range(1, histogram.GetNbinsX() + 1):
        err = histogram.GetBinError(b)
        content = histogram.GetBinContent(b)

        # Safety
        if (content <= 0) or (err / content < 0.001):
            logger.critical(f"Stat. unc. undefined in bin {b} of hist '{histogram.GetName()}': content = {content}, error = {err}.", exception_cls=ValueError)

        # Careful: The bin count "b" in this loop starts at 1
        # In the combine model, we want it to start from 0!
        up = histogram.Clone(f"{proc}_weights_{cid}_{cid}_stat_error_{region}_bin{b-1}_Up")
        down = histogram.Clone(f"{proc}_weights_{cid}_{cid}_stat_error_{region}_bin{b-1}_Down")
        up.SetBinContent(b, content + err)
        down.SetBinContent(b, content - err)
        outfile.WriteTObject(up)
        outfile.WriteTObject(down)

        logger.info(f"Adding statistical variation with absolute error = {err:.4f}, relative error = {err / content:.4f}: {up.GetName()}")
        CR.add_nuisance_shape(f"{cid}_stat_error_{region}_bin{b-1}", outfile, functype="lognorm")


def add_variation(
    nominal: ROOT.TH1,
    unc_file: ROOT.TFile,
    unc_name: str,
    new_name: str,
    outfile: ROOT.TFile,
) -> None:
    # TODO: remove
    unc_name = unc_name.replace("znunu_over_", "signal_qcdzjets_over_").replace("zmumu_qcd", "Zmm_qcdzll").replace("zee_qcd", "Zee_qcdzll")
    unc_name = unc_name.replace("zmumu_zjets_", "Zmm_qcdzll_").replace("zee_zjets_", "Zee_qcdzll_")
    unc_name = unc_name.replace("Zmm_qcdzll_zjets", "Zmm_qcdzll").replace("Zee_qcdzll_zjets", "Zee_qcdzll")
    unc_name = unc_name.replace("wlnu_qcd", "signal_qcdwjets").replace("gjets_qcd", "gjets_qcdgjets")
    unc_name = unc_name.replace("wlnu_over_", "signal_qcdwjets_over_").replace("wmunu_qcd", "Wmn_qcdwjets").replace("wenu_qcd", "Wen_qcdwjets")
    unc_name = unc_name.replace("wlnu_zjets", "signal_qcdwjets_zjets")
    unc_name = unc_name.replace("qcdwjets_zjets", "qcdwjets").replace("gjets_zjets", "gjets_qcdgjets")
    unc_name = unc_name.replace("wmunu_wjets", "Wmn_qcdwjets").replace("wenu_wjets", "Wen_qcdwjets")
    if "_ewk" in unc_name:
        unc_name = unc_name.replace("signal_qcdzjets", "signal_ewkzjets")
        unc_name = unc_name.replace("zmumu_ewk", "Zmm_ewkzll").replace("zee_ewk", "Zee_ewkzll")
        unc_name = unc_name.replace("wlnu_ewk", "signal_ewkwjets").replace("wmunu_ewk", "Wmn_ewkwjets").replace("wenu_ewk", "Wen_ewkwjets")
        unc_name = unc_name.replace("gjets_ewk", "gjets_ewkgjets")

    unc_name = unc_name.replace("zoverw_over_qcd_w_qcd", "signal_qcdzjets_over_signal_qcdwjets")
    unc_name = unc_name.replace("zoverw_over_ewk_w_ewk", "signal_ewkzjets_over_signal_ewkwjets")
    unc_name = unc_name.replace("goverz_over_qcd_photon_qcd", "signal_qcdzjets_over_gjets_qcdgjets")
    unc_name = unc_name.replace("goverz_over_ewk_photon_ewk", "signal_ewkzjets_over_gjets_ewkgjets")

    unc_name = unc_name.replace("zoverw_over_qcd_w_zjets", "signal_qcdzjets_over_signal_qcdwjets")
    unc_name = unc_name.replace("goverz_over_qcd_photon_zjets", "signal_qcdzjets_over_gjets_qcdgjets")

    # TODO: why are these ratios missing?
    unc_name = unc_name.replace("signal_qcdwjets_over_Wmn_ewkwjets", "signal_qcdzjets_over_signal_ewkzjets")
    unc_name = unc_name.replace("signal_qcdwjets_over_Wen_ewkwjets", "signal_qcdzjets_over_signal_ewkzjets")
    factor = unc_file.Get(unc_name)
    variation = nominal.Clone(new_name)
    if factor.GetNbinsX() == 1:
        factor_value = factor.GetBinContent(1)
        variation.Scale(factor_value)
    else:
        assert variation.Multiply(factor)
        # TODO
        if "_pdf" in unc_name:
            for bin_idx in range(0, factor.GetNbinsX() + 2):  # includes underflow (0) and overflow (nbins+1)
                content = factor.GetBinContent(bin_idx)
                new_content = (content - 1) / 10 + 1
                factor.SetBinContent(bin_idx, new_content)
        variation.Multiply(factor)
    outfile.WriteTObject(variation)
