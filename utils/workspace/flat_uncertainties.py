from typing import Any


def get_processes(analysis: str, region: str, type: str) -> list[str]:
    """Return the list of processes for a given analysis, region, and type."""
    processes = {
        "vbf": {
            "signal": {
                "signals": ["zh", "wh", "vbf", "ggh"],
                "models": ["qcd_zjets", "qcd_wjets", "ewk_zjets", "ewk_wjets"],
                "backgrounds": ["qcdzll", "ewkzll", "top", "diboson"],
            },
            "dimuon": {
                "models": ["qcd_zll", "ewk_zll"],
                "backgrounds": ["top", "diboson"],
            },
            "dielec": {
                "models": ["qcd_zll", "ewk_zll"],
                "backgrounds": ["top", "diboson"],
            },
            "singlemu": {
                "models": ["qcd_wjets", "ewk_wjets"],
                "backgrounds": ["qcdzll", "ewkzll", "top", "diboson"],
            },
            "singleel": {
                "models": ["qcd_wjets", "ewk_wjets"],
                "backgrounds": ["qcdzll", "qcdgjets", "ewkzll", "top", "diboson"],
            },
            "photon": {
                "models": ["qcd_gjets", "ewk_gjets"],
            },
        },
        # TODO: "monojet": {...}
    }

    return processes.get(analysis, {}).get(region, {}).get(type, [])


def get_region_label_map() -> list[tuple[str, str]]:
    return [
        ("dielec", "Zee"),
        ("dimuon", "Zmm"),
        ("signal", "signal"),
        ("singleel", "Wen"),
        ("singlemu", "Wmn"),
        ("photon", "gjets"),
    ]


def get_process_model_map(region: str) -> dict[str, dict[str, str]]:
    return {
        "dielec": {
            "ewk_zll": "ewk_dielectron_ewk_zjets",
            "qcd_zll": "qcd_dielectron_qcd_zjets",
        },
        "dimuon": {
            "ewk_zll": "ewk_dimuon_ewk_zjets",
            "qcd_zll": "qcd_dimuon_qcd_zjets",
        },
        "signal": {
            "ewk_wjets": "ewk_wjetssignal_ewk_zjets",
            "ewk_zjets": "ewkqcd_signal_qcd_zjets",
            "qcd_wjets": "qcd_wjetssignal_qcd_zjets",
            "qcd_zjets": "signal_qcd_zjets",
        },
        "singleel": {
            "ewk_wjets": "ewk_singleelectron_ewk_wjets",
            "qcd_wjets": "qcd_singleelectron_qcd_wjets",
        },
        "singlemu": {
            "ewk_wjets": "ewk_singlemuon_ewk_wjets",
            "qcd_wjets": "qcd_singlemuon_qcd_wjets",
        },
        "photon": {
            "ewk_gjets": "ewk_photon_ewk_zjets",
            "qcd_gjets": "qcd_photon_qcd_zjets",
        },
    }[region]


def get_flat_unc(process: str) -> dict[str, float]:
    return {
        "zmm": {"pu": 0.01, "id": 0.01, "trigger": 0.01},
        "zee": {"pu": 0.01, "id": 0.03, "trigger": 0.005},
        "photon": {"pu": 0.01, "id": 0.015, "trigger": 0.01},
        "w": {"pu": 0.01},
        "wen": {"pu": 0.01, "id": 0.015, "trigger": 0.01},
        "wmn": {"pu": 0.01, "id": 0.005, "trigger": 0.01},
    }[process]


def get_veto_unc(model: str) -> dict[str, float]:
    # TODO split or with the same name?
    return {
        # VBF
        "ewk_wjets": {"t": 0.01, "m": 0.02, "e": 0.03},
        "qcd_wjets": {"t": 0.01, "m": 0.015, "e": 0.03},
        "ewk_zjets": {"t": -0.01, "m": -0.02, "e": -0.03},
        "qcd_zjets": {"t": -0.01, "m": -0.015, "e": -0.03},
        # Mono
        "wjets": {"t": 0.01, "m": 0.015, "e": 0.03},
        "zjets": {"t": -0.01, "m": -0.015, "e": -0.03},
    }[model]


def get_processes_by_type(analysis: str, types: list[str] = ["signals", "backgrounds"]) -> set[str]:
    """Return the set of all processes of the given types used in all regions of the given analysis."""
    regions = [region for region, _ in get_region_label_map()]
    return {proc for region in regions for category in types for proc in get_processes(analysis=analysis, region=region, type=category)}


def get_processes_by_region(analysis: str, region: str, types: list[str] = ["backgrounds", "models"]) -> set[str]:
    """Return the set of all processes of the given types used in all regions of the given analysis."""
    return {proc for category in types for proc in get_processes(analysis=analysis, region=region, type=category)}


def get_lumi_unc(year: str, analysis: str) -> dict[str, Any]:
    """Return luminosity lnN systematics for a given year and analysis."""
    _ = analysis  # Currently unused

    proc_list = get_processes_by_type(analysis=analysis)
    return {
        "Run3": {
            f"lumi_13p6TeV_{year}": {"value": 1.014, "processes": proc_list},
        },
    }[year]


def get_lepton_eff_unc(year: str, analysis: str) -> dict[str, Any]:
    """Return lepton and photon efficiency lnN systematics for a given year and analysis."""
    _ = analysis  # Currently unused
    procs_by_region = {
        "dielec": get_processes_by_region(analysis=analysis, region="dielec"),
        "singleel": get_processes_by_region(analysis=analysis, region="singleel"),
        "dimuon": get_processes_by_region(analysis=analysis, region="dimuon"),
        "singlemu": get_processes_by_region(analysis=analysis, region="singlemu"),
        "photon": get_processes_by_region(analysis=analysis, region="photon"),
        "signal": get_processes_by_type(analysis=analysis),
    }

    # Run-dependent efficiencies
    m_id_eff = {"Run3": 0.004}[year]
    m_iso_eff = {"Run3": 0.005}[year]
    e_id_eff = {"Run3": 0.014}[year]
    e_reco_eff = {"Run3": 0.010}[year]  # TODO update for run3
    g_id_eff = {"Run3": 0.014}[year]
    results = {
        f"CMS_eff_b_{year}": {"value": 1.03, "processes": ["top"]},
        f"CMS_fake_b_{year}": {"value": 1.01, "processes": procs_by_region["signal"] - {"top"}},
    }
    lepton_unc = {
        f"CMS_eff_e_id_{year}": {"dielec": 2 * e_id_eff, "singleel": e_id_eff},
        f"CMS_eff_e_reco_{year}": {"dielec": 2 * e_reco_eff, "singleel": e_reco_eff},
        f"CMS_eff_m_id_{year}": {"dimuon": 2 * m_id_eff, "singlemu": m_id_eff},
        f"CMS_eff_m_iso_{year}": {"dimuon": 2 * m_iso_eff, "singlemu": m_iso_eff},
        f"CMS_eff_g_id_{year}": {"photon": g_id_eff},
    }

    for name, entries in lepton_unc.items():
        results[name] = {region: {"value": 1 + value, "processes": procs_by_region[region]} for region, value in entries.items()}

    return results


def get_trigger_unc(year: str, analysis: str) -> dict[str, Any]:
    """Return trigger-related lnN systematics for a given year and analysis."""
    _ = analysis  # Currently unused
    # All processes
    proc_list = get_processes_by_type(analysis=analysis, types=["signals", "backgrounds", "models"])
    return {
        "Run3": {
            f"CMS_eff_g_trigger_{year}_13p6TeV": {
                "photon": {"value": 1.01, "processes": proc_list},
            },
            f"CMS_eff_e_trigger_{year}_13p6TeV": {
                "dielec": {"value": 1.01, "processes": proc_list},
                "singleel": {"value": 1.01, "processes": proc_list},
            },
            f"CMS_eff_met_trigger_{year}_13p6TeV": {
                "signal": {"value": 1.02, "processes": proc_list},
                "dimuon": {"value": 1.02, "processes": proc_list},
                "singlemu": {"value": 1.02, "processes": proc_list},
            },
        },
    }[year]


def get_qcd_unc(year: str, analysis: str) -> dict[str, Any]:
    """Return QCD scale lnN systematics for a given year and analysis."""
    _ = analysis  # Currently unused
    return {
        "Run3": {
            "QCDscale_VV": {"value": 1.15, "processes": ["diboson"]},
            "QCDscale_VV_ACCEPT": {"value": 1.15, "processes": ["diboson"]},
            "QCDscale_ttbar": {"value": 1.1, "processes": ["top"]},
            "QCDscale_ttbar_ACCEPT": {"value": 1.1, "processes": ["top"]},
            "QCDscale_ggH2in": {"value": (0.933, 1.046), "processes": ["ggh"]},
            "QCDscale_ggH2in_ACCEPT": {"value": 1.4, "processes": ["ggh"]},
            "QCDscale_qqH": {"value": (0.997, 1.004), "processes": ["vbf"]},
            "QCDscale_qqH_ACCEPT": {"value": 1.02, "processes": ["vbf"]},
        },
    }[year]


def get_pdf_unc(year: str, analysis: str) -> dict[str, Any]:
    """Return PDF systematics for a given year and analysis."""
    _ = analysis  # Currently unused
    return {
        "Run3": {
            "pdf_Higgs_gg": {"value": 1.032, "processes": ["ggh"]},
            "pdf_Higgs_qq": {"value": 1.021, "processes": ["vbf"]},
            "pdf_Higgs_qq_ACCEPT": {"value": 1.01, "processes": ["vbf"]},
        },
    }[year]


def get_misc_unc(year: str, analysis: str) -> dict[str, Any]:
    """Return miscellaneous lnN systematics for a given year and analysis."""
    _ = analysis  # Currently unused
    return {
        "Run3": {
            "Top_Reweight13TeV": {"value": 1.1, "processes": ["top"]},
            "UEPS": {"value": 1.168, "processes": ["ggh"]},
            "ZJets_Norm13TeV": {"value": 1.2, "processes": ["qcdzll", "ewkzll"]},
        },
    }[year]


def get_jec_shape(year: str, analysis: str) -> dict[str, dict[str, Any]]:
    """Return JER and JES shape systematics for a given year and analysis."""
    _ = year  # Currently unused
    _ = analysis  # Currently unused

    proc_list = get_processes_by_type(analysis=analysis)
    # TODO need to change to CMS_scale_j_Absolute
    return {
        f"jer_{year}": {"value": 1.0, "processes": proc_list},
        "jesAbsolute": {"value": 1.0, "processes": proc_list},
        f"jesAbsolute_{year}": {"value": 1.0, "processes": proc_list},
        "jesBBEC1": {"value": 1.0, "processes": proc_list},
        f"jesBBEC1_{year}": {"value": 1.0, "processes": proc_list},
        "jesEC2": {"value": 1.0, "processes": proc_list},
        f"jesEC2_{year}": {"value": 1.0, "processes": proc_list},
        "jesFlavorQCD": {"value": 1.0, "processes": proc_list},
        "jesHF": {"value": 1.0, "processes": proc_list},
        f"jesHF_{year}": {"value": 1.0, "processes": proc_list},
        "jesRelativeBal": {"value": 1.0, "processes": proc_list},
        f"jesRelativeSample_{year}": {"value": 1.0, "processes": proc_list},
    }
