import matplotlib.pyplot as plt
import numpy as np
import json
import argparse
import subprocess


def plot_combined(title: str, names: list[str], all_data: list[dict], max_impact: float = 1.0, labels=[]) -> None:

    impacts = [data["impacts"] for data in all_data]
    values = [data["values"] for data in all_data]
    errors = [data["errors"] for data in all_data]

    markers = ["o", "s", "D", "^", "v", "x", "p", "*"]
    colors = ["blue", "orange", "green", "red", "purple", "brown", "pink"]

    # X axis, based on how many data points we have
    x = range(len(names))

    # Create two subplots, one for impacts and one for postfit values
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(12, 10), gridspec_kw={"height_ratios": [1, 1]})
    fig.subplots_adjust(hspace=0.1)

    # --- Plot impacts (top)
    impact_ylabel = "Impact"
    # Normalize impacts if the maximum impact is not 1.0
    if max_impact != 1.0:
        # Normalize impacts to the maximum impact value
        impacts = [[100.0 * imp / max_impact for imp in impact_list] for impact_list in impacts]
        ax1.set_ylim(-5, 105)
        ax1.plot([], [], " ", color="k", label=f"Max impact: {max_impact:.3f}")
        impact_ylabel = "Impact / max impact (%)"

    # Plot data points for the impacts in the first file
    for idx, impact_list in enumerate(impacts):
        ax1.plot(
            x,
            impact_list,
            markers[idx % len(markers)],
            color=colors[idx % len(colors)],
            markersize=5,
            linestyle="None",
            label=labels[idx] if labels else None,
        )

    # Format plot
    ax1.axhline(0, color="grey", linestyle="--", linewidth=0.7)
    ax1.set_ylabel(impact_ylabel)
    ax1.grid(True)
    ax1.legend()
    ax1.set_title(title)

    # --- Plot postfit values (bottom)
    # Plot error bars for the postfit values in the first file
    for idx, (value_list, error_list) in enumerate(zip(values, errors)):
        # list of offset x values for better visibility
        off_x = [xi + 0.2 * idx for xi in x]

        ax2.errorbar(
            off_x,
            value_list,
            yerr=error_list,
            fmt=markers[idx % len(markers)],
            color=colors[idx % len(colors)],
            markersize=5,
            linestyle="None",
            label=labels[idx] if labels else None,
        )

    # Format plot
    ax2.axhline(0, color="grey", linestyle="--", linewidth=0.7)
    ax2.set_ylabel("Postfit value")
    ax2.set_xticks(x)
    ax2.set_xticklabels(names, rotation=45, ha="right")
    ax2.grid(True)
    ax2.legend()

    plt.tight_layout()
    plt.savefig(f"{title}_combined.pdf")
    plt.close()


def filter_data(data, condition):
    filtered = [
        (el["name"].replace(condition + "_", ""), el["fit"][1], abs(el["fit"][1] - el["fit"][0]), el["impact_r"]) for el in data if condition in el["name"]
    ]
    filtered = sorted(filtered, key=lambda x: "NoTopTag" not in x[0])
    filtered = sorted(filtered, key=lambda x: "CR" not in x[0])

    names, values, errors, impacts = zip(*filtered) if filtered else ([], [], [], [])
    filt_data = {"names": list(names), "values": list(values), "errors": list(errors), "impacts": list(impacts)}
    return filt_data


parser = argparse.ArgumentParser(description="Arguments for workspace creation.")
parser.add_argument("--files", type=str, nargs="+")
args = parser.parse_args()

json_paths = args.files

json_data = [json.load(open(path))["params"] for path in json_paths if path]

max_impact = max(param["impact_r"] for param in json_data[0]) if json_data else 1.0

plotlist = [
    "CMS",
    "QCDscale",
    "Photon",
    "ZnunuWJets",
    "jes",
    "lumi",
    "pdf_Higgs",
    "dielectronCR",
    "dimuonCR",
    "photonCR",
    "singleelectron",
    "singlemuon",
    "wzCR",
    "ewkqcdzCR",
    "ewk_ewk_vbf",
    "ewkphoton_ewk_vbf",
    "qcd_ewk_vbf",
    "qcd_photon_ewk_vbf",
    "model_mu_cat_vbf_Run3_qcd_zjets",
    "misc",
]

for title in plotlist:
    if title == "jes":
        jer_data = [filter_data(data, condition="jer") for data in json_data]
        filt_data = [{key: val + jer_data[idx][key] for key, val in filter_data(data, condition=title).items()} for idx, data in enumerate(json_data)]

    elif title == "misc":
        Top_data = [filter_data(data, condition="Top_Reweight13TeV") for data in json_data]
        UEPS_data = [filter_data(data, condition="UEPS") for data in json_data]
        filt_data = [
            {key: val + Top_data[idx][key] + UEPS_data[idx][key] for key, val in filter_data(data, condition="ZJets_Norm13TeV").items()}
            for idx, data in enumerate(json_data)
        ]
    else:
        filt_data = [filter_data(data, condition=title) for data in json_data]

    # Remove entries in subsequent data that are not in the reference names
    ref_names = filt_data[0]["names"]
    for extra_data in filt_data[1:]:
        idx, names = zip(*[(idx, name) for (idx, name) in enumerate(extra_data["names"]) if name in ref_names])
        extra_data["names"] = list(names)
        for key in ["values", "errors", "impacts"]:
            extra_data[key] = [extra_data[key][i] for i in idx]

    # Pad entries in names not in names2
    for extra_data in filt_data[1:]:
        to_pad = list(set(ref_names) - set(extra_data["names"]))
        idx_to_pad = [ref_names.index(name) for name in to_pad]
        idx_to_pad.sort()
        for name_idx, idx in enumerate(reversed(idx_to_pad)):
            extra_data["names"].insert(idx, to_pad[name_idx])
            for key in ["values", "errors", "impacts"]:
                extra_data[key].insert(idx, 0.0)

    plot_combined(title=title, names=ref_names, all_data=filt_data, max_impact=max_impact, labels=[path.split("/")[0] for path in json_paths])

# Merge all generated PDFs into a single file
subprocess.run(["pdfunite"] + [f"{title}_combined.pdf" for title in plotlist] + ["impacts.pdf"])
subprocess.run(["rm"] + [f"{title}_combined.pdf" for title in plotlist])
