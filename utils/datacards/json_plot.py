import matplotlib.pyplot as plt
import numpy as np
import json
import argparse
import subprocess


def plot_combined(infos, infos2=None, max_impact=1.0):
    # Infos from first file
    title = infos["title"]
    names = infos["names"]
    impacts = infos["impacts"]
    values = infos["values"]
    errors = infos["errors"]
    file1 = infos["file"]

    # Infos from second file (if available)
    if infos2:
        impacts2 = infos2["impacts"]
        values2 = infos2["values"]
        errors2 = infos2["errors"]
        file2 = infos2["file"]

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
        impacts = [100.0 * imp / max_impact for imp in impacts]
        ax1.set_ylim(-5, 105)
        ax1.plot([], [], " ", color="k", label=f"Max impact: {max_impact:.3f}")
        impact_ylabel = "Impact / max impact (%)"
        if infos2:
            impacts2 = [100.0 * imp / max_impact for imp in impacts2]
    # Plot data points for the impacts in the first file
    ax1.plot(x, impacts, "o", markersize=5, linestyle="None", label=file1.split("/")[0])

    # Plot data points for the impacts in the secound file (if available)
    if infos2:
        ax1.plot(x, impacts2, "s", markersize=5, linestyle="None", label=file2.split("/")[0], color="orange")

    # Format plot
    ax1.axhline(0, color="grey", linestyle="--", linewidth=0.7)
    ax1.set_ylabel(impact_ylabel)
    ax1.grid(True)
    ax1.legend()
    ax1.set_title(title)

    # --- Plot postfit values (bottom)
    # Plot error bars for the postfit values in the first file
    ax2.errorbar(x, values, yerr=errors, fmt="o", markersize=5, linestyle="None", label=file1.split("/")[0])
    # Plot error bars for the postfit values in the second file (if available)
    if infos2:
        # Add a small offset to the x values for the second file to avoid overlap
        x2 = [xi + 0.2 for xi in x]
        ax2.errorbar(x2, values2, yerr=errors2, fmt="s", markersize=5, linestyle="None", label=file2.split("/")[0], color="orange")

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
    return filtered


parser = argparse.ArgumentParser(description="Arguments for workspace creation.")
parser.add_argument("--file", type=str)
parser.add_argument("--file2", type=str)
args = parser.parse_args()

path1 = args.file  # , args.file2
path2 = args.file2 if args.file2 else None

data = json.load(open(path1))["params"]
data2 = json.load(open(path2))["params"] if path2 else None

max_impact = max(param["impact_r"] for param in data)

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
        filt_data = filter_data(data, condition=title) + filter_data(data, condition="jer")
        filt_data2 = filter_data(data2, condition=title) + filter_data(data2, condition="jer") if data2 else []
    elif title == "misc":
        filt_data = filter_data(data, condition="Top_Reweight13TeV") + filter_data(data, condition="UEPS") + filter_data(data, condition="ZJets_Norm13TeV")
        filt_data2 = (
            filter_data(data2, condition="Top_Reweight13TeV") + filter_data(data2, condition="UEPS") + filter_data(data2, condition="ZJets_Norm13TeV")
            if data2
            else []
        )
    else:
        filt_data = filter_data(data, condition=title)
        filt_data2 = filter_data(data2, condition=title) if data2 else []

    names, values, errors, impacts = zip(*filt_data)
    names2, values2, errors2, impacts2 = zip(*filt_data2) if data2 else ([], [], [])

    # Remove entries from names2 that are not in names
    idx, names2 = zip(*[(idx, name) for (idx, name) in enumerate(names2) if name in names])
    names2 = list(names2)
    values2 = [values2[i] for i in idx]
    errors2 = [errors2[i] for i in idx]
    impacts2 = [impacts2[i] for i in idx]

    # Pad entries in names not in names2
    to_pad = set(names) - set(names2)
    idx_to_pad = [names.index(name) for name in to_pad]
    idx_to_pad.sort()
    for idx in reversed(idx_to_pad):
        names2.insert(idx, names[idx])
        values2.insert(idx, 0.0)
        errors2.insert(idx, 0.0)
        impacts2.insert(idx, 0.0)

    names2, values2, errors2, impacts2 = tuple(names2), tuple(values2), tuple(errors2), tuple(impacts2)

    info_1 = {"title": title, "names": names, "values": values, "errors": errors, "file": path1, "impacts": impacts}
    info_2 = {"title": title, "names": names2, "values": values2, "errors": errors2, "file": path2, "impacts": impacts2} if data2 else None

    # plot_postfit(info_1, info_2)
    # plot_impacts(info_1, info_2)
    plot_combined(info_1, info_2, max_impact=max_impact)


subprocess.run(["pdfunite"] + [f"{title}_combined.pdf" for title in plotlist] + ["impacts.pdf"])
subprocess.run(["rm"] + [f"{title}_combined.pdf" for title in plotlist])
