import matplotlib.pyplot as plt
import numpy as np
import json
import argparse
import subprocess


def plot_impacts(infos, infos2=None):
    title, names, impacts, file = infos["title"], infos["names"], infos["impacts"], infos["file"]
    plt.figure(figsize=(12, 6))

    # X positions for the first set
    x = range(len(impacts))
    # import pdb

    # pdb.set_trace()
    plt.plot(x, impacts, "o", markersize=5, linestyle="None", label=file.split("/")[0])

    # Plot second set if provided
    if infos2:
        names2, impacts2, file2 = infos2["names"], infos2["impacts"], infos2["file"]
        # Offset x positions slightly
        # x2 = [xi + 0.2 for xi in x]
        plt.plot(x, impacts2, "s", markersize=5, linestyle="None", label=file2.split("/")[0], color="orange")

    plt.axhline(0, color="grey", linestyle="--", linewidth=0.7)
    plt.xticks(x, names, rotation=45, ha="right")
    plt.title(title)
    plt.ylabel("Impacts")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"{title}_impacts.pdf")
    plt.close()


def plot_postfit(infos, infos2=None):
    title, names, values, errors, file = infos["title"], infos["names"], infos["values"], infos["errors"], infos["file"]
    plt.figure(figsize=(12, 6))

    # X positions for the first set
    x = range(len(values))
    plt.errorbar(x, values, yerr=errors, fmt="o", markersize=5, linestyle="None", label=file.split("/")[0])

    # Plot second set if provided
    if infos2:
        names2, values2, errors2, file2 = infos2["names"], infos2["values"], infos2["errors"], infos2["file"]
        # Offset x positions slightly
        x2 = [xi + 0.2 for xi in x]
        plt.errorbar(x2, values2, yerr=errors2, fmt="s", markersize=5, linestyle="None", label=file2.split("/")[0], color="orange")

    plt.axhline(0, color="grey", linestyle="--", linewidth=0.7)
    plt.xticks(x, names, rotation=45, ha="right")
    plt.title(title)
    plt.ylabel("Postfit value")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"{title}_postfit.pdf")
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

    plot_postfit(info_1, info_2)
    plot_impacts(info_1, info_2)


# Merge all plots into a single PDF
subprocess.run(["pdfunite"] + [f"{title}_postfit.pdf" for title in plotlist] + ["postfit_np.pdf"])
subprocess.run(["rm"] + [f"{title}_postfit.pdf" for title in plotlist])

subprocess.run(["pdfunite"] + [f"{title}_impacts.pdf" for title in plotlist] + ["impacts_np.pdf"])
subprocess.run(["rm"] + [f"{title}_impacts.pdf" for title in plotlist])
