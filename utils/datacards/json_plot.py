import matplotlib.pyplot as plt
import numpy as np
import json
import argparse
import subprocess

# Load the two JSON files
parser = argparse.ArgumentParser(description="Arguments for workspace creation.")
parser.add_argument("--file", type=str)
parser.add_argument("--file2", type=str)
args = parser.parse_args()

path1 = args.file  # , args.file2
path2 = args.file2 if args.file2 else None


def plot(title, names, values, errors):
    plt.figure(figsize=(12, 6))
    plt.errorbar(range(len(values)), values, yerr=errors, fmt="o", capsize=5, capthick=2, markersize=5, linestyle="None")
    plt.axhline(0, color="grey", linestyle="--", linewidth=0.7)
    plt.xticks(range(len(values)), names, rotation=45, ha="right")
    plt.title(title)
    plt.ylabel("Postfit value")
    plt.grid(True)
    plt.tight_layout()
    # plt.show()
    plt.savefig(f"{title}.pdf")


def filter_data(data, condition):
    filtered = [(el["name"].replace(condition + "_", ""), el["fit"][1], abs(el["fit"][1] - el["fit"][0])) for el in data if condition in el["name"]]
    filtered = sorted(filtered, key=lambda x: "NoTopTag" not in x[0])
    filtered = sorted(filtered, key=lambda x: "CR" not in x[0])
    return filtered


data = json.load(open(path1))["params"]
data2 = json.load(open(path2))["params"] if path2 else None

# title = "QCDscale"
# qcdscale_data = filter_data(data, condition=title)
# names, values, errors = zip(*qcdscale_data)
# plot(title, names, values, errors)

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
    elif title == "misc":
        filt_data = filter_data(data, condition="Top_Reweight13TeV") + filter_data(data, condition="UEPS") + filter_data(data, condition="ZJets_Norm13TeV")
    else:
        filt_data = filter_data(data, condition=title)

    names, values, errors = zip(*filt_data)
    plot(title, names, values, errors)

subprocess.run(["pdfunite"] + [f"{title}.pdf" for title in plotlist] + ["out.pdf"])
subprocess.run(["rm"] + [f"{title}.pdf" for title in plotlist])


# title = "rate_st"
# rate_st = filter_data(data, condition=title)
# names, values, errors = zip(*rate_st)
# plot(title, names, values, errors)

# title = "rate_others2"
# rate_others2 = filter_data(data, condition=title)
# names, values, errors = zip(*rate_others2)
# plot(title, names, values, errors)

# title = "prop"
# prop = filter_data(data, condition=title)
# names, values, errors = zip(*prop)
# values, errors = np.array(values), np.array(errors)
# mask = abs(values / errors) > 3
# names = [name for name, m in zip(names, mask) if m]
# values = values[mask]
# errors = errors[mask]

# plot(title, names, values, errors)
