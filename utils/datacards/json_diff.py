import json
import argparse

# Load the two JSON files
parser = argparse.ArgumentParser(description="Arguments for workspace creation.")
parser.add_argument("--file1", type=str)
parser.add_argument("--file2", type=str)
args = parser.parse_args()

path1, path2 = args.file1, args.file2

with open(path1) as f1, open(path2) as f2:
    data1 = json.load(f1)
    data2 = json.load(f2)

# Build dictionaries for fast lookup by parameter name
impacts1 = {param["name"]: param["impact_r"] for param in data1["params"]}
impacts2 = {param["name"]: param["impact_r"] for param in data2["params"]}

# Combine all unique parameter names
all_params = sorted(set(impacts1) | set(impacts2))

# Compare and print differences
for name in all_params:
    val1 = impacts1.get(name, None)
    val2 = impacts2.get(name, None)

    if val1 is None:
        print(f"{name}: only in second file (impact_r = {val2})")
    elif val2 is None:
        print(f"{name}: only in first file (impact_r = {val1})")
    else:
        diff = val2 - val1
        print(f"{name}: Δimpact_r = {diff:+.6f} (file1: {val1:.6f}, file2: {val2:.6f})")
