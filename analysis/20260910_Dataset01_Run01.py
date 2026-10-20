# from shared import tem_chromatin_analyzer as TCA

# TCA.process_image(
#     "../data/Sample #1/1.jpg",
#     "../QuPath/Sample 1/export/1-labels.png",
#     "../processed/2026-09-24 Dev",
# )

# TCA.process_folder("../QuPath/Sample 1", "../processed/2026-10-09/Sample 1")
# TCA.process_folder("../QuPath/Sample 2", "../processed/2026-10-09/Sample 2")

# TCA.process_folder("../QuPath/Sample 3", "../processed/2026-10-09/Sample 3")
# TCA.process_folder("../QuPath/Sample 4", "../processed/2026-10-09/Sample 4")

# TCA.process_folder("../QuPath/Sample 5", "../processed/2026-10-09/Sample 5")
# TCA.process_folder("../QuPath/Sample 6", "../processed/2026-10-09/Sample 6")


# Make plots

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# 1. Define the root folder containing the sample subdirectories ("Sample 1", "Sample 2", etc.)
root_dir = Path("../processed/2026-10-09")  # Replace with your actual directory path

# List of metric properties to plot
properties = [
    "average_hchromatin_area",
    "nuclear_area",
    "nuclear_equiv_diameter_area",
    "nuclear_eccentricity",
    "nuclear_circularity",
    "nuclear_solidity",
    "nuclear_tortuosity",
]

# 2. Find and read all CSV files across subfolders
df_list = []

for csv_path in root_dir.rglob("summary_data.csv"):
    # Extract the subfolder name as the sample label (e.g., "Sample 1")
    sample_name = csv_path.parent.name

    # Read CSV
    temp_df = pd.read_csv(csv_path)

    # Assign the sample identifier column
    temp_df["sample"] = sample_name
    df_list.append(temp_df)

# Combine all data into a single DataFrame
if not df_list:
    raise FileNotFoundError(f"No CSV files found in {root_dir}")

df_all = pd.concat(df_list, ignore_index=True)

# Sort samples naturally if desired (e.g., "Sample 1", "Sample 2")
df_all = df_all.sort_values(by="sample")

# 3. Create Box-and-Whisker Plots
sns.set_theme(style="whitegrid")

# Option A: Plot all properties side-by-side in a multi-panel figure
fig, axes = plt.subplots(nrows=2, ncols=4, figsize=(20, 10))
axes = axes.flatten()

for idx, prop in enumerate(properties):
    sns.boxplot(data=df_all, x="sample", y=prop, ax=axes[idx], palette="Set2")
    axes[idx].set_title(prop.replace("_", " ").title(), fontsize=12, fontweight="bold")
    axes[idx].set_xlabel("Sample", fontsize=10)
    axes[idx].set_ylabel("Value", fontsize=10)
    axes[idx].tick_params(axis="x", rotation=45)

# Hide any remaining unused subplot axes (if 7 metrics in an 8-panel grid)
if len(properties) < len(axes):
    for unused_ax in axes[len(properties) :]:
        unused_ax.set_visible(False)

plt.tight_layout()
plt.show()

# # Option B (Optional): Save individual high-resolution plots per property
# output_dir = Path("./boxplots_output")
# output_dir.mkdir(exist_ok=True)

# for prop in properties:
#     plt.figure(figsize=(8, 6))
#     sns.boxplot(data=df_all, x="sample", y=prop, palette="Set2")
#     # Overlay jittered data points for better visibility of distributions
#     sns.stripplot(
#         data=df_all, x="sample", y=prop, color="black", alpha=0.3, jitter=0.2, size=4
#     )
#     plt.title(f"{prop.replace('_', ' ').title()} by Sample", fontsize=14)
#     plt.xlabel("Sample", fontsize=12)
#     plt.ylabel(prop, fontsize=12)
#     plt.xticks(rotation=45)
#     plt.tight_layout()
#     plt.savefig(output_dir / f"boxplot_{prop}.png", dpi=300)
#     plt.close()
