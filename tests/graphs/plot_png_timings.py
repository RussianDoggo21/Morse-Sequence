#Launch this file from the graphs repertory

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

# === File paths configuration ===
csv_cpp = "timing_results_cpp.csv"
csv_py = "timing_results_python.csv"

# Check if CSV files exist
if not os.path.exists(csv_cpp) or not os.path.exists(csv_py):
    raise FileNotFoundError("Missing CSV files. Run timer_m_sequence() and timer_comparison() first.")

# === Load data ===
df_cpp = pd.read_csv(csv_cpp)
df_py = pd.read_csv(csv_py)
df = pd.concat([df_cpp, df_py], ignore_index=True)

# Clean up the data
expected_columns = ["size", "operation", "language", "time_seconds"]
df = df.loc[:, expected_columns]
df = df[~df.apply(lambda row: any(val == col for val, col in zip(row, df.columns)), axis=1)]
df["size"] = df["size"].astype(int)
df["time_seconds"] = df["time_seconds"].astype(float)
df = df[df["operation"] != "total"]

# === Create output folders ===
output_dir_ops = "plots/by_grid_sizes"
output_dir_sizes = "plots/by_operation"
os.makedirs(output_dir_ops, exist_ok=True)
os.makedirs(output_dir_sizes, exist_ok=True)

# === Fixed color palette ===
palette = {
    "cpp": "#1f77b4",
    "python": "#2ca02c",
    "cpp/python": "#ff7f0e"
}

# === Seaborn style ===
sns.set_theme(style="whitegrid")

# === 1. Line plots per operation (one plot per operation, X axis = grid size) ===
for op in df["operation"].unique():
    plt.figure(figsize=(8, 5))
    subset = df[df["operation"] == op]

    ax = sns.lineplot(
        data=subset,
        x="size",
        y="time_seconds",
        hue="language",
        palette=palette,
        marker="o"
    )
    ax.set_title(f"Execution time for '{op}' by grid size")
    ax.set_xlabel("Grid size (k × k)")
    ax.set_ylabel("Time (seconds)")
    ax.legend(title="Language")
    plt.tight_layout()
    output_path = f"{output_dir_ops}/{op}.png"
    plt.savefig(output_path)
    plt.close()
    print(f"✅ Plot saved: {output_path}")

# === 2. Bar plots per grid size (one plot per grid size, X axis = operation) ===
for size in sorted(df["size"].unique()):
    df_size = df[df["size"] == size]

    plt.figure(figsize=(10, 6))
    sns.barplot(
        data=df_size,
        x="operation",
        y="time_seconds",
        hue="language",
        palette="Set2"
    )

    plt.title(f"Execution Time by Operation — Grid {size}×{size}")
    plt.xlabel("Operation")
    plt.ylabel("Time (seconds)")
    plt.xticks(rotation=30)
    plt.tight_layout()
    plt.legend(title="Implementation")
    plt.grid(axis='y', linestyle='--', linewidth=0.5)

    output_path = f"{output_dir_sizes}/grid_{size}x{size}.png"
    plt.savefig(output_path)
    plt.close()
    print(f"✅ Plot saved: {output_path}")

print("\n🎉 All plots generated successfully.")
