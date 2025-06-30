import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

# Load CSV files
csv_cpp = "timing_results_cpp.csv"
csv_py = "timing_results_python.csv"

if not os.path.exists(csv_cpp) or not os.path.exists(csv_py):
    raise FileNotFoundError("Missing CSV files. Run timer_m_sequence() and timer_comparison() first.")

# Read and merge the data
df_cpp = pd.read_csv(csv_cpp)
df_py = pd.read_csv(csv_py)

# Combine both DataFrames
df = pd.concat([df_cpp, df_py], ignore_index=True)

# Keep only the relevant columns (and drop any extra ones)
expected_columns = ["size", "operation", "language", "time_seconds"]
df = df.loc[:, expected_columns]

# Remove rows where any column is literally the name of the column (i.e., repeated headers)
df = df[~df.apply(lambda row: any(val == col for val, col in zip(row, df.columns)), axis=1)]

# Ensure 'size' is integer and 'time_seconds' is float
df["size"] = df["size"].astype(int)
df["time_seconds"] = df["time_seconds"].astype(float)

# Convert 'size' to integer (in case it's a float)
df["size"] = df["size"].astype(int)

# Remove total timing if present
df = df[df["operation"] != "total"]

# Create output directory for plots
os.makedirs("plots/by_grid_sizes", exist_ok=True)

# Set seaborn style
sns.set_theme(style="whitegrid")

# List of operations to plot
operations = df["operation"].unique()

# Fixed color palette
palette = {
    "cpp": "#1f77b4",
    "python": "#2ca02c",
    "cpp/python": "#ff7f0e"
}

# Generate one plot per operation
for op in operations:
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
    plt.savefig(f"plots/by_grid_sizes/{op}.png")
    plt.close()
    print(f"✅ Plot saved: plots/{op}.png")
