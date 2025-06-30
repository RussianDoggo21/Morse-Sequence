import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

# File paths
CSV_CPP = "timing_results_cpp.csv"
CSV_PY = "timing_results_python.csv"

# Load data
df_cpp = pd.read_csv(CSV_CPP)
df_py = pd.read_csv(CSV_PY)

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


# Get all unique grid sizes
grid_sizes = sorted(df["size"].unique())

# Create output directory for plots
output_dir = "plots/by_operation"
os.makedirs(output_dir, exist_ok=True)

# Generate one plot per grid size
for size in grid_sizes:
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

    # Save the plot
    plt.savefig(f"{output_dir}/grid_{size}x{size}.png")
    plt.close()

print("✅ Plots saved in:", output_dir)
