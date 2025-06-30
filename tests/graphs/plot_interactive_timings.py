#Launch this file from the graphs repertory

import pandas as pd
import plotly.express as px
import os

# === File paths ===
csv_cpp = "timing_results_cpp.csv"
csv_py = "timing_results_python.csv"

# === Check if CSV files exist ===
if not os.path.exists(csv_cpp) or not os.path.exists(csv_py):
    raise FileNotFoundError("Missing CSV files. Run timer_m_sequence() and timer_comparison() first.")

# === Load data ===
df_cpp = pd.read_csv(csv_cpp)
df_py = pd.read_csv(csv_py)
df = pd.concat([df_cpp, df_py], ignore_index=True)

# === Clean up the data ===
expected_columns = ["size", "operation", "language", "time_seconds"]
df = df.loc[:, expected_columns]
df = df[~df.apply(lambda row: any(val == col for val, col in zip(row, df.columns)), axis=1)]
df["size"] = df["size"].astype(int)
df["time_seconds"] = df["time_seconds"].astype(float)
df = df[df["operation"] != "total"]

# === Create output directories ===
output_dir_ops = "plots_interactive/by_grid_sizes"
output_dir_sizes = "plots_interactive/by_operation"
os.makedirs(output_dir_ops, exist_ok=True)
os.makedirs(output_dir_sizes, exist_ok=True)

# === 1. One interactive plot per operation (X = grid size, Y = time) ===
for op in df["operation"].unique():
    subset = df[df["operation"] == op]

    fig = px.line(
        subset,
        x="size",
        y="time_seconds",
        color="language",
        markers=True,
        title=f"Execution time for '{op}' by grid size",
        labels={
            "size": "Grid size (k × k)",
            "time_seconds": "Time (seconds)",
            "language": "Language"
        }
    )
    fig.update_layout(template="plotly_white")

    output_path = f"{output_dir_ops}/{op}.html"
    fig.write_html(output_path)
    print(f"✅ Interactive plot saved: {output_path}")

# === 2. One interactive bar plot per grid size (X = operation, Y = time) ===
for size in sorted(df["size"].unique()):
    df_size = df[df["size"] == size]

    fig = px.bar(
        df_size,
        x="operation",
        y="time_seconds",
        color="language",
        barmode="group",
        title=f"Execution Time by Operation — Grid {size}×{size}",
        labels={
            "operation": "Operation",
            "time_seconds": "Time (seconds)",
            "language": "Language"
        }
    )
    fig.update_layout(template="plotly_white", xaxis_tickangle=30)

    output_path = f"{output_dir_sizes}/grid_{size}x{size}.html"
    fig.write_html(output_path)
    print(f"✅ Interactive plot saved: {output_path}")

print("\n🎉 All interactive plots generated successfully.")
