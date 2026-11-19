import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import re

CSV_FILES = [
    "Temperature_Bi_0.001000.csv",
    "Temperature_Bi_0.010000.csv",
    "Temperature_Bi_0.100000.csv",
    "Temperature_Bi_1.000000.csv",
    "Temperature_Bi_10.000000.csv",
    "Temperature_Bi_100.000000.csv",
]

# Stepped-solid boundary in dimensionless coordinates.
# OUTLINE = np.array([
#     [0.0, 0.0],
#     [2.0, 0.0],
#     [2.0, 0.5],
#     [3.0, 0.5],
#     [3.0, 1.5],
#     [2.0, 1.5],
#     [2.0, 2.0],
#     [0.0, 2.0],
#     [0.0, 0.0]
# ])

fig, axes = plt.subplots(
    2, 3,
    figsize=(16, 9),
    constrained_layout=True
)

for ax, filename in zip(axes.flat, CSV_FILES):

    df = pd.read_csv(filename)

    required = {"X", "Y", "Theta"}
    if not required.issubset(df.columns):
        raise ValueError(
            f"{filename} must contain columns X, Y and Theta. "
            f"Found columns: {list(df.columns)}"
        )

    if df.duplicated(["X", "Y"]).any():
        raise ValueError(f"Duplicate X,Y coordinates in {filename}")

    # Reconstruct the 2D field from coordinate data.
    table = (
        df.pivot(index="Y", columns="X", values="Theta")
          .sort_index()
          .sort_index(axis=1)
    )

    x = table.columns.to_numpy(dtype=float)
    y = table.index.to_numpy(dtype=float)
    theta = np.ma.masked_invalid(table.to_numpy(dtype=float))
    X, Y = np.meshgrid(x, y)
# Thus each Biot number has an independent colour scale.
    vmin = float(theta.min())
    vmax = float(theta.max())

    mesh = ax.pcolormesh(
        X, Y, theta,
        shading="nearest",
        cmap="inferno",
        vmin=vmin,
        vmax=vmax
    )

    # Draw the stepped geometry; missing outside cells stay unfilled.
    # ax.plot(
    #     # OUTLINE[:, 0],
    #     # OUTLINE[:, 1],
    #     color="cyan",
    #     linewidth=1.5
    # )

    match = re.search(r"Bi_([0-9.]+)", Path(filename).name)
    bi = match.group(1) if match else Path(filename).stem

    ax.set_title(
        f"Bi = {bi}\nTemperature range: {vmin:.4f} to {vmax:.4f}"
    )
    ax.set_xlabel("X = x/L")
    ax.set_ylabel("Y = y/L")
    ax.set_xlim(0, 3)
    ax.set_ylim(0, 2)
    ax.set_aspect("equal")
    ax.set_xticks(np.arange(0, 3.01, 0.5))
    ax.set_yticks(np.arange(0, 2.01, 0.5))

    # A separate colour bar for this individual heatmap.
    cbar = fig.colorbar(mesh, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Dimensionless temperature, θ")

fig.suptitle(
    "2D Steady-State Conduction: Temperature Heatmaps \n"
    "Independent Colour Scale for Each Biot Number",
    fontsize=15
)

plt.savefig(
    "CFD_temperature_heatmaps_individual_scales.png",
    dpi=300,
    bbox_inches="tight"
)
plt.show()
