#!/usr/bin/env python3

import argparse
import math
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def parse_arguments():
    script_dir = Path(__file__).resolve().parent
    default_postprocessed_dir = (
        script_dir.parent / "pipeline" / "output" / "postprocessed"
    )

    parser = argparse.ArgumentParser(
        description=(
            "Plot the per-residue mean and variance of Q over the final "
            "time window for one mutation."
        )
    )
    parser.add_argument("mutation", help="Mutation directory name, e.g. WT or A45G")
    parser.add_argument("window", type=float, help="Length of the final time window")
    parser.add_argument("time_unit", choices=("ps", "ns"), help="Unit of the time window")
    parser.add_argument(
        "--postprocessed-dir",
        type=Path,
        default=default_postprocessed_dir,
        help=(
            "Directory containing the mutation directories "
            f"(default: {default_postprocessed_dir})"
        ),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Output directory (default: <mutation>/plots)",
    )

    args = parser.parse_args()

    if not math.isfinite(args.window) or args.window <= 0:
        parser.error("window must be a finite number greater than zero")
    if Path(args.mutation).name != args.mutation or args.mutation in {".", ".."}:
        parser.error("mutation must be a directory name, not a path")

    return args


def get_residue_number(path):
    match = re.fullmatch(r"residue_(-?\d+)\.dat", path.name)
    if match is None:
        raise ValueError(f"Unexpected per-residue Q filename: {path.name}")
    return int(match.group(1))


def read_residue_file(path):
    with path.open() as file:
        header = [field.strip().lstrip("#").strip()
                  for field in file.readline().rstrip("\r\n").split("\t")]

    if len(header) < 3:
        raise ValueError(f"Invalid header in {path}")

    time_match = re.fullmatch(r"Time\s*\((ps|ns)\)", header[1], re.IGNORECASE)
    if header[0] != "Frame" or time_match is None or header[2] != "Q":
        raise ValueError(
            f"Expected tab-separated columns 'Frame', 'Time (ps|ns)', 'Q' "
            f"in {path}; found: {header}"
        )

    data = np.genfromtxt(
        path,
        delimiter="\t",
        skip_header=1,
        usecols=(1, 2),
        dtype=float,
        missing_values=("NA", ""),
        filling_values=np.nan,
        ndmin=2,
    )
    if data.size == 0:
        raise ValueError(f"No data rows in {path}")

    times = data[:, 0]
    q_values = data[:, 1]
    if not np.all(np.isfinite(times)):
        raise ValueError(f"Invalid time values in {path}")
    if np.any(np.diff(times) <= 0):
        raise ValueError(f"Time values must be strictly increasing in {path}")

    return times, q_values, time_match.group(1).lower()


def main():
    args = parse_arguments()

    unit_to_ps = {"ps": 1.0, "ns": 1000.0}
    mutation_dir = args.postprocessed_dir / args.mutation
    residue_dir = mutation_dir / "analysis" / "qvalue_per_residue"

    if not residue_dir.is_dir():
        raise FileNotFoundError(f"Per-residue Q directory not found: {residue_dir}")

    files = sorted(residue_dir.glob("residue_*.dat"), key=get_residue_number)
    if not files:
        raise FileNotFoundError(f"No residue_*.dat files found in {residue_dir}")

    residues = []
    mean_values = []
    variance_values = []
    reference_times = None
    data_unit = None

    window_in_ps = args.window * unit_to_ps[args.time_unit]

    for residue_file in files:
        residue = get_residue_number(residue_file)
        times, q_values, file_unit = read_residue_file(residue_file)

        if data_unit is None:
            data_unit = file_unit
        elif file_unit != data_unit:
            raise ValueError(f"Inconsistent time units in {residue_dir}")

        if reference_times is None:
            reference_times = times
        elif (
            times.shape != reference_times.shape
            or not np.allclose(times, reference_times, rtol=1e-7, atol=1e-10)
        ):
            raise ValueError(f"Inconsistent time columns in {residue_dir}")

        window_in_data_unit = window_in_ps / unit_to_ps[data_unit]
        final_time = reference_times[-1]
        in_window = times >= final_time - window_in_data_unit
        values = q_values[in_window & np.isfinite(q_values)]

        if np.any(values < -1e-6) or np.any(values > 1.0 + 1e-6):
            raise ValueError(f"Q values outside [0, 1] in {residue_file}")

        if values.size:
            mean_q = float(np.mean(values))
            variance_q = float(np.var(values, ddof=0))
            residues.append(residue)
            mean_values.append(mean_q)
            variance_values.append(variance_q)
        else:
            mean_q = np.nan
            variance_q = np.nan

    if not residues:
        raise ValueError(
            f"No valid Q values found in the final {args.window:g} {args.time_unit} "
            f"for mutation {args.mutation}"
        )

    output_dir = args.output_dir or mutation_dir / "plots"
    output_dir.mkdir(parents=True, exist_ok=True)

    window_label = f"{args.window:g}{args.time_unit}"
    plot_path = output_dir / f"qvalue_per_residue_mean_variance_last_{window_label}.png"

    fig_width = max(10, 0.32 * len(residues))
    fig, mean_ax = plt.subplots(
        figsize=(fig_width, 5.5),
        layout="constrained",
    )

    mean_ax.errorbar(
        residues,
        mean_values,
        yerr=np.sqrt(variance_values),
        fmt="o",
        markersize=3.5,
        capsize=2,
        linewidth=1,
        label="Mean Q ± standard deviation",
    )
    mean_ax.set(
        xlabel="Residue number",
        ylabel="Per-residue Q",
        title=(
            f"{args.mutation}: per-residue Q over the final "
            f"{args.window:g} {args.time_unit}"
        ),
        ylim=(0, 1),
    )
    mean_ax.set_xticks(residues)
    mean_ax.tick_params(axis="x", labelrotation=90)
    mean_ax.grid(axis="y", alpha=0.25)
    mean_ax.legend()

    fig.savefig(plot_path, dpi=180)
    plt.close(fig)

    print(f"Mutation: {args.mutation}")
    print(
        f"Window: final {args.window:g} {args.time_unit} "
        f"(trajectory time unit: {data_unit})"
    )
    print(f"Residues plotted: {len(residues)}")
    print(f"Saved plot: {plot_path}")


if __name__ == "__main__":
    main()
