#!/usr/bin/env python3

import argparse
import math
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


SCRIPT_DIR = Path(__file__).resolve().parent
RESIDUE_SUBDIR = Path("analysis/qvalue_per_residue")


def parse_arguments():
    default_postprocessed_dir = SCRIPT_DIR.parent / "pipeline/output/postprocessed"

    parser = argparse.ArgumentParser(
        usage="%(prog)s [options] [MUTATION ...] [WINDOW] [{ps,ns}]",
        description="Plot per-residue Q boxplots or compare mean Q across mutations.",
    )
    parser.add_argument(
        "selection", nargs="*", metavar="MUTATION",
        help="Mutation names or 'all' (default), followed by a window (default: 2500) "
             "and its unit (default: ps). Example: WT A45G 2.5 ns.",
    )
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
        help="Output directory (default: <mutation>/plots or eda/output for comparisons)",
    )

    args = parser.parse_intermixed_args()
    mutations = args.selection
    args.time_unit = mutations.pop() if mutations and mutations[-1] in {"ps", "ns"} else "ps"
    args.window = 2500.0
    if mutations:
        try:
            window = float(mutations[-1])
        except ValueError:
            pass
        else:
            args.window = window
            mutations.pop()

    if not math.isfinite(args.window) or args.window <= 0:
        parser.error("window must be a finite number greater than zero")
    if any(Path(name).name != name or name in {"", ".", ".."} for name in mutations):
        parser.error("mutations must be directory names, not paths")
    if "all" in mutations and mutations != ["all"]:
        parser.error("use 'all' by itself or specify individual mutations")

    args.mutations = list(dict.fromkeys(mutations)) or ["all"]
    return args


def get_mutations(args):
    if args.mutations != ["all"]:
        return args.mutations

    mutations = [
        path.name for path in args.postprocessed_dir.iterdir()
        if (path / RESIDUE_SUBDIR).is_dir()
    ]
    if not mutations:
        raise FileNotFoundError(f"No per-residue Q output found in {args.postprocessed_dir}")
    return sorted(mutations, key=lambda name: (name != "WT", name.count("_"), name))


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


def load_residue_samples(args, mutation, *, required=True):
    unit_to_ps = {"ps": 1.0, "ns": 1000.0}
    residue_dir = args.postprocessed_dir / mutation / RESIDUE_SUBDIR

    if not residue_dir.is_dir():
        raise FileNotFoundError(f"Per-residue Q directory not found: {residue_dir}")

    files = sorted(residue_dir.glob("residue_*.dat"), key=get_residue_number)
    if not files:
        raise FileNotFoundError(f"No residue_*.dat files found in {residue_dir}")

    residues = []
    samples = []
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
            residues.append(residue)
            samples.append(values)

    if required and not residues:
        raise ValueError(
            f"No valid Q values found in the final {args.window:g} {args.time_unit} "
            f"for mutation {mutation}"
        )

    return residues, samples, data_unit


def plot_boxplots(mutation, residues, samples, window_title):
    fig_width = max(8, 0.53 * len(residues))
    fig, ax = plt.subplots(
        figsize=(fig_width, 5.5),
        layout="constrained",
    )

    ax.boxplot(
        samples, tick_labels=residues, showmeans=True, meanline=True,
        meanprops={"linestyle": "-", "color": "black"}, medianprops={"visible": False},
    )
    ax.set(
        xlabel="Residue number",
        ylabel="Per-residue Q",
        title=f"{mutation}: per-residue Q over the final {window_title}",
        ylim=(0, 1),
    )
    ax.tick_params(axis="x", labelrotation=90)
    return fig


def plot_comparison(means_by_mutation, window_title):
    residues = sorted({residue for means in means_by_mutation.values() for residue in means})
    positions = {residue: index for index, residue in enumerate(residues)}
    series = [(name, means) for name, means in means_by_mutation.items() if name != "WT"]
    markers = ("o", "s", "^", "D", "v", "P", "X", "h")
    colors = plt.get_cmap("tab10").colors
    per_panel = len(markers)
    panel_count = math.ceil(len(series) / per_panel)
    fig, axes = plt.subplots(
        panel_count, 1, sharex=True, sharey=True, squeeze=False,
        figsize=(max(8, 0.53 * len(residues)), max(4.5, 3.0 * panel_count)),
        layout="constrained",
    )
    bottom_ax = axes[-1, 0]
    bottom_ax.set_xticks(range(len(residues)), labels=residues)
    bottom_ax.set(xlim=(-0.5, len(residues) - 0.5), ylim=(-0.025, 1.025))

    for panel_index, ax in enumerate(axes[:, 0]):
        start = panel_index * per_panel
        group = [("WT", means_by_mutation.get("WT", {})), *series[start:start + per_panel]]
        # Give each mutation its own position within each residue's tick interval.
        spacing = 0.8 / len(group)
        for index, (mutation, means) in enumerate(group):
            offset = (index - (len(group) - 1) / 2) * spacing
            ax.plot(
                [positions[residue] + offset for residue in means], list(means.values()),
                linestyle="none", marker="*" if mutation == "WT" else markers[(index + panel_index) % per_panel],
                color="black" if mutation == "WT" else colors[index - 1], markersize=6 if mutation == "WT" else 3.6,
                markeredgewidth=0.35, label=mutation,
            )
        ax.grid(axis="x", alpha=0.18)
        ax.set_axisbelow(True)
        ax.tick_params(axis="x", labelbottom=True, labelrotation=90)
        for residue, tick in zip(residues, ax.get_xticklabels(), strict=True):
            if residue not in group[0][1]:
                ax.annotate("×", (0.5, 0), xycoords=tick, color="red", ha="center", va="top", annotation_clip=False)
        ax.legend(loc="lower left", bbox_to_anchor=(0, 1), ncols=min(4, len(group)),
                  frameon=False, borderaxespad=0)

    fig.supxlabel("Residue number")
    fig.supylabel("Mean per-residue Q")
    fig.suptitle(f"Mean per-residue Q over the final {window_title}")
    return fig, residues


def main():
    args = parse_arguments()
    mutations = get_mutations(args)
    window_title = f"{args.window:g} {args.time_unit}"
    window_label = f"{args.window:g}{args.time_unit}"

    if len(mutations) == 1:
        mutation = mutations[0]
        residues, samples, data_unit = load_residue_samples(args, mutation)
        fig = plot_boxplots(mutation, residues, samples, window_title)
        output_dir = args.output_dir or args.postprocessed_dir / mutation / "plots"
        filename = f"qvalue_per_residue_mean_variance_last_{window_label}.png"
        print(f"Mutation: {mutation}")
        print(f"Window: final {window_title} (trajectory time unit: {data_unit})")
    else:
        means_by_mutation = {}
        wt = ["WT"] if any((args.postprocessed_dir / "WT" / RESIDUE_SUBDIR).glob("residue_*.dat")) else []
        for mutation in [*(name for name in mutations if name != "WT"), *wt]:
            residues, samples, _ = load_residue_samples(args, mutation, required=mutation != "WT")
            means_by_mutation[mutation] = {
                residue: values.mean() for residue, values in zip(residues, samples, strict=True)
            }
        fig, residues = plot_comparison(means_by_mutation, window_title)
        output_dir = args.output_dir or SCRIPT_DIR / "output"
        filename = f"qvalue_per_residue_mean_comparison_last_{window_label}.png"
        print(f"Mutations compared: {len(mutations)}")
        print(f"Window: final {window_title}")

    output_dir.mkdir(parents=True, exist_ok=True)
    plot_path = output_dir / filename
    try:
        fig.savefig(plot_path, dpi=180)
    finally:
        plt.close(fig)

    print(f"Residues plotted: {len(residues)}")
    print(f"Saved plot: {plot_path}")


if __name__ == "__main__":
    main()
