import re
import sys
from collections import defaultdict
from heapq import nsmallest
from math import isfinite
from pathlib import Path

import pandas as pd
from matplotlib import pyplot as plt


PIPELINE_DIR = Path(__file__).resolve().parents[1]
MUTATION_ID = re.compile(r"[A-Za-z][A-Za-z0-9]([0-9]+)[A-Za-z]")
SMOOTHING_WINDOW = 5

# Helices and strands from Fig. 1 (p. 366) of yeast_cofilin.pdf.
SECONDARY_STRUCTURE = {
    "α1": range(10, 22),
    "α2": range(50, 56),
    "α3": range(95, 112),
    "α4": range(129, 138),
    "β1": range(7, 10),
    "β2": range(27, 34),
    "β3": range(38, 46),
    "β4": range(63, 74),
    "β5": range(78, 90),
    "β6": range(117, 122),
}

# Conserved residues are listed on p. 367 of yeast_cofilin.pdf.
RESIDUE_GROUPS = (
    (
        "Hydrophobic cores",
        "hydrophobic_cores.png",
        (8, 16, 19, 28, 30, 32, 39, 53, 57, 83, 86, 108, 119, 127, 132),
    ),
    (
        "Highly conserved basic residues",
        "conserved_basic_residues.png",
        (20, 80, 96, 98, 109),
    ),
    ("β4", "beta4.png", SECONDARY_STRUCTURE["β4"]),
    ("β5", "beta5.png", SECONDARY_STRUCTURE["β5"]),
)


def score_positions(input_file: Path) -> dict[int, float]:
    """Sum up to three negative energies per residue and reverse the sign."""
    energies_by_residue = defaultdict(list)
    with input_file.open(encoding="utf-8") as source:
        for line_number, line in enumerate(source, start=1):
            if not line.strip():
                continue
            try:
                mutation, energy_text = line.split()
                match = MUTATION_ID.fullmatch(mutation)
                energy = float(energy_text)
                if match is None or not isfinite(energy):
                    raise ValueError
            except ValueError as exc:
                raise ValueError(
                    f"{input_file}:{line_number}: invalid mutation row"
                ) from exc

            negative_energies = energies_by_residue[int(match.group(1))]
            if energy < 0:
                negative_energies.append(energy)

    if not energies_by_residue:
        raise ValueError(f"{input_file}: no mutation rows found")

    return {
        residue: -sum(nsmallest(3, energies))
        for residue, energies in sorted(energies_by_residue.items())
    }


def minmax_normalize(values: pd.Series) -> pd.Series:
    """Scale values to [0, 1], mapping a constant series to zero."""
    minimum = values.min()
    span = values.max() - minimum
    return (values - minimum) / (span if span else 1)


def smooth(values: pd.Series) -> pd.Series:
    return values.rolling(SMOOTHING_WINDOW, center=True, min_periods=1).mean()


def plot_residue_groups(scores: pd.Series, rmsf: pd.Series, output_dir: Path):
    from common.plot import finish_graph

    for title, filename, residues in RESIDUE_GROUPS:
        group_scores = scores.loc[list(residues)]
        group_rmsf = rmsf.loc[list(residues)]
        positions = range(len(group_scores))
        width = 0.38
        fig, ax = plt.subplots(figsize=(8, 5))
        try:
            ax.bar(
                [x - width / 2 for x in positions],
                group_scores,
                width,
                label="Mutation score (−Σ ΔΔG)",
            )
            ax.bar(
                [x + width / 2 for x in positions],
                group_rmsf,
                width,
                label="WT Cα RMSF",
            )
            ax.set_xticks(positions, labels=group_scores.index)
            ax.set_xlabel("Residue number")
            ax.set_ylabel("Normalized value")
            ax.set_ylim(0, 1.05)
            ax.set_title(title)
            ax.set_axisbelow(True)
            ax.grid(axis="y", alpha=0.3)
            ax.legend()
            finish_graph(fig, output_dir / filename)
        finally:
            plt.close(fig)


def plot_wt_rmsf_by_structure(rmsf: pd.Series, output_dir: Path):
    from common.plot import finish_graph

    rmsf = rmsf.sort_index()
    structured_residues = {
        residue for segment in SECONDARY_STRUCTURE.values() for residue in segment
    }
    in_structure = rmsf.index.isin(structured_residues)
    groups = (
        (in_structure, "α-helices and β-strands", "wt_rmsf_secondary_structure.png"),
        (
            ~in_structure,
            "Outside α-helices and β-strands",
            "wt_rmsf_outside_secondary_structure.png",
        ),
    )
    for mask, label, filename in groups:
        fig, ax = plt.subplots(figsize=(8, 5))
        try:
            ax.plot(rmsf.index, rmsf.where(mask), linewidth=1.5, marker="o", markersize=3)
            ax.set_xlim(rmsf.index.min() - 1, rmsf.index.max() + 1)
            ax.set_ylim(0, rmsf.max() * 1.05 or 1)
            ax.set_xlabel("Residue number")
            ax.set_ylabel("RMSF (Å)")
            ax.set_title(f"WT Cα RMSF · {label}")
            ax.grid(alpha=0.3)
            finish_graph(fig, output_dir / filename)
        finally:
            plt.close(fig)


def main():
    sys.path.insert(0, str(PIPELINE_DIR.parent))
    from common.plot import create_graph, read_data

    output_dir = PIPELINE_DIR / "output"
    residue_label = "Residue number"
    score_label = "−Σ ΔE (up to 3 lowest negative values)"
    scores = score_positions(PIPELINE_DIR / "sorted_positions_single.txt")
    mutation_scores = pd.Series(scores, name=score_label).rename_axis(residue_label)
    create_graph(
        smooth(mutation_scores).reset_index(),
        f"Promising positions ({SMOOTHING_WINDOW}-residue moving average)",
        output_dir / "promising_positions.png",
    )

    rmsf_path = output_dir / "postprocessed" / "WT" / "analysis" / "rmsfCA.dat"
    rmsf = read_data(rmsf_path).set_index("Residue")["RMSF (Å)"]
    if not rmsf.index.is_unique:
        raise ValueError(f"{rmsf_path}: duplicate residue numbers")
    normalized_scores, normalized_rmsf = minmax_normalize(mutation_scores).align(
        minmax_normalize(rmsf), join="inner"
    )
    if normalized_scores.empty:
        raise ValueError("Mutation scores and WT RMSF have no residues in common")

    delta = (normalized_scores - normalized_rmsf).abs()
    delta = delta.rename_axis(residue_label).rename("Absolute difference (0–1)")
    create_graph(
        smooth(delta).reset_index(),
        f"Normalized mutation score vs. WT Cα RMSF ({SMOOTHING_WINDOW}-residue moving average)",
        output_dir / "delta.png",
    )
    plot_residue_groups(normalized_scores, normalized_rmsf, output_dir)
    plot_wt_rmsf_by_structure(rmsf, output_dir)


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError) as exc:
        raise SystemExit(f"Error: {exc}")
