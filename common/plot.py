import argparse
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def read_data(filename: Path) -> pd.DataFrame:
    """Read a tab-separated analysisMD table with a header."""
    data = pd.read_csv(filename, sep="\t", comment="#", encoding="utf-8")
    if data.shape[1] < 2 or data.empty:
        raise ValueError(
            f"{filename}: expected a header and at least two columns of data."
        )
    try:
        data = data.apply(pd.to_numeric, errors="raise")
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{filename}: contains non-numeric data.") from exc
    if not np.isfinite(data.to_numpy(dtype=float)).all():
        raise ValueError(f"{filename}: contains missing or non-finite values.")
    return data


def finish_graph(fig, output_file=None):
    """Display or save a figure."""
    fig.tight_layout()
    if output_file is None:
        plt.show()
    else:
        output_file = Path(output_file)
        output_file.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_file, dpi=300, bbox_inches="tight")
        print(f"Saved: {output_file}")


def create_graph(data, title=None, output_file=None):
    """Plot one dependent column against the first column."""

    x_label, y_label = data.columns
    fig, ax = plt.subplots(figsize=(8, 5))

    try:
        ax.plot(
            data[x_label],
            data[y_label],
            linewidth=1.5,
            marker="o" if len(data) == 1 else None,
        )
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)
        ax.set_title(title if title is not None else y_label)
        ax.grid(alpha=0.3)
        finish_graph(fig, output_file)
    finally:
        plt.close(fig)


def compare_graph(datasets, labels, title=None, output_file=None):
    """Overlay compatible datasets, preserving each dataset's x coordinates."""

    if len(datasets) != len(labels) or not datasets:
        raise ValueError("Provide exactly one legend label per dataset.")

    x_label, y_label = datasets[0].columns
    fig, ax = plt.subplots(figsize=(8, 5))

    try:
        for data, label in zip(datasets, labels):
            if list(data.columns) != [x_label, y_label]:
                raise ValueError("Comparison requires matching column names and units.")
            ax.plot(
                data[x_label],
                data[y_label],
                label=label,
                linewidth=1.5,
                marker="o" if len(data) == 1 else None,
            )
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)
        ax.set_title(title if title is not None else f"{y_label} comparison")
        ax.legend()
        ax.grid(alpha=0.3)
        finish_graph(fig, output_file)
    finally:
        plt.close(fig)


def existing_file(path_str: str) -> Path:
    """Validate that the path exists and is a regular file."""
    path = Path(path_str)
    if not path.exists():
        raise argparse.ArgumentTypeError(f"File '{path_str}' does not exist.")
    if not path.is_file():
        raise argparse.ArgumentTypeError(f"Path '{path_str}' is not a regular file.")
    return path


class VerboseHelpFormatter(
    argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter
):
    pass


def parse_args():
    """Configure and parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Plot analysisMD tables independently or compare matching datasets.",
        formatter_class=VerboseHelpFormatter,
        epilog="""Examples:
  python plot.py rmsd.dat
  python plot.py rmsd.dat -o rmsd.png
  python plot.py analysis/equilibration/*.dat -o plots/equilibration
  python plot.py alpha.dat --title "Cofilin at 300 K" -o plots
  python plot.py 300K/rmsd.dat 400K/rmsd.dat --compare --labels "300 K" "400 K" -o plots
""",
    )

    parser.add_argument(
        "files",
        type=existing_file,
        metavar="DAT_FILE",
        nargs="+",
        help="One or more tab-separated .dat files.",
    )

    parser.add_argument(
        "--title",
        type=str,
        help="Custom title (prefix when generating multiple plots).",
    )

    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        help="Image path for one plot, or directory for separate PNGs. "
        "Without this option, display plots.",
    )

    parser.add_argument(
        "--compare",
        action="store_true",
        help="Overlay input files with matching column names, one figure per metric.",
    )

    parser.add_argument(
        "--labels",
        nargs="+",
        metavar="LABEL",
        help="Legend labels in input-file order (comparison mode only).",
    )

    args = parser.parse_args()

    if args.labels is not None and not args.compare:
        parser.error("--labels requires --compare.")

    if args.compare and len(args.files) < 2:
        parser.error("--compare requires at least two input files.")

    if args.labels is not None and len(args.labels) != len(args.files):
        parser.error("Provide exactly one --labels value per input file.")

    return args


def safe_name(label):
    """Convert a column name into a filename-safe component."""

    return re.sub(r"[^\w.-]+", "_", label).strip("_.") or "metric"


def main():
    args = parse_args()

    # Read and validate every input before generating any figures.
    tables = [read_data(filename) for filename in args.files]
    jobs = []
    labels = None

    if args.compare:
        columns = list(tables[0].columns)

        for filename, data in zip(args.files, tables):
            if list(data.columns) != columns:
                raise ValueError(
                    f"{filename}: comparison requires matching headers and units. "
                    f"Expected {columns}, got {list(data.columns)}. "
                    "Compare one type of analysis at a time."
                )

        labels = args.labels

        if labels is None:
            labels = [p.resolve().parent.name for p in args.files]
            if len(set(labels)) != len(labels):
                labels = [str(p) for p in args.files]

        if len(set(labels)) != len(labels):
            raise ValueError("Legend labels must be distinct; use --labels.")

        x_label = columns[0]

        for y_label in columns[1:]:
            datasets = [data[[x_label, y_label]] for data in tables]
            jobs.append((args.files[0], datasets, y_label, len(columns) > 2))

    else:
        for filename, data in zip(args.files, tables):
            x_label = data.columns[0]
            for y_label in data.columns[1:]:
                jobs.append(
                    (filename, data[[x_label, y_label]], y_label, data.shape[1] > 2)
                )

    # A path ending in a known image extension represents a single figure.
    # Other output paths are treated as directories.
    image_extensions = {
        ".png",
        ".pdf",
        ".svg",
        ".jpg",
        ".jpeg",
        ".tif",
        ".tiff",
        ".webp",
        ".eps",
        ".ps",
    }

    output = args.output

    single_image = (
        output is not None
        and not output.is_dir()
        and output.suffix.lower() in image_extensions
    )

    if single_image and len(jobs) != 1:
        raise ValueError("Use an output directory when generating multiple plots.")

    if output is not None and not single_image:
        output.mkdir(parents=True, exist_ok=True)

    used_paths = set()

    for filename, data, y_label, multiple_columns in jobs:
        title = args.title
        if title is not None and len(jobs) > 1:
            title = f"{title} — {y_label}"

        output_file = None
        if single_image:
            output_file = output
        elif output is not None:
            stem = filename.stem

            if multiple_columns:
                stem = f"{stem}_{safe_name(y_label)}"

            if args.compare:
                stem += "_comparison"

            # Avoid overwriting another plot in this invocation when input
            # files from different directories share the same basename.
            output_file = output / f"{stem}.png"
            index = 2

            while output_file in used_paths:
                output_file = output / f"{stem}_{index}.png"
                index += 1

            used_paths.add(output_file)

        if args.compare:
            compare_graph(data, labels, title, output_file)
        else:
            create_graph(data, title, output_file)


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError) as exc:
        raise SystemExit(f"Error: {exc}")
