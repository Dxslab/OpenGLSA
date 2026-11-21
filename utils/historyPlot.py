#!/usr/bin/env python3
# Copyright (c) 2026 Ruiqing Yang and Baofeng Ma
# SPDX-License-Identifier: Apache-2.0

"""Monitor resid.history with NumPy and Matplotlib.

Run from the solver's working directory:
    python3 /path/to/utils/historyPlot.py
    python3 /path/to/utils/historyPlot.py other.history --refresh 5
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import FuncFormatter


def read_history(path):
    """Read whitespace-separated data, tolerating an unfinished final row."""
    rows = []
    columns = None
    try:
        with path.open() as stream:
            for line in stream:
                fields = line.split("#", 1)[0].split()
                if len(fields) < 2:
                    continue
                try:
                    row = [float(field.replace("D", "E").replace("d", "e"))
                           for field in fields]
                except ValueError:
                    continue
                if columns is None:
                    # Wait for a complete first row while the solver writes it.
                    if not line.endswith("\n"):
                        continue
                    columns = len(row)
                if len(row) == columns and line.endswith("\n"):
                    rows.append(row)
    except FileNotFoundError:
        pass  # The solver may not have created the file yet.
    return np.asarray(rows, dtype=float)


def residual_titles(columns):
    if columns == 6:
        return ["Continuity", "Momentum-X", "Momentum-Y", "Energy", "Total"]
    if columns == 7:
        return ["Continuity", "Momentum-X", "Momentum-Y", "Momentum-Z",
                "Energy", "Total"]
    return [f"Residual-{index}" for index in range(1, columns)]


def configure_axes(ax):
    ax.set_title("Residual History (auto-refresh)")
    ax.set_xlabel("Step")
    ax.set_ylabel("Residual")
    ax.set_yscale("log")
    ax.yaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:.0e}"))
    ax.grid(True)


def update_plot(ax, path, linewidth):
    data = read_history(path)
    if data.size == 0:
        return
    ax.clear()
    configure_axes(ax)
    for column, title in enumerate(residual_titles(data.shape[1]), start=1):
        residual = np.ma.masked_where(
            ~np.isfinite(data[:, column]) | (data[:, column] <= 0)
            | ~np.isfinite(data[:, 0]), data[:, column]
        )
        ax.plot(data[:, 0], residual, linewidth=linewidth, label=title)
    ax.legend(loc="upper right")


def positive_float(value):
    number = float(value)
    if not np.isfinite(number) or number <= 0:
        raise argparse.ArgumentTypeError("must be a finite positive number")
    return number


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("resfile", nargs="?", type=Path, default=Path("resid.history"))
    parser.add_argument("--refresh", type=positive_float, default=10,
                        help="refresh interval in seconds (default: 10)")
    parser.add_argument("--linewidth", type=positive_float, default=2,
                        help="line width (default: 2)")
    args = parser.parse_args()

    fig, ax = plt.subplots(figsize=(10, 6), dpi=100)
    configure_axes(ax)
    fig.tight_layout()
    update_plot(ax, args.resfile, args.linewidth)

    def refresh():
        update_plot(ax, args.resfile, args.linewidth)
        fig.canvas.draw_idle()

    timer = fig.canvas.new_timer(interval=max(1, round(args.refresh * 1000)))
    timer.add_callback(refresh)
    timer.start()
    fig.canvas.mpl_connect("close_event", lambda event: timer.stop())
    try:
        plt.show()
    finally:
        timer.stop()


if __name__ == "__main__":
    main()
