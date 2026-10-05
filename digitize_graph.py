#!/usr/bin/env python3
"""Digitize the Lotka–Volterra graph in a PNG and write a CSV.

The default calibration matches the included data.png:

    python digitize_graph.py data.png data.csv

The output columns are t,x,y:

* t: year
* x: hares (red curve)
* y: lynx (black curve)

For a different image with the same graph layout, adjust the pixel
calibration with --x-start, --x-end, --y-zero and --pixels-per-unit.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import numpy as np
from PIL import Image


def parse_years(value: str) -> np.ndarray:
    """Parse a year range such as 1900:1920."""
    try:
        start_text, end_text = value.split(":", 1)
        start, end = int(start_text), int(end_text)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(
            "years must use START:END, for example 1900:1920"
        ) from exc

    if end < start:
        raise argparse.ArgumentTypeError("year END must be greater than or equal to START")
    return np.arange(start, end + 1, dtype=int)


def load_rgb(path: Path) -> np.ndarray:
    """Load an image as an RGB NumPy array."""
    try:
        return np.asarray(Image.open(path).convert("RGB"))
    except FileNotFoundError as exc:
        raise SystemExit(f"Input image not found: {path}") from exc


def color_mask(image: np.ndarray, series: str) -> np.ndarray:
    """Return a mask for the red hare curve or black lynx curve."""
    red, green, blue = image[..., 0], image[..., 1], image[..., 2]

    if series == "hares":
        # Red curve and markers, tolerant of anti-aliased pixels.
        return (
            (red > 170)
            & (green < 120)
            & (blue < 120)
            & (red > green * 1.45)
            & (red > blue * 1.45)
        )

    # The lynx curve uses a dark neutral RGB color around (35, 31, 32).
    # Excluding pure black prevents the axes and text from being selected.
    channel_spread = image.max(axis=2) - image.min(axis=2)
    return (
        (red >= 20)
        & (red <= 70)
        & (green >= 15)
        & (green <= 60)
        & (blue >= 15)
        & (blue <= 60)
        & (channel_spread < 25)
    )


def extract_series(
    image: np.ndarray,
    years: np.ndarray,
    series: str,
    *,
    x_start: float,
    x_end: float,
    y_zero: float,
    pixels_per_unit: float,
    plot_top: int,
    plot_bottom: int,
    x_window: int,
    marker_radius: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Extract one curve by finding the densest marker neighborhood."""
    mask = color_mask(image, series)

    # The source graph contains a legend in this area. It has the same
    # colors as the curves, so remove it before sampling the data region.
    legend_x0, legend_x1 = 980, min(mask.shape[1], 1135)
    legend_y0, legend_y1 = 40, min(mask.shape[0], 160)
    mask[legend_y0:legend_y1, legend_x0:legend_x1] = False

    y0 = max(0, plot_top)
    y1 = min(mask.shape[0] - 1, plot_bottom)
    pixel_x = np.linspace(x_start, x_end, len(years))
    pixel_y = []
    scores = []
    integral = np.pad(mask.astype(np.int32), ((1, 0), (1, 0))).cumsum(0).cumsum(1)

    for x in pixel_x:
        center = int(round(x))
        candidate_x = range(
            max(marker_radius, center - x_window),
            min(mask.shape[1] - marker_radius, center + x_window) + 1,
        )
        candidate_y = np.arange(y0 + marker_radius, y1 - marker_radius + 1)
        best_score = -1
        detected_y = None

        for candidate_center_x in candidate_x:
            left = candidate_center_x - marker_radius
            right = candidate_center_x + marker_radius
            top = candidate_y - marker_radius
            bottom = candidate_y + marker_radius
            scores_for_y = (
                integral[bottom + 1, right + 1]
                - integral[top, right + 1]
                - integral[bottom + 1, left]
                + integral[top, left]
            )
            candidate_index = int(np.argmax(scores_for_y))
            candidate_score = int(scores_for_y[candidate_index])
            if candidate_score > best_score:
                best_score = candidate_score
                detected_y = float(candidate_y[candidate_index])

        if best_score < 10 or detected_y is None:
            raise RuntimeError(
                f"Could not detect the {series} curve near pixel x={center}. "
                "Adjust the calibration or color thresholds."
            )

        pixel_y.append(detected_y)
        scores.append(best_score)

    values = (y_zero - np.asarray(pixel_y)) / pixels_per_unit
    return values, np.asarray(scores)


def write_csv(path: Path, years: np.ndarray, hares: np.ndarray, lynx: np.ndarray) -> None:
    """Write the extracted observations using the project CSV format."""
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(("t", "x", "y"))
        for year, hare, lynx_value in zip(years, hares, lynx):
            writer.writerow((int(year), f"{hare:.2f}", f"{lynx_value:.2f}"))


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Digitize the hare and lynx curves from a graph PNG."
    )
    parser.add_argument("input_png", type=Path, help="PNG containing the plotted curves")
    parser.add_argument("output_csv", type=Path, help="CSV file to create")
    parser.add_argument(
        "--years",
        type=parse_years,
        default=parse_years("1900:1920"),
        metavar="START:END",
        help="year range (default: 1900:1920)",
    )
    parser.add_argument("--x-start", type=float, default=138.0, help="pixel x for first year")
    parser.add_argument("--x-end", type=float, default=1094.0, help="pixel x for last year")
    parser.add_argument(
        "--y-zero", type=float, default=452.0, help="pixel y corresponding to population 0"
    )
    parser.add_argument(
        "--pixels-per-unit",
        type=float,
        default=5.34,
        help="vertical pixels per population unit (default: 5.34)",
    )
    parser.add_argument("--plot-top", type=int, default=35, help="top pixel of plot area")
    parser.add_argument("--plot-bottom", type=int, default=454, help="bottom pixel of plot area")
    parser.add_argument(
        "--x-window",
        type=int,
        default=4,
        help="half-width of marker search window (default: 4)",
    )
    parser.add_argument(
        "--marker-radius",
        type=int,
        default=8,
        help="half-size of the marker search box (default: 8)",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="print detected pixel positions and values",
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    image = load_rgb(args.input_png)

    if args.pixels_per_unit <= 0:
        raise SystemExit("--pixels-per-unit must be greater than zero")
    if args.x_window < 1:
        raise SystemExit("--x-window must be at least 1")
    if args.marker_radius < 1:
        raise SystemExit("--marker-radius must be at least 1")

    hares, hare_scores = extract_series(
        image,
        args.years,
        "hares",
        x_start=args.x_start,
        x_end=args.x_end,
        y_zero=args.y_zero,
        pixels_per_unit=args.pixels_per_unit,
        plot_top=args.plot_top,
        plot_bottom=args.plot_bottom,
        x_window=args.x_window,
        marker_radius=args.marker_radius,
    )
    lynx, lynx_scores = extract_series(
        image,
        args.years,
        "lynx",
        x_start=args.x_start,
        x_end=args.x_end,
        y_zero=args.y_zero,
        pixels_per_unit=args.pixels_per_unit,
        plot_top=args.plot_top,
        plot_bottom=args.plot_bottom,
        x_window=args.x_window,
        marker_radius=args.marker_radius,
    )

    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    write_csv(args.output_csv, args.years, hares, lynx)

    if args.verbose:
        for year, hare, lynx_value, hare_score, lynx_score in zip(
            args.years, hares, lynx, hare_scores, lynx_scores
        ):
            print(
                f"{year}: hares={hare:.2f}, lynx={lynx_value:.2f} "
                f"(scores={hare_score}/{lynx_score})"
            )
    else:
        print(f"Wrote {len(args.years)} observations to {args.output_csv}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
