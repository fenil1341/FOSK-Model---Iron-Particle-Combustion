import re
import sys
import argparse
import matplotlib.pyplot as plt
import numpy as np
import os
from scipy.interpolate import interp1d
from pathlib import Path

def parse_file(path):
    with open(path, "r") as f:
        text = f.read()

    title_match = re.search(r'\(title\s+"([^"]*)"\)', text)
    title = title_match.group(1) if title_match else ""
    # title = "Simulation"

    labels_match = re.search(r'\(labels\s+"([^"]*)"\s+"([^"]*)"\)', text)
    xlabel, ylabel = labels_match.groups() if labels_match else ("X", "Y")

    series = {}
    for block_match in re.finditer(
        r'\(\(xy/key/label\s+"([^"]*)"\)(.*?)\)', text, re.DOTALL
    ):
        name = block_match.group(1)
        body = block_match.group(2)
        nums = [float(n) for n in re.findall(r'[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?', body)]
        xs = nums[0::2]
        ys = nums[1::2]
        series[name] = (xs, ys)

    return title, xlabel, ylabel, series


def shift_time_to_temperature(xs, ys, target_temperature=950.0):
    xs = np.asarray(xs, dtype=float)
    ys = np.asarray(ys, dtype=float)

    exact = np.flatnonzero(ys == target_temperature)
    if exact.size:
        crossing_index = exact[0]
        crossing_x = xs[crossing_index]
    else:
        crossings = np.flatnonzero(
            (ys[:-1] < target_temperature) & (ys[1:] >= target_temperature)
        )
        if not crossings.size:
            raise ValueError(
                f"Temperature never reaches {target_temperature:g} K."
            )

        index = crossings[0]
        crossing_x = np.interp(
            target_temperature,
            [ys[index], ys[index + 1]],
            [xs[index], xs[index + 1]],
        )

    return xs - crossing_x


def main():
    input_file_path = r"C:\Users\fp25tevo\WORK\Simulations\udfTest\UDF\New\Plot_data\first_plateau.xy"
    
    # 1. Load Fluent Data
    title, xlabel, ylabel, series = parse_file(input_file_path)

    if not series:
        print("No xy/key/label series found in file.", file=sys.stderr)
        sys.exit(1)

    # 2. Setup Paper Data (Converted ms to seconds to match Fluent)
    paper_time_s = np.array([-1.7, -1.0, 0.0, 0.5, 1.0, 1.4, 1.55, 1.6, 2.0, 2.1, 3.0, 4.0, 5.0, 6.0, 7.0, 7.8, 8.0, 10.0, 12.0, 15.0]) / 1000.0
    paper_temp = np.array([300, 480, 850, 950, 1080, 1400, 1642, 1642, 1811, 1811, 2180, 2600, 2950, 3180, 3280, 3350, 3100, 2480, 2180, 1900])
    
    f_interp = interp1d(paper_time_s, paper_temp, kind='linear')
    smooth_paper_time = np.linspace(paper_time_s.min(), paper_time_s.max(), 500)
    smooth_paper_temp = f_interp(smooth_paper_time)
    
    # Shift paper data using the exact same function so 900 K sits perfectly at t=0
    shifted_paper_time = shift_time_to_temperature(smooth_paper_time, smooth_paper_temp, target_temperature=1000.0)

    # 3. Setup Plot
    fig, ax = plt.subplots(figsize=(10, 5)) # Slightly taller to accommodate the inset clearly
    axins = ax.inset_axes([0.45, 0.2, 0.4, 0.4]) 

    # 4. Shift and Plot Fluent Series
    shifted_series = {}
    for name, (xs, ys) in series.items():
        try:
            shifted_series[name] = (shift_time_to_temperature(xs, ys), ys)
        except ValueError as error:
            print(f"Skipping {name}: {error}", file=sys.stderr)

    if not shifted_series:
        print("No series reaches 1000 K.", file=sys.stderr)
        sys.exit(1)

    for name, (xs, ys) in shifted_series.items():
        ax.plot(xs, ys, linewidth=1.5, alpha=1, label=name, color='b')
        axins.plot(xs, ys, linewidth=1.5, alpha=1, color='b')

    # 5. Plot the Paper Profile
    ax.plot(shifted_paper_time, smooth_paper_temp, linewidth=2, color='r', label='FOSK (Paper)')
    axins.plot(shifted_paper_time, smooth_paper_temp, linewidth=2, color='r', label='Simulation')

    # 6. Format Main Plot
    ax.set_title(title if title else "Temperature vs. Time", fontsize=18)
    ax.set_xlabel(xlabel + " (t) [s]", fontsize=14)
    ax.set_ylabel(ylabel + ' $(T_{p})$ [K]', fontsize=14)
    
    all_shifted_xs = [x for xs, _ in shifted_series.values() for x in xs] + list(shifted_paper_time)
    ax.set_xlim(min(all_shifted_xs), 0.015) 
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=10, loc="upper left")

    # 7. Format Inset Plot
    axins.set_ylim(1500, 1900)
    axins.set_xlim(7e-4, 2.1e-3) # Adjusted slightly to ensure both plateaus remain visible
    axins.grid(True, alpha=0.3)
    
    ax.indicate_inset_zoom(axins, edgecolor="black")

    fig.tight_layout()

    output_path = os.path.dirname(input_file_path)
    file_stem = Path(input_file_path).stem
    file_name = f"{file_stem}.png"
    # final_save_path = os.path.join(output_path, "Temp_vs_time_plateau_both_radON_moreFeO_correct_Hf.png")
    final_save_path = r"C:\Users\fp25tevo\WORK\Simulations\udfTest\UDF\New\Plot_data"
    
    print(f"Saving plot to: {final_save_path}")
    final_save_path = os.path.join(final_save_path, file_stem)

    plt.savefig(final_save_path, bbox_inches="tight", dpi=300)
    # plt.show()

if __name__ == "__main__":
    main()