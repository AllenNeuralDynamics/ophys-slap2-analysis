"""
Script to generate a StimInfo.csv table from photodiode, frameClk, and orientation data files.

This script searches for the necessary CSV files in a specified directory, find trial start and stop times,
and outputs a table containing trial information, including orientation, stimulus onset, and offset times.

Usage:
    python StimInfo.py -dr /path/to/your/behavior/directory

Arguments:
    -dr, --directory  Path to the main directory containing the CSV files.
"""

import numpy as np
import pandas as pd
import os
from pathlib import Path
import argparse

def generate_table(main_path):
    """
    Generate a StimInfo.csv file with trial information.
    
    Args:
        main_path (Path): The main directory containing photodiode, clock, and orientation CSV files.
    """
    NSTIM = 8  # Number of stimuli per trial
    
    # Find the relevant CSV files
    for i in range(100):
        photodiode_fp = main_path / f'photodiode_{i}.csv'
        clock_fp = main_path / f'frameClk_{i}.csv'
        orientation_fp = main_path / f'orientations_{i}.csv'
        if os.path.exists(photodiode_fp) and os.path.exists(clock_fp) and os.path.exists(orientation_fp):
            break
    else:
        raise FileNotFoundError("Required CSV files not found in the specified directory.")

    print('Found the following files:\n', photodiode_fp, '\n', clock_fp, '\n', orientation_fp)   

    # Read CSV files into pandas DataFrames and rename columns
    photodiode = pd.read_csv(photodiode_fp)
    clock = pd.read_csv(clock_fp, header=None)
    orientation = pd.read_csv(orientation_fp, header=None)
    clock = clock.rename(columns={0: 'value', 1: 'time'})
    photodiode = photodiode.rename(columns={'Item2': 'value', 'Item1': 'time'})
    stim_rad = list(orientation.iloc[:, 0])

    # Convert radians to degrees and round to the nearest multiple of 45
    stim_deg = np.round(np.degrees(stim_rad) / 45) * 45
    stim_deg = [int(deg) for deg in stim_deg]

    # Process clock and photodiode data to determine trial and stimulus times
    if clock.shape[0] > 1:
        # Calculate trial start times
        diff = np.diff(clock.time)
        stop_ind = np.where(diff > 0.05)[0]
        start_ind = stop_ind + 1
        start_t = np.array(clock.time)[start_ind]

        # Calculate stimulus onset and offset times
        diode_values = np.array(photodiode.value, dtype=float)
        diode_values[diode_values < 400] = 0
        diode_values[diode_values >= 400] = 1
        diode_diff = np.diff(diode_values)
        stim_onset_ind = np.where(diode_diff > 0.5)[0]
        stim_offset_ind = np.where(diode_diff < -0.5)[0]
        stim_onset_t = np.array(photodiode.time)[stim_onset_ind]
        stim_offset_t = np.array(photodiode.time)[stim_offset_ind]
    else:
        print('SLAP2 clock data is missing, assuming trial structure of 1s OFF-2s ON after start of each trial')

    # Create the StimInfo table
    stim_df = pd.DataFrame()
    for trial_idx in range(1, 1 + len(stim_deg) // NSTIM):
        trial_degrees = stim_deg[NSTIM * (trial_idx - 1):NSTIM * trial_idx]

        if clock.shape[0] > 1:
            trial_start_time = start_t[trial_idx]
            trial_stim_onset_times = stim_onset_t[NSTIM * (trial_idx - 1):NSTIM * trial_idx] - trial_start_time
            trial_stim_offset_times = stim_offset_t[NSTIM * (trial_idx - 1):NSTIM * trial_idx] - trial_start_time
        else:
            trial_start_time = 0
            trial_stim_onset_times = [1, 4, 7, 10, 13, 16, 19, 22]
            trial_stim_offset_times = 2 + np.array(trial_stim_onset_times)

        _df = pd.DataFrame({
            'trial': NSTIM * [trial_idx + 1],
            'orientation': trial_degrees,
            'onset(s)': trial_stim_onset_times,
            'offset(s)': trial_stim_offset_times
        })
        stim_df = pd.concat([stim_df, _df], ignore_index=True)

    # Display the DataFrame and save to CSV
    print(stim_df)
    stim_df.to_csv(main_path / 'StimInfo.csv', index=False)
    print(f"Table saved in {main_path / 'StimInfo.csv'}")

def main():
    parser = argparse.ArgumentParser(description="Generate StimInfo.csv from photodiode, clock, and orientation data.")
    parser.add_argument('-dr', '--directory', type=str, required=True, help="Path to the main directory containing the behavior CSV files.")
    args = parser.parse_args()

    main_path = Path(args.directory)
    generate_table(main_path)

if __name__ == "__main__":
    main()

