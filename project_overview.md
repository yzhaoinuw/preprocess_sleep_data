# Project Overview

This document orients a new agent or collaborator to the active MATLAB preprocessing path. Keep it current when the main entry point, output contract, or alignment assumptions change.

## What This Repo Is

This repository preprocesses sleep-study signals for downstream use in the Sleep Scoring App. The active MATLAB pipeline extracts EEG/EMG from Viewpoint `.exp` files or TDT recordings, optionally extracts fiber photometry from TDT recordings, optionally imports manual sleep scores, aligns signals when needed, and writes Sleep Scoring App-ready `.mat` files.

The project currently targets MATLAB R2025a.

## Active Runtime Path

### 1. Entrypoint

[`preprocess_sleep_data.m`](preprocess_sleep_data.m)

- Main user-facing function for Viewpoint and/or TDT preprocessing.
- Accepts a required EEG/EMG path plus name-value options such as `ne_dir`, TDT stream/channel names, `sleep_score_file`, `save_path`, and `show_figure`. The legacy `video_alignment` option is still accepted for compatibility.
- Saves `eeg`, `emg`, `ne`, `sleep_scores`, `start_time`, `video_start_time`, `num_class`, `eeg_frequency`, `ne_frequency`, `video_name`, and `video_path`.

### 2. Viewpoint Loading

[`EEGtoolbox/loadEXP.m`](EEGtoolbox/loadEXP.m)

- Parses Viewpoint `.exp` metadata, including bin files, video files, start times, durations, and video frame/timestamp metadata when available.

[`ExtractContinuousData`](EEGtoolbox/)

- Used by `preprocess_sleep_data.m` to read EEG/EMG traces from Viewpoint bin data.

### 3. TDT Loading

`TDTbin2mat`

- Used to load TDT EEG/EMG streams and/or fiber photometry streams.
- For photometry, the active code expects 465 and 405 streams and optionally a TTL epoc channel when synchronizing separate Viewpoint EEG/EMG with TDT photometry.

### 4. Optional Sleep Scores

`readmatrix`

- Reads manual sleep score spreadsheets.
- The current code expects wake, SWS, and REM onset/duration columns in the layout documented by the existing examples.

## Repo Structure Map

```text
project_root/
|- AGENTS.md
|- README.md
|- project_overview.md
|- next_steps.md
|- work_log.md
|- preprocess_sleep_data.m
|- preprocess_sirenia.m
|- EEGtoolbox/
|- *_sketch.m and test_*.m exploratory scripts
```

## What Looks Active vs. Legacy

### Active / relevant now

- [`preprocess_sleep_data.m`](preprocess_sleep_data.m) - main preprocessing function for Viewpoint/TDT sleep data.
- [`preprocess_sirenia.m`](preprocess_sirenia.m) - separate active path for Pinnacle Sirenia plus TDT preprocessing.
- [`README.md`](README.md) - user-facing usage examples and output contract.
- [`next_steps.md`](next_steps.md) - active follow-up list, currently including video-alignment checks.
- [`work_log.md`](work_log.md) - newest session notes and verification breadcrumbs.

### Likely older or secondary

- `*_sketch.m` files - exploratory scripts and implementation notes; useful for context, but do not treat them as the active pipeline without checking.
- `LoadFPandEEG_Klaudia.m`, `Sleep_comp_3_groups.m`, and small `test_*.m` scripts - reference or ad hoc analysis scripts unless a task specifically names them.
- `EEGtoolbox/` - large helper toolbox. Read targeted files only when the active function calls them or metadata behavior needs confirmation.

## Video Alignment Mental Model

`video_start_time` is a signed offset in seconds that maps EEG/EMG time to video time:

```matlab
video_time = eeg_time + video_start_time;
```

As of `preprocess_sleep_data.m` version 0.2.8, single-bin, single-video Viewpoint recordings with no TDT photometry input compute the offset from metadata starts:

```matlab
video_start_time = (eeg_start_time - video_start_time_metadata) * 24 * 3600;
```

Negative values mean the video starts after the EEG/EMG trace, so early EEG/EMG samples have no matching video. Positive values mean the video starts before the EEG/EMG trace, so EEG/EMG time zero maps later into the source video. Multi-bin Viewpoint and TTL-aligned layouts still need targeted follow-up before changing their behavior.

## Tests And Fixtures

There is not yet a canonical automated test suite or fixture set documented for this repo. For now, use MATLAB `checkcode` as a lightweight syntax/style check after code edits, and record any data-specific smoke test commands in `work_log.md`.

## User Data Expectations

- Viewpoint EEG/EMG input: path to a `.exp` file.
- TDT EEG/EMG input: path to a TDT recording folder plus EEG/EMG stream names.
- Optional TDT photometry input: `ne_dir`, `chan_465`, `chan_405`, and sometimes `chan_ttl_pulse`.
- Optional manual scores: spreadsheet path passed as `sleep_score_file`.
- Output: one `.mat` file for single-bin or short recordings, or segmented `.mat` files for multi-bin Viewpoint data and long TDT recordings.

## Practical Mental Model

If you only want to understand the current product, read files in this order:

1. [`AGENTS.md`](AGENTS.md)
2. [`README.md`](README.md)
3. [`preprocess_sleep_data.m`](preprocess_sleep_data.m)
4. [`next_steps.md`](next_steps.md)
5. Targeted helper files under [`EEGtoolbox/`](EEGtoolbox/) only when needed.

## Questions Worth Clarifying Later

- What sample recordings should serve as canonical smoke-test fixtures?
- Should multi-bin Viewpoint exports use per-bin signed metadata offsets?
- For Viewpoint EEG/EMG plus separate TDT photometry, should `video_start_time` combine TTL-based EEG trimming with the Viewpoint video offset?
- Does the Sleep Scoring App correctly handle signed fractional-second `video_start_time` offsets and no-video regions?
