# Next Steps

Use this checklist alongside `work_log.md`.

## Currently Hot

Active threads - read these first to know what work is in flight:

- Video alignment follow-ups - single-bin Viewpoint metadata offset is implemented; confirm behavior for TTL-aligned and multi-bin layouts.

## Video Alignment Follow-Ups

Status: in progress

The current patch uses signed metadata offsets for single-bin, single-video Viewpoint exports with no TDT photometry input. Follow-up inspection showed that the end-aligned offsets from 2026-05-26 were artifacts of extraction truncation, not true video-vs-EEG start offsets.

Findings from `20260505_box3_M64_LongBaseline` and `20260505_box4_M65_LongBaseline`:

- `Info.BinFiles(1).TStart` and `Info.VideosFiles(1).Files(1).TStart` expose plausible absolute starts in the `.exp` metadata.
- The video starts after the EEG/EMG bin start: about `2.941 s` later for box3/M64 and `6.141 s` later for box4/M65.
- `Info.BinFiles(1).Duration` is computed from the actual `.bin` file size and sampling rate; for these files it is `11617.000 s`.
- The raw `.exp` acquisition duration is `11618.030 s`, so it should not be treated as more authoritative than the actual `.bin` length.
- `Info.VideosFiles(1).Files(1).Duration` comes from the `.exp` video duration in milliseconds. MATLAB `VideoReader` gave durations about one frame longer, so video duration metadata is close but not exact.
- `ExtractContinuousData(..., TimeRelEndSec = Inf, ...)` rounds the extraction end down to the last full minute. For a `11617 s` bin, that exports `11580 s`, dropping `37 s`. This is why end alignment produced `video_start_time = 35.088 s` and `31.889 s`; those offsets compare video duration to a truncated EEG trace.

Implemented single-bin Viewpoint approach:

- Stop using end alignment as the default alignment rule for single-bin Viewpoint data. It is brittle because exported EEG duration can be shortened by extraction behavior unrelated to video timing.
- Build alignment from absolute metadata: `Info.BinFiles(1).TStart` for EEG/EMG start and `Info.VideosFiles(1).Files(1).TStart` for video start.
- Keep EEG/EMG as the canonical downstream sleep-scoring timeline.
- Store `video_start_time = (eeg_start_time - video_start_time_metadata) * 24 * 3600`, so downstream can locate matching video time by `video_time = eeg_time + video_start_time`.
- Preserve signed fractional offsets. Negative values mean video starts after EEG/EMG and early EEG/EMG samples have no matching video; positive values mean video starts before EEG/EMG and EEG/EMG time zero maps later into the source video.
- For full single-bin exports, avoid the `ExtractContinuousData(..., Inf, ...)` minute truncation by passing an explicit end time based on actual `.bin` duration.

Implemented supported follow-up cases:

- For currently supported multi-bin Viewpoint exports with one Viewpoint video per bin, save a per-bin metadata-derived `video_start_time` instead of only setting the first segment.
- For the existing Viewpoint EEG/EMG plus TDT photometry TTL-sync path, add the rounded EEG TTL trim to the first saved segment's metadata-derived Viewpoint video offset.
- Keep other layouts out of scope unless a real supported fixture requires them; do not add speculative support for mismatched video/bin counts or separate Viewpoint video metadata for TDT-only EEG/EMG.

Remaining work:

- Implement a diagnostic summary for single-bin Viewpoint `.exp` files that reports bin/video start times, durations, actual extracted EEG duration, signed start offset, and no-video regions at the beginning or end of the EEG/EMG timeline.
- Confirm how the Sleep Scoring App interprets signed fractional-second `video_start_time` and displays no-video regions.
- If a real supported fixture appears with mismatched Viewpoint video/bin counts or multiple video groups, inspect it first and add only the minimum needed handling.

## Background / Paused

Sections below this line are older threads kept for context. They're not the current focus, but recording the state they were left in saves the next agent from re-deriving it.
