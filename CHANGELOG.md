# Changelog

This file summarizes user-facing changes in each tagged release.

## v0.2.9 - 2026-07-16

### Added

- MAT outputs now include `recording_start_time`, an ISO-8601 character string for the first EEG/EMG sample actually saved in each file.
- Absolute start timestamps are read from Viewpoint bin metadata, TDT UTC block metadata, or Sirenia EDF headers and adjusted for synchronization, trimming, and output segmentation.
- `preprocess_sirenia` can process a Sirenia EDF file without requiring TDT photometry data.
- `preprocess_sirenia` supports optional interval extraction.

### Changed

- Supported Viewpoint exports use signed, metadata-derived `video_start_time` offsets, including per-bin offsets and the existing Viewpoint-plus-TDT synchronization path.
- Sirenia outputs consistently include empty NE fields when no TDT block is supplied, keeping the MAT-file schema compatible with downstream consumers.
- Output documentation now describes absolute recording timestamps, signed video offsets, EDF-only preprocessing, and source timezone behavior.

### Timestamp format

- TDT timestamps are UTC and end in `Z`, for example `2023-09-28T08:38:30.000Z`.
- Viewpoint and Sirenia timestamps preserve their source clock without inventing a timezone, for example `2026-05-17T10:18:23.000`.

## v0.2.7 - 2025-02-27

- Added video metadata to MAT outputs, including `video_name`, `video_path`, and `video_start_time`.
