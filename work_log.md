# Work Log

Prepend new session notes to the top of this file.

## 2026-07-16

### Released absolute recording start timestamps as v0.2.9

- Added `recording_start_time` to both preprocessing outputs as an ISO-8601 character string identifying the first EEG/EMG sample saved in each MAT file.
- Derived the timestamp from Viewpoint bin metadata, TDT UTC block metadata, or the Sirenia EDF header, including per-segment, TTL-trim, and Sirenia interval offsets.
- Kept source timezone semantics explicit: TDT timestamps use a trailing `Z`; Viewpoint and Sirenia source clocks remain timezone-unspecified.
- Added a user-facing `CHANGELOG.md` covering the changes since `v0.2.7` and prepared the repository for the `v0.2.9` tag and GitHub Release.
- Verification:
  - `matlab -batch "issues1=checkcode('preprocess_sleep_data.m','-id'); issues2=checkcode('preprocess_sirenia.m','-id'); ..."` confirmed `checkcode_main=0` and `checkcode_sirenia=0`.
  - A real Sirenia EDF-only export with `interval = [10 20]` saved `recording_start_time = '2026-05-17T10:18:23.000'` as a character array, correctly advancing the EDF header time by 10 seconds.
  - A real single-bin Viewpoint export saved `recording_start_time = '2026-05-05T09:01:48.997'`, exactly matching `Info.BinFiles(1).TStart` formatted as ISO-8601.
  - `TDTbin2mat(..., 'T2', 1, 'NODATA', true)` on a real TDT block confirmed its date and UTC start metadata convert to `2023-09-28T08:38:30.000Z`.
  - `git diff --check`
  - `Get-Content README.md -Encoding Unicode | Select-String -Pattern 'recording_start_time' -Context 1,1`, plus a byte check confirming the UTF-16 LE BOM remains `FF FE`.

## 2026-07-02

### Added Sirenia EDF-only preprocessing

- Updated `preprocess_sirenia.m` so the first argument can be a Sirenia `.edf` file, and name-value-only calls using `edf_file` also parse.
- Detects whether a supplied folder contains TDT `.tsq` files. If not, the Sirenia EEG/EMG path skips TTL/TDT NE synchronization and saves `ne = []` with `ne_frequency = NaN`.
- Kept the existing TDT+EDF path guarded behind actual TDT block detection rather than adding TDT-only support or new segmentation behavior.
- Updated `README.md` to document the EDF-only Sirenia call form while preserving its UTF-16 LE BOM encoding.
- Verification:
  - `git diff --check -- preprocess_sirenia.m`
  - `git diff --check -- README.md work_log.md preprocess_sirenia.m`
  - `Get-Content -Path README.md -Encoding Unicode -TotalCount 48`
  - `matlab -batch "checkcode('preprocess_sirenia.m')"`
  - `matlab -batch "addpath(pwd); out=fullfile(pwd,'.codex_tmp','sirenia_edf_only_smoke.mat'); if exist(out,'file'), delete(out); end; preprocess_sirenia('C:\Users\yzhao\matlab_projects\sleep_data\Default_2026-05-17_10_18_13_export.edf','save_path',out,'show_figure',false); vars=whos('-file',out); fprintf('vars=%s\n', strjoin(string({vars.name}), ',')); load(out,'eeg','emg','ne','eeg_frequency','ne_frequency','video_start_time','video_name','video_path'); fprintf('eeg=%d emg=%d ne=%d eeg_frequency=%.6g ne_frequency=%.6g video_start_time=%.6g video_name_len=%d video_path_len=%d\n', numel(eeg), numel(emg), numel(ne), eeg_frequency, ne_frequency, video_start_time, strlength(string(video_name)), strlength(string(video_path))); delete(out);"` confirmed `eeg=3744000`, `emg=3744000`, `ne=0`, `eeg_frequency=400`, `ne_frequency=NaN`, `video_start_time=0`.
  - `matlab -batch "addpath(pwd); out=fullfile(pwd,'.codex_tmp','sirenia_edf_only_namevalue_smoke.mat'); if exist(out,'file'), delete(out); end; preprocess_sirenia('edf_file','C:\Users\yzhao\matlab_projects\sleep_data\Default_2026-05-17_10_18_13_export.edf','save_path',out,'show_figure',false,'interval',[0 10]); load(out,'eeg','emg','ne','eeg_frequency','ne_frequency','video_start_time'); fprintf('namevalue_eeg=%d namevalue_emg=%d namevalue_ne=%d eeg_frequency=%.6g ne_frequency=%.6g video_start_time=%.6g\n', numel(eeg), numel(emg), numel(ne), eeg_frequency, ne_frequency, video_start_time); delete(out);"` confirmed `namevalue_eeg=4001`, `namevalue_emg=4001`, `namevalue_ne=0`, `eeg_frequency=400`, `ne_frequency=NaN`, `video_start_time=0`.

<!--
Each session entry follows this shape:

## YYYY-MM-DD

### Short title for what was done

- bullet describing what was added or changed
- another bullet - keep them high-level and user/agent-facing, not implementation play-by-play
- if relevant, intended profiling signal or measurement:
  - what to look for in logs / output
  - what numbers were observed
- Verification:
  - the exact command(s) that were actually run
  - what passed / what was confirmed

Newest entry goes on top. If the session did multiple distinct pieces of work, use multiple `###` subsections under one `##` date header.
-->
