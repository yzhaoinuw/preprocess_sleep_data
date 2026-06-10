# Work Log

Prepend new session notes to the top of this file.

## 2026-06-10

### Prepared treaty docs and release handoff

- Added the official tri-color Agent Collab Treaty adoption badge to `README.md`.
- Confirmed the treaty documentation set is present in this repo: `AGENTS.md`, `project_overview.md`, `next_steps.md`, `work_log.md`, and `work_log_archive/README.md`.
- Verification:
  - `Get-Content -Path README.md -Encoding Unicode -TotalCount 12`
  - `Get-Content -Path work_log_archive/README.md`
  - `git diff --check`
  - `matlab -batch "addpath(pwd); m=checkcode('preprocess_sleep_data.m','-id'); if isempty(m), disp('preprocess_sleep_data checkcode clean'); else, for k=1:numel(m), fprintf('%d:%s:%s\n',m(k).line,m(k).id,m(k).message); end; end"`
  - `matlab -batch "m=checkcode('test_preprocess_func.m','-id'); if isempty(m), disp('test_preprocess_func checkcode clean'); else, for k=1:numel(m), fprintf('%d:%s:%s\n',m(k).line,m(k).id,m(k).message); end; end"`

## 2026-06-01

### Exported all box3/box4 Viewpoint files

- Updated `test_preprocess_func.m` to batch-export all `box3` and `box4` folders under `C:\Users\yzhao\matlab_projects\sleep_data_extraction`.
- Wrote one output `.mat` per folder directly under `sleep_data_extraction`, named `<folder>.mat`; existing files with those names were overwritten.
- Verified the eight requested outputs have signed `video_start_time` values and full actual EEG/EMG durations.
- Verification:
  - `matlab -batch "run('C:\Users\yzhao\matlab_projects\preprocess_sleep_data\test_preprocess_func.m')"`
  - `matlab -batch "base='C:\Users\yzhao\matlab_projects\sleep_data_extraction'; files=dir(fullfile(base,'*.mat')); names={files.name}; keep=~cellfun('isempty', regexp(names,'box[34]','once')) & ~contains(names,{'_end_alignment','_start_alignment'}); files=files(keep); [~,idx]=sort({files.name}); files=files(idx); fprintf('file\tvideo_start_time\teeg_duration_s\teeg_len\teeg_frequency\tvideo_name\n'); for i=1:numel(files), p=fullfile(files(i).folder,files(i).name); s=load(p); fprintf('%s\t%.6f\t%.6f\t%d\t%.6f\t%s\n', files(i).name, s.video_start_time, numel(s.eeg)/s.eeg_frequency, numel(s.eeg), s.eeg_frequency, string(s.video_name)); end"`
  - `matlab -batch "m=checkcode('test_preprocess_func.m','-id'); if isempty(m), disp('test_preprocess_func checkcode clean'); else, for k=1:numel(m), fprintf('%d:%s:%s\n',m(k).line,m(k).id,m(k).message); end; end"`

### Implemented signed Viewpoint video offsets

- Updated `preprocess_sleep_data.m` so single-bin, single-video Viewpoint-only exports preserve EEG/EMG as the canonical timeline and save `video_start_time` as a signed offset where `video_time = eeg_time + video_start_time`.
- Removed the end-alignment behavior from the single-bin Viewpoint path while keeping the legacy `video_alignment` argument accepted for compatibility.
- Avoided the `ExtractContinuousData(..., Inf, ...)` last-full-minute truncation for single-bin Viewpoint exports by passing the actual `.bin` duration as the extraction end time.
- Updated `README.md`, `AGENTS.md`, `project_overview.md`, and `next_steps.md` to document the signed-offset convention and leave TTL/multi-bin alignment as follow-up work.
- Verification:
  - `matlab -batch "addpath(pwd); m=checkcode('preprocess_sleep_data.m','-id'); if isempty(m), disp('checkcode clean'); else, for k=1:numel(m), fprintf('%d:%s:%s\n',m(k).line,m(k).id,m(k).message); end; end"`
  - `matlab -batch "addpath('C:\Users\yzhao\matlab_projects\preprocess_sleep_data'); addpath('C:\Users\yzhao\matlab_projects\preprocess_sleep_data\EEGtoolbox'); exps={'C:\Users\yzhao\matlab_projects\sleep_data_extraction\20260505_box3_M64_LongBaseline\20260505_box3_M64_LongBaseline_2026-05-05_09-01-47-463.exp','C:\Users\yzhao\matlab_projects\sleep_data_extraction\20260505_box4_M65_LongBaseline\20260505_box4_M65_LongBaseline_2026-05-05_09-01-47-463.exp'}; expected=[-2.9411244,-6.1406189]; for i=1:numel(exps), out=[tempname '.mat']; preprocess_sleep_data(exps{i},'save_path',out); s=load(out); fprintf('\n%s\n', exps{i}); fprintf('video_start_time=%.6f expected_about=%.6f\n', s.video_start_time, expected(i)); fprintf('eeg_duration=%.6f eeg_len=%d fs=%.6f\n', numel(s.eeg)/s.eeg_frequency, numel(s.eeg), s.eeg_frequency); delete(out); end"`
  - Smoke outputs: box3 `video_start_time = -2.940992`, box4 `video_start_time = -6.140995`; both exported `11617.000000` seconds of EEG/EMG.

### Documented single-bin Viewpoint alignment diagnosis

- Updated `next_steps.md` with the finding that the previous 31-35 second end-alignment offsets came from `ExtractContinuousData(..., Inf, ...)` truncating a 11617 second bin down to 11580 seconds.
- Recorded the single-bin Viewpoint recommendation: align from absolute `Info.BinFiles(1).TStart` and `Info.VideosFiles(1).Files(1).TStart`, use actual file-duration coverage checks, avoid minute truncation for full-bin exports, and trim to the overlapping EEG/video window when video starts after EEG.
- Left TTL-pulse and multi-bin Viewpoint alignment as explicit follow-up work rather than extending the single-bin rule prematurely.
- Verification:
  - `Get-Content -Path next_steps.md`

## 2026-05-26

### Extracted Box4/M65 alignment comparison via test driver

- Updated `test_preprocess_func.m` to prepend a user-style comparison run for `C:\Users\yzhao\matlab_projects\sleep_data_extraction\20260505_box4_M65_LongBaseline\20260505_box4_M65_LongBaseline_2026-05-05_09-01-47-463.exp`.
- Commented out the previously active SOD1 example block so the test driver only runs the Box4/M65 end/start alignment comparison.
- Saved an end-aligned export to `C:\Users\yzhao\matlab_projects\sleep_data_extraction\20260505_box4_M65_LongBaseline_end_alignment.mat` and a legacy start-aligned export to `C:\Users\yzhao\matlab_projects\sleep_data_extraction\20260505_box4_M65_LongBaseline_start_alignment.mat`.
- Verified both outputs have matching EEG/EMG lengths and video path metadata; the end-aligned output has `video_start_time = 31.889` seconds, while the start-aligned output has `video_start_time = 0`.
- Verification:
  - `matlab -batch "run('C:\Users\yzhao\matlab_projects\preprocess_sleep_data\test_preprocess_func.m')"`
  - `matlab -batch "files={'C:\Users\yzhao\matlab_projects\sleep_data_extraction\20260505_box4_M65_LongBaseline_end_alignment.mat','C:\Users\yzhao\matlab_projects\sleep_data_extraction\20260505_box4_M65_LongBaseline_start_alignment.mat'}; for i=1:numel(files), s=load(files{i}); [~,name]=fileparts(files{i}); fprintf('%s\n',name); fprintf('  video_start_time: %.6f\n', s.video_start_time); fprintf('  eeg_len: %d, emg_len: %d, eeg_frequency: %.6f\n', numel(s.eeg), numel(s.emg), s.eeg_frequency); fprintf('  duration_sec: %.6f\n', numel(s.eeg)/s.eeg_frequency); fprintf('  video_name: %s\n', string(s.video_name)); fprintf('  video_path: %s\n', string(s.video_path)); end"`

### Extracted Viewpoint alignment comparison files

- Generated two `.mat` exports from `C:\Users\yzhao\matlab_projects\sleep_data_extraction\20260505_box3_M64_LongBaseline\20260505_box3_M64_LongBaseline_2026-05-05_09-01-47-463.exp`.
- Saved an end-aligned export to `C:\Users\yzhao\matlab_projects\sleep_data_extraction\20260505_box3_M64_LongBaseline_end_alignment.mat` and a legacy start-aligned export to `C:\Users\yzhao\matlab_projects\sleep_data_extraction\20260505_box3_M64_LongBaseline_start_alignment.mat`.
- Verified both outputs have matching EEG/EMG lengths and video path metadata; the end-aligned output has `video_start_time = 35.088` seconds, while the start-aligned output has `video_start_time = 0`.
- Verification:
  - `matlab -batch "addpath('C:\Users\yzhao\matlab_projects\preprocess_sleep_data'); addpath('C:\Users\yzhao\matlab_projects\preprocess_sleep_data\EEGtoolbox'); expPath='C:\Users\yzhao\matlab_projects\sleep_data_extraction\20260505_box3_M64_LongBaseline\20260505_box3_M64_LongBaseline_2026-05-05_09-01-47-463.exp'; outEnd='C:\Users\yzhao\matlab_projects\sleep_data_extraction\20260505_box3_M64_LongBaseline_end_alignment.mat'; outStart='C:\Users\yzhao\matlab_projects\sleep_data_extraction\20260505_box3_M64_LongBaseline_start_alignment.mat'; preprocess_sleep_data(expPath,'save_path',outEnd,'video_alignment','end'); preprocess_sleep_data(expPath,'save_path',outStart,'video_alignment','start'); fprintf('Wrote %s\nWrote %s\n', outEnd, outStart);"`
  - `matlab -batch "files={'C:\Users\yzhao\matlab_projects\sleep_data_extraction\20260505_box3_M64_LongBaseline_end_alignment.mat','C:\Users\yzhao\matlab_projects\sleep_data_extraction\20260505_box3_M64_LongBaseline_start_alignment.mat'}; for i=1:numel(files), s=load(files{i}); [~,name]=fileparts(files{i}); fprintf('%s\n',name); fprintf('  video_start_time: %.6f\n', s.video_start_time); fprintf('  eeg_len: %d, emg_len: %d, eeg_frequency: %.6f\n', numel(s.eeg), numel(s.emg), s.eeg_frequency); fprintf('  duration_sec: %.6f\n', numel(s.eeg)/s.eeg_frequency); fprintf('  video_name: %s\n', string(s.video_name)); fprintf('  video_path: %s\n', string(s.video_path)); end"`

### Updated agent-facing project docs

- Added project-specific reminders to `AGENTS.md` for the active preprocessing entry point, scoped video-alignment behavior, and `README.md` encoding.
- Replaced the placeholder `project_overview.md` with a concise map of the active MATLAB runtime path, relevant helper layers, active vs. secondary files, video-alignment mental model, data expectations, and open questions.
- Verification:
  - `Get-Content -Path AGENTS.md`
  - `Get-Content -Path project_overview.md`

### Defaulted scoped Viewpoint video end alignment

- Added a `video_alignment` name-value parameter to `preprocess_sleep_data.m`, defaulting to `end`.
- Scoped the new end-alignment calculation to single-bin, single-video Viewpoint exports with no TDT photometry input; other acquisition layouts keep their current `video_start_time` behavior for now.
- Updated `README.md` to document the new parameter and clarify that `video_start_time` is the source-video seek offset at exported EEG/EMG time zero.
- Replaced the placeholder `next_steps.md` content with concrete follow-ups for multi-bin Viewpoint, Viewpoint plus separate TDT photometry, TDT-only no-video outputs, and Sleep Scoring App offset handling.
- Verification:
  - `matlab -batch "addpath(pwd); m=checkcode('preprocess_sleep_data.m','-id'); if isempty(m), disp('checkcode clean'); else, for k=1:numel(m), fprintf('%d:%s:%s\n',m(k).line,m(k).id,m(k).message); end; end"`
  - MATLAB `checkcode` completed; it reported only existing unused-variable style warnings in `preprocess_sleep_data.m`.

## 2026-05-25

### Onboarded to preprocessing pipeline

- Read `AGENTS.md`, `README.md`, and `preprocess_sleep_data.m` to identify the active entry point and repo workflow expectations.
- Confirmed `preprocess_sleep_data.m` is the main MATLAB function for exporting Sleep Scoring App-ready `.mat` files from Viewpoint `.exp` EEG/EMG data, TDT EEG/EMG data, and optional TDT fiber photometry data.
- Key flow understood: parse name-value inputs, load EEG/EMG from Viewpoint or TDT, optionally load and normalize 465/405 photometry channels, optionally import sleep scores, align separate EEG/EMG and photometry recordings with TTL pulses, interpolate missing samples, plot optional diagnostics, then save one output file or segmented 12-hour/bin-specific files.
- Noted open project hygiene question: `AGENTS.md` still contains placeholder runtime, test, formatter, and project-specific reminder values, so future verification recipes are not yet authoritative.
- Verification:
  - `Get-Content -Path AGENTS.md`
  - `Get-Content -Path README.md`
  - `Get-Content -Path preprocess_sleep_data.m`
  - `rg -n '^(#|##|###)|preprocess_sleep_data|Output|Note|More Examples' README.md`
  - `rg -n '^function|^%%|TDTbin2mat|loadEXP|ExtractContinuousData|readmatrix|polyfit|downsample|sleep_scores|save\(' preprocess_sleep_data.m`
  - `git config --global --add safe.directory C:/Users/yzhao/matlab_projects/preprocess_sleep_data`
  - `git status --short --branch`

### Documented MATLAB runtime

- Updated `AGENTS.md` to record MATLAB R2025a as the runtime for this project.
- Left unresolved template items such as formatter, tests, and project-specific reminders for future sessions when those choices become concrete.
- Verification:
  - Not run; documentation-only update.

Rotation policy: the live log holds at most the **5 most recent unique calendar dates**. When a new date would push the file past 5 unique dates, move the oldest 5 dates as a chunk into a new file at `work_log_archive/work_log_<earliest>_to_<latest>.md`. The live file always holds at most 5 unique dates; each archive file always holds exactly 5.

If today's date already has a `## YYYY-MM-DD` header at the top, add a new `###` session subsection under it rather than starting a second `## YYYY-MM-DD` header for the same date.

<!--
Each session entry follows this shape:

## YYYY-MM-DD

### Short title for what was done

- bullet describing what was added or changed
- another bullet — keep them high-level and user/agent-facing, not implementation play-by-play
- if relevant, intended profiling signal or measurement:
  - what to look for in logs / output
  - what numbers were observed
- Verification:
  - the exact command(s) that were actually run
  - what passed / what was confirmed

Newest entry goes on top. If the session did multiple distinct pieces of work, use multiple `###` subsections under one `##` date header.
-->
