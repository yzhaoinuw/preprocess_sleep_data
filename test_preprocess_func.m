addpath 'C:\Users\yzhao\matlab_projects\preprocess_sleep_data\'
addpath 'C:\Users\yzhao\matlab_projects\preprocess_sleep_data\EEGtoolbox'


%%
% Batch-export all current single-bin Viewpoint box3/box4 recordings.
base_dir = 'C:\Users\yzhao\matlab_projects\sleep_data_extraction';
folders = [dir(fullfile(base_dir, '*box3*')); dir(fullfile(base_dir, '*box4*'))];
folders = folders([folders.isdir]);
[~, sort_idx] = sort({folders.name});
folders = folders(sort_idx);

for i = 1:numel(folders)
    exp_files = dir(fullfile(folders(i).folder, folders(i).name, '*.exp'));
    if isempty(exp_files)
        warning('No .exp file found in %s; skipping.', folders(i).name);
        continue
    elseif numel(exp_files) > 1
        warning('Multiple .exp files found in %s; using %s.', folders(i).name, exp_files(1).name);
    end

    exp_path = fullfile(exp_files(1).folder, exp_files(1).name);
    save_path = fullfile(base_dir, [folders(i).name '.mat']);
    fprintf('Exporting %s -> %s\n', exp_path, save_path);
    preprocess_sleep_data(exp_path, ...
        'save_path', save_path, ...
        'show_figure', false ...
    );
end
