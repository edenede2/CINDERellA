
% CINDERellA Analysis Script v2
disp('--- CINDERellA Analysis Started ---');

% Add paths
addpath('C:/Users/PsyLab-6028/Documents/CINDERellA');
addpath(genpath('C:/Users/PsyLab-6028/Documents/CINDERellA'));

% Parameters for CINDERellA
max_parents = 3;
runtime_minutes = 30;  % Reduced for initial testing
num_samples = 200;     % Reduced for initial testing
edge_threshold = 0.3;

% Read data matrix directly (no headers, just numeric data)
disp('Reading expression data...');
data_matrix = readmatrix('C:/Users/PsyLab-6028/Documents/CINDERellA/notebooks/CIND_ME_plus_phenos/expression_data_for_cinderella.txt');

disp(['Data matrix size: ', num2str(size(data_matrix))]);
disp(['Samples (rows): ', num2str(size(data_matrix, 1))]);
disp(['Variables (cols): ', num2str(size(data_matrix, 2))]);

% Check for NaN or Inf values
nan_count = sum(isnan(data_matrix(:)));
inf_count = sum(isinf(data_matrix(:)));
disp(['NaN values: ', num2str(nan_count)]);
disp(['Inf values: ', num2str(inf_count)]);

if nan_count > 0 || inf_count > 0
    disp('Warning: Data contains NaN or Inf values. Replacing with median...');
    for j = 1:size(data_matrix, 2)
        col = data_matrix(:, j);
        med_val = median(col(~isnan(col) & ~isinf(col)));
        col(isnan(col) | isinf(col)) = med_val;
        data_matrix(:, j) = col;
    end
    disp('Data cleaning completed.');
end

% Run CINDERellA
disp('Starting CINDERellA analysis...');
try
    CINDERellA(data_matrix, ...
        'output_dir', 'C:/Users/PsyLab-6028/Documents/CINDERellA/notebooks/CIND_ME_plus_phenos', ...
        'sampler', 'M.c2PB', ...
        'max_parents', max_parents, ...
        'runtime_minutes', runtime_minutes, ...
        'num_samples', num_samples, ...
        'edge_threshold', edge_threshold);
    
    disp('--- CINDERellA Analysis Completed Successfully ---');
catch ME
    disp('--- CINDERellA Analysis Failed ---');
    disp(['Error: ', ME.message]);
    rethrow(ME);
end
