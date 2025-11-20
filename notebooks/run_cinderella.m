
% CINDERellA Analysis Script
disp('--- CINDERellA Analysis Started ---');

% Add paths
addpath('C:/Users/PsyLab-6028/Documents/CINDERellA');
addpath(genpath('C:/Users/PsyLab-6028/Documents/CINDERellA'));

% Parameters for CINDERellA
max_parents = 3;
runtime_minutes = 60;
num_samples = 500;
edge_threshold = 0.3;

% Read expression data
disp('Reading expression data...');
expdata = read_exp('C:/Users/PsyLab-6028/Documents/CINDERellA/notebooks/cinderella_input_ME_mods_plus_phenotypes.tsv');
data_matrix = expdata.data;

disp(['Data matrix size: ', num2str(size(data_matrix))]);

% Run CINDERellA
disp('Starting CINDERellA analysis...');
CINDERellA(data_matrix, ...
    'output_dir', 'C:/Users/PsyLab-6028/Documents/CINDERellA/notebooks/CIND_ME_plus_phenos', ...
    'sampler', 'M.c2PB', ...
    'max_parents', max_parents, ...
    'runtime_minutes', runtime_minutes, ...
    'num_samples', num_samples, ...
    'edge_threshold', edge_threshold);

disp('--- CINDERellA Analysis Completed ---');
