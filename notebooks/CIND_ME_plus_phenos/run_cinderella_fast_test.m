
% CINDERellA FAST TEST - 30 second runtime
% This tests if the setup works before running longer analysis

clear all; close all; clc;

fprintf('\n=== CINDERELLA FAST TEST (30 SECONDS) ===\n');

% Add CINDERellA to path
addpath('C:/Users/PsyLab-6028/Documents/CINDERellA');
addpath('C:/Users/PsyLab-6028/Documents/CINDERellA/functions');

% Set paths
input_file = 'C:/Users/PsyLab-6028/Documents/CINDERellA/notebooks/cinderella_input_ME_mods_plus_phenotypes.tsv';
output_dir = 'C:/Users/PsyLab-6028/Documents/CINDERellA/notebooks/CIND_ME_plus_phenos/cinderella_fast_test';

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

try
    fprintf('\n=== LOADING DATA ===\n');
    data = readtable(input_file, 'FileType', 'text', 'Delimiter', '\t');
    fprintf('✓ Loaded: %d variables x %d samples\n', width(data)-1, height(data));
    
    data_matrix = table2array(data(:, 2:end))';
    var_names = data.Properties.VariableNames(2:end);
    [nGenes, nSamples] = size(data_matrix);
    
    fprintf('✓ Matrix: %d genes x %d samples\n', nGenes, nSamples);
    
    %% CREATE PRIORS
    fprintf('\n=== CREATING PRIORS ===\n');
    
    module_indices = [];
    phenotype_indices = [];
    
    for i = 1:length(var_names)
        if contains(var_names{i}, 'ME_')
            module_indices = [module_indices, i];
        elseif contains(var_names{i}, 'PHENO_')
            phenotype_indices = [phenotype_indices, i];
        end
    end
    
    fprintf('Modules: %d, Phenotypes: %d\n', length(module_indices), length(phenotype_indices));
    
    % Prior matrix
    prior_matrix = ones(nGenes, nGenes);
    prior_matrix = setdiag(prior_matrix, 0);
    
    % Block Phenotype → Module
    for i = phenotype_indices
        for j = module_indices
            prior_matrix(i, j) = 0;
        end
    end
    
    % Block Phenotype → Phenotype
    for i = phenotype_indices
        for j = phenotype_indices
            if i ~= j
                prior_matrix(i, j) = 0;
            end
        end
    end
    
    fprintf('✓ Prior matrix created\n');
    fprintf('  Allowed edges: %d\n', sum(prior_matrix(:)));
    fprintf('  Blocked edges: %d\n', sum(prior_matrix(:)==0) - nGenes);
    
    %% RUN FAST TEST
    fprintf('\n=== RUNNING CINDERELLA (30 SECONDS) ===\n');
    
    % MINIMAL PARAMETERS FOR TESTING
    runtime_minutes = 0.2;    % 30 seconds only!
    max_parents = 2;          % Fewer parents
    num_samples = 10;         % Fewer samples
    edge_threshold = 0.2;     % Lower threshold
    
    fprintf('Test parameters:\n');
    fprintf('  Runtime: %.1f minutes (30 seconds)\n', runtime_minutes);
    fprintf('  Max parents: %d\n', max_parents);
    fprintf('  Samples: %d\n', num_samples);
    fprintf('  Threshold: %.2f\n', edge_threshold);
    
    fprintf('\n🚀 Starting test...\n');
    
    tic;
    CINDERellA(data_matrix, ...
               'output_dir', output_dir, ...
               'sampler', 'M.c2PB', ...
               'max_parents', max_parents, ...
               'runtime_minutes', runtime_minutes, ...
               'num_samples', num_samples, ...
               'edge_threshold', edge_threshold, ...
               'prior_matrix', prior_matrix, ...
               'layout', 'layered');
    elapsed = toc;
    
    fprintf('\n✅ Test completed in %.1f seconds!\n', elapsed);
    
    %% CHECK RESULTS
    fprintf('\n=== CHECKING RESULTS ===\n');
    
    edgefrq_file = fullfile(output_dir, 'edgefrq.txt');
    if exist(edgefrq_file, 'file')
        edgefrq_data = dlmread(edgefrq_file);
        fprintf('✓ Edge frequency file found: %d edges\n', size(edgefrq_data, 1));
        
        if ~isempty(edgefrq_data)
            edgefrq_matrix = sparse(edgefrq_data(:,1), edgefrq_data(:,2), edgefrq_data(:,3), nGenes, nGenes);
            edgefrq_full = full(edgefrq_matrix);
            
            % Quick analysis
            n_mod_to_pheno = 0;
            n_pheno_to_mod = 0;
            
            for i = 1:nGenes
                for j = 1:nGenes
                    if edgefrq_full(i, j) > 0.1
                        if ismember(i, module_indices) && ismember(j, phenotype_indices)
                            n_mod_to_pheno = n_mod_to_pheno + 1;
                        elseif ismember(i, phenotype_indices) && ismember(j, module_indices)
                            n_pheno_to_mod = n_pheno_to_mod + 1;
                        end
                    end
                end
            end
            
            fprintf('\nEdge counts (>0.1 frequency):\n');
            fprintf('  Module → Phenotype: %d\n', n_mod_to_pheno);
            fprintf('  Phenotype → Module: %d\n', n_pheno_to_mod);
            
            if n_pheno_to_mod == 0
                fprintf('\n✅ SUCCESS! Priors working - no Pheno→Module edges!\n');
            else
                fprintf('\n⚠️  Found %d Pheno→Module edges despite priors\n', n_pheno_to_mod);
            end
        else
            fprintf('⚠️  Edge frequency file is empty\n');
        end
    else
        fprintf('❌ Edge frequency file not found\n');
    end
    
    fprintf('\n=== FAST TEST COMPLETE ===\n');
    fprintf('If this worked, you can run longer analysis!\n');
    
catch ME
    fprintf('\n❌ Error: %s\n', ME.message);
    if ~isempty(ME.stack)
        fprintf('Location: %s (line %d)\n', ME.stack(1).file, ME.stack(1).line);
    end
    rethrow(ME);
end
