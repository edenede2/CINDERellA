
% CINDERellA High-Level Analysis with Module→Phenotype Priors
% FIXED: Explicitly specify FileType for MATLAB R2023a compatibility

clear all; close all; clc;

fprintf('\n=== CINDERELLA WITH DIRECTIONAL PRIORS (FIXED) ===\n');
fprintf('Using high-level CINDERellA() function\n');

% Add CINDERellA to path
addpath('C:/Users/PsyLab-6028/Documents/CINDERellA');
addpath('C:/Users/PsyLab-6028/Documents/CINDERellA/functions');

% Set paths
input_file = 'C:/Users/PsyLab-6028/Documents/CINDERellA/notebooks/cinderella_input_ME_mods_plus_phenotypes.tsv';
output_dir = 'C:/Users/PsyLab-6028/Documents/CINDERellA/notebooks/CIND_ME_plus_phenos/cinderella_with_priors';

% Create output directory
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

try
    fprintf('\n=== STEP 1: LOAD DATA ===\n');
    fprintf('Loading: %s\n', input_file);
    
    % FIXED: Explicitly specify FileType='text' for .tsv files
    data = readtable(input_file, 'FileType', 'text', 'Delimiter', '\t');
    fprintf('✓ Data loaded: %d variables x %d samples\n', width(data)-1, height(data));
    
    % Extract data matrix and variable names
    data_matrix = table2array(data(:, 2:end))';  % Transpose: genes=rows, samples=cols
    var_names = data.Properties.VariableNames(2:end);
    
    [nGenes, nSamples] = size(data_matrix);
    fprintf('✓ Matrix size: %d genes (variables) x %d samples\n', nGenes, nSamples);
    
    % Display variables
    fprintf('\nVariables:\n');
    for i = 1:length(var_names)
        fprintf('  %d. %s\n', i, var_names{i});
    end
    
    %% === STEP 2: CREATE PRIOR MATRIX ===
    fprintf('\n=== STEP 2: CREATE DIRECTIONAL PRIOR MATRIX ===\n');
    
    % Identify modules and phenotypes
    module_indices = [];
    phenotype_indices = [];
    
    for i = 1:length(var_names)
        if contains(var_names{i}, 'ME_')
            module_indices = [module_indices, i];
        elseif contains(var_names{i}, 'PHENO_')
            phenotype_indices = [phenotype_indices, i];
        end
    end
    
    fprintf('\nIdentified:\n');
    fprintf('  Modules (%d): ', length(module_indices));
    for idx = module_indices
        fprintf('%s ', var_names{idx});
    end
    fprintf('\n');
    
    fprintf('  Phenotypes (%d): ', length(phenotype_indices));
    for idx = phenotype_indices
        fprintf('%s ', var_names{idx});
    end
    fprintf('\n');
    
    % Initialize prior matrix
    % In CINDERellA: 1 = edge ALLOWED, 0 = edge DISALLOWED
    prior_matrix = ones(nGenes, nGenes);  % Start with all edges allowed
    prior_matrix = setdiag(prior_matrix, 0);  % No self-loops
    
    fprintf('\n=== APPLYING DIRECTIONAL CONSTRAINTS ===\n');
    
    % CONSTRAINT 1: DISALLOW Phenotype → Module edges
    for i = phenotype_indices
        for j = module_indices
            prior_matrix(i, j) = 0;  % Phenotype CANNOT cause Module
        end
    end
    fprintf('✓ Blocked all Phenotype → Module edges\n');
    
    % CONSTRAINT 2: ALLOW Module → Phenotype edges (already 1, but make explicit)
    for i = module_indices
        for j = phenotype_indices
            prior_matrix(i, j) = 1;  % Module CAN cause Phenotype
        end
    end
    fprintf('✓ Allowed all Module → Phenotype edges\n');
    
    % CONSTRAINT 3: Allow Module → Module edges (for regulatory networks)
    fprintf('✓ Module → Module edges allowed (default)\n');
    
    % CONSTRAINT 4: Block Phenotype → Phenotype edges (they are both outcomes)
    for i = phenotype_indices
        for j = phenotype_indices
            if i ~= j
                prior_matrix(i, j) = 0;
            end
        end
    end
    fprintf('✓ Blocked Phenotype → Phenotype edges\n');
    
    % Count constraints
    n_disallowed = sum(prior_matrix(:) == 0) - nGenes;  % Subtract diagonal
    n_allowed = sum(prior_matrix(:) == 1);
    fprintf('\nPrior matrix summary:\n');
    fprintf('  Allowed edges: %d\n', n_allowed);
    fprintf('  Disallowed edges: %d\n', n_disallowed);
    fprintf('  Total possible edges: %d\n', nGenes*(nGenes-1));
    
    %% === STEP 3: RUN CINDERELLA ===
    fprintf('\n=== STEP 3: RUN CINDERELLA WITH PRIORS ===\n');
    
    % CINDERellA parameters
    runtime_minutes = 2;      % 2 minutes for reasonable convergence
    max_parents = 3;           % Max 3 parents per node
    num_samples = 200;         % 200 network samples
    edge_threshold = 0.3;      % Show edges with >30%% frequency
    
    fprintf('\nCINDERellA Parameters:\n');
    fprintf('  Runtime: %d minutes\n', runtime_minutes);
    fprintf('  Max parents: %d\n', max_parents);
    fprintf('  Network samples: %d\n', num_samples);
    fprintf('  Edge threshold: %.2f\n', edge_threshold);
    fprintf('  Sampler: M.REV50 (multi-chain with reversible jump)\n');
    
    fprintf('\n🚀 Starting CINDERellA analysis...\n');
    fprintf('This will take approximately %d minutes...\n', runtime_minutes);
    fprintf('Progress will be shown during execution.\n\n');
    
    % Run CINDERellA with priors
    CINDERellA(data_matrix, ...
               'output_dir', output_dir, ...
               'sampler', 'M.REV50', ...
               'max_parents', max_parents, ...
               'runtime_minutes', runtime_minutes, ...
               'num_samples', num_samples, ...
               'edge_threshold', edge_threshold, ...
               'prior_matrix', prior_matrix, ...
               'layout', 'layered');
    
    fprintf('\n✅ CINDERellA analysis completed!\n');
    
    %% === STEP 4: ANALYZE RESULTS ===
    fprintf('\n=== STEP 4: ANALYZE DIRECTIONAL RESULTS ===\n');
    
    % Load edge frequencies
    edgefrq_file = fullfile(output_dir, 'edgefrq.txt');
    if exist(edgefrq_file, 'file')
        edgefrq_data = dlmread(edgefrq_file);
        edgefrq_matrix = sparse(edgefrq_data(:,1), edgefrq_data(:,2), edgefrq_data(:,3), nGenes, nGenes);
        
        % Convert to full matrix for analysis
        edgefrq_full = full(edgefrq_matrix);
        
        % Analyze directionality
        module_to_phenotype = 0;
        phenotype_to_module = 0;
        module_to_module = 0;
        phenotype_to_phenotype = 0;
        
        total_edges = 0;
        
        fprintf('\nStrong directional edges (>0.3):\n');
        fprintf('----------------------------------------\n');
        
        for i = 1:nGenes
            for j = 1:nGenes
                if edgefrq_full(i, j) > 0.3  % Strong edges only
                    total_edges = total_edges + 1;
                    
                    if ismember(i, module_indices) && ismember(j, phenotype_indices)
                        module_to_phenotype = module_to_phenotype + 1;
                        fprintf('  ✓ Module → Phenotype: %s → %s (%.3f)\n', ...
                            var_names{i}, var_names{j}, edgefrq_full(i, j));
                    elseif ismember(i, phenotype_indices) && ismember(j, module_indices)
                        phenotype_to_module = phenotype_to_module + 1;
                        fprintf('  ⚠️  Phenotype → Module: %s → %s (%.3f) [SHOULD BE BLOCKED]\n', ...
                            var_names{i}, var_names{j}, edgefrq_full(i, j));
                    elseif ismember(i, module_indices) && ismember(j, module_indices)
                        module_to_module = module_to_module + 1;
                        fprintf('  • Module → Module: %s → %s (%.3f)\n', ...
                            var_names{i}, var_names{j}, edgefrq_full(i, j));
                    elseif ismember(i, phenotype_indices) && ismember(j, phenotype_indices)
                        phenotype_to_phenotype = phenotype_to_phenotype + 1;
                    end
                end
            end
        end
        
        fprintf('\n=== DIRECTIONALITY SUMMARY (edges > 0.3) ===\n');
        fprintf('Total strong edges: %d\n', total_edges);
        fprintf('  Module → Phenotype: %d (%.1f%%)\n', module_to_phenotype, 100*module_to_phenotype/max(total_edges,1));
        fprintf('  Phenotype → Module: %d (%.1f%%)\n', phenotype_to_module, 100*phenotype_to_module/max(total_edges,1));
        fprintf('  Module → Module: %d (%.1f%%)\n', module_to_module, 100*module_to_module/max(total_edges,1));
        fprintf('  Phenotype → Phenotype: %d (%.1f%%)\n', phenotype_to_phenotype, 100*phenotype_to_phenotype/max(total_edges,1));
        
        % Validation check
        if phenotype_to_module == 0
            fprintf('\n✅ SUCCESS: No Phenotype→Module edges detected!\n');
            fprintf('   Priors successfully enforced directional causality.\n');
        else
            fprintf('\n⚠️  WARNING: %d Phenotype→Module edges found despite priors.\n', phenotype_to_module);
        end
        
        % Save detailed results with variable names
        results_file = fullfile(output_dir, 'directed_edges_with_names.txt');
        fid = fopen(results_file, 'w');
        fprintf(fid, 'from_idx\tto_idx\tfrequency\tfrom_name\tto_name\tedge_type\n');
        
        for i = 1:nGenes
            for j = 1:nGenes
                if edgefrq_full(i, j) > 0.1
                    % Determine edge type
                    if ismember(i, module_indices) && ismember(j, phenotype_indices)
                        edge_type = 'Module→Phenotype';
                    elseif ismember(i, phenotype_indices) && ismember(j, module_indices)
                        edge_type = 'Phenotype→Module';
                    elseif ismember(i, module_indices) && ismember(j, module_indices)
                        edge_type = 'Module→Module';
                    elseif ismember(i, phenotype_indices) && ismember(j, phenotype_indices)
                        edge_type = 'Phenotype→Phenotype';
                    else
                        edge_type = 'Other';
                    end
                    
                    fprintf(fid, '%d\t%d\t%.6f\t%s\t%s\t%s\n', ...
                        i, j, edgefrq_full(i, j), var_names{i}, var_names{j}, edge_type);
                end
            end
        end
        fclose(fid);
        fprintf('\nDetailed results saved to: %s\n', results_file);
        
        fprintf('\n🎉 Analysis complete! Check output directory: %s\n', output_dir);
        fprintf('Key files:\n');
        fprintf('  • edgefrq.txt - Edge frequencies\n');
        fprintf('  • directed_edges_with_names.txt - Human-readable results\n');
        fprintf('  • network_visualization.png - Network plot\n');
        fprintf('  • mcmc_diagnostics.png - Convergence diagnostics\n');
        
    else
        fprintf('⚠️  Edge frequency file not found: %s\n', edgefrq_file);
    end
    
catch ME
    fprintf('\n❌ Error: %s\n', ME.message);
    fprintf('Location: %s (line %d)\n', ME.stack(1).file, ME.stack(1).line);
    rethrow(ME);
end

fprintf('\n=== CINDERELLA ANALYSIS COMPLETE ===\n');
