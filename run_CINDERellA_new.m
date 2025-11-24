%% ============================
%  CINDERellA on selected modules + phenotypes
%  Input: cinderella_input_ME_mods_plus_phenotypes.tsv
%  Rows  = modules (ME_M*) + phenotypes (PHENO_*)
%  Cols  = donors
% =============================

%% --- User paths ---
% Path to the CINDERellA repo (folder that contains "functions" dir)
CINDERellA_ROOT = 'C:\Users\PsyLab-6028\Documents\CINDERellA';  % <-- CHANGE THIS
addpath(fullfile(CINDERellA_ROOT, 'functions'));

% Check available memory before starting (optional but recommended for large networks)
if exist('check_memory.m', 'file')
    check_memory();
end

% Path to your TSV file (put it in the same folder as this script or give full path)
inputFile = 'C:\Users\PsyLab-6028\Documents\CINDERellA\notebooks\cinderella_input_ME_mods_plus_phenotypes.tsv';

% Output directory for CINDERellA results
outputDir = 'CINDERellA_results_selectedMods03';

%% --- Load TSV ---
fprintf('Loading TSV: %s\n', inputFile);
T = readtable(inputFile, 'FileType', 'text', 'Delimiter', '\t');

% First column = feature names (modules / phenotypes)
featureNames = T.Gene;                       % column name in your file
sampleNames  = T.Properties.VariableNames(2:end);  % all donor IDs
dataMatrix   = T{:, 2:end};                  % numeric matrix (features x donors)

fprintf('Loaded %d features (rows) x %d donors (columns).\n', ...
    size(dataMatrix,1), size(dataMatrix,2));

%% --- Define selection: modules + phenotypes ---

% These are the module numbers you gave (WITHOUT the "M" prefix)
selectedModules = [2,4,5,7,15,18,19,26,31,32,33,36,37,40,43,44,47,49,54,55,57,58,59,62,63,65,66,67,69,70,72,76,79,81,82,84,85,87,88,89,90,91,92,93,98,101,102,103,104,105,107,109,110,112,113,114,116,118,119,124,127,128,130,135,136,137,138,140,144,145,147,148,150,152,153,154,155,157,160,161,162,163,164,167,168,170,171,172,173,175,177,181,182,185,186,187,188,190,191,192,196,197,199,200,203,207,208,210,212,214,215,216,217,220,221,223,226,227,230,232,233,234,235,236,241,244,245,247,251,252,253,257,261,262,267,268,269,271,272,273,277,278,279,280,283,284,288,289,290,292,295,296,297,299,306,307,309,312,314,317,319,322,323,328,330,335,336,337,338,341,343,346,348,350,353,356,358,359,360,362,373,375,376,377,378,379,380,381,384,387,388,395,401,403,404,405,408,409,411,412,413,416,417,418,421,424,425,427,428,430,431,437,439,440,443,446,450,452,453,456,457,459,461,463,470,473,476,477,488,490,495,496,498,499,500,502];

% Your note: "add M to the start of each module number"
% In THIS TSV, module rows are named like "ME_M2", "ME_M3", ...,
% so we create target names "ME_M<number>".
selectedModuleNames = arrayfun(@(m) sprintf('ME_M%d', m), ...
                               selectedModules, 'UniformOutput', false);

% Phenotypes to keep
selectedPhenotypes = {'PHENO_pseudotime', 'PHENO_cogn_global'};

% Logical mask of rows to keep
keepMask = ismember(featureNames, selectedModuleNames) | ...
           ismember(featureNames, selectedPhenotypes);

fprintf('Keeping %d / %d rows (modules + phenotypes).\n', ...
    sum(keepMask), numel(featureNames));

dataSel    = dataMatrix(keepMask, :);    % [nNodes x nDonors]
nodeNames  = featureNames(keepMask);     % cell array of names

% Sanity check: how many modules vs phenotypes?
isPheno    = ismember(nodeNames, selectedPhenotypes);
nPheno     = sum(isPheno);
nModules   = numel(nodeNames) - nPheno;
fprintf('  -> %d modules, %d phenotypes.\n', nModules, nPheno);

%% --- (Optional) Create "pretty" labels for modules (M<number>) ---
prettyNames = nodeNames;  % default
for i = 1:numel(prettyNames)
    if startsWith(prettyNames{i}, 'ME_M')
        % Extract numeric part and rename as "M<number>"
        numVal = sscanf(prettyNames{i}, 'ME_M%d');
        prettyNames{i} = sprintf('M%d', numVal);
    end
end

%% --- Save mapping so you can interpret nodes later ---
nodeTable = table((1:numel(nodeNames))', nodeNames, prettyNames, isPheno, ...
    'VariableNames', {'node_id', 'raw_name', 'pretty_name', 'is_phenotype'});
writetable(nodeTable, 'cinderella_selected_nodes.tsv', ...
           'FileType', 'text', 'Delimiter', '\t');

%% --- Create prior matrix: modules -> phenotypes only ---
% Prior matrix controls allowed edges:
%   priorMatrix(i,j) = 1  means edge i -> j is ALLOWED
%   priorMatrix(i,j) = 0  means edge i -> j is FORBIDDEN
%
% We want: modules can point to phenotypes, but not vice versa
% - No phenotype -> module edges
% - No phenotype -> phenotype edges  
% - Modules -> modules allowed (set to 0 to also forbid these)
% - Modules -> phenotypes allowed

nNodes = numel(nodeNames);
priorMatrix = ones(nNodes, nNodes);  % Start with all edges allowed

% Find indices of modules and phenotypes
moduleIdx = find(~isPheno);   % modules are NOT phenotypes
phenoIdx  = find(isPheno);    % phenotypes

% Forbid edges FROM phenotypes TO anything (phenotypes cannot be parents)
priorMatrix(phenoIdx, :) = 0;

% Forbid self-loops (optional, but recommended)
priorMatrix(1:nNodes+1:end) = 0;

fprintf('\nPrior matrix created:\n');
fprintf('  %d modules can be parents of any node\n', numel(moduleIdx));
fprintf('  %d phenotypes CANNOT be parents (no outgoing edges)\n', numel(phenoIdx));
fprintf('  Allowed edges: modules -> modules, modules -> phenotypes\n');
fprintf('  Forbidden edges: phenotypes -> anything, self-loops\n');
fprintf('  Search space reduced by prior: %.1f%%\n\n', 100 * sum(priorMatrix(:)==0) / numel(priorMatrix));

% Memory optimization: clear large temporary variables
clear dataMatrix T keepMask selectedModuleNames selectedModules;

%% --- Run CINDERellA (high-resolution settings) ---
% CINDERellA expects: rows = "genes/nodes", columns = samples.
% We use multi-chain sampler with optimized settings for large networks

% MEMORY-OPTIMIZED SETTINGS FOR 238 NODES:
% =========================================
% For 238 nodes, memory usage is very high. If you get out-of-memory errors:
%   1. Reduce num_samples (try 20, 30, or 50 instead of 100)
%   2. Reduce runtime_minutes (0.25 to 0.5 minutes for testing)
%   3. Use 'M.REV50' sampler (more memory efficient than 'M.c2PB')
%   4. Reduce max_parents from 3 to 2
%   5. Close other programs to free RAM
%   6. Consider running on a machine with more memory (16GB+ recommended)

runtime_minutes = 0.5;   % Short runtime for memory efficiency (increase for better results)
num_samples     = 50;    % Reduced samples for memory (increase to 50-100 when stable)

fprintf('\nRunning CINDERellA on %d nodes x %d donors...\n', ...
    size(dataSel,1), size(dataSel,2));
fprintf('  runtime_minutes = %d, num_samples = %d\n', ...
    runtime_minutes, num_samples);

% Make sure output directory exists
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

% Run native CINDERellA
CINDERellA(dataSel, ...
    'output_dir',      outputDir, ...
    'sampler',         'M.REV50', ...      % Memory-efficient sampler (use M.c2PB for better mixing if enough RAM)
    'max_parents',     3, ...              % Reduced from 3 to save memory (increase if no memory issues)
    'runtime_minutes', runtime_minutes, ...
    'num_samples',     num_samples, ...
    'force_recompute',           true, ...
    'prior_matrix',    priorMatrix, ...   % Apply directional constraint
    'edge_threshold',  0.30, ...          % for internal PNG it generates
    'layout',          'layered');        % a DAG-like visualization

fprintf('CINDERellA finished. Results in folder: %s\n', outputDir);

%% --- Build directed graph from edge frequencies ---
% edgefrq.txt has 3 columns: [source_idx, target_idx, frequency] :contentReference[oaicite:4]{index=4}
edgeFile = fullfile(outputDir, 'edgefrq.txt');
if ~exist(edgeFile, 'file')
    error('edgefrq.txt not found at %s. Did CINDERellA finish successfully?', edgeFile);
end

edgefrq_data = dlmread(edgeFile);
nNodes = size(dataSel, 1);

% Sparse matrix of edge frequencies
edgeFrq = sparse(edgefrq_data(:,1), edgefrq_data(:,2), edgefrq_data(:,3), ...
                 nNodes, nNodes);

% Choose threshold for graph (this is your "high resolution" view:
% low threshold = dense, high threshold = sparse).
graphEdgeThresh = 0.50;   % you can tune this (e.g. 0.25, 0.35, 0.5, ...)

[srcIdx, dstIdx, w] = find(edgeFrq >= graphEdgeThresh);

fprintf('Graph threshold %.2f keeps %d directed edges.\n', ...
    graphEdgeThresh, numel(w));

srcNames = prettyNames(srcIdx);
dstNames = prettyNames(dstIdx);

%% --- Create digraph object & plot ---
G = digraph(srcNames, dstNames, full(w));  % Convert sparse to full for digraph

figure('Color', 'w');
p = plot(G, 'Layout', 'layered', 'EdgeLabel', round(w, 2));
title(sprintf('CINDERellA causal network (thr=%.2f)', graphEdgeThresh), ...
      'Interpreter', 'none');

% Highlight phenotypes (e.g. red) vs modules (default)
hold on;
isPhenoNode = ismember(G.Nodes.Name, prettyNames(isPheno));
highlight(p, find(isPhenoNode), 'NodeColor', 'r', 'MarkerSize', 8);
hold off;

%% --- Export edges & nodes for Cytoscape (or other tools) ---

edgeTable = table(srcNames, dstNames, w, ...
    'VariableNames', {'source', 'target', 'edge_weight'});
writetable(edgeTable, 'cinderella_selected_edges.tsv', ...
           'FileType', 'text', 'Delimiter', '\t');

writetable(G.Nodes, 'cinderella_selected_nodes_forCytoscape.tsv', ...
           'FileType', 'text', 'Delimiter', '\t');

fprintf('Exported:\n');
fprintf('  cinderella_selected_edges.tsv\n');
fprintf('  cinderella_selected_nodes_forCytoscape.tsv\n');
fprintf('  cinderella_selected_nodes.tsv (id/name mapping)\n');
