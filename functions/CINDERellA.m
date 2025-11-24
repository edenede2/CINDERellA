function CINDERellA(data_matrix, varargin)
% CINDERellA - Easy-to-use Bayesian Network Learning Tool
% 
% This function learns causal networks from gene expression data using
% Markov Chain Monte Carlo (MCMC) methods.
%
% USAGE EXAMPLES:
%   % Load your data first
%   expdata = read_exp('my_expression_data.txt');
%   
%   % Basic usage with default parameters
%   CINDERellA(expdata.data);
%
%   % With custom parameters
%   CINDERellA(expdata.data, 'output_dir', 'my_results', ...
%                   'max_parents', 3, 'runtime_minutes', 30);
%
%   % With prior knowledge to constrain search space
%   prior = ones(nGenes, nGenes);  % Start with all edges allowed
%   prior(1,2) = 0;  % Disallow edge from gene 1 to gene 2
%   CINDERellA(expdata.data, 'prior_matrix', prior);
%
%   % With layout and visualization options
%   CINDERellA(expdata.data, 'layout', 'force', 'edge_threshold', 0.4);
%
% INPUT PARAMETERS:
%   data_matrix        - Expression data matrix (required)
%                       Format: rows=genes, columns=samples
%   'output_dir'       - Output directory (default: './CINDERellA_results/')
%   'sampler'          - MCMC sampler method (default: 'M.c2PB', or 'M.REV50' if hprior used)
%                       Single chain samplers: 'STR','c2PB','c3PB','c4PB','1PB','2PB','3PB','4PB','REV50'
%                       Multi-chain samplers (add 'M.' prefix): 'M.STR','M.c2PB','M.c3PB','M.c4PB','M.1PB','M.2PB','M.3PB','M.4PB','M.REV50'
%                       Single chain: samples networks every (runtime/num_samples) seconds
%                       Multi-chain: collects final network state from each chain after (runtime/num_samples) seconds
%   'max_parents'      - Maximum parents per gene (default: 3)
%   'runtime_minutes'  - Total runtime in minutes (default: 30)
%   'num_samples'      - Number of network samples to collect (default: 100)
%   'edge_threshold'   - Edge frequency threshold for visualization (default: 0.3)
%   'prior_matrix'     - nGenes x nGenes binary matrix (1=allowed, 0=disallowed edges)
%   'force_recompute'  - Force recomputation even if results exist (default: false)
%   'layout'           - Network layout algorithm (default: 'force')
%                       Options: 'circle', 'force', 'layered', 'subspace'

%
% OUTPUTS:
%   - edgefrq.txt: Edge frequencies from sampled networks
%   - Mcmc.mat: All sampled networks
%   - Param.mat: Parameters used
%   - LS.mat: Local scores
%   - mcmc_diagnostics.png: Log likelihood trace plot
%   - network_visualization.png: Network plot with edges > threshold
%
% DATA PREPARATION:
%   Before using this function, load your data:
%   expdata = read_exp('your_data_file.txt');
%   Then pass expdata.data to this function
%
%   Data matrix should have:
%   - Rows: genes/variables
%   - Columns: samples/conditions
%
% Author: Shinya Tasaki, Ph.D. (stasaki@gmail.com)
% Modified for easier use with diagnostics and visualization
% License: 3-clause BSD License

%% Parse input arguments
p = inputParser;
addParameter(p, 'data', '', @(x) ischar(x) && exist(x, 'file'));
addParameter(p, 'output_dir', './CINDERellA_results/', @ischar);
addParameter(p, 'sampler', '', @ischar); % Empty default - will be set based on hprior
addParameter(p, 'max_parents', 3, @(x) isnumeric(x) && x > 0);
addParameter(p, 'runtime_minutes', 0.5, @(x) isnumeric(x) && x > 0);
addParameter(p, 'num_samples', 100, @(x) isnumeric(x) && x > 0);
addParameter(p, 'edge_threshold', 0.3, @(x) isnumeric(x) && x >= 0 && x <= 1);
addParameter(p, 'prior_matrix', [], @(x) isempty(x) || (isnumeric(x) && all(ismember(x(:), [0, 1]))));
addParameter(p, 'force_recompute', false, @islogical);
addParameter(p, 'layout', 'force', @ischar);

parse(p, varargin{:});

% Extract parameters
out_dir = p.Results.output_dir;
sampler = p.Results.sampler;
pa_limit = p.Results.max_parents;
runt = p.Results.runtime_minutes * 60; % Convert to seconds for internal use
nSampGraph = p.Results.num_samples;
edge_threshold = p.Results.edge_threshold;
prior_matrix = p.Results.prior_matrix;
force_recompute = p.Results.force_recompute;
layout = p.Results.layout;

%% Input validation
if isempty(data_matrix)
    error('Data matrix cannot be empty');
end

[nGenes, nSamples] = size(data_matrix);
if nGenes < 2 || nSamples < 2
    error('Data matrix must have at least 2 genes and 2 samples');
end

% Set default sampler based on whether hprior is used
if isempty(sampler)
    if isempty(prior_matrix)
        sampler = 'M.c2PB';  % Default for no prior
    else
        sampler = 'M.REV50';  % Default for with prior (hierarchical prior)
    end
end

% Validate sampler
valid_samplers = {'STR','1PB','2PB','3PB','4PB','c2PB','c3PB','c4PB','REV50',...
                 'M.STR','M.1PB','M.2PB','M.3PB','M.4PB','M.c2PB','M.c3PB','M.c4PB','M.REV50'};
if ~ismember(sampler, valid_samplers)
    warning('Unknown sampler: %s. Using default based on prior usage', sampler);
    if isempty(prior_matrix)
        sampler = 'REV50';
    else
        sampler = 'M.REV50';
    end
end

% Validate prior matrix if provided
if ~isempty(prior_matrix)
    if size(prior_matrix, 1) ~= nGenes || size(prior_matrix, 2) ~= nGenes
        error('Prior matrix must be %d x %d (nGenes x nGenes)', nGenes, nGenes);
    end
    if ~all(ismember(prior_matrix(:), [0, 1]))
        error('Prior matrix must contain only 0s and 1s (0=disallowed, 1=allowed)');
    end
    fprintf('Using provided prior matrix to constrain edge search space\n');
    fprintf('Allowed edges: %d out of %d possible\n', sum(prior_matrix(:)), nGenes^2);
end

%% Display settings
fprintf('\n=== CINDERellA Bayesian Network Learning ===\n');
fprintf('Data dimensions: %d genes x %d samples\n', nGenes, nSamples);
fprintf('Output directory: %s\n', out_dir);
fprintf('MCMC sampler: %s\n', sampler);

if strncmp(sampler, 'M.', 2)
    fprintf('Sampling strategy: Multi-chain (%d chains)\n', nSampGraph);
    fprintf('Each chain runs for %.1f seconds\n', runt/nSampGraph);
    fprintf('Final network state collected from each chain\n');
else
    fprintf('Sampling strategy: Single chain\n');
    fprintf('Networks sampled every %.1f seconds\n', runt/nSampGraph);
    fprintf('Total samples collected: %d\n', nSampGraph);
end

fprintf('Max parents per gene: %d\n', pa_limit);
fprintf('Total runtime: %.1f minutes\n', runt/60);
fprintf('Edge threshold for visualization: %.2f\n', edge_threshold);
fprintf('===============================================\n\n');



%% Setup parameters
Param.pa_limit = pa_limit;
Param.runt = runt;  % Now in seconds
Param.nSampGraph = nSampGraph;
Param.mcmethod = sampler;
Param.edge_threshold = edge_threshold;

% Configure hierarchical prior based on user input
if isempty(prior_matrix)
    Param.hprior = 'no';
    Param.cmat = 0;
else
    Param.hprior = 'yes';
    Param.cmat = 0;  % Set to 0 so flearnstruct loads LS from file
    Param.savepath.hprior = [out_dir, '/Hprior.mat'];
    % Save the prior matrix for flearnstruct to load
    Hprior.skeleton = prior_matrix;
    save(Param.savepath.hprior, 'Hprior');
    fprintf('Prior matrix applied: constraining %d edges out of %d total\n', ...
            sum(prior_matrix(:) == 0), nGenes^2);
end

Param.fixedge = 'no';
Param.eQTLnet = 0;

%% Create output directory
Param.out_dir = [out_dir, '/'];
if ~exist(Param.out_dir, 'dir')
    mkdir(Param.out_dir);
    fprintf('Created output directory: %s\n', Param.out_dir);
end

%% Read expression data
fprintf('Using provided data matrix...\n');
fprintf('Data dimensions: %d genes, %d samples\n', nGenes, nSamples);

Param.nSample = nSamples;
Param.nGene = nGenes;

%% Set output file paths
Param.savepath.param = [Param.out_dir, 'Param.mat'];
Param.savepath.ls = [Param.out_dir, 'LS.mat'];
Param.savepath.mcmc = [Param.out_dir, 'Mcmc.mat'];
Param.savepath.edgefrq = [Param.out_dir, 'edgefrq.txt'];
Param.savepath.diagnostics = [Param.out_dir, 'mcmc_diagnostics.png'];
Param.savepath.network_viz = [Param.out_dir, 'network_visualization.png'];

% Save parameters
save(Param.savepath.param, 'Param');

fprintf('param.savepath.ls: %s\n', Param.savepath.ls);
fprintf('force_recompute: %d\n', force_recompute);
%% Calculate Local Scores
if ~force_recompute && exist(Param.savepath.ls, 'file')
    fprintf('Loading existing local scores from: %s\n', Param.savepath.ls);
    load(Param.savepath.ls, 'LS');
else
    fprintf('Calculating local scores...\n');
    try
        if isempty(prior_matrix)
            LS = CcalcLS(data_matrix, 'pa_limit', Param.pa_limit);
        else
            fprintf('Applying prior matrix constraints to local score calculation...\n');
            LS = CcalcLS(data_matrix, 'pa_limit', Param.pa_limit, 'hprior', prior_matrix);
        end
        save(Param.savepath.ls, 'LS', '-v7.3');
        fprintf('Local scores calculated and saved.\n');
    catch ME
        error('Error calculating local scores: %s', ME.message);
    end
end

%% Learn network structure
if ~force_recompute && exist(Param.savepath.mcmc, 'file')
    fprintf('Loading existing MCMC results from: %s\n', Param.savepath.mcmc);
    load(Param.savepath.mcmc, 'Mcmc');
else
    fprintf('Learning network structure (this may take a while)...\n');
    
    % Memory optimization: check available memory
    try
        if ispc
            [~, sys] = memory;
            avail_mem_gb = sys.PhysicalMemory.Available / 1e9;
            fprintf('Available memory: %.2f GB\n', avail_mem_gb);
            if avail_mem_gb < 2
                warning('Low memory detected (%.2f GB). Consider reducing nodes or runtime.', avail_mem_gb);
            end
        end
    catch
        % Memory check not available on this system
    end
    
    % Clear any large temporary variables before learning
    clear data_matrix;
    
    rng(1234); % Set random seed for reproducibility

    try
        tic;
        Mcmc = flearnstruct(Param, 1);
        elapsed_time = toc;
        fprintf('Network learning completed in %.1f minutes!\n', elapsed_time/60);
        
        % Save immediately to free memory
        save(Param.savepath.mcmc, 'Mcmc', '-v7.3');
        fprintf('MCMC results saved successfully.\n');
    catch ME
        % Enhanced error reporting for memory issues
        if contains(ME.message, 'memory') || contains(ME.message, 'Out of memory')
            error(['Out of memory during network learning. Suggestions:\n' ...
                   '1. Reduce number of nodes/features\n' ...
                   '2. Reduce runtime_minutes or num_samples\n' ...
                   '3. Use a machine with more RAM\n' ...
                   '4. Reduce max_parents parameter\n' ...
                   'Original error: %s'], ME.message);
        else
            error('Error during network learning: %s', ME.message);
        end
    end
end

%% Calculate edge frequencies
if ~force_recompute && exist(Param.savepath.edgefrq, 'file')
    fprintf('Loading existing edge frequencies from: %s\n', Param.savepath.edgefrq);
    edgefrq_data = dlmread(Param.savepath.edgefrq, '\t');
    row = edgefrq_data(:,1);
    col = edgefrq_data(:,2);
    v = edgefrq_data(:,3);
    edgefrq = sparse(row, col, v, nGenes, nGenes);
else
    fprintf('Calculating edge frequencies...\n');
    edgefrq = countedgefreq(Mcmc.sampled_graphs);
    [row, col, v] = find(edgefrq);

    % Save edge frequencies
    dlmwrite(Param.savepath.edgefrq, [row col v], 'delimiter', '\t');
    fprintf('Edge frequencies saved to: %s\n', Param.savepath.edgefrq);
end

%% Generate MCMC Diagnostics Plot
if ~force_recompute && exist(Param.savepath.diagnostics, 'file')
    fprintf('MCMC diagnostics plot already exists: %s\n', Param.savepath.diagnostics);
else
    fprintf('Generating MCMC diagnostics plot...\n');
    try
        plot_mcmc_diagnostics(Mcmc, Param.savepath.diagnostics);
        fprintf('MCMC diagnostics plot saved to: %s\n', Param.savepath.diagnostics);
    catch ME
        warning('Could not generate diagnostics plot: %s', ME.message);
    end
end

%% Generate Network Visualization
if ~force_recompute && exist(Param.savepath.network_viz, 'file')
    fprintf('Network visualization already exists: %s\n', Param.savepath.network_viz);
else
    fprintf('Generating network visualization...\n');
    try
        plot_network(edgefrq, edge_threshold, Param.savepath.network_viz, nGenes, layout);
        fprintf('Network visualization saved to: %s\n', Param.savepath.network_viz);
    catch ME
        warning('Could not generate network visualization: %s', ME.message);
    end
end

%% Summary
fprintf('\n=== Results Summary ===\n');
fprintf('Total edges found: %d\n', length(row));
fprintf('Edges above threshold (%.2f): %d\n', edge_threshold, sum(v >= edge_threshold));
fprintf('Average edge frequency: %.3f\n', mean(v));
fprintf('Sampling method used: %s\n', sampler);
if strncmp(sampler, 'M.', 2)
    fprintf('Multi-chain sampling: %d chains × %.1fs each\n', nSampGraph, runt/nSampGraph);
else
    fprintf('Single-chain sampling: %d samples over %.1f minutes\n', nSampGraph, runt/60);
end
fprintf('Results saved in: %s\n', Param.out_dir);
fprintf('Key output files:\n');
fprintf('  - edgefrq.txt: Edge frequencies (use for evaluation)\n');
fprintf('  - Mcmc.mat: Sampled networks\n');
fprintf('  - Param.mat: Analysis parameters\n');
fprintf('  - mcmc_diagnostics.png: Log likelihood trace plot\n');
fprintf('  - network_visualization.png: Network with edges > %.2f\n', edge_threshold);
fprintf('\nTo evaluate results:\n');
fprintf('  Use the edgefrq.txt file with your evaluation script\n');
fprintf('  Load true network and compare with learned edge frequencies\n');
fprintf('======================\n\n');

end

%% Helper function to plot MCMC diagnostics
function plot_mcmc_diagnostics(Mcmc, save_path)
    if ~isfield(Mcmc, 'LL') || isempty(Mcmc.LL)
        warning('No log likelihood data found in MCMC results');
        return;
    end
    
    figure('Position', [100, 100, 800, 400], 'Visible', 'off');
    
    % Plot 1: Log likelihood trace
    subplot(1, 2, 1);
    plot(Mcmc.LL, 'b-', 'LineWidth', 1.5);
    xlabel('MCMC Iteration');
    ylabel('Log Likelihood');
    title('MCMC Trace Plot');
    grid on;
    
    % Plot 2: Log likelihood histogram
    subplot(1, 2, 2);
    histogram(Mcmc.LL, 30, 'FaceColor', [0.5 0.8 1], 'EdgeColor', 'black');
    xlabel('Log Likelihood');
    ylabel('Frequency');
    title('Log Likelihood Distribution');
    grid on;
    
    % Add super title (compatible with older MATLAB versions)
    if exist('sgtitle', 'file')
        sgtitle('MCMC Diagnostics', 'FontSize', 14, 'FontWeight', 'bold');
    else
        % Alternative for older MATLAB versions
        annotation('textbox', [0.35 0.95 0.3 0.05], 'String', 'MCMC Diagnostics', ...
                   'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', ...
                   'EdgeColor', 'none');
    end
    
    % Save the plot
    print(save_path, '-dpng', '-r300');
    close;
end

%% Helper function to visualize network
function plot_network(edgefrq, threshold, save_path, nGenes, layout)
    % Find edges above threshold
    [row, col, v] = find(edgefrq >= threshold);
    
    if isempty(row)
        warning('No edges found above threshold %.2f', threshold);
        return;
    end
    
    figure('Position', [100, 100, 800, 800], 'Visible', 'off');
    
    % Create adjacency matrix for visualization
    adj_matrix = sparse(row, col, v, nGenes, nGenes);
    
    % Create graph object
    G = digraph(adj_matrix);
    
    % Validate layout option
    valid_layouts = {'circle', 'force', 'layered', 'subspace'};
    if ~ismember(layout, valid_layouts)
        warning('Unknown layout "%s". Using "force" layout.', layout);
        layout = 'force';
    end
    
    % Plot the network with specified layout
    h = plot(G, 'Layout', layout, 'MarkerSize', 8, 'NodeColor', [0.7 0.9 1], ...
             'EdgeColor', 'black', 'ArrowSize', 10);
    
    % Color edges by frequency (compatible approach)
    try
        % Try modern approach first (R2016a+)
        edge_weights = G.Edges.Weight;
        colormap(hot);
        h.EdgeCData = edge_weights;
        h.LineWidth = 2 * edge_weights; % Scale line width by frequency
    catch
        % Fallback for older MATLAB versions
        colormap(hot);
        % Set uniform edge properties
        h.LineWidth = 2;
        % Manual edge coloring would require more complex code for older versions
        fprintf('Note: Advanced edge coloring not available in this MATLAB version\n');
    end
    
    % Add colorbar (only if edge coloring worked)
    try
        cb = colorbar;
        cb.Label.String = 'Edge Frequency';
    catch
        % Fallback: create manual colorbar or skip
        fprintf('Note: Colorbar not available for this configuration\n');
    end
    
    % Add labels
    node_labels = cellstr(num2str((1:nGenes)'));
    h.NodeLabel = node_labels;
    
    title(sprintf('Bayesian Network (%s layout, edges ≥ %.2f)', layout, threshold), ...
          'FontSize', 14, 'FontWeight', 'bold');
    
    % Add legend/text box with statistics
    stats_text = sprintf('Nodes: %d\nEdges: %d\nLayout: %s\nThreshold: %.2f\nMax frequency: %.3f', ...
                        nGenes, length(row), layout, threshold, max(v));
    annotation('textbox', [0.02, 0.02, 0.25, 0.18], 'String', stats_text, ...
               'BackgroundColor', 'white', 'EdgeColor', 'black', 'FontSize', 10);
    
    % Save the plot
    print(save_path, '-dpng', '-r300');
    close;
end

%% Helper function to display usage
function display_help()
    fprintf('\nCINDERellA - Bayesian Network Learning\n');
    fprintf('======================================\n\n');
    fprintf('REQUIRED:\n');
    fprintf('  data - Path to expression data file\n\n');
    fprintf('OPTIONAL PARAMETERS:\n');
    fprintf('  output_dir       - Output directory (default: ./CINDERellA_results/)\n');
    fprintf('  sampler          - MCMC sampler (default: REV50, or M.REV50 if hprior used)\n');
    fprintf('                    Single-chain: STR, c2PB, c3PB, c4PB, 1PB, 2PB, 3PB, 4PB, REV50\n');
    fprintf('                    Multi-chain: M.STR, M.c2PB, M.c3PB, M.c4PB, M.1PB, M.2PB, M.3PB, M.4PB, M.REV50\n');
    fprintf('  max_parents      - Max parents per gene (default: 3)\n');
    fprintf('  runtime_minutes  - Total runtime in minutes (default: 30)\n');
    fprintf('  num_samples      - Number of networks to sample (default: 100)\n');
    fprintf('                    Single-chain: samples every (runtime/num_samples) time\n');
    fprintf('                    Multi-chain: runs num_samples chains for (runtime/num_samples) each\n');
    fprintf('  edge_threshold   - Edge frequency threshold for visualization (default: 0.3)\n');
    fprintf('  force_recompute  - Force recomputation even if results exist (default: false)\n');
    fprintf('  layout           - Network layout algorithm: circle, force, layered, subspace (default: force)\n');
    fprintf('  prior_matrix     - nGenes x nGenes binary matrix to constrain edges (optional)\n\n');
    fprintf('EXAMPLES:\n');
    fprintf('  CINDERellA(''data'', ''expression_data.txt'');\n');
    fprintf('  CINDERellA(''data'', ''my_data.txt'', ''max_parents'', 5);\n');
    fprintf('  CINDERellA(expdata.data, ''edge_threshold'', 0.5);\n');
    fprintf('  \n');
    fprintf('  %% With prior knowledge:\n');
    fprintf('  prior = ones(nGenes, nGenes); prior(1,2) = 0;  %% Disallow specific edge\n');
    fprintf('  CINDERellA(expdata.data, ''prior_matrix'', prior);\n\n');
    fprintf('DATA FORMAT:\n');
    fprintf('  Tab-separated file with genes as rows, samples as columns\n\n');
end