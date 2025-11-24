function check_memory()
% CHECK_MEMORY - Display current memory usage and availability
%
% This function helps diagnose memory issues before running CINDERellA

fprintf('\n=== MEMORY STATUS ===\n');

try
    if ispc
        % Windows system
        [~, sys] = memory;
        
        total_gb = sys.PhysicalMemory.Total / 1e9;
        avail_gb = sys.PhysicalMemory.Available / 1e9;
        used_gb = total_gb - avail_gb;
        used_pct = 100 * used_gb / total_gb;
        
        fprintf('Total Physical Memory:     %.2f GB\n', total_gb);
        fprintf('Available Memory:          %.2f GB\n', avail_gb);
        fprintf('Used Memory:               %.2f GB (%.1f%%)\n', used_gb, used_pct);
        
        % Memory recommendations
        fprintf('\n=== RECOMMENDATIONS ===\n');
        if avail_gb < 2
            fprintf('⚠ WARNING: Very low memory (%.2f GB available)\n', avail_gb);
            fprintf('   Strongly recommend:\n');
            fprintf('   - Close other programs\n');
            fprintf('   - Reduce num_samples to 20-30\n');
            fprintf('   - Reduce max_parents to 2\n');
            fprintf('   - Reduce runtime_minutes to 0.25-0.5\n');
        elseif avail_gb < 4
            fprintf('⚠ CAUTION: Limited memory (%.2f GB available)\n', avail_gb);
            fprintf('   Recommended settings:\n');
            fprintf('   - num_samples: 30-50\n');
            fprintf('   - max_parents: 2\n');
            fprintf('   - runtime_minutes: 0.5-1.0\n');
        elseif avail_gb < 8
            fprintf('✓ Moderate memory (%.2f GB available)\n', avail_gb);
            fprintf('   Recommended settings:\n');
            fprintf('   - num_samples: 50-100\n');
            fprintf('   - max_parents: 2-3\n');
            fprintf('   - runtime_minutes: 1-5\n');
        else
            fprintf('✓ Good memory (%.2f GB available)\n', avail_gb);
            fprintf('   You can use standard settings:\n');
            fprintf('   - num_samples: 100+\n');
            fprintf('   - max_parents: 3\n');
            fprintf('   - runtime_minutes: 5-60\n');
        end
        
    elseif ismac || isunix
        % Mac/Linux system
        if ismac
            [~, result] = system('vm_stat | grep "Pages free" | awk ''{print $3}'' | sed ''s/\.//''');
            free_pages = str2double(result);
            avail_gb = (free_pages * 4096) / 1e9; % Page size is typically 4KB
            fprintf('Available Memory (estimate): %.2f GB\n', avail_gb);
        else
            [~, result] = system('free -g | grep Mem | awk ''{print $7}''');
            avail_gb = str2double(result);
            fprintf('Available Memory: %.2f GB\n', avail_gb);
        end
        
        if avail_gb < 2
            fprintf('⚠ WARNING: Low memory detected\n');
        elseif avail_gb < 4
            fprintf('⚠ CAUTION: Limited memory\n');
        else
            fprintf('✓ Sufficient memory available\n');
        end
    end
    
catch ME
    fprintf('Could not determine memory status: %s\n', ME.message);
    fprintf('Manual check recommended before running large analyses.\n');
end

fprintf('\n=== MATLAB MEMORY ===\n');
try
    % Show MATLAB's current memory usage
    whos_output = evalc('whos');
    fprintf('Current workspace variables:\n');
    disp(whos_output);
catch
    fprintf('Could not display workspace memory.\n');
end

fprintf('\n=== TIPS FOR 238 NODES ===\n');
fprintf('For your network with 238 nodes:\n');
fprintf('• Expected peak memory usage: 4-8 GB\n');
fprintf('• Start with conservative settings and increase gradually\n');
fprintf('• Use force_recompute=false to avoid recalculating LS\n');
fprintf('• Save intermediate results frequently\n');
fprintf('• Close MATLAB and restart if you encounter memory fragmentation\n');
fprintf('=======================\n\n');

end
