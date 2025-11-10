function A_reshaped = reshape_for_multilayer(A_g, n_subjects, T_g, group_name)
    total_windows = size(A_g, 1);
    wins_per_subj = total_windows / n_subjects;
    
    if mod(total_windows, n_subjects) ~= 0
        error('Total windows (%d) not evenly divisible by subjects (%d) for %s', ...
              total_windows, n_subjects, group_name);
    end
    
    fprintf('--- %s: %d windows → %d subjects × %d windows\n', ...
            group_name, total_windows, n_subjects, wins_per_subj);
    
    A_reshaped = cell(n_subjects, wins_per_subj);
    for s = 1:n_subjects
        start_idx = (s-1)*wins_per_subj + 1;
        end_idx = s*wins_per_subj;
        A_reshaped(s, :) = A_g(start_idx:end_idx);
    end
end