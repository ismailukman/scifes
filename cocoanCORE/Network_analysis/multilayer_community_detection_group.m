function N_all = multilayer_community_detection_group(N, varargin)
    %%1. Default variables
    n_repeat = 100;
    thresh_type = 'max';
    gamma = 1;
    allegiance_mat_out = 0;

    %2. Process optional input arguments (varargin) and update parameters if provided
    for i = 1:length(varargin)
        if ischar(varargin{i})
            switch varargin{i}
                case 'n_repeat'
                    n_repeat = varargin{i+1};
                case 'thresh_type'
                    thresh_type = varargin{i+1};
                case 'gamma'
                    gamma = varargin{i+1};
            end
        end
    end

    %3. Extract key matrix dimensions
    n_subj = numel(N);
    n_run = size(N{1}.multi_module_consensus, 2);
    n_node = size(N{1}.multi_module_consensus, 1);

    fprintf('n_subj is %d, n_run is %d\n and n_node is %d \n', n_subj, n_run, n_node);

    for run_i = 1:n_run
        for cons_iter = 1:10
            if cons_iter == 1
                cons_comm = NaN(n_node, n_subj);
                for subj_i = 1:n_subj
                    cons_comm(:, subj_i) = N{subj_i}.multi_module_consensus(:, run_i);
                end
            end

            % Compute allegiance matrix
            [allegiance_mat, perm_max(cons_iter, 1)] = calc_allegiance(cons_comm, 'threshold', thresh_type);

            % Filter allegiance matrix
            wh_conn = any(allegiance_mat ~= 0);
            allegiance_mat = allegiance_mat(wh_conn, wh_conn);

            % Initialize consensus community
            cons_comm = NaN(n_node, n_repeat);

            %run checks on A and g length
            % fprintf('A, g length checksB.');
            % disp(size(allegiance_mat, 1));
            % disp(length(gamma));

            % fprintf('check [A(i, j)] %d\n', allegiance_mat_out); % ?

            % fprintf('[A(1, 1)] %d\n', allegiance_mat(1, 1) ); % ?

            [B, ~] = modularity(allegiance_mat, gamma);
            allegiance_mat = [];
            for rep_i = 1:n_repeat
                [S, ~] = genlouvain(B, 10000, 0);
                cons_comm(wh_conn, rep_i) = S;
            end
            B = [];

            % Check for consensus
            if all(all(diff(cons_comm(wh_conn, :), [], 2) == 0))
                N_all.multi_module_consensus(:, run_i) = cons_comm(:, 1);
                N_all.multi_niter_consensus(:, run_i) = cons_iter;
                N_all.multi_permmax_consensus{run_i} = perm_max;
                perm_max = [];
                break;
            end
        end
    end

    % Match community affiliation
    N_all.multi_module_consensus_sorted = match_community_affiliation(N_all.multi_module_consensus);
end
