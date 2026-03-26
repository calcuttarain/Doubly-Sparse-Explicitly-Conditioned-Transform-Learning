clear; clc; close all;
addpath('TLAlgorithms/');
addpath('utils/');

%%%%%%%%%%%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%

T0_list = [6, 8, 10];                                        % sparsity level for each representation

patch_size_list = [64, 121, 196, 256];                       % patch size, transform size

bresler_lambdas = linspace(2.1e-1, 2.1e-13, 20);             % input \lambda parameter for Bresler method

n_pareto_points = 50;                                        % number of points/trials for a pareto plot

ell_1_lambdas = linspace(2.1e-1, 2.1e-13, n_pareto_points);  % input \lambda parameter for ell_1 penalization in the Structured Conditioned Method
T1_bresler = linspace(0.05, 1.0, n_pareto_points);           % input sparsity for Bresler method

numiter = 3000;                                              % number of iterations for the Alternating Minimization algorithms

%%% datasets
input_folders_datasets = {'data', 'DIV2K_valid_HR_data'};
train_sets = { {'Barbara', 'Couple', 'Lena', 'Cameraman'}, ...
               {'0801', '0802', '0803', '0804', '0805', '0806', '0807'} };
test_sets  = { {'Hill', 'Man', 'Baboon', 'MRI'}, ...
               {'0808', '0809', '0810'} };

% generate parameters grid and datasets configuration
[parameters_settings, datasets_by_n] = generate_settings_pareto(T0_list, patch_size_list, input_folders_datasets, train_sets, test_sets);

% save global settings
timestamp = datestr(now, 'yyyy-mm-dd_HHMM');
folder = sprintf('results_%s', timestamp);

if ~exist(folder, 'dir')
    mkdir(folder);
end

file_name = 'global_settings.mat';
path = fullfile(folder, file_name);
save(path);

%%%%%%%%%%%%%%%%%%%%%%%% For Loop %%%%%%%%%%%%%%%%%%%%%%%%

for bl_idx = 1:length(bresler_lambdas)
    for ps_idx = 1:length(parameters_settings)

        setting = parameters_settings(ps_idx);

        T0 = setting.T0;
        n = setting.patch_size;
        idx_n = setting.patch_size_idx;
        dataset_idx = setting.dataset_idx;

        bresler_lambda = bresler_lambdas(bl_idx);

        YH_train = datasets_by_n(idx_n).YH_train{dataset_idx};
        YH_test  = datasets_by_n(idx_n).YH_test{dataset_idx};

        W0 = kron(dctmtx(sqrt(n)), dctmtx(sqrt(n)));
        B0 = eye(n);

        STY_tr = T0 * ones(1, size(YH_train, 2)); STY_te = T0 * ones(1, size(YH_test, 2));

        l2_bresler = bresler_lambda * (norm(YH_train, 'fro'))^2;
        l3_bresler = l2_bresler;
        
        %%% run Bresler with 100% sty
        [W_bresler, ~, ~] = StructuredBreslerCF(B0, YH_train, YH_test, numiter, l2_bresler, l3_bresler, n^2, STY_tr, STY_te);

        %%% get rho and tau
        rho = cond(W_bresler);
        tau = norm(W_bresler, 'fro');

        setting.rho = rho;
        setting.tau = tau;

        T1_list = round(T1_bresler * n^2);

        results = repmat(struct(), 1, n_pareto_points);

        for idx = 1:n_pareto_points
            %%% Bresler Method
            tic;
            [W_bresler, ~, error_bresler] = StructuredBreslerCF(B0, YH_train, YH_test, numiter, l2_bresler, l3_bresler, T1_list(idx), STY_tr, STY_te);
            time_bresler = toc

            total = numel(W_bresler);           
            curr_sty = nnz(W_bresler(:) ~= 0);  
            sty_pct = 100 * curr_sty / total;

            results(idx).bresler_method.transform = sparse(W_bresler);
            results(idx).bresler_method.error = error_bresler;
            results(idx).bresler_method.sty_pct = sty_pct;
            results(idx).bresler_method.time = time_bresler;

            %%% Structured Conditioned Method
            tic;
            [T_sc, ~, error_sc, sty_pct_sc, sty_vec_sc] = StructuredConditioned(W0, YH_train, YH_test, numiter, STY_tr, STY_te, rho, tau, ell_1_lambdas(idx));
            time_sc = toc

            results(idx).structured_conditioned.transform = sparse(T_sc);
            results(idx).structured_conditioned.error = error_sc;
            results(idx).structured_conditioned.sty_pct = sty_pct_sc;
            results(idx).structured_conditioned.sty_vec = sty_vec_sc;
            results(idx).structured_conditioned.time = time_sc;
        end

        file = sprintf('results_blidx_%07d_psidx_%07d.mat', bl_idx, ps_idx);
        path = fullfile(folder, file);
        save(path, 'setting', 'results');
    end
end