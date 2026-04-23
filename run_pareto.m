clear; clc; close all;
addpath('TLAlgorithms/');
addpath('utils/');

%%%%%%%%%%%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%

% T0_list = [6, 8, 10];                                        % sparsity level for each representation
T0_list = [6, 8];                                        % sparsity level for each representation

% patch_size_list = [64, 121, 196, 256];                       % patch size, transform size

% patch_size_list = [196, 256];                                % patch size, transform size

patch_size_list = [64];                                 % patch size, transform size

bresler_lambdas = linspace(2.1e-4, 2.1e-14, 10);             % input \lambda parameter for Bresler method

n_pareto_points = 20;                                        % number of points/trials for a pareto plot

ell_1_lambda = 1e-2;                                         % input \lambda parameter for ell_1 penalization in the Structured Conditioned Method

homotopy_steps = 100;

T1_bresler = linspace(1.0, 0.05, n_pareto_points);           % input sparsity for Bresler method

numiter = 2000;                                               % number of iterations for the Alternating Minimization algorithms with warm start

%%% datasets
% input_folders_datasets = {'data', 'DIV2K_valid_HR_data'};
% train_sets = { {'Barbara', 'Couple', 'Lena', 'Cameraman'}, ...
%                {'0801', '0802', '0803', '0804', '0805', '0806', '0807'} };
% test_sets  = { {'Hill', 'Man', 'Baboon', 'MRI'}, ...
%                {'0808', '0809', '0810'} };

% input_folders_datasets = {'DIV2K_valid_HR_data'};
% train_sets = { {'0801', '0802', '0803', '0804', '0805', '0806', '0807'} };
% test_sets  = { {'0808', '0809'} };
% test_sets  = { {'0808', '0809', '0810'} };

input_folders_datasets = {'data'};
train_sets = { {'Barbara', 'Couple', 'Lena', 'Cameraman'} };
test_sets  = { {'Hill', 'Man', 'Baboon', 'MRI'}};

% generate parameters grid and datasets configuration
[parameters_settings, datasets_by_n] = generate_settings_pareto(bresler_lambdas, T0_list, patch_size_list, input_folders_datasets, train_sets, test_sets);

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

for ps_idx = 1:length(parameters_settings)

    time_config = tic;

    fprintf('\n-------------- Starting Configuration %d out of %d --------------\n', ps_idx, length(parameters_settings));

    setting = parameters_settings(ps_idx);

    T0 = setting.T0;
    n = setting.patch_size;
    idx_n = setting.patch_size_idx;
    dataset_idx = setting.dataset_idx;

    bresler_lambda = setting.bresler_lambda;

    YH_train = datasets_by_n(idx_n).YH_train{dataset_idx};
    YH_test  = datasets_by_n(idx_n).YH_test{dataset_idx};

    W0 = kron(dctmtx(sqrt(n)), dctmtx(sqrt(n)));
    B0 = eye(n);

    YH2_train = W0 * YH_train; YH2_test = W0 * YH_test;

    STY_tr = T0 * ones(1, size(YH_train, 2)); STY_te = T0 * ones(1, size(YH_test, 2));

    l2_bresler = bresler_lambda * (norm(YH_train, 'fro'))^2;
    l3_bresler = l2_bresler;

    T1_list = round(T1_bresler * n^2);

    results = repmat(struct(), 1, n_pareto_points);

    % run methods
    for idx = 1:n_pareto_points

        %%% Bresler Method
        tic;
        [W_bresler, ~, error_bresler] = StructuredBreslerCF(B0, YH2_train, YH2_test, numiter, l2_bresler, l3_bresler, T1_list(idx), STY_tr, STY_te);
        time_bresler = toc;

        total = numel(W_bresler);           
        curr_sty = nnz(W_bresler(:) ~= 0);  
        sty_pct = 100 * curr_sty / total;

        fprintf('-> Done Bresler with sparsity sty_pct %.2f - iteration %d out of %d - %.2f seconds\n', sty_pct, idx, n_pareto_points, time_bresler);

        results(idx).bresler_method.transform = sparse(W_bresler);
        results(idx).bresler_method.error = error_bresler;
        results(idx).bresler_method.sty_pct = sty_pct;
        results(idx).bresler_method.time = time_bresler;

        %%% get rho and tau
        rho = cond(W_bresler);
        tau = norm(W_bresler, 'fro');

        results(idx).rho = rho;
        results(idx).tau = tau;

        %%% Structured Conditioned Method
        tic;
        [T_sc, ~, error_sc, sty_pct_sc, sty_vec_sc] = StructuredConditioned(W0, YH2_train, YH2_test, numiter, STY_tr, STY_te, rho, tau, ell_1_lambda, homotopy_steps);
        time_sc = toc;

        fprintf('-> Done Structured Conditioned with sparsity %.2f - iteration %d out of %d - %.2f seconds\n', sty_pct_sc, idx, n_pareto_points, time_sc);

        results(idx).structured_conditioned.transform = sparse(T_sc);
        results(idx).structured_conditioned.error = error_sc;
        results(idx).structured_conditioned.sty_pct = sty_pct_sc;
        results(idx).structured_conditioned.sty_vec = sty_vec_sc;
        results(idx).structured_conditioned.time = time_sc;
    end

    time_config = toc(time_config);

    time_config;

    file = sprintf('results_%d.mat', ps_idx);
    path = fullfile(folder, file);
    save(path, 'setting', 'results', 'time_config');
end