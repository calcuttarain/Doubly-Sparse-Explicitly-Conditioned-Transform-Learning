clear; clc; close all;
addpath('utils/');

folder_path = 'results_2026-03-28_1920'; 
file_name = 'results_1.mat';
full_path = fullfile(folder_path, file_name);

target_idx = 1; 

show_plots = true;
output_folder = folder_path; 

if ~exist(full_path, 'file')
    error('File not found at: %s', full_path);
end

load(full_path, 'setting', 'results');

kappa = setting.rho;
tau = setting.tau;
T = setting.T0;

res_bresler = results(target_idx).bresler_method;
res_sc = results(target_idx).structured_conditioned;

n_pareto_points = length(results);
ell_1_lambdas = linspace(0, 2.1e-1, n_pareto_points);
lambda_sc = ell_1_lambdas(target_idx);

errors_tr = {res_bresler.error.m1.tr, res_sc.error.m1.tr};
errors2_tr = {res_bresler.error.m2.tr, res_sc.error.m2.tr};

errors_te = {res_bresler.error.m1.te, res_sc.error.m1.te};
errors2_te = {res_bresler.error.m2.te, res_sc.error.m2.te};

labels = {'Bresler', 'Structured Conditioned'};
numiter = length(errors_tr{1}); 
sty_pct_sc = res_sc.sty_pct;

base_name_tr = sprintf("plot_conv_train_idx%02d", target_idx);
plot_convergence(numiter, errors_tr, errors2_tr, labels, kappa, tau, ...
                 sty_pct_sc, T, 'Training', output_folder, base_name_tr, show_plots);

base_name_te = sprintf("plot_conv_test_idx%02d", target_idx);
plot_convergence(numiter, errors_te, errors2_te, labels, kappa, tau, ...
                 sty_pct_sc, T, 'Testing', output_folder, base_name_te, show_plots);

base_name_spy = sprintf("plot_sparsity_idx%02d", target_idx);
plot_sparsity(res_sc.transform, numiter, res_sc.error.m1.tr, res_sc.sty_vec, ...
              kappa, tau, lambda_sc, T, output_folder, base_name_spy, show_plots);