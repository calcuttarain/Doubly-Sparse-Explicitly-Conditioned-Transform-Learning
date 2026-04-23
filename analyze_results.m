clear; clc; close all;
addpath('utils/');

chosen_timestamp = "2026-04-21_1231";
target_result = "17";
target_idx = 8; 

folder_path = "results_" + chosen_timestamp; 

target_dir = fullfile("plots_" + chosen_timestamp, 'general', "results_" + target_result, "idx_" + string(target_idx));

if ~exist(char(target_dir), "dir")
    mkdir(char(target_dir));
end

file_name = "results_" + target_result + ".mat";

full_path = fullfile(folder_path, file_name);

show_plots = false;
output_folder = target_dir; 

if ~exist(full_path, 'file')
    error('File not found at: %s', full_path);
end

load(full_path, 'setting', 'results');

T = setting.T0;

res_bresler = results(target_idx).bresler_method;
res_sc = results(target_idx).structured_conditioned;

kappa = results(target_idx).rho;
tau = results(target_idx).tau;

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

base_name_tr = "plot_conv_train";
plot_convergence(numiter, errors_tr, errors2_tr, labels, kappa, tau, ...
                 sty_pct_sc, T, 'Training', output_folder, base_name_tr, show_plots);

base_name_te = "plot_conv_test";
plot_convergence(numiter, errors_te, errors2_te, labels, kappa, tau, ...
                 sty_pct_sc, T, 'Testing', output_folder, base_name_te, show_plots);

base_name_spy = "plot_sparsity";
plot_sparsity(res_sc.transform, numiter, res_sc.error.m1.tr, res_sc.sty_vec, ...
              kappa, tau, lambda_sc, T, output_folder, base_name_spy, show_plots);
