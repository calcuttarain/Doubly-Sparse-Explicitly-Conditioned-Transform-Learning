clear; clc; close all;
addpath('utils/');

folder = 'results';
files = dir(fullfile(folder, '*.mat'));

[~, idx] = sort([files.datenum], 'descend');
file = files(idx(1)).name;

path = fullfile(folder, file);

load_data(path);

proposed_method = methods{end}; % 'structured_conditioned'

n_methods = length(methods);

errors_train = cell(1, n_methods); errors2_train = cell(1, n_methods);
errors_test = cell(1, n_methods); errors2_test = cell(1, n_methods);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Train / Test Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

labels = replace(methods, '_', ' ');
labels = regexprep(labels, '(^|[ ]).', '${upper($0)}');

Transform_proposed = paramsout.(proposed_method).transform;
error_proposed = paramsout.(proposed_method).error.m1.tr;
sty_pct_proposed = paramsout.(proposed_method).sty_pct;
sty_vec_proposed = paramsout.(proposed_method).sty_vec;
lambda_proposed = paramsout.(proposed_method).lambda;

for i = 1:n_methods
    errors_train{i} = paramsout.(methods{i}).error.m1.tr; 
    errors_test{i} = paramsout.(methods{i}).error.m1.te;

    errors2_train{i} = paramsout.(methods{i}).error.m2.tr; 
    errors2_test{i} = paramsout.(methods{i}).error.m2.te;
end

title_text = 'Train Data';
file_name = "convergence_" + proposed_method + "_train";
plot_convergence(numiter, errors_train, errors2_train, labels, rho, tau, sty_pct_proposed, T0, title_text, file_name);

file_name = "convergence_" + proposed_method + "_test";
title_text = 'Test Data';
plot_convergence(numiter, errors_test, errors2_test, labels, rho, tau, sty_pct_proposed, T0, title_text, file_name);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Sparsity %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

file_name = proposed_method;
plot_sparsity(Transform_proposed, numiter, error_proposed, sty_vec_proposed, rho, tau, lambda_proposed, T0, file_name);
