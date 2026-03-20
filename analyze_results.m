clear; clc; close all;
addpath('utils/');

%%% test parameters
timestamp = '2026-03-20_1422';     % for example
parameters_settings_idx = 1;       % for example

input_folder = "results_" + timestamp; 
output_folder = "plots_" + timestamp;

%%% load results from batches
load(fullfile(input_folder,"global_settings.mat"));

files = dir(fullfile(input_folder, 'results_batch_*_to_*.mat'));
target_file = '';
start_idx_of_batch = 0;

for k = 1:length(files)
    tokens = regexp(files(k).name, 'results_batch_(\d+)_to_(\d+)\.mat', 'tokens');
    if ~isempty(tokens)
        s_idx = str2double(tokens{1}{1});
        e_idx = str2double(tokens{1}{2});
        
        if parameters_settings_idx >= s_idx && parameters_settings_idx <= e_idx
            target_file = files(k).name;
            start_idx_of_batch = s_idx;
            break;
        end
    end
end

if isempty(target_file)
    error('Batch with parameters_settings_idx %d not found!', parameters_settings_idx);
end

load(fullfile(input_folder, target_file));

local_idx = parameters_settings_idx - start_idx_of_batch + 1;
results = batch_results{local_idx}.result;

subfolder = sprintf('plots_P%d', parameters_settings_idx);
target_folder = fullfile(output_folder, subfolder);

% Bresler methods for loop
for iter = 1:length(results)
    status = results{iter}.status;

    %%% handle error
    if strcmp(status, 'Failed') || strcmp(status, 'Skipped_BadCondition')
        fprintf('--- Error Report ---\n');
        fprintf('Status: %s\n', status);
        fprintf('Method: %s\n', bresler_methods{iter});
        
        if isfield(results{iter}, 'ME') && ~isempty(results{iter}.ME)
            ME = results{iter}.ME;
            fprintf('Message: %s\n', ME.message);
            fprintf('Identifier: %s\n', ME.identifier);
            if ~isempty(ME.stack)
                st = ME.stack(1);
                fprintf('Location: File "%s", Line %d, Function "%s"\n', ...
                    st.file, st.line, st.name);
            end
        end
        fprintf('--------------------\n');

        continue;
    end

    %%% extract parameters
    paramsin = results{iter}.paramsin;
    paramsout = results{iter}.paramsout;
    global_params = parameters_settings(paramsin.parameters_settings_idx);

    methods = fieldnames(paramsout);

    n = global_params.patch_size;

    T0 = global_params.T0;
    T1 = global_params.T1;

    rho = paramsin.rho;
    tau = paramsin.tau;

    proposed_method = methods{end}; % 'StructuredConditioned'
    bresler_method = paramsin.bresler_method;

    final_folder = fullfile(target_folder, bresler_method);

    if ~exist(final_folder, 'dir')
        mkdir(final_folder);
    end

    n_methods = length(methods);

    %%% plots
    errors_train = cell(1, n_methods); errors2_train = cell(1, n_methods);
    errors_test = cell(1, n_methods); errors2_test = cell(1, n_methods);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Train / Test Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

    show_plots = false;

    title_text = 'Train Data';
    file_name = "convergence_train";
    plot_convergence(numiter, errors_train, errors2_train, methods, rho, tau, sty_pct_proposed, T0, title_text, final_folder, file_name, show_plots);

    file_name = "convergence_test";
    title_text = 'Test Data';
    plot_convergence(numiter, errors_test, errors2_test, methods, rho, tau, sty_pct_proposed, T0, title_text, final_folder, file_name, show_plots);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Sparsity %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    file_name = proposed_method;
    plot_sparsity(Transform_proposed, numiter, error_proposed, sty_vec_proposed, rho, tau, lambda_proposed, T0, final_folder, file_name, show_plots);
end