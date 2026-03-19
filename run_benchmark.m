warning('off','MATLAB:rankDeficientMatrix');
clear; clc; close all;
addpath('TLAlgorithms/');
addpath('utils/');

%%%%%%%%%%%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%

% T0_list = [6, 8, 10, 15];
% T1_list = [0.15, 0.20, 0.25, 0.30];
% 
% patch_size_list = [64, 121, 196, 256];
T0_list = [6, 8];   % sparsity level for each representation
T1_list = [20];     % Structured Transforms sparsity percent

patch_size_list = [64, 256]; % patch size, transform size

% bresler_lambdas = linspace(2.1e-1, 2.1e-13, 2);
bresler_lambdas = [2.1e-9];
numiter = 200;

%bresler_methods = {'UnstructuredBresler', 'StructuredBresler', 'StructuredBreslerCF'};
bresler_methods = {'UnstructuredBresler'};

%%% datasets
input_folders_datasets = {'data', 'DIV2K_valid_HR_data'};
train_sets = { {'Barbara', 'Couple', 'Lena', 'Cameraman'}, ...
               {'0801', '0802', '0803', '0804', '0805', '0806', '0807'} };
test_sets  = { {'Hill', 'Man', 'Baboon', 'MRI'}, ...
               {'0808', '0809', '0810'} };

% generate parameters grid and datasets configuration
[parameters_settings, datasets_by_n] = generate_settings(T0_list, T1_list, patch_size_list, input_folders_datasets, train_sets, test_sets);

% datasets are frozen for read-only 
SharedDatasets = parallel.pool.Constant(datasets_by_n);

% lambda search loop parameters for Doubly Sparse Conditioned Transform
max_iter = 3;                                      % maximum number of iterations
tol_sty_pct = 1;                                   % sparsity percent tolerance
global_lambda_min = 1e-10;                         % left endpoint of search interval
global_lambda_max = 1e3;                           % right endpoint of search interval

timestamp = datestr(now, 'yyyy-mm-dd_HHMM');
folder = sprintf('results_%s', timestamp);
        
if ~exist(folder, 'dir')
    mkdir(folder);
end

file_name = 'global_settings.mat';
path = fullfile(folder, file_name);
save(path);

%%%%%%%%%%%%%%%%%%%%%%%% Big For Loop %%%%%%%%%%%%%%%%%%%%%%%%
for iter_bresler_lambda = 1:length(bresler_lambdas)
    bresler_lambda = bresler_lambdas(iter_bresler_lambda);

    % paralelize on parameters grid
    parfor iter_param_setting = 1:length(parameters_settings)
        % t = getCurrentTask();
        % if isempty(t)
        %     worker_id = 0; 
        % else
        %     worker_id = t.ID; 
        % end
        % fprintf('Worker %d started iteration %d\n', worker_id, iter_param_setting);
        paramsin = struct(); 
        results = cell(1, length(bresler_methods));

        %%% extract parameters

        paramsin.T0 = parameters_settings(iter_param_setting).T0;
        paramsin.T1 = parameters_settings(iter_param_setting).T1;
        n = parameters_settings(iter_param_setting).patch_size;
        idx_n = parameters_settings(iter_param_setting).patch_size_idx;
        dataset_idx = parameters_settings(iter_param_setting).dataset_idx;

        paramsin.YH_train = SharedDatasets.Value(idx_n).YH_train{dataset_idx};
        paramsin.YH_test  = SharedDatasets.Value(idx_n).YH_test{dataset_idx};

        %%% compute the rest of the parameters

        % initializations
        paramsin.W0 = kron(dctmtx(sqrt(n)), dctmtx(sqrt(n)));
        paramsin.B0 = eye(n);

        % set the sparsity levels
        paramsin.STY_tr = paramsin.T0 * ones(1, size(paramsin.YH_train, 2)); paramsin.STY_te = paramsin.T0 * ones(1, size(paramsin.YH_test, 2));

        paramsin.l2_bresler = bresler_lambda * (norm(paramsin.YH_train, 'fro'))^2;
        paramsin.l3_bresler = paramsin.l2_bresler;

        paramsin.numiter = numiter;

        %%% methods for
        for iter_bresler_method = 1:length(bresler_methods)
            bresler_method = bresler_methods{iter_bresler_method};

            local_paramsin = paramsin;
            rho_max = 1e12;
            paramsout = struct(); 
            result = struct();

            try
                %%%%%%%%%%% Run Bresler Method %%%%%%%%%%%
    
                lastwarn('');
                tic;
                [paramsout.(bresler_method).transform, paramsout.(bresler_method).representation, paramsout.(bresler_method).error] = feval(bresler_method, local_paramsin);
                paramsout.(bresler_method).time = toc;
    
                local_paramsin.rho = cond(paramsout.(bresler_method).transform);
                local_paramsin.tau = norm(paramsout.(bresler_method).transform, 'fro');

                % if the condition number yielded by the bresler method is
                % unstable, throw an error
                if local_paramsin.rho > rho_max || isnan(local_paramsin.rho)
                    error('TRANSFORM:PoorConditioningBreslerMethod', ...
                          'Condition number too high (%.2e).', local_paramsin.rho);
                end

                if ~strcmp(bresler_method, 'UnstructuredBresler')
                    paramsout.(bresler_method).transform = sparse(paramsout.(bresler_method).transform);
                end
                paramsout.(bresler_method).representation = sparse(paramsout.(bresler_method).representation);
    
                [msg, id] = lastwarn;
                if ~isempty(msg)
                    paramsout.(bresler_method).warning.msg = msg;
                    paramsout.(bresler_method).warning.id = id;
                end
    
                %%%%%%%%%%% Unstructured Conditioned Method %%%%%%%%%%%
    
                method = 'UnstructuredConditioned';
                lastwarn('');
                tic;
                [paramsout.(method).transform, paramsout.(method).representation, paramsout.(method).error] = UnstructuredConditioned(local_paramsin);
                paramsout.(method).time = toc;
    
                paramsout.(method).representation = sparse(paramsout.(method).representation);
    
                [msg, id] = lastwarn;
                if ~isempty(msg)
                    paramsout.(method).warning.msg = msg;
                    paramsout.(method).warning.id = id;
                end

    
                %%%%%%%%%%% Structured Conditioned Method %%%%%%%%%%%
    
                lastwarn('');
    
                sty_pct_lambda_search = NaN(1, max_iter);
                lambdas_lambda_search = NaN(1, max_iter);
                lambda_min = global_lambda_min;
                lambda_max = global_lambda_max;
                local_paramsin.relaxation = 0.01;

                lambda_mid = NaN;
                T_proposed = [];
                X_proposed = [];
                error_proposed = [];
                sty_pct = NaN;
                sty_vec = [];
                time_doubly_cond = 0;
    
                % secant method could be used instead
                for iter = 1:max_iter
                    lambda_mid = exp((log(lambda_min) + log(lambda_max)) / 2);
    
                    local_paramsin.lambda = lambda_mid;
                    
                    tic;
                    [T_proposed, X_proposed, error_proposed, sty_pct, sty_vec] = StructuredConditioned(local_paramsin);
                    time_doubly_cond = toc;

                    % 'restart' lambda search loop with a tighter relaxation parameter 
                    % if the sparsity percentage stagnates
                    if (iter > 1) && (abs(sty_pct - sty_pct_lambda_search(iter - 1)) <= 0.00001)
                        local_paramsin.relaxation = local_paramsin.relaxation / 10;

                        lambdas_lambda_search(iter) = local_paramsin.lambda;
                        sty_pct_lambda_search(iter) = sty_pct;

                        lambda_min = global_lambda_min;
                        lambda_max = global_lambda_max;

                        continue;
                    end
            
                    lambdas_lambda_search(iter) = local_paramsin.lambda;
                    sty_pct_lambda_search(iter) = sty_pct;
    
                    if abs(sty_pct - local_paramsin.T1) <= tol_sty_pct
                        break;
                    end
                    
                    if sty_pct > local_paramsin.T1
                        lambda_min = lambda_mid;
                    else
                        lambda_max = lambda_mid;
                    end
                
                end
    
                method = 'StructuredConditioned';
                
                paramsout.(method).transform = sparse(T_proposed);
                paramsout.(method).sparse_representation = sparse(X_proposed);
                paramsout.(method).error = error_proposed;
                paramsout.(method).lambda = lambda_mid;
                paramsout.(method).sty_pct = sty_pct;
                paramsout.(method).sty_vec = sty_vec;
                paramsout.(method).lambda_search_iterations = iter;
                paramsout.(method).sty_pct_lambda_search = sty_pct_lambda_search(1:iter);
                paramsout.(method).lambdas_lambda_search = lambdas_lambda_search(1:iter);
                paramsout.(method).time = time_doubly_cond;
            
                [msg, id] = lastwarn;
                if ~isempty(msg)
                    paramsout.(method).warning.msg = msg;
                    paramsout.(method).warning.id = id;
                end

                %%% save result
    
                result.paramsin.bresler_lambda = bresler_lambda;
                result.paramsin.bresler_method = bresler_method;
                result.paramsin.parameters_settings_idx = iter_param_setting;
                result.paramsin.rho = local_paramsin.rho;
                result.paramsin.tau = local_paramsin.tau;
                result.paramsin.l2_bresler = local_paramsin.l2_bresler;
                result.paramsin.l3_bresler = local_paramsin.l3_bresler;
                result.paramsout = paramsout;
                result.status = 'Success';

            catch ME
                result.ME = ME;

                result.paramsin.bresler_lambda = bresler_lambda;
                result.paramsin.parameters_settings_idx = iter_param_setting;
                if strcmp(ME.identifier, 'TRANSFORM:PoorConditioningBreslerMethod')
                    result.status = 'Skipped_BadCondition';
                else
                    result.status = 'Failed';
                end
            end

            results{iter_bresler_method} = result;
        end

        %%% save to disk
        file = sprintf('results_L%d_P%d.mat', iter_bresler_lambda, iter_param_setting);
        
        path = fullfile(folder, file);

        parsave_results(path, results);
    end
end

% parsave_results func is used by workers
function parsave_results(filepath, results)
    save(filepath, 'results');
end
