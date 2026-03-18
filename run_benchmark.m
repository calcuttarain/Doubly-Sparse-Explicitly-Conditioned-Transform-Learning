clear; clc; close all;
addpath('TLAlgorithms/');
addpath('utils/');

%%%%%%%%%%%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%

% T0_list = [6, 8, 10, 15];
% T1_list = [0.15, 0.20, 0.25, 0.30];
% 
% patch_size_list = [64, 121, 196, 256];
T0_list = [6];
T1_list = [20];

patch_size_list = [64, 256];

input_folders_datasets = {'data', 'DIV2K_valid_HR_data'};
train_sets = { {'Barbara', 'Couple', 'Lena', 'Cameraman'}, ...
               {'0801', '0802', '0803', '0804', '0805', '0806', '0807'} };
test_sets  = { {'Hill', 'Man', 'Baboon', 'MRI'}, ...
               {'0808', '0809', '0810'} };

[parameters_settings, datasets_by_n] = generate_settings(T0_list, T1_list, patch_size_list, input_folders_datasets, train_sets, test_sets);

bresler_methods = {'UnstructuredBresler', 'StructuredBresler', 'StructuredBreslerCF'};
% bresler_lambdas = linspace(2.1e-1, 2.1e-13, 2);
bresler_lambdas = [2.1e-6, 2.1e-9];
numiter = 1000;

% lambda search loop parameters for Doubly Sparse Conditioned Transform
max_iter = 30;                                     % maximum number of iterations
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

    for iter_param_setting = 1:length(parameters_settings)
        results = cell(1, length(bresler_methods));

        %%% extract parameters

        paramsin.T0 = parameters_settings(iter_param_setting).T0;
        paramsin.T1 = parameters_settings(iter_param_setting).T1;
        n = parameters_settings(iter_param_setting).patch_size;
        dataset_idx = parameters_settings(iter_param_setting).dataset_idx;

        paramsin.YH_train = datasets_by_n([datasets_by_n.n] == n).YH_train{dataset_idx};
        paramsin.YH_test = datasets_by_n([datasets_by_n.n] == n).YH_test{dataset_idx};

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

            try
                %%%%%%%%%%% Run Bresler Method %%%%%%%%%%%
    
                lastwarn('');
                tic;
                [paramsout.(bresler_method).transform, paramsout.(bresler_method).representation, paramsout.(bresler_method).error] = feval(bresler_method, paramsin);
                paramsout.(bresler_method).time = toc;
    
                paramsin.rho = cond(paramsout.(bresler_method).transform);
                paramsin.tau = norm(paramsout.(bresler_method).transform, 'fro');

                if ~strcmp(bresler_method, 'UnstructuredBresler')
                    paramsout.(bresler_method).transform = sparse(paramsout.(bresler_method).transform);
                end
                paramsout.(bresler_method).representation = sparse(paramsout.(bresler_method).representation);
    
                [msg, id] = lastwarn;
                if ~isempty(msg)
                    paramsout.(bresler_method).warning.msg = msg;
                    paramsout.(bresler_method).warning.id = id;
                end
    
                bresler_method
    
                %%%%%%%%%%% Unstructured Conditioned Method %%%%%%%%%%%
    
                method = 'UnstructuredConditioned';
                lastwarn('');
                tic;
                [paramsout.(method).transform, paramsout.(method).representation, paramsout.(method).error] = UnstructuredConditioned(paramsin);
                paramsout.(method).time = toc;
    
                paramsout.(method).representation = sparse(paramsout.(method).representation);
    
                [msg, id] = lastwarn;
                if ~isempty(msg)
                    paramsout.(method).warning.msg = msg;
                    paramsout.(method).warning.id = id;
                end
    
                method
                
    
                %%%%%%%%%%% Structured Conditioned Method %%%%%%%%%%%
    
                lastwarn('');
    
                sty_pct_lambda_search = NaN(1, max_iter);
                lambda_min = global_lambda_min;
                lambda_max = global_lambda_max;
                paramsin.relaxation = 0.01;
    
                for iter = 1:max_iter
                    lambda_mid = exp((log(lambda_min) + log(lambda_max)) / 2);
    
                    paramsin.lambda = lambda_mid;
                    
                    tic;
                    [T, X, error, sty_pct, sty_vec] = StructuredConditioned(paramsin);
                    time_doubly_cond = toc;

                    if (iter > 1) && (abs(sty_pct - sty_pct_lambda_search(iter - 1)) <= 0.00001)
                        paramsin.relaxation = paramsin.relaxation / 10
                    end
            
                    sty_pct_lambda_search(iter) = sty_pct;
    
                    sty_pct
                    
                    if abs(sty_pct - paramsin.T1) <= tol_sty_pct
                        break;
                    end
                    
                    if sty_pct > paramsin.T1
                        lambda_min = lambda_mid;
                    else
                        lambda_max = lambda_mid;
                    end
                
                end
    
                method = 'StructuredConditioned';
                
                paramsout.(method).transform = sparse(T);
                paramsout.(method).sparse_representation = sparse(X);
                paramsout.(method).error = error;
                paramsout.(method).lambda = lambda_mid;
                paramsout.(method).sty_pct = sty_pct;
                paramsout.(method).sty_vec = sty_vec;
                paramsout.(method).lambda_search_iterations = iter;
                paramsout.(method).sty_pct_lambda_search = sty_pct_lambda_search(1:iter);
                paramsout.(method).time = time_doubly_cond;
            
                [msg, id] = lastwarn;
                if ~isempty(msg)
                    paramsout.(method).warning.msg = msg;
                    paramsout.(method).warning.id = id;
                end
                method
                
                %%% save result
    
                result.paramsin.bresler_lambda = bresler_lambda;
                result.paramsin.bresler_method = bresler_method;
                result.paramsin.parameters_settings_idx = iter_param_setting;
                result.paramsin.rho = paramsin.rho;
                result.paramsin.tau = paramsin.tau;
                result.paramsin.l2_bresler = paramsin.l2_bresler;
                result.paramsin.l3_bresler = paramsin.l3_bresler;
                result.paramsout = paramsout;
                result.status = 'Success';

            catch ME
                result.ME = ME;

                result.paramsin.bresler_lambda = bresler_lambda;
                result.paramsin.parameters_settings_idx = iter_param_setting;
                result.status = 'Failed';
            end

            results{iter_bresler_method} = result;

            clear paramsout;
        end

        %%% save to disk
        file = sprintf('results_L%d_P%d.mat', iter_bresler_lambda, iter_param_setting);
        
        path = fullfile(folder, file);
        
        save(path, 'results');
    
        clear results;
    end
end
