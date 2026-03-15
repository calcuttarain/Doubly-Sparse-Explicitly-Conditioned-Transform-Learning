clear; clc; close all;
addpath('TLAlgorithms/');
addpath('utils/');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Set parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bresler_methods = {'unstructured_bresler', 'structured_bresler', 'structured_bresler_cf'};
conditioned_methods = {'unstructured_conditioned', 'structured_conditioned'};

methods = {'structured_bresler_cf', 'unstructured_conditioned', 'structured_conditioned'};

% train / test sets
train_set = {'Barbara', 'Couple', 'Lena'};
test_set  = {'Hill', 'Man'};

n = 64;                                            % patch size 

T0 = 6;                                            % sparsity level for each representation

numiter = 300;                                     % number of iterations for AM algorithm

W0 = kron(dctmtx(sqrt(n)), dctmtx(sqrt(n)));       % 2D DCT initialization, canonical transform factor

lambda0 = 2.1e-7;                                  % Bresler Method parameter

% Bresler Doubly Sparse Transform parameters
T1 = round((0.25)*(n^2));                          % Bresler Doubly Sparse Transform sparsity
B0 = eye(n);                                       % Bresler Doubly Sparse Transform initialization
mu = 2e-9;
numg = 30;
cbb = 0; stopcn = 0; stopth = 0;

% lambda search loop parameters
max_iter = 30;                                     % maximum number of iterations
tol_sty_pct = 1;                                   % sparsity percent tolerance
lambda_min = 1e-10;                                % left endpoint of search interval
lambda_max = 1e3;                                  % right endpoint of search interval

% save
paramsin.n = n;
paramsin.T0 = T0;
paramsin.numiter = numiter;
paramsin.W0 = W0;
paramsin.lambda0 = lambda0;
paramsin.T1 = T1;
paramsin.B0 = B0; 
paramsin.mu = mu; 
paramsin.numg = numg; 
paramsin.cbb = cbb; 
paramsin.stopcn = stopcn; 
paramsin.stopth = stopth;
paramsin.max_iter = max_iter;
paramsin.lambda_min = lambda_min;
paramsin.lambda_max = lambda_max;
paramsin.tol_sty_pct = tol_sty_pct;
paramsin.train_set = train_set;
paramsin.test_set = test_set;
paramsin.rng_state = rng;
paramsin.methods = methods;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Data Loading and Preparation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

blocks_tr = [];
blocks_te = [];

for i = 1:length(train_set)
    img = train_set{i};

    filepath = fullfile('data', [img, '.mat']);

    data = struct2cell(load(filepath));
    img = data{1};

    blocks_tr = [blocks_tr, my_im2col(img, [sqrt(n), sqrt(n)], sqrt(n))];
end

for i = 1:length(test_set)
    img = test_set{i};

    filepath = fullfile('data', [img, '.mat']);

    data = struct2cell(load(filepath));
    img = data{1};

    % concatenate
    blocks_te = [blocks_te, my_im2col(img, [sqrt(n), sqrt(n)], sqrt(n))];
end

% subtract the means
br_tr = mean(blocks_tr); br_te = mean(blocks_te);
TE_tr = blocks_tr - (ones(n, 1) * br_tr); TE_te = blocks_te - (ones(n, 1) * br_te);
YH_train = TE_tr; YH_test = TE_te;

% data in analytical transform domain
YH2_train = W0 * YH_train; YH2_test = W0 * YH_test;

% set the sparsity levels
STY_tr = T0 * ones(1, size(YH_train, 2)); STY_te = T0 * ones(1, size(YH_test, 2)); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Run Transforms %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
%%%%%%%%%%%%%%%%%%%%%%%%%% Structured Bresler Method %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    l2_bresler = lambda0 * (norm(YH_train, 'fro'))^2;
    l3_bresler = l2_bresler;
    
    %[W_bresler_doubly, X_bresler_doubly, ~, error_bresler_doubly] = StructuredBresler(B0, YH2_train, YH2_test, numiter, l2_bresler, l3_bresler, T1, mu, numg, STY_tr, STY_te, cbb, stopcn, stopth);
    lastwarn('');
    tic;
    [W_bresler_doubly_cf, X_bresler_doubly_cf, error_bresler_doubly_cf] = StructuredBreslerCF(B0, YH2_train, YH2_test, numiter, l2_bresler, l3_bresler, T1, STY_tr, STY_te, cbb);
    time_bresler_doubly_cf = toc;
    
    % set rho and tau based on Bresler Learnt Transform for Explicitly Conditioned Methods
    rho = cond(W_bresler_doubly_cf);
    tau = norm(W_bresler_doubly_cf, 'fro');
    
    % set target sparsity percent for Structured Explicitly Conditioned Method
    total = numel(W_bresler_doubly_cf);           
    curr_sty = nnz(W_bresler_doubly_cf(:) ~= 0);
    target_sty_pct = 100 * curr_sty / total;
    
    % save
    paramsin.l2_bresler = l2_bresler;
    paramsin.l3_bresler = l3_bresler;
    paramsin.rho = rho;
    paramsin.tau = tau;
    
    paramsout.(methods{1}).transform = sparse(W_bresler_doubly_cf);
    paramsout.(methods{1}).sparse_representation = sparse(X_bresler_doubly_cf);
    paramsout.(methods{1}).error = error_bresler_doubly_cf;
    paramsout.(methods{1}).time = time_bresler_doubly_cf;

    [msg, id] = lastwarn;
    if ~isempty(msg)
        paramsout.(methods{1}).warning.msg = msg;
        paramsout.(methods{1}).warning.id = id;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%% Unstructured Conditioned Method %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    lastwarn('');
    tic;
    [W_cond, X_cond, error_cond] = UnstructuredConditioned(W0, YH_train, YH_test, numiter, STY_tr, STY_te, rho, tau);
    time_cond = toc;
    
    paramsout.(methods{2}).transform = W_cond;
    paramsout.(methods{2}).sparse_representation = sparse(X_cond);
    paramsout.(methods{2}).error = error_cond;
    paramsout.(methods{2}).time = time_cond;

    [msg, id] = lastwarn;
    if ~isempty(msg)
        paramsout.(methods{2}).warning.msg = msg;
        paramsout.(methods{2}).warning.id = id;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%% Structured Conditioned Method %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    lastwarn('');
    for iter = 1:max_iter
        lambda_mid = exp((log(lambda_min) + log(lambda_max)) / 2);
        
        tic;
        [T_doubly_cond, X_doubly_cond, error_doubly_cond, sty_pct, sty_vec] = StructuredConditioned(W0, YH2_train, YH2_test, numiter, STY_tr, STY_te, rho, tau, lambda_mid);
        time_doubly_cond = toc;
        
        if abs(sty_pct - target_sty_pct) <= tol_sty_pct
            best_lambda = lambda_mid;
            lambda_search_iterations = iter;
            break;
        end
        
        if sty_pct > target_sty_pct
            lambda_min = lambda_mid;
        else
            lambda_max = lambda_mid;
        end
    
    end
    
    paramsout.(methods{3}).transform = sparse(T_doubly_cond);
    paramsout.(methods{3}).sparse_representation = sparse(X_doubly_cond);
    paramsout.(methods{3}).error = error_doubly_cond;
    paramsout.(methods{3}).lambda = best_lambda;
    paramsout.(methods{3}).sty_pct = sty_pct;
    paramsout.(methods{3}).sty_vec = sty_vec;
    paramsout.(methods{3}).lambda_search_iterations = lambda_search_iterations;
    paramsout.(methods{3}).time = time_doubly_cond;

    [msg, id] = lastwarn;
    if ~isempty(msg)
        paramsout.(methods{3}).warning.msg = msg;
        paramsout.(methods{3}).warning.id = id;
    end

    paramsout.status = 'Success';

    %%% save to disk

    folder = 'results';
    
    if ~exist(folder, 'dir')
        mkdir(folder);
    end
    
    timestamp = datestr(now, 'yyyy-mm-dd_HHMMSS');
    file = sprintf('%s.mat', timestamp);
    
    path = fullfile(folder, file);
    
    save(path, 'paramsin', 'paramsout');

    %%% clear output
    clear paramsout;

catch ME
    % beep;

    paramsout.status = 'Failed';

    folder = 'debug_logs';
    
    if ~exist(folder, 'dir')
        mkdir(folder);
    end

    timestamp = datestr(now, 'yyyy-mm-dd_HHMMSS');

    path = fullfile('debug_logs', ['CRASH_', '_', timestamp, '.mat']);

    save(path, 'ME', 'paramsin', 'paramsout');

end
