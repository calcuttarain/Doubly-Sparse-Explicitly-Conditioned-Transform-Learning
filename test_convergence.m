warning('off','MATLAB:rankDeficientMatrix');
clear; clc; %close all;
rng(0);
addpath('TLAlgorithms/');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Set parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = 64;                                       % patch size 

T0 = 6;                                       % sparsity level for each representation

numiter = 1000;                               % Number of iterations for AM algorithm

W0 = kron(dctmtx(sqrt(n)), dctmtx(sqrt(n)));  % 2D DCT initialization, canonical transform factor                      

%lambda0 = 2.1e-7;                            % Bresler Method parameter
lambda0 = 2.1e-12;


% DoublySparseConditionedTL parameters
lambda = 10;                               % parameter for \ell_1 regularization term

homotopy_steps = 100;                          % number of steps with bigger \lambda

debias_start = numiter - 100;                 % iteration from which the sparse structure is fixed
clipping_eps = 10e-14;                         % post-projection clipping tolerance


% Bresler Doubly Sparse Transform parameters
T1 = round((0.25)*(n^2));                     % Bresler Doubly Sparse Transform sparsity percent
B0 = eye(n); 
mu = 2e-9; 
numg = 30; 
cbb = 0; stopcn = 0; stopth = 0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Data Loading and Preparation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load training set
barbara = struct2cell(load('img/Barbara/sigma5/I7.mat')); barbara = barbara{1}; 
couple = struct2cell(load('img/Couple/sigma5/I7.mat')); couple = couple{1}; 
lena = struct2cell(load('img/Lena/sigma5/I7.mat')); lena = lena{1}; 
cameraman = struct2cell(load('img/Cameraman/sigma5/I7.mat')); cameraman = cameraman{1}; 
hill = struct2cell(load('img/Hill/sigma5/I7.mat')); hill = hill{1}; 
man = struct2cell(load('img/Man/sigma5/I7.mat')); man = man{1}; 
mri = struct2cell(load('img/MRI/sigma5/I7.mat')); mri = mri{1}; 
baboon = struct2cell(load('img/Baboon/sigma5/I7.mat')); baboon = baboon{1};

% vectorize
[blocks_barbara] = my_im2col(barbara, [sqrt(n), sqrt(n)], sqrt(n));
[blocks_couple] = my_im2col(couple, [sqrt(n), sqrt(n)], sqrt(n));
[blocks_lena] = my_im2col(lena, [sqrt(n), sqrt(n)], sqrt(n));
[blocks_cameraman] = my_im2col(cameraman, [sqrt(n), sqrt(n)], sqrt(n));
[blocks_hill] = my_im2col(hill, [sqrt(n), sqrt(n)], sqrt(n));
[blocks_man] = my_im2col(man, [sqrt(n), sqrt(n)], sqrt(n));
[blocks_mri] = my_im2col(mri, [sqrt(n), sqrt(n)], sqrt(n));
[blocks_baboon] = my_im2col(baboon, [sqrt(n), sqrt(n)], sqrt(n));

% concatenate
[blocks_tr] = [blocks_barbara, blocks_couple, blocks_lena];
[blocks_te] = [blocks_hill, blocks_man];
% blocks_tr = blocks_tr(1:16, :);
% blocks_te = blocks_te(1:16, :);

% subtract the means
br_tr = mean(blocks_tr); br_te = mean(blocks_te);
TE_tr = blocks_tr - (ones(n, 1) * br_tr); TE_te = blocks_te - (ones(n, 1) * br_te);
YH_train = TE_tr; YH_test = TE_te;

% data in analytical transform domain 
YH2_train = W0 * YH_train; YH2_test = W0 * YH_test;

% set the sparsity levels
STY_tr = T0 * ones(1, size(YH_train, 2)); STY_te = T0 * ones(1, size(YH_test, 2)); 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Run Transforms %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DCT Method %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[s_tr]=sort(abs(YH2_train),'descend'); [s_te]=sort(abs(YH2_test),'descend'); 
X_tr = YH2_train.*(bsxfun(@ge,abs(YH2_train),s_tr(STY_tr))); X_te = YH2_test.*(bsxfun(@ge,abs(YH2_test),s_te(STY_te))); 

error_dct.m1.tr = ones(1, numiter) * norm(X_tr - YH2_train, 'fro'); error_dct.m1.te = ones(1, numiter) * norm(X_te - YH2_test, 'fro');
error_dct.m2.tr = ones(1, numiter) * (norm(X_tr - YH2_train, 'fro') / norm(YH2_train, 'fro')); error_dct.m2.te = ones(1, numiter) * (norm(X_te - YH2_test, 'fro') / norm(YH2_test, 'fro'));

errors(1) = error_dct;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Bresler Method %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

l2_bresler = lambda0 * (norm(YH_train, 'fro'))^2;
l3_bresler = l2_bresler;

[B_bresler, X_bresler, error_bresler]= TLclosedformmethod(W0, YH_train, YH_test, numiter, l2_bresler, l3_bresler, STY_tr, STY_te);
fprintf('Bresler Method Done\n');

errors(2) = error_bresler;

% set rho and tau based on Bresler Learnt Transform for Explicitly Conditioned Methods
rho = cond(B_bresler);
tau = norm(B_bresler, 'fro');


%%%%%%%%%%%%%%%%%%%%%%%%%% Structured Bresler Method %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%[W_bresler_doubly, X_bresler_doubly, ~, error_bresler_doubly] = BreslerDoublySparseTL(B0, YH2_train, YH2_test, numiter, l2_bresler, l3_bresler, T1, mu, numg, STY_tr, STY_te, cbb, stopcn, stopth);
% [W_bresler_doubly, X_bresler_doubly, error_bresler_doubly] = ClosedFormBreslerDoublySparse(B0, YH2_train, YH2_test, numiter, l2_bresler, l3_bresler, T1, STY_tr, STY_te, cbb);
% fprintf('Structured Bresler Method Done\n');
%
% errors(2) = error_bresler_doubly;

% set rho and tau based on Bresler Learnt Transform for Explicitly Conditioned Methods
% rho = cond(W_bresler_doubly);
% tau = norm(W_bresler_doubly, 'fro');


%%%%%%%%%%%%%%%%%%%%%%% Unstructured Conditioned Method %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[W_cond, X_cond, error_cond] = ConditionedTransformLearning(W0, YH_train, YH_test, numiter, STY_tr, STY_te, rho, tau);
fprintf('Unstructured Conditioned Method Done\n');

errors(3) = error_cond;


%%%%%%%%%%%%%%%%%%%%%%% Structured Conditioned Method %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[T_doubly_cond, X_doubly_cond, error_doubly_cond, sty_pct, sty_vec] = DoublySparseConditionedTL(W0, YH2_train, YH2_test, numiter, STY_tr, STY_te, rho, tau, lambda, homotopy_steps, debias_start, clipping_eps);
fprintf('Structured Conditioned Method with %.2f%% Sparsity Done\n', sty_pct);

errors(4) = error_doubly_cond;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Save Results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

labels = {'DCT', "Bresler Method", "Conditioned Unstructured", "Conditioned Structured"};
errors_train = cell(1, 4); errors2_train = cell(1, 4);
errors_test = cell(1, 4); errors2_test = cell(1, 4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Train / Test Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:4
    errors_train{i} = errors(i).m1.tr; errors_test{i} = errors(i).m1.te;
    errors2_train{i} = errors(i).m2.tr; errors2_test{i} = errors(i).m2.te;
end

title_text = 'Train Data';
plot_convergence(numiter, errors_train, errors2_train, labels, rho, tau, sty_pct, T0, title_text);

% title_text = 'Test Data';
% plot_convergence(numiter, errors_test, errors2_test, labels, rho, tau, sty_pct, T0, title_text);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Sparsity %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%plot_sparsity(T_doubly_cond, numiter, errors(4).m1.tr, sty_vec, rho, tau, lambda, T0);