warning('off','MATLAB:rankDeficientMatrix');
clear; clc;
rng(0);
addpath('TLAlgorithms/');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Set parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = 16;                                       % patch size 

T0 = 6;                                       % sparsity level for each representation

numiter = 2000;                               % Number of iterations for AM algorithm

W0 = kron(dctmtx(sqrt(n)), dctmtx(sqrt(n)));  % 2D DCT initialization, canonical transform factor

lambda = 0.00001;                             % penalization parameter for \ell_1 norm on T

lambda0 = 1.1e-8;                             % Bresler Method parameter

% Bresler Doubly Sparse Transform parameters
% T1 = round((0.25)*(n^2));                     % Bresler Doubly Sparse Transform sparsity percent
% B0 = eye(n); 
% mu = 2e-9; 
% numg = 30; 
% cbb = 0; stopcn = 0; stopth = 0;


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
[blocks] = [blocks_barbara, blocks_couple, blocks_lena];
blocks = blocks(1:16, :);

% subtract the means
br = mean(blocks);
TE = blocks - (ones(n, 1) * br);
YH = TE; 

% data in analytical transform domain 
YH2 = W0 * YH;

% set the sparsity levels
STY = T0 * ones(1, size(YH, 2)); 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Run Transforms %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DCT Method %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[s]=sort(abs(YH2),'descend'); 
X = YH2.*(bsxfun(@ge,abs(YH2),s(STY))); 

error_dct = ones(1, numiter) * norm(X - YH2, 'fro');
error2_dct = ones(1, numiter) * (norm(X - YH2, 'fro') / norm(YH2, 'fro'));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Bresler Method %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
l2_bresler = lambda0 * (norm(YH, 'fro'))^2;
l3_bresler = l2_bresler;
[B_bresler, X_bresler, error_bresler,error2_bresler]= TLclosedformmethod(W0, YH, numiter, l2_bresler, l3_bresler, STY);
fprintf('Bresler Method Done\n');

% set rho and tau based on Bresler Learnt Transform for Explicitly Conditioned Methods
rho = cond(B_bresler);
tau = norm(B_bresler,'fro');


%%%%%%%%%%%%%%%%%%%%%%%%%% Structured Bresler Method %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% [W_bresler_doubly, X_bresler_doubly, error_bresler_doubly, error2_bresler_doubly] = BreslerDoublySparseTL(B0, YH2, numiter, l2_bresler, l3_bresler, T1, mu, numg, STY, cbb, stopcn, stopth);
% fprintf('Structured Bresler Method Done\n');


%%%%%%%%%%%%%%%%%%%%%%% Unstructured Conditioned Method %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[W_cond, X_cond, error_cond, error2_cond] = ConditionedTransformLearning(W0, YH, numiter, STY, rho, tau);
fprintf('Unstructured Conditioned Method Done\n');


%%%%%%%%%%%%%%%%%%%%%%% Structured Conditioned Method %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[T_doubly_cond, X_doubly_cond, error_doubly_cond, error2_doubly_cond, sty_pct, sty_vec] = DoublySparseConditionedTL(W0, YH2 ,numiter, STY, rho, tau, lambda);
fprintf('Structured Conditioned Method with %.2f%% Sparsity Done\n', sty_pct);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Save Results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

errors = {error_dct, error_bresler, error_cond, error_doubly_cond};
errors2 = {error2_dct, error2_bresler, error2_cond, error2_doubly_cond};
labels = {'DCT', "Bresler Method", "Conditioned Unstructured", "Conditioned Structured"};

plot_convergence(numiter, errors, errors2, labels, rho, tau, T0);
plot_sparsity(T_doubly_cond, numiter, error_doubly_cond, sty_vec, rho, tau, lambda, T0);
