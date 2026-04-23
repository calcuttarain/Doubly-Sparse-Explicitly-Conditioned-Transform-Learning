function [T, XT, error, sty_pct, sty_vec] = StructuredConditioned(T, Y, Y_test, numiter, STY, STY_te, rho, tau, lambda, homotopy_steps)
    addpath('TLAlgorithms/StructuredConditionedTLRoutines/');

    rng(0);

    [K, n] = size(T); 
    XT = zeros(K, size(Y, 2)); 

    ix = find(STY > 0); ix_te = find(STY_te > 0);
    q = Y(:, ix); q_te = Y_test(:, ix_te);
    STY = STY(:, ix); STY_te = STY_te(:, ix_te);
    N = size(q, 2); N_te = size(q_te, 2); 

    ez = K * (0:(N-1)); ez_te = K * (0:(N_te-1)); 
    STY = STY + ez; STY_te = STY_te + ez_te;
    Y = Y(:, ix); Y_test = Y_test(:, ix_te);

    error.m1.tr = zeros(1, numiter); error.m1.te = zeros(1, numiter);
    error.m2.tr = zeros(1, numiter); error.m2.te = zeros(1, numiter);
    sty_vec = zeros(1, numiter);

    D_ant = zeros(n, n);
    T_ant = zeros(n, n);
    D_curr = zeros(n, n);
    T_curr = T;
    t_curr = 1;
    t_ant = 1;

    for i = 1:numiter + 1
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Sparse Representation Update Step %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        X1 = T_curr * q; X1_te = T_curr * q_te;
        [s] = sort(abs(X1), 'descend'); [s_te] = sort(abs(X1_te), 'descend');

        X = X1 .* (bsxfun(@ge, abs(X1), s(STY))); X_test = X1_te .* (bsxfun(@ge, abs(X1_te), s_te(STY_te)));


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Transform Update Step %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % initialise T
        if i == 1
            if rho <= 2
                T_curr = Y*X';
            else
                T_curr = (Y'\X')';
            end

            % fix non-singularity
            [U, Sigma, V] = svd(T_curr, 'econ');

            sigmas = diag(Sigma);
            tol = 1e-12 * sigmas(1);
            sigmas = max(sigmas, tol);

            T_curr = U * diag(sigmas) * V';

            % one unconstrained gradient step
            D_curr = (T_curr * Y - X) * Y';

            % get alpha
            alpha = get_alpha(T_curr, D_curr, T_ant, D_ant, Y, "exact_line_search");

            % compute the transform
            T_ant = T_curr;
            T_curr = T_curr - alpha * D_curr;

            % generate lambdas list

            T_temp = (T_curr + T_curr') / 2; 

            [Q, L] = eig(T_temp);
            lambdas = get_spectrum_doubly(L, rho, tau);
            T_temp = Q * diag(lambdas) * Q';

            lambda_start = max(abs(T_temp(:))) / alpha;
            
            decreasing_lambdas = logspace(log10(lambda_start), log10(lambda), homotopy_steps);

            continue;
        end

        % nesterov acceleration
        beta = (t_ant - 1) / t_curr;

        T_curr = T_curr + beta * (T_curr - T_ant);

        t_ant = t_curr;
        t_curr = (1 + sqrt(1 + 4 * t_ant^2)) / 2;

        % compute the gradient at the current step
        D_ant = D_curr;
        D_curr = (T_curr * Y - X) * Y';

        % get alpha
        alpha = get_alpha(T_curr, D_curr, T_ant, D_ant, Y, "barzilai_borwein");

        % compute the transform
        T_ant = T_curr;
        T_curr = T_curr - alpha * D_curr;

        % soft-thersholding
        if i - 1 <= length(decreasing_lambdas)
            chosen_lambda = decreasing_lambdas(i - 1);
        else
            chosen_lambda = lambda;
        end

        T_curr = sign(T_curr) .* max(abs(T_curr) - alpha * chosen_lambda, 0);

        % projection onto feasible space
        T_curr = (T_curr + T_curr') / 2; % symmetry

        [Q, L] = eig(T_curr);
        lambdas = get_spectrum_doubly(L, rho, tau);
        T_curr = Q * diag(lambdas) * Q';

        % clipping
        indices = abs(T_curr) <= 10e-7;
        T_curr(indices) = 0;

        T = T_curr;

        total = numel(T);           
        curr_sty = nnz(T(:) ~= 0);
        sty_vec(i - 1) = 100 * curr_sty / total;

        TY = T * Y;
        TY_test = T * Y_test;

        error.m1.tr(i - 1) = norm(X - TY, 'fro');
        error.m2.tr(i - 1) = error.m1.tr(i - 1) / norm(TY, 'fro');

        error.m1.te(i - 1) = norm(X_test - TY_test, 'fro');
        error.m2.te(i - 1) = error.m1.te(i - 1) / norm(TY_test, 'fro');
    end

    total = numel(T);           
    curr_sty = nnz(T(:) ~= 0);  
    sty_pct = 100 * curr_sty / total;   

    XT(:, ix) = X;
end
