function [T, XT, error, sty_pct, sty_vec] = DoublySparseConditionedTL(T, Y, Y_test, numiter, STY, STY_te, rho, tau, lambda, homotopy_steps, debias_start, clipping_eps)
    addpath('TLAlgorithms/DoublyConditionedTLRoutines/');

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
    T_curr = eye(n);
    t_curr = 1;

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
            T_ant_Y = T_curr;
            T_ant = T_curr;
            T_curr = T_curr - alpha * D_curr;

            % generate lambdas list
            [U, S, V] = svd(T_curr);
            sigmas = get_spectrum_doubly(S, rho, tau);
            T_temp = U * diag(sigmas) * V';

            %T_curr = T_temp;

            lambda_start = max(abs(T_temp(:))) / alpha * 10e-2;
            
            decreasing_lambdas = logspace(log10(lambda_start), log10(lambda), homotopy_steps);

            continue;
        end

        % nesterov acceleration
        t_next = (1 + sqrt(1 + 4 * t_curr^2)) / 2;
        beta = (t_curr - 1) / t_next;

        Y_k = T_curr + beta * (T_curr - T_ant);

        t_curr = t_next;

        % compute the gradient at the current step
        D_ant = D_curr;
        D_curr = (Y_k * Y - X) * Y';

        % get alpha
        alpha = get_alpha(Y_k, D_curr, T_ant_Y, D_ant, Y, "exact_line_search");
        T_ant_Y = Y_k;

        % compute the transform
        T_ant = T_curr;
        T_curr = Y_k - alpha * D_curr;
        %T_curr = X * Y_pinv; % exact solution

        % soft-thersholding
        if i <= debias_start
            if i - 1 <= length(decreasing_lambdas)
                chosen_lambda = decreasing_lambdas(i - 1);
            else
                chosen_lambda = lambda;
            end
            T_curr = sign(T_curr) .* max(abs(T_curr) - alpha * chosen_lambda, 0);
        end

        % projection onto feasible space
        [U, S, V] = svd(T_curr);
        sigmas = get_spectrum_doubly(S, rho, tau);
        T_curr = U * diag(sigmas) * V';

        % clipping / hard-thresholding
        if i < debias_start
            indices = abs(T_curr) <= clipping_eps;
            T_curr(indices) = 0;
        elseif i == debias_start
            indices = abs(T_curr) <= clipping_eps;
            T_curr(indices) = 0;
            support_mask = (T_curr ~= 0);
        else
            T_curr = T_curr .* support_mask;
        end
        

        T = T_curr;

        total = numel(T);           
        curr_sty = nnz(T(:) ~= 0);
        sty_vec(i - 1) = 100 * curr_sty / total;

        error.m1.tr(i - 1) = norm(X - T_curr * Y, 'fro');
        error.m2.tr(i - 1) = norm(X - T_curr * Y, 'fro') / norm(T_curr * Y, 'fro');

        error.m1.te(i - 1) = norm(X_test - T_curr * Y_test, 'fro');
        error.m2.te(i - 1) = norm(X_test - T_curr * Y_test, 'fro') / norm(T_curr * Y_test, 'fro');
    end

    total = numel(T);           
    curr_sty = nnz(T(:) ~= 0);  
    sty_pct = 100 * curr_sty / total;   

    XT(:, ix) = X;
end
