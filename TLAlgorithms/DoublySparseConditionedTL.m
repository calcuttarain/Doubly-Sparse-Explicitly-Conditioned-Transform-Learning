function [T, XT, error, error2, sty_pct] = DoublySparseConditionedTL(T, Y, numiter, STY, rho, tau, lambda)
    addpath('TLAlgorithms/DoublyConditionedTLRoutines/');

    rng(0);

    [K, n] = size(T); 
    XT = zeros(K, size(Y, 2)); 

    ix = find(STY > 0);
    q = Y(:, ix);
    STY = STY(:, ix); 
    N = size(q, 2); 

    ez = K * (0:(N-1)); 
    STY = STY + ez; 
    Y = Y(:, ix);

    error = zeros(1, numiter);
    error2 = zeros(1, numiter);

    D_ant = zeros(n, n);
    T_ant = zeros(n, n);
    D_curr = zeros(n, n);
    T_curr = T;
    t_curr = 1;
    t_ant = 1;

    for i = 1:numiter + 1
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Sparse Representation Update Step %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        X1 = T_curr * q;
        [s] = sort(abs(X1), 'descend');

        X = X1 .* (bsxfun(@ge, abs(X1), s(STY)));


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
        T_curr = sign(T_curr) .* max(abs(T_curr) - lambda, 0);

        % projection onto feasible space
        T_curr = (T_curr + T_curr') / 2; % symmetry

        [Q, L] = eig(T_curr);
        lambdas = get_spectrum_doubly(L, rho, tau);
        T_curr = Q * diag(lambdas) * Q';

        T = T_curr;

        error(i - 1) = norm(X - T_curr * Y, 'fro');
        error2(i - 1) = norm(X - T_curr * Y, 'fro') / norm(T_curr * Y, 'fro');
    end

    total = numel(T);           
    num_zero = nnz(T(:) == 0);  
    sty_pct = 100 * num_zero / total;   

    XT(:, ix) = X;
end
