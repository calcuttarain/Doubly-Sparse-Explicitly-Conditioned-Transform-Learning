function [solution] = get_spectrum_doubly(L, rho, tau)
    d = diag(L);

    %%% positive eigenvalues
    eps_min = 1e-250;
    d(d < eps_min) = eps_min;

    [d_sorted, d_indices] = sort(d, 'ascend');

    %%% conditioning constraint
    [k1, k2, l_star] = kappa_search_doubly(d_sorted, rho);

    if l_star ~= 0
        solution = d_sorted;
        solution(1:k1) = l_star;
        solution(k2:end) = rho * l_star;
    else
        solution = d_sorted;
    end

    %%% frobenius norm constraint
    solution = solution * (tau / norm(solution, 2));

    new_indices(d_indices) = 1:length(solution);
    solution = solution(new_indices);
