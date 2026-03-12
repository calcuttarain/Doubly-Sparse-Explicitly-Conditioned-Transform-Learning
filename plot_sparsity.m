function plot_sparsity(T_doubly_cond, numiter, errors, sty_vec, kappa, tau, lambda, T)
    f = figure('Position', [100, 100, 1000, 400]);
    f.Theme = 'light';

    subplot(1, 2, 1);
    plot(1:numiter, errors, 'b-', 'LineWidth', 1.5);
    set(gca, 'YScale', 'log');
    xlabel('Iterations', 'Interpreter', 'latex');
    ylabel('$\|X - W Y\|_F$', 'Interpreter', 'latex');
    grid on;

    subplot(1, 2, 2);
    plot(1:numiter, sty_vec, 'r-', 'LineWidth', 1.5);
    xlabel('Iterations', 'Interpreter', 'latex');
    ylabel('$T$ Learnt Transform Sparsity Percent', 'Interpreter', 'latex');
    grid on;

    sgtitle(['$ r = ', num2str(T), ',\ \rho = ', num2str(kappa, '%.2f'), ...
             ',\ \tau = ', num2str(tau, '%.2f'), ',\ \lambda = ', num2str(lambda, '%.2e'), '$'], ...
             'Interpreter', 'latex');

    % exportgraphics(gcf, 'plots/sparsity_pct.pdf', 'ContentType', 'vector');

    f = figure;
    f.Theme = 'light';
    
    spy(sparse(T_doubly_cond));
    xlabel('');

    sgtitle('$T$ Learnt Transform Sparsity Structure', 'Interpreter', 'latex');

    % exportgraphics(gcf, 'plots/sparsity_struct.pdf', 'ContentType', 'vector');
end
