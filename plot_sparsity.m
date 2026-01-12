function plot_sparsity(numiter, errors, sty_vec, kappa, tau, lambda, T)
    figure('Position', [100, 100, 1000, 400]);

    subplot(1, 2, 1);
    plot(1:numiter, errors, 'b-', 'LineWidth', 1.5);
    set(gca, 'YScale', 'log');
    xlabel('Iterations', 'Interpreter', 'latex');
    ylabel('$\|X - W Y\|_F$', 'Interpreter', 'latex');
    grid on;

    subplot(1, 2, 2);
    plot(1:numiter, sty_vec, 'r-', 'LineWidth', 1.5);
    xlabel('Iterations', 'Interpreter', 'latex');
    ylabel('$T$ Transform sparsity percent', 'Interpreter', 'latex');
    grid on;

    sgtitle(['$ r = ', num2str(T), ',\ \rho = ', num2str(kappa, '%.2f'), ...
             ',\ \tau = ', num2str(tau, '%.2f'), ',\ \lambda = ', num2str(lambda, '%.4f'), '$'], ...
             'Interpreter', 'latex');
end
