function plot_convergence(numiter, errors, errors2, labels, kappa, tau, sty_pct, T, title_text)
    f = figure('Position', [100, 100, 1000, 400]);
    f.Theme = 'light';

    colors = {'k-', 'c-', 'b-', 'r-'}; 
    n_methods = length(labels);

    % regular error
    subplot(1, 2, 1);
    hold on;
    set(gca, 'YScale', 'log');
    for i = 1:n_methods
        plot(1:numiter, errors{i}, colors{i}, 'LineWidth', 1.5);
    end
    xlabel('Iterations', 'Interpreter', 'latex');
    ylabel('$\|X - W Y\|_F$', 'Interpreter', 'latex');
    legend(labels, 'Location', 'northeast');
    grid on;

    % normalized error
    subplot(1, 2, 2);
    hold on;
    set(gca, 'YScale', 'log');
    for i = 1:n_methods
        plot(1:numiter, errors2{i}, colors{i}, 'LineWidth', 1.5);
    end
    xlabel('Iterations', 'Interpreter', 'latex');
    ylabel('$\|X - W Y\|_F / \|W Y\|_F$', 'Interpreter', 'latex');
    legend(labels, 'Location', 'northeast');
    grid on;

    % include r-sparsity for each representation, sparsity percent for doubly sparse conditioned transform and the parameters rho and tau in the title
    sgtitle({ ['\bf ', title_text], ...
        ['$r = ', num2str(T), ...
        ',\ \rho = ', num2str(kappa,'%.2f'), ...
        ',\ \tau = ', num2str(tau,'%.2f'), ...
        ',\ sparsity = ', num2str(sty_pct,'%.2f'), '\%$']}, ...
        'Interpreter','latex');
    
    exportgraphics(gcf, 'convergence.pdf', 'ContentType', 'vector');

end
