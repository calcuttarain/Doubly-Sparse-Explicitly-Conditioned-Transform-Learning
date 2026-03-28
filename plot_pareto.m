clear; clc; close all;

folder_path = 'results_2026-03-28_1931'; 
file_name = 'results_1.mat';
full_path = fullfile(folder_path, file_name);

load(full_path, 'setting', 'results');

n_pts = length(results);
err_b_tr = zeros(1, n_pts);
err_b_te = zeros(1, n_pts);
sty_b = zeros(1, n_pts);
err_sc_tr = zeros(1, n_pts);
err_sc_te = zeros(1, n_pts);
sty_sc = zeros(1, n_pts);

for i = 2:n_pts
    err_b_tr(i) = results(i).bresler_method.error.m1.tr(end);
    err_b_te(i) = results(i).bresler_method.error.m1.te(end);
    sty_b(i) = results(i).bresler_method.sty_pct;
    
    err_sc_tr(i) = results(i).structured_conditioned.error.m1.tr(end);
    err_sc_te(i) = results(i).structured_conditioned.error.m1.te(end);
    sty_sc(i) = results(i).structured_conditioned.sty_pct;
end

% [sty_b, idx_b] = sort(sty_b, 'ascend');
% err_b_tr = err_b_tr(idx_b);
% err_b_te = err_b_te(idx_b);
% 
% [sty_sc, idx_sc] = sort(sty_sc, 'ascend');
% err_sc_tr = err_sc_tr(idx_sc);
% err_sc_te = err_sc_te(idx_sc);

plot_test(err_b_tr, err_sc_tr, err_b_te, err_sc_te, sty_b, sty_sc, folder_path, file_name);

function plot_test(err_b_tr, err_sc_tr, err_b_te, err_sc_te, sty_b, sty_sc, folder, file_name)
    f = figure('Position', [100, 100, 1000, 400], 'Visible', 'on');
    f.Theme = 'light';
    
    subplot(1, 2, 1);
    hold on;
    plot(err_b_tr, sty_b, 'k-o', 'LineWidth', 1.5);
    plot(err_sc_tr, sty_sc, 'r-s', 'LineWidth', 1.5);

    set(gca, 'XScale', 'log');

    xlabel('$\|X - W Y_{train}\|_F$', 'Interpreter', 'latex');
    ylabel('Sparsity (\%)', 'Interpreter', 'latex');
    title('Training Pareto Front', 'Interpreter', 'latex');
    legend({'Bresler', 'Structured Conditioned'}, 'Location', 'northeast');
    grid on;
    
    subplot(1, 2, 2);
    hold on;
    plot(err_b_te, sty_b, 'k-o', 'LineWidth', 1.5);
    plot(err_sc_te, sty_sc, 'r-s', 'LineWidth', 1.5);

    set(gca, 'XScale', 'log');

    xlabel('$\|X - W Y_{test}\|_F$', 'Interpreter', 'latex');
    ylabel('Sparsity (\%)', 'Interpreter', 'latex');
    title('Testing Pareto Front', 'Interpreter', 'latex');
    legend({'Bresler', 'Structured Conditioned'}, 'Location', 'northeast');
    grid on;
    
    sgtitle('\bf Pareto Front: Sparsity vs. Error', 'Interpreter', 'latex');
    
    [~, name, ~] = fileparts(file_name);
    path = fullfile(folder, string(name) + "_pareto.pdf");
    exportgraphics(gcf, path, 'ContentType', 'vector');
end