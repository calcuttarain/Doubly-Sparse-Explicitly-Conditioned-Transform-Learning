clear; clc; close all;

chosen_timestamp = "2026-04-21_1231";

input_folder = "results_" + chosen_timestamp; 

target_dir = fullfile("plots_" + chosen_timestamp, 'pareto');

if ~exist(char(target_dir), "dir")
    mkdir(char(target_dir));
end

files = dir(fullfile(input_folder, "results_*"));

for k = 1:length(files)
    current_file = files(k).name;
    full_path = fullfile(input_folder, current_file);
    
    load(full_path, 'setting', 'results');
    
    n_pts = length(results);
    err_b_tr = zeros(1, n_pts);
    err_b_te = zeros(1, n_pts);
    sty_b = zeros(1, n_pts);
    err_sc_tr = zeros(1, n_pts);
    err_sc_te = zeros(1, n_pts);
    sty_sc = zeros(1, n_pts);
    
    for i = 1:n_pts
        err_b_tr(i) = results(i).bresler_method.error.m1.tr(end);
        err_b_te(i) = results(i).bresler_method.error.m1.te(end);
        sty_b(i) = results(i).bresler_method.sty_pct;
        
        err_sc_tr(i) = results(i).structured_conditioned.error.m1.tr(end);
        err_sc_te(i) = results(i).structured_conditioned.error.m1.te(end);
        sty_sc(i) = results(i).structured_conditioned.sty_pct;
    end
    
    f = figure('Position', [100, 100, 1000, 400], 'Visible', 'off');
    f.Theme = 'light';
    
    colors = lines(n_pts); 
    
    % training plot
    subplot(1, 2, 1);
    hold on;
    for i = 1:n_pts
        plot([err_b_tr(i), err_sc_tr(i)], [sty_b(i), sty_sc(i)], '-', 'Color', colors(i,:), 'LineWidth', 1);
    end
    plot(err_b_tr, sty_b, 'ko', 'MarkerFaceColor', 'k', 'DisplayName', 'Bresler');
    plot(err_sc_tr, sty_sc, 'rs', 'MarkerFaceColor', 'r', 'DisplayName', 'Structured Conditioned');
    
    set(gca, 'XScale', 'log');
    xlabel('$\|X - W Y_{train}\|_F$', 'Interpreter', 'latex');
    ylabel('Sparsity (\%)', 'Interpreter', 'latex');
    title('Training Pareto Front', 'Interpreter', 'latex');
    grid on;
    
    % testing plot
    subplot(1, 2, 2);
    hold on;
    for i = 1:n_pts
        plot([err_b_te(i), err_sc_te(i)], [sty_b(i), sty_sc(i)], '-', 'Color', colors(i,:), 'LineWidth', 1);
    end
    plot(err_b_te, sty_b, 'ko', 'MarkerFaceColor', 'k');
    plot(err_sc_te, sty_sc, 'rs', 'MarkerFaceColor', 'r');
    
    set(gca, 'XScale', 'log');
    xlabel('$\|X - W Y_{test}\|_F$', 'Interpreter', 'latex');
    ylabel('Sparsity (\%)', 'Interpreter', 'latex');
    title('Testing Pareto Front', 'Interpreter', 'latex');
    grid on;
    
    sgtitle('\bf Pareto Front: Sparsity vs. Error', 'Interpreter', 'latex');
    
    [~, current_file_name, ~] = fileparts(current_file);
    path = fullfile(target_dir, string(current_file_name) + ".pdf");
    exportgraphics(gcf, path, 'ContentType', 'vector');

    close(f);
end
