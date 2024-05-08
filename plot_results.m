% Plotting the results

clc; clear; close all;

fc = 0;

% Here, the K and eps1coeff values should match that of the results file
datanametoload= ['K_vec_' sprintf('%d_', [20]) sprintf('eps1coeff_%02d_', 100*0.80) ...
    'methods' sprintf('_%d', [0 0 0 0 0 0 1 2 3 4 5 6])]; 

load(datanametoload);

% Color settings for each approach
color_mtx = [repmat((0:0.1:0.5)', 1, 3);[[0, 1, 0.5]; [0, 0, 1]; [1, 0, 0]; [1, 0.5, 0]; [0, 0, 0.6]; [0, 0.5, 1]]];

for i1 = 1:L_K % 1:L_K
    K = K_vec(i1);
    for i2 = 1:L_r
        rho_coeff = rho_coeff_vec(i2);
        fc = fc + 1; f = figure(fc);
        f.Position = [100 100 300 150];
        for i3 = 1:L_m
            loss_type = loss_type_vec(i3);
            boxplot(squeeze(TV{i1, i2}(:, i3, :))', eps_DP_vec, 'positions',...
                (1:L_e)+(i3-6.5)*0.07, 'PlotStyle', 'compact', 'color', ...
                color_mtx(i3, :), 'Symbol', '.');
            hold on;
        end
        hold off;
        set(gca, 'xtick', 1:3, 'xticklabel', eps_DP_vec);
        set(findobj(gca,'type','line'),'linew',2);
        h = findobj(gca,'Tag','Box');

        xlabel('$\epsilon$', 'Interpreter', 'Latex');
        set(gca, 'ylim', [0, 0.6], 'xlim', [0.5, 3.5]);
        title(sprintf('$K =$ %d, $\\rho$ = %.2f', K, rho_coeff), 'Interpreter', 'Latex');
        grid on;
        exportgraphics(gca, [datanametoload sprintf('rhocoeff_%03d', 100*rho_coeff), '.pdf']);

        tv_array = mean(K_selected{i1, i2}, 3);

        method_names = ["non-adapt ($\alpha = 1$)", ...
                "semi adapt $\alpha = 0.95$", "semi adapt $\alpha =0.90$",...
                "semi adapt $\alpha = 0.80$", "semi adapt $\alpha =0.60$" , "$\alpha =0.20$",...
                "FIM", "Entropy", "TV1", "TV2", "MSE", "$P(X = Y)$"];

        eps_names = ["$\epsilon = 0.5$", "$\epsilon = 1$", "$\epsilon = 5$"];
        fc = fc + 1; f = figure(fc);
        f.Position = [100 100 400 200];
        heatmap(tv_array, 'XData', method_names, 'YData', eps_names, 'Interpreter', 'Latex');
        title(sprintf('$K =$ %d, $\\rho$ = %.2f', K, rho_coeff));
        exportgraphics(gca, [datanametoload sprintf('heatmap_rhocoeff_%03d', 100*rho_coeff), '.pdf']);

    end
end