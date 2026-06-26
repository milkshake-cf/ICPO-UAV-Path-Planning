%_________________________________________________________________________%
%  5-Algorithm Convergence Curve + Wilcoxon Test                           %
%  Uses existing 20-run batch results                                      %
%_________________________________________________________________________%

clc; clear; close all;

%% Load data
load('results/batch_comparison_results.mat');  % SPSO, CPO, ICPO
spsocpoicpo = results;
load('results/batch_gwo_woa_results.mat');      % GWO, WOA

%% Merge into unified structure
algs = {'SPSO', 'GWO', 'ICPO', 'WOA', 'CPO'};
allResults = struct();
allResults.SPSO = spsocpoicpo.SPSO;
allResults.CPO  = spsocpoicpo.CPO;
allResults.ICPO = spsocpoicpo.ICPO;
allResults.GWO  = results.GWO;
allResults.WOA  = results.WOA;

N_RUNS = 20;
MaxIt = 200;

%% ===== FIGURE 1: Convergence Curves =====
figure(1); clf; hold on;
colors = {[0.0 0.4 0.8], [0.2 0.7 0.2], [0.8 0.2 0.2], [0.4 0.8 0.8], [0.8 0.5 0.2]};
styles = {'-', '-', '-', '--', '--'};
lw = [2.5, 2.0, 2.0, 1.5, 1.5];

for a = 1:5
    alg = algs{a};
    mc = mean(allResults.(alg).bestCosts, 1);
    plot(1:MaxIt, mc, styles{a}, 'LineWidth', lw(a), 'Color', colors{a});
end

legend(algs, 'Location', 'northeast', 'FontSize', 11);
xlabel('Iteration', 'FontSize', 12);
ylabel('Mean Best Cost', 'FontSize', 12);
title(sprintf('Convergence Curves — 5 Algorithms (%d runs each)', N_RUNS), 'FontSize', 13);
grid on;
set(gca, 'FontSize', 11);
saveas(gcf, 'figures/convergence_5alg.png');
fprintf('Figure 1 saved: figures/convergence_5alg.png\n');

%% ===== FIGURE 2: Final Cost Boxplot =====
figure(2); clf;
boxData = zeros(N_RUNS, 5);
for a = 1:5
    boxData(:,a) = allResults.(algs{a}).finalCosts';
end
boxplot(boxData, algs);
ylabel('Final Best Cost', 'FontSize', 12);
title(sprintf('Distribution of Final Costs (%d runs)', N_RUNS), 'FontSize', 13);
set(gca, 'FontSize', 11);
grid on;
saveas(gcf, 'figures/boxplot_5alg.png');
fprintf('Figure 2 saved: figures/boxplot_5alg.png\n');

%% ===== FIGURE 3: Bar Chart (Mean + Std) =====
figure(3); clf;
means = zeros(1,5); stds = zeros(1,5);
for a = 1:5
    means(a) = mean(allResults.(algs{a}).finalCosts);
    stds(a) = std(allResults.(algs{a}).finalCosts);
end
bar_colors = {[0.0 0.4 0.8], [0.2 0.7 0.2], [0.8 0.2 0.2], [0.4 0.8 0.8], [0.8 0.5 0.2]};
b = bar(means, 'FaceColor', 'flat');
for a = 1:5; b.CData(a,:) = colors{a}; end
hold on;
errorbar(1:5, means, stds, 'k.', 'LineWidth', 1.5);
set(gca, 'XTickLabel', algs, 'FontSize', 11);
ylabel('Mean Final Cost', 'FontSize', 12);
title(sprintf('Mean ± Std Final Cost (%d runs)', N_RUNS), 'FontSize', 13);
grid on;
% Add value labels
for a = 1:5
    text(a, means(a)+stds(a)+80, sprintf('%.0f', means(a)), ...
        'HorizontalAlignment', 'center', 'FontSize', 10, 'FontWeight', 'bold');
end
saveas(gcf, 'figures/barchart_5alg.png');
fprintf('Figure 3 saved: figures/barchart_5alg.png\n');

%% ===== Wilcoxon Rank-Sum Test =====
fprintf('\n========== WILCOXON RANK-SUM TESTS ==========\n');
fprintf('Comparing each algorithm against ICPO\n');
fprintf('H0: distributions are equal (p > 0.05 = not significantly different)\n\n');

fprintf('%-12s %-12s %12s %12s %s\n', 'Alg 1', 'Alg 2', 'p-value', 'Significant?', 'Interpretation');
fprintf('%-12s %-12s %12s %12s %s\n', '-----', '-----', '-------', '-----------', '-------------');

icpo_fc = allResults.ICPO.finalCosts;
for a = [1,2,4,5]  % Compare SPSO, GWO, WOA, CPO vs ICPO
    alg = algs{a};
    other_fc = allResults.(alg).finalCosts;
    p = ranksum(icpo_fc, other_fc);
    sig = 'NO';
    interp = 'Not significantly different';
    if p < 0.01
        sig = 'YES (p<0.01)';
        if mean(icpo_fc) < mean(other_fc)
            interp = 'ICPO is significantly BETTER';
        else
            interp = sprintf('%s is significantly BETTER', alg);
        end
    elseif p < 0.05
        sig = 'YES (p<0.05)';
        if mean(icpo_fc) < mean(other_fc)
            interp = 'ICPO is significantly BETTER';
        else
            interp = sprintf('%s is significantly BETTER', alg);
        end
    end
    fprintf('%-12s %-12s %12.6f %12s %s\n', 'ICPO', alg, p, sig, interp);
end

% Also compare within ablation variants if data exists
fprintf('\n--- Key pairwise comparisons ---\n');
fprintf('SPSO vs GWO:  p = %.6f\n', ranksum(allResults.SPSO.finalCosts, allResults.GWO.finalCosts));
fprintf('CPO vs ICPO:  p = %.6f\n', ranksum(allResults.CPO.finalCosts, allResults.ICPO.finalCosts));
fprintf('CPO vs WOA:   p = %.6f\n', ranksum(allResults.CPO.finalCosts, allResults.WOA.finalCosts));

fprintf('\nAll figures and tests complete.\n');
