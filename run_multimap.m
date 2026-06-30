% Multi-map batch: 5 algorithms × 3 maps × 10 runs
maps = {'mountain','urban','coastal'};
algs = {'SPSO','GWO','AGWO','WOA','CPO'};
N_RUNS = 10;
MAX_IT = 200;

Results = struct();
for mi = 1:3
    load(['model_' maps{mi} '.mat'], 'model');
    fprintf('\n=== MAP %d: %s ===\n', mi, maps{mi});

    nVar = model.n; VarSize = [1 nVar];
    VarMin.x=model.xmin; VarMax.x=model.xmax;
    VarMin.y=model.ymin; VarMax.y=model.ymax;
    VarMin.z=model.zmin; VarMax.z=model.zmax;
    VarMax.r=2*norm(model.start-model.end)/nVar; VarMin.r=0;
    AngleRange=pi/4; VarMin.psi=-AngleRange; VarMax.psi=AngleRange;
    dirVector=model.end-model.start;
    phi0=atan2(dirVector(2),dirVector(1));
    VarMin.phi=phi0-AngleRange; VarMax.phi=phi0+AngleRange;

    for ai = 1:5
        Costs=zeros(1,N_RUNS); Times=zeros(1,N_RUNS);
        for run=1:N_RUNS
            tic;
            switch algs{ai}
                case 'SPSO', cost = runSPSO_mm(model,500,MAX_IT);
                case 'GWO',  cost = runGWO_mm(model,150,MAX_IT);
                case 'AGWO', cost = runAGWO_mm(model,150,MAX_IT);
                case 'WOA',  cost = runWOA_mm(model,150,MAX_IT);
                case 'CPO',  cost = runCPO_mm(model,150,MAX_IT);
            end
            Costs(run)=cost; Times(run)=toc;
            fprintf('.');
        end
        Results(mi).algs(ai).name = algs{ai};
        Results(mi).algs(ai).costs = Costs;
        Results(mi).algs(ai).times = Times;
        Results(mi).algs(ai).mean_cost = mean(Costs);
        Results(mi).algs(ai).std_cost = std(Costs);
        Results(mi).algs(ai).min_cost = min(Costs);
        fprintf(' Mean=%.1f +/- %.1f (%.1fs)\n', mean(Costs), std(Costs), mean(Times));
    end
end
mkdir('results');
save('results/multi_map_results.mat','Results','maps','algs','N_RUNS','MAX_IT');

fprintf('\n=== MULTI-MAP SUMMARY ===\n');
for mi=1:3
    fprintf('\n--- %s ---\n', maps{mi});
    for ai=1:5
        a=Results(mi).algs(ai);
        fprintf('  %-6s: %8.1f +/- %7.1f\n', a.name, a.mean_cost, a.std_cost);
    end
end

%% Generate multi-map bar chart
figure('Position',[100 100 1200 450]);
colors = {[0 0.45 0.74],[0.85 0.33 0.10],[0.93 0.69 0.13],[0.49 0.18 0.56],[0.64 0.08 0.18]};
for mi=1:3
    subplot(1,3,mi);
    means=zeros(1,5); stds=zeros(1,5);
    for ai=1:5
        means(ai)=Results(mi).algs(ai).mean_cost;
        stds(ai)=Results(mi).algs(ai).std_cost;
    end
    b=bar(means); hold on;
    for i=1:5, b.FaceColor='flat'; b.CData(i,:)=colors{i}; end
    x_ends = zeros(5,1);
    for i=1:5, x_ends(i)=b.XEndPoints(i); end
    errorbar(x_ends, means, stds, 'k','LineStyle','none','LineWidth',1.5);
    set(gca,'XTickLabel',algs);
    title(sprintf('%s', maps{mi})); ylabel('Mean Best Cost'); grid on;
end
sgtitle('Algorithm Performance Across Different Terrains (10 runs)');
saveas(gcf,'figures/multi_map_comparison.png');
fprintf('Saved: figures/multi_map_comparison.png\n');
