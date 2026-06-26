%_________________________________________________________________________%
%  Ablation Batch: ICPO vs no-pBest vs no-Adaptive                         %
%  20 runs each, isolates each improvement's contribution                  %
%_________________________________________________________________________%

clc; clear; close all;
N_RUNS = 20; MaxIt = 200; nPop = 150;

algs = {'ICPO_full', 'ICPO_noPbest', 'ICPO_noAdapt'};
results = struct();
for a = 1:3; alg = algs{a};
    results.(alg).bestCosts = zeros(N_RUNS, MaxIt);
    results.(alg).finalCosts = zeros(1, N_RUNS);
    results.(alg).times = zeros(1, N_RUNS);
end

%% Run each variant
for variant = 1:3
    vname = algs{variant};
    disp(['========== Running ' vname ' ==========']);
    
    for run = 1:N_RUNS
        tic;
        model = CreateModel();
        CostFunction=@(x) MyCost(x,model);
        nVar=model.n; VarSize=[1 nVar];
        VarMin.x=model.xmin; VarMax.x=model.xmax;
        VarMin.y=model.ymin; VarMax.y=model.ymax;
        VarMin.z=model.zmin; VarMax.z=model.zmax;
        VarMax.r=2*norm(model.start-model.end)/nVar; VarMin.r=0;
        AngleRange=pi/4; VarMin.psi=-AngleRange; VarMax.psi=AngleRange;
        dirVector=model.end-model.start; phi0=atan2(dirVector(2),dirVector(1));
        VarMin.phi=phi0-AngleRange; VarMax.phi=phi0+AngleRange;
        alpha=0.2; Tf=0.8;
        
        use_pBest = ~contains(vname, 'noPbest');
        use_Adapt = ~contains(vname, 'noAdapt');
        
        empty_agent.Position=[]; empty_agent.Cost=[];
        if use_pBest; empty_agent.pBest.Position=[]; empty_agent.pBest.Cost=[]; end
        GlobalBest.Cost=inf; GlobalBest.Position=[];
        pop=repmat(empty_agent,nPop,1);
        prev_pos=cell(nPop,1);
        
        isInit=false;
        while ~isInit
            for i=1:nPop
                pop(i).Position=CreateRandomSolution(VarSize,VarMin,VarMax);
                cp=SphericalToCart(pop(i).Position,model);
                if any(isnan(cp.x))||any(isnan(cp.y))||any(isnan(cp.z)); pop(i).Cost=inf;
                else; try; pop(i).Cost=CostFunction(cp); catch; pop(i).Cost=inf; end; end
                if use_pBest; pop(i).pBest.Position=pop(i).Position; pop(i).pBest.Cost=pop(i).Cost; end
                prev_pos{i}=pop(i).Position;
                if pop(i).Cost<GlobalBest.Cost; GlobalBest.Position=pop(i).Position; GlobalBest.Cost=pop(i).Cost; isInit=true; end
            end
        end
        BestCost=zeros(MaxIt,1);
        
        for t=1:MaxIt
            BestCost(t)=GlobalBest.Cost;
            if use_Adapt; explRatio=0.7*(1-t/MaxIt)^0.5+0.3; else; explRatio=0.5; end
            
            for i=1:nPop
                U1=rand(VarSize)>rand();
                if rand()<explRatio
                    if rand()<rand()
                        k=randi(nPop); m=randi(nPop);
                        if use_pBest
                            gr=(pop(k).pBest.Position.r+pop(m).pBest.Position.r)/2;
                            gpsi=(pop(k).pBest.Position.psi+pop(m).pBest.Position.psi)/2;
                            gphi=(pop(k).pBest.Position.phi+pop(m).pBest.Position.phi)/2;
                        else
                            gr=(pop(k).Position.r+pop(m).Position.r)/2;
                            gpsi=(pop(k).Position.psi+pop(m).Position.psi)/2;
                            gphi=(pop(k).Position.phi+pop(m).Position.phi)/2;
                        end
                        pop(i).Position.r=pop(i).Position.r+randn(VarSize).*abs(2*rand()*GlobalBest.Position.r-gr);
                        pop(i).Position.psi=pop(i).Position.psi+randn(VarSize).*abs(2*rand()*GlobalBest.Position.psi-gpsi);
                        pop(i).Position.phi=pop(i).Position.phi+randn(VarSize).*abs(2*rand()*GlobalBest.Position.phi-gphi);
                    else
                        k=randi(nPop); m=randi(nPop);
                        if use_pBest
                            yr=(pop(i).Position.r+pop(k).pBest.Position.r)/2; dr=pop(m).pBest.Position.r-pop(k).pBest.Position.r;
                            ypsi=(pop(i).Position.psi+pop(k).pBest.Position.psi)/2; dpsi=pop(m).pBest.Position.psi-pop(k).pBest.Position.psi;
                            yphi=(pop(i).Position.phi+pop(k).pBest.Position.phi)/2; dphi=pop(m).pBest.Position.phi-pop(k).pBest.Position.phi;
                        else
                            yr=(pop(i).Position.r+pop(k).Position.r)/2; dr=pop(m).Position.r-pop(k).Position.r;
                            ypsi=(pop(i).Position.psi+pop(k).Position.psi)/2; dpsi=pop(m).Position.psi-pop(k).Position.psi;
                            yphi=(pop(i).Position.phi+pop(k).Position.phi)/2; dphi=pop(m).Position.phi-pop(k).Position.phi;
                        end
                        pop(i).Position.r=U1.*pop(i).Position.r+(1-U1).*(yr+rand()*dr);
                        pop(i).Position.psi=U1.*pop(i).Position.psi+(1-U1).*(ypsi+rand()*dpsi);
                        pop(i).Position.phi=U1.*pop(i).Position.phi+(1-U1).*(yphi+rand()*dphi);
                    end
                else
                    Yt=2*rand()*(1-t/MaxIt)^(t/MaxIt); U2=(rand(VarSize)<0.5)*2-1; S=rand()*U2;
                    sc=0; for j=1:nPop
                        if use_pBest; sc=sc+pop(j).pBest.Cost; else; sc=sc+pop(j).Cost; end
                    end; sf=sc+eps;
                    if rand()<Tf
                        if use_pBest; cst=pop(i).pBest.Cost; else; cst=pop(i).Cost; end
                        St=exp(cst/sf); S=S.*Yt.*St; k=randi(nPop); m=randi(nPop);
                        if use_pBest
                            pop(i).Position.r=(1-U1).*pop(i).Position.r+U1.*(pop(k).pBest.Position.r+St*(pop(m).pBest.Position.r-pop(k).pBest.Position.r)-S);
                            pop(i).Position.psi=(1-U1).*pop(i).Position.psi+U1.*(pop(k).pBest.Position.psi+St*(pop(m).pBest.Position.psi-pop(k).pBest.Position.psi)-S);
                            pop(i).Position.phi=(1-U1).*pop(i).Position.phi+U1.*(pop(k).pBest.Position.phi+St*(pop(m).pBest.Position.phi-pop(k).pBest.Position.phi)-S);
                        else
                            pop(i).Position.r=(1-U1).*pop(i).Position.r+U1.*(pop(k).Position.r+St*(pop(m).Position.r-pop(k).Position.r)-S);
                            pop(i).Position.psi=(1-U1).*pop(i).Position.psi+U1.*(pop(k).Position.psi+St*(pop(m).Position.psi-pop(k).Position.psi)-S);
                            pop(i).Position.phi=(1-U1).*pop(i).Position.phi+U1.*(pop(k).Position.phi+St*(pop(m).Position.phi-pop(k).Position.phi)-S);
                        end
                    else
                        if use_pBest; cst=pop(i).pBest.Cost; else; cst=pop(i).Cost; end
                        Mt=exp(cst/sf); k=randi(nPop); r2_p=rand();
                        if use_pBest
                            Ft_r=rand(VarSize).*(Mt*(-pop(i).Position.r+pop(k).pBest.Position.r)); S_r=S.*Yt.*Ft_r;
                            Ft_psi=rand(VarSize).*(Mt*(-pop(i).Position.psi+pop(k).pBest.Position.psi)); S_psi=S.*Yt.*Ft_psi;
                            Ft_phi=rand(VarSize).*(Mt*(-pop(i).Position.phi+pop(k).pBest.Position.phi)); S_phi=S.*Yt.*Ft_phi;
                        else
                            Ft_r=rand(VarSize).*(Mt*(-pop(i).Position.r+pop(k).Position.r)); S_r=S.*Yt.*Ft_r;
                            Ft_psi=rand(VarSize).*(Mt*(-pop(i).Position.psi+pop(k).Position.psi)); S_psi=S.*Yt.*Ft_psi;
                            Ft_phi=rand(VarSize).*(Mt*(-pop(i).Position.phi+pop(k).Position.phi)); S_phi=S.*Yt.*Ft_phi;
                        end
                        pop(i).Position.r=GlobalBest.Position.r+(alpha*(1-r2_p)+r2_p)*(U2.*GlobalBest.Position.r-pop(i).Position.r)-S_r;
                        pop(i).Position.psi=GlobalBest.Position.psi+(alpha*(1-r2_p)+r2_p)*(U2.*GlobalBest.Position.psi-pop(i).Position.psi)-S_psi;
                        pop(i).Position.phi=GlobalBest.Position.phi+(alpha*(1-r2_p)+r2_p)*(U2.*GlobalBest.Position.phi-pop(i).Position.phi)-S_phi;
                    end
                end
                pop(i).Position.r=max(pop(i).Position.r,VarMin.r); pop(i).Position.r=min(pop(i).Position.r,VarMax.r);
                pop(i).Position.psi=max(pop(i).Position.psi,VarMin.psi); pop(i).Position.psi=min(pop(i).Position.psi,VarMax.psi);
                pop(i).Position.phi=max(pop(i).Position.phi,VarMin.phi); pop(i).Position.phi=min(pop(i).Position.phi,VarMax.phi);
                cp=SphericalToCart(pop(i).Position,model);
                if any(isnan(cp.x))||any(isnan(cp.y))||any(isnan(cp.z)); nc=inf;
                else; try; nc=CostFunction(cp); catch; nc=inf; end; end
                if pop(i).Cost<nc; pop(i).Position=prev_pos{i};
                else; prev_pos{i}=pop(i).Position; pop(i).Cost=nc;
                    if use_pBest&&nc<pop(i).pBest.Cost; pop(i).pBest.Position=pop(i).Position; pop(i).pBest.Cost=nc; end
                    if nc<GlobalBest.Cost; GlobalBest.Position=pop(i).Position; GlobalBest.Cost=nc; end
                end
            end
        end
        results.(vname).bestCosts(run,:)=BestCost;
        results.(vname).finalCosts(run)=GlobalBest.Cost;
        results.(vname).times(run)=toc;
        fprintf('%s Run %d/%d: Best=%.2f, Time=%.1fs\n',vname,run,N_RUNS,GlobalBest.Cost,results.(vname).times(run));
    end
end

%% Stats
fprintf('\n========== ABLATION STATISTICS (N_RUNS=%d, MaxIt=%d) ==========\n', N_RUNS, MaxIt);
fprintf('%-18s %10s %10s %10s %10s\n', 'Variant', 'Best', 'Worst', 'Mean', 'Std');
for a = 1:3
    alg = algs{a}; fc = results.(alg).finalCosts;
    fprintf('%-18s %10.2f %10.2f %10.2f %10.2f\n', alg, min(fc), max(fc), mean(fc), std(fc));
end

%% Convergence plot
figure(1); clf; hold on;
colors = {[0.8 0.2 0.2], [0.8 0.5 0.2], [0.5 0.5 0.5]};
styles = {'-', '--', ':'};
legends = {'ICPO (full)', 'ICPO w/o pBest', 'ICPO w/o Adaptive'};
for a = 1:3
    mc = mean(results.(algs{a}).bestCosts, 1);
    plot(1:MaxIt, mc, styles{a}, 'LineWidth', 2, 'Color', colors{a});
end
legend(legends); xlabel('Iteration'); ylabel('Mean Best Cost');
title(sprintf('Ablation Study (%d runs)', N_RUNS)); grid on;
saveas(gcf, 'figures/ablation_convergence.png');
disp('Plot saved to figures/ablation_convergence.png');

save('results/ablation_results.mat', 'results', 'algs', 'N_RUNS', 'MaxIt');
disp('Results saved to results/ablation_results.mat');
