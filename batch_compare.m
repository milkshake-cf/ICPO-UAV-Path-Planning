%_________________________________________________________________________%
%  Batch Comparison: SPSO vs CPO vs ICPO for UAV Path Planning            %
%  Runs each algorithm for N_RUNS independent times, saves statistics     %
%_________________________________________________________________________%

clc; clear; close all;

%% Settings
N_RUNS = 5;        % Number of independent runs per algorithm (test run)
MaxIt = 200;        % Iterations (consistent comparison)

% Results storage
results = struct();
results.SPSO.bestCosts = zeros(N_RUNS, MaxIt);
results.SPSO.finalCosts = zeros(1, N_RUNS);
results.SPSO.times = zeros(1, N_RUNS);

results.CPO.bestCosts = zeros(N_RUNS, MaxIt);
results.CPO.finalCosts = zeros(1, N_RUNS);
results.CPO.times = zeros(1, N_RUNS);

results.ICPO.bestCosts = zeros(N_RUNS, MaxIt);
results.ICPO.finalCosts = zeros(1, N_RUNS);
results.ICPO.times = zeros(1, N_RUNS);

%% Run SPSO
disp('========== Running SPSO ==========');
for run = 1:N_RUNS
    tic;
    % ---- SPSO Core ----
    model = CreateModel();
    CostFunction=@(x) MyCost(x,model);
    nVar=model.n;
    VarSize=[1 nVar];
    
    VarMin.x=model.xmin; VarMax.x=model.xmax;
    VarMin.y=model.ymin; VarMax.y=model.ymax;
    VarMin.z=model.zmin; VarMax.z=model.zmax;
    VarMax.r=2*norm(model.start-model.end)/nVar; VarMin.r=0;
    AngleRange=pi/4; VarMin.psi=-AngleRange; VarMax.psi=AngleRange;
    dirVector=model.end-model.start; phi0=atan2(dirVector(2),dirVector(1));
    VarMin.phi=phi0-AngleRange; VarMax.phi=phi0+AngleRange;
    
    alpha_v=0.5; VelMax.r=alpha_v*(VarMax.r-VarMin.r); VelMin.r=-VelMax.r;
    VelMax.psi=alpha_v*(VarMax.psi-VarMin.psi); VelMin.psi=-VelMax.psi;
    VelMax.phi=alpha_v*(VarMax.phi-VarMin.phi); VelMin.phi=-VelMax.phi;
    
    nPop=500; w=1; wdamp=0.98; c1=1.5; c2=1.5;
    
    empty_particle.Position=[]; empty_particle.Velocity=[]; empty_particle.Cost=[];
    empty_particle.Best.Position=[]; empty_particle.Best.Cost=[];
    GlobalBest.Cost=inf;
    particle=repmat(empty_particle,nPop,1);
    
    isInit=false;
    while ~isInit
        for i=1:nPop
            particle(i).Position=CreateRandomSolution(VarSize,VarMin,VarMax);
            particle(i).Velocity.r=zeros(VarSize);
            particle(i).Velocity.psi=zeros(VarSize);
            particle(i).Velocity.phi=zeros(VarSize);
            particle(i).Cost=CostFunction(SphericalToCart(particle(i).Position,model));
            particle(i).Best.Position=particle(i).Position;
            particle(i).Best.Cost=particle(i).Cost;
            if particle(i).Best.Cost < GlobalBest.Cost
                GlobalBest=particle(i).Best; isInit=true;
            end
        end
    end
    BestCostSPSO=zeros(MaxIt,1);
    for it=1:MaxIt
        BestCostSPSO(it)=GlobalBest.Cost;
        for i=1:nPop
            particle(i).Velocity.r=w*particle(i).Velocity.r+c1*rand(VarSize).*(particle(i).Best.Position.r-particle(i).Position.r)+c2*rand(VarSize).*(GlobalBest.Position.r-particle(i).Position.r);
            particle(i).Velocity.r=max(particle(i).Velocity.r,VelMin.r); particle(i).Velocity.r=min(particle(i).Velocity.r,VelMax.r);
            particle(i).Position.r=particle(i).Position.r+particle(i).Velocity.r;
            OutOfRange=(particle(i).Position.r<VarMin.r|particle(i).Position.r>VarMax.r);
            particle(i).Velocity.r(OutOfRange)=-particle(i).Velocity.r(OutOfRange);
            particle(i).Position.r=max(particle(i).Position.r,VarMin.r); particle(i).Position.r=min(particle(i).Position.r,VarMax.r);
            
            particle(i).Velocity.psi=w*particle(i).Velocity.psi+c1*rand(VarSize).*(particle(i).Best.Position.psi-particle(i).Position.psi)+c2*rand(VarSize).*(GlobalBest.Position.psi-particle(i).Position.psi);
            particle(i).Velocity.psi=max(particle(i).Velocity.psi,VelMin.psi); particle(i).Velocity.psi=min(particle(i).Velocity.psi,VelMax.psi);
            particle(i).Position.psi=particle(i).Position.psi+particle(i).Velocity.psi;
            OutOfRange=(particle(i).Position.psi<VarMin.psi|particle(i).Position.psi>VarMax.psi);
            particle(i).Velocity.psi(OutOfRange)=-particle(i).Velocity.psi(OutOfRange);
            particle(i).Position.psi=max(particle(i).Position.psi,VarMin.psi); particle(i).Position.psi=min(particle(i).Position.psi,VarMax.psi);
            
            particle(i).Velocity.phi=w*particle(i).Velocity.phi+c1*rand(VarSize).*(particle(i).Best.Position.phi-particle(i).Position.phi)+c2*rand(VarSize).*(GlobalBest.Position.phi-particle(i).Position.phi);
            particle(i).Velocity.phi=max(particle(i).Velocity.phi,VelMin.phi); particle(i).Velocity.phi=min(particle(i).Velocity.phi,VelMax.phi);
            particle(i).Position.phi=particle(i).Position.phi+particle(i).Velocity.phi;
            OutOfRange=(particle(i).Position.phi<VarMin.phi|particle(i).Position.phi>VarMax.phi);
            particle(i).Velocity.phi(OutOfRange)=-particle(i).Velocity.phi(OutOfRange);
            particle(i).Position.phi=max(particle(i).Position.phi,VarMin.phi); particle(i).Position.phi=min(particle(i).Position.phi,VarMax.phi);
            
            particle(i).Cost=CostFunction(SphericalToCart(particle(i).Position,model));
            if particle(i).Cost < particle(i).Best.Cost
                particle(i).Best.Position=particle(i).Position; particle(i).Best.Cost=particle(i).Cost;
                if particle(i).Best.Cost < GlobalBest.Cost
                    GlobalBest=particle(i).Best;
                end
            end
        end
        w=w*wdamp;
    end
    results.SPSO.bestCosts(run,:) = BestCostSPSO;
    results.SPSO.finalCosts(run) = GlobalBest.Cost;
    results.SPSO.times(run) = toc;
    fprintf('SPSO Run %d/%d: Best=%8.2f, Time=%.1fs\n', run, N_RUNS, GlobalBest.Cost, results.SPSO.times(run));
end

%% Run CPO
disp('========== Running CPO ==========');
for run = 1:N_RUNS
    tic;
    % ---- CPO Core (simplified inline) ----
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
    
    nPop=150; alpha=0.2; Tf=0.8;
    empty_agent.Position=[]; empty_agent.Cost=[];
    GlobalBest.Cost=inf;
    pop=repmat(empty_agent,nPop,1);
    prev_pos=cell(nPop,1);
    
    isInit=false;
    while ~isInit
        for i=1:nPop
            pop(i).Position=CreateRandomSolution(VarSize,VarMin,VarMax);
            cp=SphericalToCart(pop(i).Position,model);
            if any(isnan(cp.x))||any(isnan(cp.y))||any(isnan(cp.z)); pop(i).Cost=inf;
            else; try; pop(i).Cost=CostFunction(cp); catch; pop(i).Cost=inf; end; end
            prev_pos{i}=pop(i).Position;
            if pop(i).Cost<GlobalBest.Cost; GlobalBest.Position=pop(i).Position; GlobalBest.Cost=pop(i).Cost; isInit=true; end
        end
    end
    BestCostCPO=zeros(MaxIt,1);
    for t=1:MaxIt
        BestCostCPO(t)=GlobalBest.Cost;
        for i=1:nPop
            U1=rand(VarSize)>rand();
            if rand()<rand()
                if rand()<rand()
                    k=randi(nPop); y_r=(pop(i).Position.r+pop(k).Position.r)/2;
                    pop(i).Position.r=pop(i).Position.r+randn(VarSize).*abs(2*rand()*GlobalBest.Position.r-y_r);
                    y_psi=(pop(i).Position.psi+pop(k).Position.psi)/2;
                    pop(i).Position.psi=pop(i).Position.psi+randn(VarSize).*abs(2*rand()*GlobalBest.Position.psi-y_psi);
                    y_phi=(pop(i).Position.phi+pop(k).Position.phi)/2;
                    pop(i).Position.phi=pop(i).Position.phi+randn(VarSize).*abs(2*rand()*GlobalBest.Position.phi-y_phi);
                else
                    k=randi(nPop); m=randi(nPop);
                    y_r=(pop(i).Position.r+pop(k).Position.r)/2;
                    pop(i).Position.r=U1.*pop(i).Position.r+(1-U1).*(y_r+rand()*(pop(m).Position.r-pop(k).Position.r));
                    y_psi=(pop(i).Position.psi+pop(k).Position.psi)/2;
                    pop(i).Position.psi=U1.*pop(i).Position.psi+(1-U1).*(y_psi+rand()*(pop(m).Position.psi-pop(k).Position.psi));
                    y_phi=(pop(i).Position.phi+pop(k).Position.phi)/2;
                    pop(i).Position.phi=U1.*pop(i).Position.phi+(1-U1).*(y_phi+rand()*(pop(m).Position.phi-pop(k).Position.phi));
                end
            else
                Yt=2*rand()*(1-t/MaxIt)^(t/MaxIt); U2=(rand(VarSize)<0.5)*2-1; S=rand()*U2;
                sc=0; for j=1:nPop; sc=sc+pop(j).Cost; end; sf=sc+eps;
                if rand()<Tf
                    St=exp(pop(i).Cost/sf); S=S.*Yt.*St; k=randi(nPop); m=randi(nPop);
                    pop(i).Position.r=(1-U1).*pop(i).Position.r+U1.*(pop(k).Position.r+St*(pop(m).Position.r-pop(k).Position.r)-S);
                    pop(i).Position.psi=(1-U1).*pop(i).Position.psi+U1.*(pop(k).Position.psi+St*(pop(m).Position.psi-pop(k).Position.psi)-S);
                    pop(i).Position.phi=(1-U1).*pop(i).Position.phi+U1.*(pop(k).Position.phi+St*(pop(m).Position.phi-pop(k).Position.phi)-S);
                else
                    Mt=exp(pop(i).Cost/sf); k=randi(nPop); r2_p=rand();
                    Ft_r=rand(VarSize).*(Mt*(-pop(i).Position.r+pop(k).Position.r)); S_r=S.*Yt.*Ft_r;
                    pop(i).Position.r=GlobalBest.Position.r+(alpha*(1-r2_p)+r2_p)*(U2.*GlobalBest.Position.r-pop(i).Position.r)-S_r;
                    Ft_psi=rand(VarSize).*(Mt*(-pop(i).Position.psi+pop(k).Position.psi)); S_psi=S.*Yt.*Ft_psi;
                    pop(i).Position.psi=GlobalBest.Position.psi+(alpha*(1-r2_p)+r2_p)*(U2.*GlobalBest.Position.psi-pop(i).Position.psi)-S_psi;
                    Ft_phi=rand(VarSize).*(Mt*(-pop(i).Position.phi+pop(k).Position.phi)); S_phi=S.*Yt.*Ft_phi;
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
            else
                prev_pos{i}=pop(i).Position; pop(i).Cost=nc;
                if nc<GlobalBest.Cost; GlobalBest.Position=pop(i).Position; GlobalBest.Cost=nc; end
            end
        end
    end
    results.CPO.bestCosts(run,:) = BestCostCPO;
    results.CPO.finalCosts(run) = GlobalBest.Cost;
    results.CPO.times(run) = toc;
    fprintf('CPO  Run %d/%d: Best=%8.2f, Time=%.1fs\n', run, N_RUNS, GlobalBest.Cost, results.CPO.times(run));
end

%% Run ICPO
disp('========== Running ICPO ==========');
for run = 1:N_RUNS
    tic;
    % ---- ICPO Core (simplified inline) ----
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
    
    nPop=150; alpha=0.2; Tf=0.8;
    empty_agent.Position=[]; empty_agent.Cost=[]; empty_agent.pBest.Position=[]; empty_agent.pBest.Cost=[];
    GlobalBest.Cost=inf;
    pop=repmat(empty_agent,nPop,1);
    prev_pos=cell(nPop,1);
    
    isInit=false;
    while ~isInit
        for i=1:nPop
            pop(i).Position=CreateRandomSolution(VarSize,VarMin,VarMax);
            cp=SphericalToCart(pop(i).Position,model);
            if any(isnan(cp.x))||any(isnan(cp.y))||any(isnan(cp.z)); pop(i).Cost=inf;
            else; try; pop(i).Cost=CostFunction(cp); catch; pop(i).Cost=inf; end; end
            pop(i).pBest.Position=pop(i).Position; pop(i).pBest.Cost=pop(i).Cost;
            prev_pos{i}=pop(i).Position;
            if pop(i).Cost<GlobalBest.Cost; GlobalBest.Position=pop(i).Position; GlobalBest.Cost=pop(i).Cost; isInit=true; end
        end
    end
    BestCostICPO=zeros(MaxIt,1);
    for t=1:MaxIt
        BestCostICPO(t)=GlobalBest.Cost;
        explRatio=0.7*(1-t/MaxIt)^0.5+0.3;
        for i=1:nPop
            U1=rand(VarSize)>rand();
            if rand()<explRatio
                if rand()<rand()
                    k=randi(nPop); m=randi(nPop);
                    gr=(pop(k).pBest.Position.r+pop(m).pBest.Position.r)/2;
                    pop(i).Position.r=pop(i).Position.r+randn(VarSize).*abs(2*rand()*GlobalBest.Position.r-gr);
                    gpsi=(pop(k).pBest.Position.psi+pop(m).pBest.Position.psi)/2;
                    pop(i).Position.psi=pop(i).Position.psi+randn(VarSize).*abs(2*rand()*GlobalBest.Position.psi-gpsi);
                    gphi=(pop(k).pBest.Position.phi+pop(m).pBest.Position.phi)/2;
                    pop(i).Position.phi=pop(i).Position.phi+randn(VarSize).*abs(2*rand()*GlobalBest.Position.phi-gphi);
                else
                    k=randi(nPop); m=randi(nPop);
                    yr=(pop(i).Position.r+pop(k).pBest.Position.r)/2; dr=pop(m).pBest.Position.r-pop(k).pBest.Position.r;
                    pop(i).Position.r=U1.*pop(i).Position.r+(1-U1).*(yr+rand()*dr);
                    ypsi=(pop(i).Position.psi+pop(k).pBest.Position.psi)/2; dpsi=pop(m).pBest.Position.psi-pop(k).pBest.Position.psi;
                    pop(i).Position.psi=U1.*pop(i).Position.psi+(1-U1).*(ypsi+rand()*dpsi);
                    yphi=(pop(i).Position.phi+pop(k).pBest.Position.phi)/2; dphi=pop(m).pBest.Position.phi-pop(k).pBest.Position.phi;
                    pop(i).Position.phi=U1.*pop(i).Position.phi+(1-U1).*(yphi+rand()*dphi);
                end
            else
                Yt=2*rand()*(1-t/MaxIt)^(t/MaxIt); U2=(rand(VarSize)<0.5)*2-1; S=rand()*U2;
                sc=0; for j=1:nPop; sc=sc+pop(j).pBest.Cost; end; sf=sc+eps;
                if rand()<Tf
                    St=exp(pop(i).pBest.Cost/sf); S=S.*Yt.*St; k=randi(nPop); m=randi(nPop);
                    pop(i).Position.r=(1-U1).*pop(i).Position.r+U1.*(pop(k).pBest.Position.r+St*(pop(m).pBest.Position.r-pop(k).pBest.Position.r)-S);
                    pop(i).Position.psi=(1-U1).*pop(i).Position.psi+U1.*(pop(k).pBest.Position.psi+St*(pop(m).pBest.Position.psi-pop(k).pBest.Position.psi)-S);
                    pop(i).Position.phi=(1-U1).*pop(i).Position.phi+U1.*(pop(k).pBest.Position.phi+St*(pop(m).pBest.Position.phi-pop(k).pBest.Position.phi)-S);
                else
                    Mt=exp(pop(i).pBest.Cost/sf); k=randi(nPop); r2_p=rand();
                    Ft_r=rand(VarSize).*(Mt*(-pop(i).Position.r+pop(k).pBest.Position.r)); S_r=S.*Yt.*Ft_r;
                    pop(i).Position.r=GlobalBest.Position.r+(alpha*(1-r2_p)+r2_p)*(U2.*GlobalBest.Position.r-pop(i).Position.r)-S_r;
                    Ft_psi=rand(VarSize).*(Mt*(-pop(i).Position.psi+pop(k).pBest.Position.psi)); S_psi=S.*Yt.*Ft_psi;
                    pop(i).Position.psi=GlobalBest.Position.psi+(alpha*(1-r2_p)+r2_p)*(U2.*GlobalBest.Position.psi-pop(i).Position.psi)-S_psi;
                    Ft_phi=rand(VarSize).*(Mt*(-pop(i).Position.phi+pop(k).pBest.Position.phi)); S_phi=S.*Yt.*Ft_phi;
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
            else
                prev_pos{i}=pop(i).Position; pop(i).Cost=nc;
                if nc<pop(i).pBest.Cost; pop(i).pBest.Position=pop(i).Position; pop(i).pBest.Cost=nc; end
                if nc<GlobalBest.Cost; GlobalBest.Position=pop(i).Position; GlobalBest.Cost=nc; end
            end
        end
    end
    results.ICPO.bestCosts(run,:) = BestCostICPO;
    results.ICPO.finalCosts(run) = GlobalBest.Cost;
    results.ICPO.times(run) = toc;
    fprintf('ICPO Run %d/%d: Best=%8.2f, Time=%.1fs\n', run, N_RUNS, GlobalBest.Cost, results.ICPO.times(run));
end

%% Compute Statistics & Display
fprintf('\n========== FINAL STATISTICS (N_RUNS=%d, MaxIt=%d) ==========\n', N_RUNS, MaxIt);
fprintf('%-10s %12s %12s %12s %12s %12s\n', 'Algorithm', 'Best', 'Worst', 'Mean', 'Std', 'Avg Time(s)');

algs = {'SPSO', 'CPO', 'ICPO'};
stat_table = zeros(3, 5);
for a = 1:3
    alg = algs{a};
    fc = results.(alg).finalCosts;
    stat_table(a,:) = [min(fc), max(fc), mean(fc), std(fc), mean(results.(alg).times)];
    fprintf('%-10s %12.2f %12.2f %12.2f %12.2f %12.1f\n', alg, stat_table(a,1), stat_table(a,2), stat_table(a,3), stat_table(a,4), stat_table(a,5));
end

%% Plot Convergence Curves (mean over runs)
figure(1); clf; hold on;
colors = {[0.0 0.4 0.8], [0.8 0.4 0.0], [0.8 0.2 0.2]};
styles = {'-', '--', '-.'};
for a = 1:3
    alg = algs{a};
    meanCurve = mean(results.(alg).bestCosts, 1);
    plot(1:MaxIt, meanCurve, styles{a}, 'LineWidth', 2, 'Color', colors{a});
end
legend('SPSO', 'CPO', 'ICPO', 'Location', 'northeast');
xlabel('Iteration'); ylabel('Mean Best Cost');
title(sprintf('Convergence Comparison (%d runs)', N_RUNS));
grid on; hold off;
saveas(gcf, 'convergence_comparison.png');

%% Boxplot
figure(2); clf;
boxplot([results.SPSO.finalCosts; results.CPO.finalCosts; results.ICPO.finalCosts]', algs);
ylabel('Final Best Cost');
title(sprintf('Distribution of Final Costs (%d runs)', N_RUNS));
grid on;
saveas(gcf, 'boxplot_comparison.png');

%% Save results
save('batch_comparison_results.mat', 'results', 'stat_table', 'algs', 'N_RUNS', 'MaxIt');
fprintf('\nResults saved to batch_comparison_results.mat\n');
fprintf('Convergence plot saved to convergence_comparison.png\n');
fprintf('Boxplot saved to boxplot_comparison.png\n');
