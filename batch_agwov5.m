%_________________________________________________________________________%
%  Batch Test: AGWO v5 (3-Leader + Fitness Weighting + Local Search)       %
%  Runs 20 independent trials                                              %
%_________________________________________________________________________%

clc; clear; close all;

N_RUNS = 20;
MaxIt = 200;

finalCosts = zeros(1, N_RUNS);
bestCosts = zeros(N_RUNS, MaxIt);
times = zeros(1, N_RUNS);

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

    nPop=150;

    empty_agent.Position=[]; empty_agent.Cost=[]; empty_agent.pBest.Position=[]; empty_agent.pBest.Cost=[];
    Alpha.Cost=inf; Alpha.Position=[]; Beta.Cost=inf; Beta.Position=[]; Delta.Cost=inf; Delta.Position=[];
    pop=repmat(empty_agent,nPop,1);

    isInit=false;
    while ~isInit
        for i=1:nPop
            pop(i).Position=CreateRandomSolution(VarSize,VarMin,VarMax);
            cp=SphericalToCart(pop(i).Position,model);
            if any(isnan(cp.x))||any(isnan(cp.y))||any(isnan(cp.z)); pop(i).Cost=inf;
            else; try; pop(i).Cost=CostFunction(cp); catch; pop(i).Cost=inf; end; end
            pop(i).pBest.Position=pop(i).Position; pop(i).pBest.Cost=pop(i).Cost;
            if pop(i).pBest.Cost<Alpha.Cost; Delta=Beta; Beta=Alpha; Alpha.Position=pop(i).pBest.Position; Alpha.Cost=pop(i).pBest.Cost; isInit=true;
            elseif pop(i).pBest.Cost<Beta.Cost; Delta=Beta; Beta.Position=pop(i).pBest.Position; Beta.Cost=pop(i).pBest.Cost;
            elseif pop(i).pBest.Cost<Delta.Cost; Delta.Position=pop(i).pBest.Position; Delta.Cost=pop(i).pBest.Cost; end
        end
    end
    if isempty(Beta.Position); Beta=Alpha; end
    if isempty(Delta.Position); Delta=Alpha; end

    BestCostRun=zeros(MaxIt,1);

    for t=1:MaxIt
        BestCostRun(t)=Alpha.Cost;
        a=2*(0.7*(1-t/MaxIt)^0.5+0.3);

        if t<MaxIt*0.3; localProb=0.1;
        elseif t<MaxIt*0.7; localProb=0.25;
        else; localProb=0.45; end

        for i=1:nPop
            % GWO target: r
            A1=2*a*rand(VarSize)-a; C1=2*rand(VarSize); X1=Alpha.Position.r-A1.*abs(C1.*Alpha.Position.r-pop(i).Position.r);
            A2=2*a*rand(VarSize)-a; C2=2*rand(VarSize); X2=Beta.Position.r-A2.*abs(C2.*Beta.Position.r-pop(i).Position.r);
            A3=2*a*rand(VarSize)-a; C3=2*rand(VarSize); X3=Delta.Position.r-A3.*abs(C3.*Delta.Position.r-pop(i).Position.r);
            wa=1/(Alpha.Cost+eps); wb=1/(Beta.Cost+eps); wd=1/(Delta.Cost+eps); ws=wa+wb+wd;
            target_r=(wa*X1+wb*X2+wd*X3)/ws;
            % GWO target: psi
            A1=2*a*rand(VarSize)-a; C1=2*rand(VarSize); X1=Alpha.Position.psi-A1.*abs(C1.*Alpha.Position.psi-pop(i).Position.psi);
            A2=2*a*rand(VarSize)-a; C2=2*rand(VarSize); X2=Beta.Position.psi-A2.*abs(C2.*Beta.Position.psi-pop(i).Position.psi);
            A3=2*a*rand(VarSize)-a; C3=2*rand(VarSize); X3=Delta.Position.psi-A3.*abs(C3.*Delta.Position.psi-pop(i).Position.psi);
            target_psi=(wa*X1+wb*X2+wd*X3)/ws;
            % GWO target: phi
            A1=2*a*rand(VarSize)-a; C1=2*rand(VarSize); X1=Alpha.Position.phi-A1.*abs(C1.*Alpha.Position.phi-pop(i).Position.phi);
            A2=2*a*rand(VarSize)-a; C2=2*rand(VarSize); X2=Beta.Position.phi-A2.*abs(C2.*Beta.Position.phi-pop(i).Position.phi);
            A3=2*a*rand(VarSize)-a; C3=2*rand(VarSize); X3=Delta.Position.phi-A3.*abs(C3.*Delta.Position.phi-pop(i).Position.phi);
            target_phi=(wa*X1+wb*X2+wd*X3)/ws;

            % Local search
            if a<1.0 && rand()<localProb
                lr=(VarMax.r-VarMin.r)*0.03*(a/2); target_r=Alpha.Position.r+lr*randn(VarSize);
                lp=(VarMax.psi-VarMin.psi)*0.03*(a/2); target_psi=Alpha.Position.psi+lp*randn(VarSize);
                lph=(VarMax.phi-VarMin.phi)*0.03*(a/2); target_phi=Alpha.Position.phi+lph*randn(VarSize);
            end

            pop(i).Position.r=target_r; pop(i).Position.psi=target_psi; pop(i).Position.phi=target_phi;
            pop(i).Position.r=max(pop(i).Position.r,VarMin.r); pop(i).Position.r=min(pop(i).Position.r,VarMax.r);
            pop(i).Position.psi=max(pop(i).Position.psi,VarMin.psi); pop(i).Position.psi=min(pop(i).Position.psi,VarMax.psi);
            pop(i).Position.phi=max(pop(i).Position.phi,VarMin.phi); pop(i).Position.phi=min(pop(i).Position.phi,VarMax.phi);

            cp=SphericalToCart(pop(i).Position,model);
            if any(isnan(cp.x))||any(isnan(cp.y))||any(isnan(cp.z))||any(isinf(cp.x))||any(isinf(cp.y))||any(isinf(cp.z)); nc=inf;
            else; try; nc=CostFunction(cp); catch; nc=inf; end; end
            pop(i).Cost=nc;
            if nc<pop(i).pBest.Cost; pop(i).pBest.Position=pop(i).Position; pop(i).pBest.Cost=nc; end
        end

        Alpha.Cost=inf; Beta.Cost=inf; Delta.Cost=inf;
        for i=1:nPop
            if pop(i).pBest.Cost<Alpha.Cost; Delta=Beta; Beta=Alpha; Alpha.Position=pop(i).pBest.Position; Alpha.Cost=pop(i).pBest.Cost;
            elseif pop(i).pBest.Cost<Beta.Cost; Delta=Beta; Beta.Position=pop(i).pBest.Position; Beta.Cost=pop(i).pBest.Cost;
            elseif pop(i).pBest.Cost<Delta.Cost; Delta.Position=pop(i).pBest.Position; Delta.Cost=pop(i).pBest.Cost; end
        end
        if isempty(Beta.Position); Beta=Alpha; end
        if isempty(Delta.Position); Delta=Alpha; end
    end

    finalCosts(run)=Alpha.Cost;
    bestCosts(run,:)=BestCostRun;
    times(run)=toc;
    fprintf('AGWOv5 Run %2d/%d: Best=%8.2f, Time=%.1fs\n', run, N_RUNS, Alpha.Cost, times(run));
end

fprintf('\n========== AGWO v5 BATCH RESULTS (N_RUNS=%d, MaxIt=%d) ==========\n', N_RUNS, MaxIt);
fprintf('Best:    %12.2f\n', min(finalCosts));
fprintf('Worst:   %12.2f\n', max(finalCosts));
fprintf('Mean:    %12.2f\n', mean(finalCosts));
fprintf('Std:     %12.2f\n', std(finalCosts));
fprintf('Avg Time: %11.1f s\n', mean(times));

fprintf('\n=== Full Comparison ===\n');
fprintf('SPSO   (500p, 200it): Mean ~ 4,802 ± 199\n');
fprintf('GWO    (150p, 200it): Mean ~ 4,932\n');
fprintf('AGWO   (150p, 200it): Mean ~ 6,372 ± 406\n');
fprintf('AGWOv3 (150p, 200it): Mean ~ 4,875 ± 297\n');
fprintf('AGWOv5 (150p, 200it): Mean ~ %.0f ± %.0f\n', mean(finalCosts), std(finalCosts));

save('results/AGWOv5_batch_results.mat', 'finalCosts', 'bestCosts', 'times', 'N_RUNS', 'MaxIt');
fprintf('\nResults saved to results/AGWOv5_batch_results.mat\n');
