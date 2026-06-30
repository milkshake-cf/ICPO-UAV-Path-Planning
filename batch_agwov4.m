%_________________________________________________________________________%
%  Batch Test: AGWO v4 (3-Leader + Velocity Momentum)                      %
%  Runs 20 independent trials, prints statistics                           %
%_________________________________________________________________________%

clc; clear; close all;

N_RUNS = 20;
MaxIt = 200;

finalCosts = zeros(1, N_RUNS);
bestCosts = zeros(N_RUNS, MaxIt);
times = zeros(1, N_RUNS);

for run = 1:N_RUNS
    tic;

    %% Problem Definition
    model = CreateModel();
    CostFunction=@(x) MyCost(x,model);
    nVar=model.n;
    VarSize=[1 nVar];

    VarMin.x=model.xmin; VarMax.x=model.xmax;
    VarMin.y=model.ymin; VarMax.y=model.ymax;
    VarMin.z=model.zmin; VarMax.z=model.zmax;
    VarMax.r=2*norm(model.start-model.end)/nVar; VarMin.r=0;
    AngleRange = pi/4;
    VarMin.psi=-AngleRange; VarMax.psi=AngleRange;
    dirVector = model.end - model.start;
    phi0 = atan2(dirVector(2),dirVector(1));
    VarMin.phi=phi0 - AngleRange; VarMax.phi=phi0 + AngleRange;

    %% Parameters
    nPop=150;
    w=1.0; wdamp=0.98; c1=1.5; c2=1.5;
    alpha_v=0.5;
    VelMax.r=alpha_v*(VarMax.r-VarMin.r); VelMin.r=-VelMax.r;
    VelMax.psi=alpha_v*(VarMax.psi-VarMin.psi); VelMin.psi=-VelMax.psi;
    VelMax.phi=alpha_v*(VarMax.phi-VarMin.phi); VelMin.phi=-VelMax.phi;

    %% Initialization
    empty_agent.Position=[]; empty_agent.Velocity=[]; empty_agent.Cost=[];
    empty_agent.pBest.Position=[]; empty_agent.pBest.Cost=[];

    Alpha.Cost=inf; Alpha.Position=[];
    Beta.Cost=inf; Beta.Position=[];
    Delta.Cost=inf; Delta.Position=[];

    pop=repmat(empty_agent,nPop,1);

    isInit = false;
    while (~isInit)
        for i=1:nPop
            pop(i).Position=CreateRandomSolution(VarSize,VarMin,VarMax);
            pop(i).Velocity.r=zeros(VarSize); pop(i).Velocity.psi=zeros(VarSize); pop(i).Velocity.phi=zeros(VarSize);
            cartPos = SphericalToCart(pop(i).Position, model);
            if any(isnan(cartPos.x)) || any(isnan(cartPos.y)) || any(isnan(cartPos.z))
                pop(i).Cost = inf;
            else
                try; pop(i).Cost = CostFunction(cartPos); catch; pop(i).Cost = inf; end
            end
            pop(i).pBest.Position = pop(i).Position;
            pop(i).pBest.Cost = pop(i).Cost;

            if pop(i).pBest.Cost < Alpha.Cost
                Delta = Beta; Beta = Alpha;
                Alpha.Position = pop(i).pBest.Position; Alpha.Cost = pop(i).pBest.Cost;
                isInit = true;
            elseif pop(i).pBest.Cost < Beta.Cost
                Delta = Beta;
                Beta.Position = pop(i).pBest.Position; Beta.Cost = pop(i).pBest.Cost;
            elseif pop(i).pBest.Cost < Delta.Cost
                Delta.Position = pop(i).pBest.Position; Delta.Cost = pop(i).pBest.Cost;
            end
        end
    end
    if isempty(Beta.Position); Beta = Alpha; end
    if isempty(Delta.Position); Delta = Alpha; end

    BestCostRun=zeros(MaxIt,1);

    %% Main Loop
    for t=1:MaxIt
        BestCostRun(t)=Alpha.Cost;
        a = 2 * (0.7 * (1 - t/MaxIt)^0.5 + 0.3);

        for i=1:nPop
            % GWO target: r
            A1=2*a*rand(VarSize)-a; C1=2*rand(VarSize); D_alpha=abs(C1.*Alpha.Position.r-pop(i).Position.r); X1=Alpha.Position.r-A1.*D_alpha;
            A2=2*a*rand(VarSize)-a; C2=2*rand(VarSize); D_beta=abs(C2.*Beta.Position.r-pop(i).Position.r); X2=Beta.Position.r-A2.*D_beta;
            A3=2*a*rand(VarSize)-a; C3=2*rand(VarSize); D_delta=abs(C3.*Delta.Position.r-pop(i).Position.r); X3=Delta.Position.r-A3.*D_delta;
            target_r=(X1+X2+X3)/3;
            % GWO target: psi
            A1=2*a*rand(VarSize)-a; C1=2*rand(VarSize); D_alpha=abs(C1.*Alpha.Position.psi-pop(i).Position.psi); X1=Alpha.Position.psi-A1.*D_alpha;
            A2=2*a*rand(VarSize)-a; C2=2*rand(VarSize); D_beta=abs(C2.*Beta.Position.psi-pop(i).Position.psi); X2=Beta.Position.psi-A2.*D_beta;
            A3=2*a*rand(VarSize)-a; C3=2*rand(VarSize); D_delta=abs(C3.*Delta.Position.psi-pop(i).Position.psi); X3=Delta.Position.psi-A3.*D_delta;
            target_psi=(X1+X2+X3)/3;
            % GWO target: phi
            A1=2*a*rand(VarSize)-a; C1=2*rand(VarSize); D_alpha=abs(C1.*Alpha.Position.phi-pop(i).Position.phi); X1=Alpha.Position.phi-A1.*D_alpha;
            A2=2*a*rand(VarSize)-a; C2=2*rand(VarSize); D_beta=abs(C2.*Beta.Position.phi-pop(i).Position.phi); X2=Beta.Position.phi-A2.*D_beta;
            A3=2*a*rand(VarSize)-a; C3=2*rand(VarSize); D_delta=abs(C3.*Delta.Position.phi-pop(i).Position.phi); X3=Delta.Position.phi-A3.*D_delta;
            target_phi=(X1+X2+X3)/3;

            % Velocity update + mirroring: r
            pop(i).Velocity.r=w*pop(i).Velocity.r+c1*rand(VarSize).*(pop(i).pBest.Position.r-pop(i).Position.r)+c2*rand(VarSize).*(target_r-pop(i).Position.r);
            pop(i).Velocity.r=max(pop(i).Velocity.r,VelMin.r); pop(i).Velocity.r=min(pop(i).Velocity.r,VelMax.r);
            pop(i).Position.r=pop(i).Position.r+pop(i).Velocity.r;
            OutOfRange=(pop(i).Position.r<VarMin.r|pop(i).Position.r>VarMax.r); pop(i).Velocity.r(OutOfRange)=-pop(i).Velocity.r(OutOfRange);
            pop(i).Position.r=max(pop(i).Position.r,VarMin.r); pop(i).Position.r=min(pop(i).Position.r,VarMax.r);

            % Velocity update + mirroring: psi
            pop(i).Velocity.psi=w*pop(i).Velocity.psi+c1*rand(VarSize).*(pop(i).pBest.Position.psi-pop(i).Position.psi)+c2*rand(VarSize).*(target_psi-pop(i).Position.psi);
            pop(i).Velocity.psi=max(pop(i).Velocity.psi,VelMin.psi); pop(i).Velocity.psi=min(pop(i).Velocity.psi,VelMax.psi);
            pop(i).Position.psi=pop(i).Position.psi+pop(i).Velocity.psi;
            OutOfRange=(pop(i).Position.psi<VarMin.psi|pop(i).Position.psi>VarMax.psi); pop(i).Velocity.psi(OutOfRange)=-pop(i).Velocity.psi(OutOfRange);
            pop(i).Position.psi=max(pop(i).Position.psi,VarMin.psi); pop(i).Position.psi=min(pop(i).Position.psi,VarMax.psi);

            % Velocity update + mirroring: phi
            pop(i).Velocity.phi=w*pop(i).Velocity.phi+c1*rand(VarSize).*(pop(i).pBest.Position.phi-pop(i).Position.phi)+c2*rand(VarSize).*(target_phi-pop(i).Position.phi);
            pop(i).Velocity.phi=max(pop(i).Velocity.phi,VelMin.phi); pop(i).Velocity.phi=min(pop(i).Velocity.phi,VelMax.phi);
            pop(i).Position.phi=pop(i).Position.phi+pop(i).Velocity.phi;
            OutOfRange=(pop(i).Position.phi<VarMin.phi|pop(i).Position.phi>VarMax.phi); pop(i).Velocity.phi(OutOfRange)=-pop(i).Velocity.phi(OutOfRange);
            pop(i).Position.phi=max(pop(i).Position.phi,VarMin.phi); pop(i).Position.phi=min(pop(i).Position.phi,VarMax.phi);

            % Evaluation
            cartPos = SphericalToCart(pop(i).Position, model);
            if any(isnan(cartPos.x))||any(isnan(cartPos.y))||any(isnan(cartPos.z))||any(isinf(cartPos.x))||any(isinf(cartPos.y))||any(isinf(cartPos.z))
                newCost=inf;
            else
                try; newCost=CostFunction(cartPos); catch; newCost=inf; end
            end
            pop(i).Cost=newCost;
            if newCost<pop(i).pBest.Cost; pop(i).pBest.Position=pop(i).Position; pop(i).pBest.Cost=newCost; end
        end

        % Re-rank leaders
        Alpha.Cost=inf; Beta.Cost=inf; Delta.Cost=inf;
        for i=1:nPop
            if pop(i).pBest.Cost<Alpha.Cost; Delta=Beta; Beta=Alpha; Alpha.Position=pop(i).pBest.Position; Alpha.Cost=pop(i).pBest.Cost;
            elseif pop(i).pBest.Cost<Beta.Cost; Delta=Beta; Beta.Position=pop(i).pBest.Position; Beta.Cost=pop(i).pBest.Cost;
            elseif pop(i).pBest.Cost<Delta.Cost; Delta.Position=pop(i).pBest.Position; Delta.Cost=pop(i).pBest.Cost; end
        end
        if isempty(Beta.Position); Beta=Alpha; end
        if isempty(Delta.Position); Delta=Alpha; end

        w=w*wdamp;
    end

    finalCosts(run)=Alpha.Cost;
    bestCosts(run,:)=BestCostRun;
    times(run)=toc;
    fprintf('AGWOv4 Run %2d/%d: Best=%8.2f, Time=%.1fs\n', run, N_RUNS, Alpha.Cost, times(run));
end

%% Statistics
fprintf('\n========== AGWO v4 BATCH RESULTS (N_RUNS=%d, MaxIt=%d) ==========\n', N_RUNS, MaxIt);
fprintf('Best:    %12.2f\n', min(finalCosts));
fprintf('Worst:   %12.2f\n', max(finalCosts));
fprintf('Mean:    %12.2f\n', mean(finalCosts));
fprintf('Std:     %12.2f\n', std(finalCosts));
fprintf('Avg Time: %11.1f s\n', mean(times));

fprintf('\n=== Comparison with baselines ===\n');
fprintf('SPSO   (500p, 200it): Mean ~ 4,802\n');
fprintf('GWO    (150p, 200it): Mean ~ 4,932\n');
fprintf('AGWOv3 (150p, 200it): Mean ~ 4,875\n');
fprintf('AGWOv4 (150p, 200it): Mean ~ %.0f\n', mean(finalCosts));

save('results/AGWOv4_batch_results.mat', 'finalCosts', 'bestCosts', 'times', 'N_RUNS', 'MaxIt');
fprintf('\nResults saved to results/AGWOv4_batch_results.mat\n');
