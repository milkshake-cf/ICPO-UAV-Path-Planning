%_________________________________________________________________________%
%  Batch Test: AGWO v3 (GWO-style 3-Leader Hierarchy)                      %
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

    %% Initialization
    empty_agent.Position=[];
    empty_agent.Cost=[];
    empty_agent.pBest.Position=[];
    empty_agent.pBest.Cost=[];

    Alpha.Cost=inf;   Alpha.Position=[];
    Beta.Cost=inf;    Beta.Position=[];
    Delta.Cost=inf;   Delta.Position=[];

    pop=repmat(empty_agent,nPop,1);

    isInit = false;
    while (~isInit)
        for i=1:nPop
            pop(i).Position=CreateRandomSolution(VarSize,VarMin,VarMax);
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
            % --- r component ---
            A1_r = 2*a*rand(VarSize) - a;  C1_r = 2*rand(VarSize);
            D_alpha_r = abs(C1_r .* Alpha.Position.r - pop(i).Position.r);
            X1_r = Alpha.Position.r - A1_r .* D_alpha_r;
            A2_r = 2*a*rand(VarSize) - a;  C2_r = 2*rand(VarSize);
            D_beta_r = abs(C2_r .* Beta.Position.r - pop(i).Position.r);
            X2_r = Beta.Position.r - A2_r .* D_beta_r;
            A3_r = 2*a*rand(VarSize) - a;  C3_r = 2*rand(VarSize);
            D_delta_r = abs(C3_r .* Delta.Position.r - pop(i).Position.r);
            X3_r = Delta.Position.r - A3_r .* D_delta_r;
            pop(i).Position.r = (X1_r + X2_r + X3_r) / 3;

            % --- psi component ---
            A1_psi = 2*a*rand(VarSize) - a;  C1_psi = 2*rand(VarSize);
            D_alpha_psi = abs(C1_psi .* Alpha.Position.psi - pop(i).Position.psi);
            X1_psi = Alpha.Position.psi - A1_psi .* D_alpha_psi;
            A2_psi = 2*a*rand(VarSize) - a;  C2_psi = 2*rand(VarSize);
            D_beta_psi = abs(C2_psi .* Beta.Position.psi - pop(i).Position.psi);
            X2_psi = Beta.Position.psi - A2_psi .* D_beta_psi;
            A3_psi = 2*a*rand(VarSize) - a;  C3_psi = 2*rand(VarSize);
            D_delta_psi = abs(C3_psi .* Delta.Position.psi - pop(i).Position.psi);
            X3_psi = Delta.Position.psi - A3_psi .* D_delta_psi;
            pop(i).Position.psi = (X1_psi + X2_psi + X3_psi) / 3;

            % --- phi component ---
            A1_phi = 2*a*rand(VarSize) - a;  C1_phi = 2*rand(VarSize);
            D_alpha_phi = abs(C1_phi .* Alpha.Position.phi - pop(i).Position.phi);
            X1_phi = Alpha.Position.phi - A1_phi .* D_alpha_phi;
            A2_phi = 2*a*rand(VarSize) - a;  C2_phi = 2*rand(VarSize);
            D_beta_phi = abs(C2_phi .* Beta.Position.phi - pop(i).Position.phi);
            X2_phi = Beta.Position.phi - A2_phi .* D_beta_phi;
            A3_phi = 2*a*rand(VarSize) - a;  C3_phi = 2*rand(VarSize);
            D_delta_phi = abs(C3_phi .* Delta.Position.phi - pop(i).Position.phi);
            X3_phi = Delta.Position.phi - A3_phi .* D_delta_phi;
            pop(i).Position.phi = (X1_phi + X2_phi + X3_phi) / 3;

            % Bounds
            pop(i).Position.r = max(pop(i).Position.r, VarMin.r);
            pop(i).Position.r = min(pop(i).Position.r, VarMax.r);
            pop(i).Position.psi = max(pop(i).Position.psi, VarMin.psi);
            pop(i).Position.psi = min(pop(i).Position.psi, VarMax.psi);
            pop(i).Position.phi = max(pop(i).Position.phi, VarMin.phi);
            pop(i).Position.phi = min(pop(i).Position.phi, VarMax.phi);

            % Evaluation
            cartPos = SphericalToCart(pop(i).Position, model);
            if any(isnan(cartPos.x)) || any(isnan(cartPos.y)) || any(isnan(cartPos.z)) || ...
               any(isinf(cartPos.x)) || any(isinf(cartPos.y)) || any(isinf(cartPos.z))
                newCost = inf;
            else
                try; newCost = CostFunction(cartPos); catch; newCost = inf; end
            end

            pop(i).Cost = newCost;
            if newCost < pop(i).pBest.Cost
                pop(i).pBest.Position = pop(i).Position;
                pop(i).pBest.Cost = newCost;
            end
        end

        % Re-rank leaders by pBest
        Alpha.Cost = inf; Beta.Cost = inf; Delta.Cost = inf;
        for i=1:nPop
            if pop(i).pBest.Cost < Alpha.Cost
                Delta = Beta; Beta = Alpha;
                Alpha.Position = pop(i).pBest.Position; Alpha.Cost = pop(i).pBest.Cost;
            elseif pop(i).pBest.Cost < Beta.Cost
                Delta = Beta;
                Beta.Position = pop(i).pBest.Position; Beta.Cost = pop(i).pBest.Cost;
            elseif pop(i).pBest.Cost < Delta.Cost
                Delta.Position = pop(i).pBest.Position; Delta.Cost = pop(i).pBest.Cost;
            end
        end
        if isempty(Beta.Position); Beta = Alpha; end
        if isempty(Delta.Position); Delta = Alpha; end
    end

    finalCosts(run) = Alpha.Cost;
    bestCosts(run,:) = BestCostRun;
    times(run) = toc;
    fprintf('AGWOv3 Run %2d/%d: Best=%8.2f, Time=%.1fs\n', run, N_RUNS, Alpha.Cost, times(run));
end

%% Statistics
fprintf('\n========== AGWO v3 BATCH RESULTS (N_RUNS=%d, MaxIt=%d) ==========\n', N_RUNS, MaxIt);
fprintf('Best:    %12.2f\n', min(finalCosts));
fprintf('Worst:   %12.2f\n', max(finalCosts));
fprintf('Mean:    %12.2f\n', mean(finalCosts));
fprintf('Std:     %12.2f\n', std(finalCosts));
fprintf('Avg Time: %11.1f s\n', mean(times));

fprintf('\n=== Comparison with baselines ===\n');
fprintf('SPSO  (500p, 200it): Mean ~ 4,802\n');
fprintf('GWO   (150p, 200it): Mean ~ 4,932\n');
fprintf('AGWO  (150p, 200it): Mean ~ 6,372\n');
fprintf('AGWOv3(150p, 200it): Mean ~ %.0f\n', mean(finalCosts));

% Convergence plot
figure(1); clf; hold on;
meanCurve = mean(bestCosts, 1);
plot(1:MaxIt, meanCurve, 'LineWidth', 2, 'Color', [0.2 0.6 0.2]);
xlabel('Iteration'); ylabel('Mean Best Cost');
title(sprintf('AGWO v3 Mean Convergence (%d runs)', N_RUNS));
grid on;

% Boxplot
figure(2); clf;
boxplot(finalCosts');
ylabel('Final Best Cost');
title(sprintf('AGWO v3 Distribution of Final Costs (%d runs)', N_RUNS));
grid on;

save('results/AGWOv3_batch_results.mat', 'finalCosts', 'bestCosts', 'times', 'N_RUNS', 'MaxIt');
fprintf('\nResults saved to results/AGWOv3_batch_results.mat\n');
