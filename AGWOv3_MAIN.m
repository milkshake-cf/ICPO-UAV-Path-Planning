%_________________________________________________________________________________%
%  AGWO v3: GWO-Style Three-Leader Hierarchy for UAV 3D Path Planning             %
%                                                                                 %
%  Improvements over AGWO:                                                        %
%  1. GWO-style 3-leader hierarchy (Alpha/Beta/Delta) replaces 4 defense strategies%
%  2. Adaptive nonlinear a-decay (AGWO's schedule mapped to [0.6, 2.0])           %
%  3. pBest memory retained for stable leader ranking                             %
%  4. Always accept new position (no greedy rejection)                            %
%_________________________________________________________________________________%

clc;
clear;
close all;

%% Problem Definition (SPSO framework)

model = CreateModel();
CostFunction=@(x) MyCost(x,model);
nVar=model.n;
VarSize=[1 nVar];

% Lower and upper Bounds
VarMin.x=model.xmin; VarMax.x=model.xmax;
VarMin.y=model.ymin; VarMax.y=model.ymax;
VarMin.z=model.zmin; VarMax.z=model.zmax;
VarMax.r=2*norm(model.start-model.end)/nVar; VarMin.r=0;
AngleRange = pi/4;
VarMin.psi=-AngleRange; VarMax.psi=AngleRange;
dirVector = model.end - model.start;
phi0 = atan2(dirVector(2),dirVector(1));
VarMin.phi=phi0 - AngleRange; VarMax.phi=phi0 + AngleRange;

%% AGWO v3 Parameters

MaxIt=200;          % Maximum Iterations (match SPSO/GWO for fair comparison)
nPop=150;           % Population Size

%% Initialization

empty_agent.Position=[];
empty_agent.Cost=[];
empty_agent.pBest.Position=[];   % Personal best position
empty_agent.pBest.Cost=[];       % Personal best cost

% Three-leader hierarchy (GWO-inspired, ranked by pBest cost)
Alpha.Cost=inf;   Alpha.Position=[];
Beta.Cost=inf;    Beta.Position=[];
Delta.Cost=inf;   Delta.Position=[];

% Keep GlobalBest for final reporting (synced with Alpha)
GlobalBest.Cost=inf;
GlobalBest.Position=[];

pop=repmat(empty_agent,nPop,1);

% Initialization with retry (same pattern as SPSO/GWO)
isInit = false;
while (~isInit)
    disp('Initializing AGWO v3 population...');
    for i=1:nPop
        pop(i).Position=CreateRandomSolution(VarSize,VarMin,VarMax);
        cartPos = SphericalToCart(pop(i).Position, model);
        if any(isnan(cartPos.x)) || any(isnan(cartPos.y)) || any(isnan(cartPos.z))
            pop(i).Cost = inf;
        else
            try; pop(i).Cost = CostFunction(cartPos); catch; pop(i).Cost = inf; end
        end

        % Initialize personal best
        pop(i).pBest.Position = pop(i).Position;
        pop(i).pBest.Cost = pop(i).Cost;

        % Rank into Alpha/Beta/Delta by pBest cost
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

% Safety: ensure Beta and Delta are valid (fallback to Alpha)
if isempty(Beta.Position); Beta = Alpha; end
if isempty(Delta.Position); Delta = Alpha; end

BestCost=zeros(MaxIt,1);

%% AGWO v3 Main Loop

disp('Starting AGWO v3 optimization...');

for t=1:MaxIt
    BestCost(t)=Alpha.Cost;

    % === Adaptive a parameter: AGWO's nonlinear schedule, scaled to [0.6, 2.0] ===
    % Early: a ~ 2.0 (strong exploration, |A| can be > 1)
    % Late:  a ~ 0.6 (mild exploitation, |A| always < 1)
    % Note: a never reaches 0, preserving some exploration via C vector
    a = 2 * (0.7 * (1 - t/MaxIt)^0.5 + 0.3);

    for i=1:nPop
        % ==================================================================
        % GWO-Style 3-Leader Position Update for r, psi, phi components
        % ==================================================================

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

        % === Enforce Bounds (clip) ===
        pop(i).Position.r = max(pop(i).Position.r, VarMin.r);
        pop(i).Position.r = min(pop(i).Position.r, VarMax.r);
        pop(i).Position.psi = max(pop(i).Position.psi, VarMin.psi);
        pop(i).Position.psi = min(pop(i).Position.psi, VarMax.psi);
        pop(i).Position.phi = max(pop(i).Position.phi, VarMin.phi);
        pop(i).Position.phi = min(pop(i).Position.phi, VarMax.phi);

        % === Evaluation ===
        cartPos = SphericalToCart(pop(i).Position, model);
        if any(isnan(cartPos.x)) || any(isnan(cartPos.y)) || any(isnan(cartPos.z)) || ...
           any(isinf(cartPos.x)) || any(isinf(cartPos.y)) || any(isinf(cartPos.z))
            newCost = inf;
        else
            try; newCost = CostFunction(cartPos); catch; newCost = inf; end
        end

        % === Always accept new position (like GWO/SPSO, NO greedy rejection) ===
        pop(i).Cost = newCost;

        % Update personal best
        if newCost < pop(i).pBest.Cost
            pop(i).pBest.Position = pop(i).Position;
            pop(i).pBest.Cost = newCost;
        end
    end

    % === Re-rank Alpha/Beta/Delta by pBest cost (stable hierarchy) ===
    % Reset leaders
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

    % Fallback
    if isempty(Beta.Position); Beta = Alpha; end
    if isempty(Delta.Position); Delta = Alpha; end

    GlobalBest = Alpha;

    if mod(t, 10) == 0 || t == 1
        disp(['Iteration ' num2str(t) ': Best Cost = ' num2str(BestCost(t))]);
    end
end

%% Plot Results
disp('==================== AGWO v3 OPTIMIZATION COMPLETE ====================');
disp(['Final Best Cost (AGWO v3) = ' num2str(GlobalBest.Cost)]);

BestPosition = SphericalToCart(GlobalBest.Position, model);
disp('Best path waypoints:');
disp(BestPosition);

smooth = 0.95;
PlotSolution(BestPosition, model, smooth);
title('AGWO v3: GWO-Style 3-Leader Hierarchy - UAV Path Planning');

figure(2);
plot(BestCost, 'LineWidth', 2, 'Color', [0.2 0.6 0.2]);
xlabel('Iteration');
ylabel('Best Cost');
title('AGWO v3 Convergence Curve');
grid on;

save('results/AGWOv3_results.mat', 'BestCost', 'GlobalBest', 'BestPosition');
disp('Results saved to results/AGWOv3_results.mat');
