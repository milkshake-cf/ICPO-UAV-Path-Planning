%_________________________________________________________________________________%
%  AGWO v6: Multi-Elite Guided CPO for UAV 3D Path Planning                       %
%                                                                                 %
%  Core innovation: Multi-Elite Guidance within CPO's 4 defense strategies        %
%  - Retains ALL 4 CPO defense strategies (visual/sound/physical/lethal)          %
%  - Introduces Alpha/Beta/Delta elite hierarchy for diversified guidance         %
%  - Elite-guided difference vectors replace random agent selection               %
%  - Rank-based fitness weighting for stability                                   %
%  - Retains pBest memory, adaptive explRatio, greedy selection, boundary clipping %
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

%% AGWO v6 Parameters

MaxIt=200;          % Maximum Iterations
nPop=150;           % Population Size

alpha_param = 0.2;  % Convergence rate (from CPO paper, used in Strategy 4)
Tf = 0.8;           % Tradeoff 3rd/4th defense

%% Initialization

empty_agent.Position=[];
empty_agent.Cost=[];
empty_agent.pBest.Position=[];   % Personal best position
empty_agent.pBest.Cost=[];       % Personal best cost

% === Multi-Elite Hierarchy (NEW) ===
Alpha.Cost=inf;   Alpha.Position=[];
Beta.Cost=inf;    Beta.Position=[];
Delta.Cost=inf;   Delta.Position=[];

% Keep GlobalBest for backward compatibility (= Alpha)
GlobalBest.Cost=inf;
GlobalBest.Position=[];

pop=repmat(empty_agent,nPop,1);
prev_positions = cell(nPop,1);

% Initialization with retry
isInit = false;
while (~isInit)
    disp('Initializing AGWO v6 population...');
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

        prev_positions{i} = pop(i).Position;

        % Rank into Alpha/Beta/Delta by pBest cost (NEW)
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

% Safety fallback
if isempty(Beta.Position); Beta = Alpha; end
if isempty(Delta.Position); Delta = Alpha; end

GlobalBest = Alpha;
BestCost=zeros(MaxIt,1);

%% Helper: tournament selection (pick best of 2 random agents by pBest cost)
function idx = tournament_select(nPop, pop)
    a = randi(nPop); b = randi(nPop);
    if pop(a).pBest.Cost < pop(b).pBest.Cost; idx = a; else; idx = b; end
end

%% AGWO v6 Main Loop

disp('Starting AGWO v6 optimization...');

for t=1:MaxIt
    BestCost(t)=Alpha.Cost;

    % === Adaptive exploration ratio (retained from AGWO) ===
    explRatio = 0.7 * (1 - t/MaxIt)^0.5 + 0.3;

    % === Rank-based fitness for stable weighting (NEW) ===
    % Sort agents by pBest cost, assign ranks (1 = best, nPop = worst)
    [~, sortIdx] = sort([pop.pBest.Cost]);
    ranks = zeros(1, nPop);
    for j = 1:nPop; ranks(sortIdx(j)) = j; end

    for i=1:nPop
        U1 = rand(VarSize) > rand();

        % === Select elite reference for this agent (NEW) ===
        % Each strategy uses a DIFFERENT elite, breaking the single-gBest bottleneck
        elite_choice = rand();

        if rand() < explRatio  % EXPLORATION phase

            if rand() < rand()  % Strategy 1: Visual Deterrence
                % === MULTI-ELITE: Use random elite instead of single GlobalBest ===
                if elite_choice < 0.5
                    elite = Alpha;
                elseif elite_choice < 0.8
                    elite = Beta;
                else
                    elite = Delta;
                end

                % Tournament selection for guide agents (NEW)
                k = tournament_select(nPop, pop);
                m = tournament_select(nPop, pop);

                % r component
                guide_r = (pop(k).pBest.Position.r + pop(m).pBest.Position.r) / 2;
                pop(i).Position.r = pop(i).Position.r + randn(VarSize) .* abs(2*rand()*elite.Position.r - guide_r);

                % psi component
                guide_psi = (pop(k).pBest.Position.psi + pop(m).pBest.Position.psi) / 2;
                pop(i).Position.psi = pop(i).Position.psi + randn(VarSize) .* abs(2*rand()*elite.Position.psi - guide_psi);

                % phi component
                guide_phi = (pop(k).pBest.Position.phi + pop(m).pBest.Position.phi) / 2;
                pop(i).Position.phi = pop(i).Position.phi + randn(VarSize) .* abs(2*rand()*elite.Position.phi - guide_phi);

            else  % Strategy 2: Sound Deterrence
                % === ELITE-GUIDED difference vector (NEW) ===
                % Use Beta-Delta difference instead of random agent difference
                % This gives structured exploration direction from elite diversity

                k = tournament_select(nPop, pop);

                % r: elite diversity as difference vector
                y_r = (pop(i).Position.r + pop(k).pBest.Position.r) / 2;
                diff_r = Beta.Position.r - Delta.Position.r;  % Elite diversity (NEW)
                pop(i).Position.r = U1.*pop(i).Position.r + (1-U1).*(y_r + rand()*diff_r);

                % psi
                y_psi = (pop(i).Position.psi + pop(k).pBest.Position.psi) / 2;
                diff_psi = Beta.Position.psi - Delta.Position.psi;
                pop(i).Position.psi = U1.*pop(i).Position.psi + (1-U1).*(y_psi + rand()*diff_psi);

                % phi
                y_phi = (pop(i).Position.phi + pop(k).pBest.Position.phi) / 2;
                diff_phi = Beta.Position.phi - Delta.Position.phi;
                pop(i).Position.phi = U1.*pop(i).Position.phi + (1-U1).*(y_phi + rand()*diff_phi);
            end

        else  % EXPLOITATION phase
            Yt = 2 * rand() * (1 - t/MaxIt)^(t/MaxIt);
            U2 = (rand(VarSize) < 0.5) * 2 - 1;
            S = rand() * U2;

            % === Rank-based fitness weighting (NEW, replaces exp(cost/sum)) ===
            % Normalized rank: 0 (best) to 1 (worst)
            rankNorm = ranks(i) / nPop;
            % Sigmoid transform for smooth weighting
            fitnessWeight = 1 / (1 + exp(5 * (rankNorm - 0.5)));  % ~0 for good, ~1 for bad

            if rand() < Tf  % Strategy 3: Physical Attack
                % === ELITE-GUIDED with tournament selection (NEW) ===
                St = fitnessWeight;  % Rank-based instead of exp(cost/sumFitness)
                S = S .* Yt .* St;

                k = tournament_select(nPop, pop);  % Tournament (NEW)
                m = tournament_select(nPop, pop);  % Tournament (NEW)

                % r: elite-guided physical attack
                pop(i).Position.r = (1-U1).*pop(i).Position.r + U1.*(pop(k).pBest.Position.r + St*(pop(m).pBest.Position.r - pop(k).pBest.Position.r) - S);
                % psi
                pop(i).Position.psi = (1-U1).*pop(i).Position.psi + U1.*(pop(k).pBest.Position.psi + St*(pop(m).pBest.Position.psi - pop(k).pBest.Position.psi) - S);
                % phi
                pop(i).Position.phi = (1-U1).*pop(i).Position.phi + U1.*(pop(k).pBest.Position.phi + St*(pop(m).pBest.Position.phi - pop(k).pBest.Position.phi) - S);

            else  % Strategy 4: Lethal Attack
                % === ALPHA-centered with Beta-Delta perturbation (NEW) ===
                Mt = fitnessWeight;  % Rank-based (NEW)
                k = tournament_select(nPop, pop);  % Tournament (NEW)
                r2_param = rand();

                % r: Alpha as center + elite diversity perturbation
                vt_r = pop(i).Position.r;
                Vtp_r = pop(k).pBest.Position.r;
                Ft_r = rand(VarSize) .* (Mt * (-vt_r + Vtp_r));
                S_r = S .* Yt .* Ft_r;
                % Use Alpha (not single GlobalBest) + Beta-Delta perturbation (NEW)
                elite_perturb_r = 0.1 * rand() * (Beta.Position.r - Delta.Position.r);
                pop(i).Position.r = Alpha.Position.r + elite_perturb_r + (alpha_param*(1-r2_param)+r2_param)*(U2.*Alpha.Position.r - pop(i).Position.r) - S_r;

                % psi
                vt_psi = pop(i).Position.psi;
                Vtp_psi = pop(k).pBest.Position.psi;
                Ft_psi = rand(VarSize) .* (Mt * (-vt_psi + Vtp_psi));
                S_psi = S .* Yt .* Ft_psi;
                elite_perturb_psi = 0.1 * rand() * (Beta.Position.psi - Delta.Position.psi);
                pop(i).Position.psi = Alpha.Position.psi + elite_perturb_psi + (alpha_param*(1-r2_param)+r2_param)*(U2.*Alpha.Position.psi - pop(i).Position.psi) - S_psi;

                % phi
                vt_phi = pop(i).Position.phi;
                Vtp_phi = pop(k).pBest.Position.phi;
                Ft_phi = rand(VarSize) .* (Mt * (-vt_phi + Vtp_phi));
                S_phi = S .* Yt .* Ft_phi;
                elite_perturb_phi = 0.1 * rand() * (Beta.Position.phi - Delta.Position.phi);
                pop(i).Position.phi = Alpha.Position.phi + elite_perturb_phi + (alpha_param*(1-r2_param)+r2_param)*(U2.*Alpha.Position.phi - pop(i).Position.phi) - S_phi;
            end
        end

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

        % === Greedy Selection (RETAINED - CPO identity) ===
        if pop(i).Cost < newCost
            % Keep old position
            pop(i).Position = prev_positions{i};
        else
            % Accept new position
            prev_positions{i} = pop(i).Position;
            pop(i).Cost = newCost;

            % Update personal best
            if newCost < pop(i).pBest.Cost
                pop(i).pBest.Position = pop(i).Position;
                pop(i).pBest.Cost = newCost;
            end

            % Update Global Best (= Alpha)
            if newCost < GlobalBest.Cost
                GlobalBest.Position = pop(i).Position;
                GlobalBest.Cost = newCost;
            end
        end
    end

    % === Re-rank Alpha/Beta/Delta by pBest cost (NEW) ===
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

    GlobalBest = Alpha;

    if mod(t, 10) == 0 || t == 1
        disp(['Iteration ' num2str(t) ': Best Cost = ' num2str(BestCost(t))]);
    end
end

%% Plot Results
disp('==================== AGWO v6 OPTIMIZATION COMPLETE ====================');
disp(['Final Best Cost (AGWO v6) = ' num2str(GlobalBest.Cost)]);

BestPosition = SphericalToCart(GlobalBest.Position, model);
disp('Best path waypoints:');
disp(BestPosition);

smooth = 0.95;
PlotSolution(BestPosition, model, smooth);
title('AGWO v6: Multi-Elite Guided CPO - UAV Path Planning');

figure(2);
plot(BestCost, 'LineWidth', 2, 'Color', [0.8 0.4 0.2]);
xlabel('Iteration');
ylabel('Best Cost');
title('AGWO v6 Convergence Curve');
grid on;

save('results/AGWOv6_results.mat', 'BestCost', 'GlobalBest', 'BestPosition');
disp('Results saved to results/AGWOv6_results.mat');
