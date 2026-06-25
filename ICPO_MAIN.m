%_________________________________________________________________________________%
%  Improved Crested Porcupine Optimizer (ICPO) for UAV 3D Path Planning           %
%                                                                                 %
%  Improvements over original CPO:                                                %
%  1. Personal Best memory (pBest) - each agent remembers its best position       %
%  2. Adaptive exploration/exploitation balance                                   %
%  3. Elite-guided position updates using pBest                                   %
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

%% ICPO Parameters

MaxIt=500;          % Maximum Iterations (extended for convergence test)
nPop=150;           % Population Size

alpha = 0.2;        % Convergence rate (from CPO paper)
Tf = 0.8;           % Tradeoff 3rd/4th defense

%% Initialization

empty_agent.Position=[];
empty_agent.Cost=[];
empty_agent.pBest.Position=[];   % Personal best position
empty_agent.pBest.Cost=[];       % Personal best cost

GlobalBest.Cost=inf;
GlobalBest.Position=[];

pop=repmat(empty_agent,nPop,1);
prev_positions = cell(nPop,1);

% Initialization with retry
isInit = false;
while (~isInit)
    disp('Initializing ICPO population...');
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
        
        if pop(i).Cost < GlobalBest.Cost
            GlobalBest.Position=pop(i).Position;
            GlobalBest.Cost=pop(i).Cost;
            isInit = true;
        end
    end
end

BestCost=zeros(MaxIt,1);

%% ICPO Main Loop

disp('Starting ICPO optimization...');

for t=1:MaxIt
    BestCost(t)=GlobalBest.Cost;
    
    % === Adaptive parameter: shifts from exploration to exploitation ===
    % Early iterations: explore more (closer to 1)
    % Late iterations: exploit more (closer to 0)
    explRatio = 0.7 * (1 - t/MaxIt)^0.5 + 0.3;
    
    for i=1:nPop
        U1 = rand(VarSize) > rand();
        
        % === Adaptive strategy selection (not purely random) ===
        if rand() < explRatio  % EXPLORATION phase (adaptive probability)
            
            if rand() < rand()  % Strategy 1: Visual Deterrence
                % Use pBest of two random agents for guided exploration
                k = randi(nPop);
                m = randi(nPop);
                
                % Blend pBest positions for smarter exploration
                % r component
                guide_r = (pop(k).pBest.Position.r + pop(m).pBest.Position.r) / 2;
                pop(i).Position.r = pop(i).Position.r + randn(VarSize) .* abs(2*rand()*GlobalBest.Position.r - guide_r);
                
                % psi component
                guide_psi = (pop(k).pBest.Position.psi + pop(m).pBest.Position.psi) / 2;
                pop(i).Position.psi = pop(i).Position.psi + randn(VarSize) .* abs(2*rand()*GlobalBest.Position.psi - guide_psi);
                
                % phi component
                guide_phi = (pop(k).pBest.Position.phi + pop(m).pBest.Position.phi) / 2;
                pop(i).Position.phi = pop(i).Position.phi + randn(VarSize) .* abs(2*rand()*GlobalBest.Position.phi - guide_phi);
                
            else  % Strategy 2: Sound Deterrence
                k = randi(nPop);
                m = randi(nPop);
                
                % Use pBest instead of current position for difference
                % r
                y_r = (pop(i).Position.r + pop(k).pBest.Position.r) / 2;
                diff_r = pop(m).pBest.Position.r - pop(k).pBest.Position.r;
                pop(i).Position.r = U1.*pop(i).Position.r + (1-U1).*(y_r + rand()*diff_r);
                
                % psi
                y_psi = (pop(i).Position.psi + pop(k).pBest.Position.psi) / 2;
                diff_psi = pop(m).pBest.Position.psi - pop(k).pBest.Position.psi;
                pop(i).Position.psi = U1.*pop(i).Position.psi + (1-U1).*(y_psi + rand()*diff_psi);
                
                % phi
                y_phi = (pop(i).Position.phi + pop(k).pBest.Position.phi) / 2;
                diff_phi = pop(m).pBest.Position.phi - pop(k).pBest.Position.phi;
                pop(i).Position.phi = U1.*pop(i).Position.phi + (1-U1).*(y_phi + rand()*diff_phi);
            end
            
        else  % EXPLOITATION phase
            Yt = 2 * rand() * (1 - t/MaxIt)^(t/MaxIt);
            U2 = (rand(VarSize) < 0.5) * 2 - 1;
            S = rand() * U2;
            
            allCosts = zeros(1, nPop);
            for j=1:nPop; allCosts(j) = pop(j).pBest.Cost; end
            sumFitness = sum(allCosts) + eps;
            
            if rand() < Tf  % Strategy 3: Physical Attack
                St = exp(pop(i).pBest.Cost / sumFitness);  % Use pBest cost
                S = S .* Yt .* St;
                
                k = randi(nPop);
                m = randi(nPop);
                
                % r: use pBest for reference
                pop(i).Position.r = (1-U1).*pop(i).Position.r + U1.*(pop(k).pBest.Position.r + St*(pop(m).pBest.Position.r - pop(k).pBest.Position.r) - S);
                % psi
                pop(i).Position.psi = (1-U1).*pop(i).Position.psi + U1.*(pop(k).pBest.Position.psi + St*(pop(m).pBest.Position.psi - pop(k).pBest.Position.psi) - S);
                % phi
                pop(i).Position.phi = (1-U1).*pop(i).Position.phi + U1.*(pop(k).pBest.Position.phi + St*(pop(m).pBest.Position.phi - pop(k).pBest.Position.phi) - S);
                
            else  % Strategy 4: Lethal Attack
                Mt = exp(pop(i).pBest.Cost / sumFitness);
                k = randi(nPop);
                r2_param = rand();
                
                % r
                vt_r = pop(i).Position.r;
                Vtp_r = pop(k).pBest.Position.r;
                Ft_r = rand(VarSize) .* (Mt * (-vt_r + Vtp_r));
                S_r = S .* Yt .* Ft_r;
                pop(i).Position.r = GlobalBest.Position.r + (alpha*(1-r2_param)+r2_param)*(U2.*GlobalBest.Position.r - pop(i).Position.r) - S_r;
                
                % psi
                vt_psi = pop(i).Position.psi;
                Vtp_psi = pop(k).pBest.Position.psi;
                Ft_psi = rand(VarSize) .* (Mt * (-vt_psi + Vtp_psi));
                S_psi = S .* Yt .* Ft_psi;
                pop(i).Position.psi = GlobalBest.Position.psi + (alpha*(1-r2_param)+r2_param)*(U2.*GlobalBest.Position.psi - pop(i).Position.psi) - S_psi;
                
                % phi
                vt_phi = pop(i).Position.phi;
                Vtp_phi = pop(k).pBest.Position.phi;
                Ft_phi = rand(VarSize) .* (Mt * (-vt_phi + Vtp_phi));
                S_phi = S .* Yt .* Ft_phi;
                pop(i).Position.phi = GlobalBest.Position.phi + (alpha*(1-r2_param)+r2_param)*(U2.*GlobalBest.Position.phi - pop(i).Position.phi) - S_phi;
            end
        end
        
        % === Enforce Bounds (clip to bounds) ===
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
        
        % === Greedy Selection with pBest update ===
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
            
            % Update Global Best
            if newCost < GlobalBest.Cost
                GlobalBest.Position = pop(i).Position;
                GlobalBest.Cost = newCost;
            end
        end
    end
    
    if mod(t, 10) == 0 || t == 1
        disp(['Iteration ' num2str(t) ': Best Cost = ' num2str(BestCost(t))]);
    end
end

%% Plot Results
disp('==================== ICPO OPTIMIZATION COMPLETE ====================');
disp(['Final Best Cost (ICPO) = ' num2str(GlobalBest.Cost)]);

BestPosition = SphericalToCart(GlobalBest.Position, model);
disp('Best path waypoints:');
disp(BestPosition);

smooth = 0.95;
PlotSolution(BestPosition, model, smooth);
title('ICPO: Improved Crested Porcupine Optimizer - UAV Path Planning');

figure(2);
plot(BestCost, 'LineWidth', 2, 'Color', [0.8 0.2 0.2]);
xlabel('Iteration');
ylabel('Best Cost');
title('ICPO Convergence Curve');
grid on;

save('results/ICPO_results.mat', 'BestCost', 'GlobalBest', 'BestPosition');
disp('Results saved to results/ICPO_results.mat');
