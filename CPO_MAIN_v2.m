%_________________________________________________________________________________%
%  Crested Porcupine Optimizer (CPO) for UAV 3D Path Planning                     %
%  Based on original CPO by Abdel-Basset & Mohamed (2024)                         %
%  Embedded in SPSO framework by Manh Duong Phung                                 %
%                                                                                 %
%  CPO Reference:                                                                 %
%    Abdel-Basset, M., Mohamed, R.                                                %
%    "Crested Porcupine Optimizer: A new nature-inspired metaheuristic"           %
%    Knowledge-Based Systems, 2024                                                %
%_________________________________________________________________________________%

clc;
clear;
close all;

%% Problem Definition (SPSO framework)

model = CreateModel(); % Create search map and parameters

CostFunction=@(x) MyCost(x,model);    % Cost Function

nVar=model.n;       % Number of Decision Variables = number of path nodes

VarSize=[1 nVar];   % Size of Decision Variables Matrix

% Lower and upper Bounds (spherical vector components)
VarMin.x=model.xmin;
VarMax.x=model.xmax;
VarMin.y=model.ymin;
VarMax.y=model.ymax;
VarMin.z=model.zmin;
VarMax.z=model.zmax;

VarMax.r=2*norm(model.start-model.end)/nVar;
VarMin.r=0;

% Inclination (elevation)
AngleRange = pi/4;
VarMin.psi=-AngleRange;
VarMax.psi=AngleRange;

% Azimuth
dirVector = model.end - model.start;
phi0 = atan2(dirVector(2),dirVector(1));
VarMin.phi=phi0 - AngleRange;
VarMax.phi=phi0 + AngleRange;

%% CPO Parameters

MaxIt=500;          % Maximum Number of Iterations
nPop=150;           % Population Size

% CPO-specific parameters (from original paper)
alpha = 0.2;        % Convergence rate for 4th defense mechanism
Tf = 0.8;           % Tradeoff between 3rd and 4th defense mechanisms

%% Initialization

% Create empty CPO agent structure
empty_agent.Position=[];
empty_agent.Cost=[];

% Initialize Global Best
GlobalBest.Cost=inf;
GlobalBest.Position=[];

% Create population
pop=repmat(empty_agent,nPop,1);

% Track previous positions for greedy selection
prev_positions = cell(nPop,1);

% Initialization Loop (with retry, same pattern as SPSO)
isInit = false;
while (~isInit)
    disp('Initializing CPO population...');
    for i=1:nPop
        % Initialize Position (same as SPSO)
        pop(i).Position=CreateRandomSolution(VarSize,VarMin,VarMax);
        
        % Evaluation (with safety)
        cartPos = SphericalToCart(pop(i).Position, model);
        if any(isnan(cartPos.x)) || any(isnan(cartPos.y)) || any(isnan(cartPos.z))
            pop(i).Cost = inf;
        else
            try
                pop(i).Cost = CostFunction(cartPos);
            catch
                pop(i).Cost = inf;
            end
        end
        
        % Save previous position
        prev_positions{i} = pop(i).Position;
        
        % Update Global Best
        if pop(i).Cost < GlobalBest.Cost
            GlobalBest.Position=pop(i).Position;
            GlobalBest.Cost=pop(i).Cost;
            isInit = true;
        end
    end
end

% Array to Hold Best Cost Values at Each Iteration
BestCost=zeros(MaxIt,1);

%% CPO Main Loop

disp('Starting CPO optimization...');

for t=1:MaxIt
    
    % Update Best Cost
    BestCost(t)=GlobalBest.Cost;
    
    for i=1:nPop
        
        % ===== Generate random vectors for position mixing =====
        U1 = rand(VarSize) > rand();  % Binary mask
        
        % ===== Exploration vs Exploitation decision =====
        if rand() < rand()  % EXPLORATION phase
            
            if rand() < rand()  % Strategy 1: Visual Deterrence
                % y_t = (X_i + X_k) / 2  where k is a random agent
                k = randi(nPop);
                
                % r component
                y_r = (pop(i).Position.r + pop(k).Position.r) / 2;
                pop(i).Position.r = pop(i).Position.r + randn(VarSize) .* abs(2*rand()*GlobalBest.Position.r - y_r);
                
                % psi component
                y_psi = (pop(i).Position.psi + pop(k).Position.psi) / 2;
                pop(i).Position.psi = pop(i).Position.psi + randn(VarSize) .* abs(2*rand()*GlobalBest.Position.psi - y_psi);
                
                % phi component
                y_phi = (pop(i).Position.phi + pop(k).Position.phi) / 2;
                pop(i).Position.phi = pop(i).Position.phi + randn(VarSize) .* abs(2*rand()*GlobalBest.Position.phi - y_phi);
                
            else  % Strategy 2: Sound Deterrence
                k = randi(nPop);
                m = randi(nPop);
                
                % r component
                y_r = (pop(i).Position.r + pop(k).Position.r) / 2;
                pop(i).Position.r = U1 .* pop(i).Position.r + (1-U1) .* (y_r + rand()*(pop(m).Position.r - pop(k).Position.r));
                
                % psi component
                y_psi = (pop(i).Position.psi + pop(k).Position.psi) / 2;
                pop(i).Position.psi = U1 .* pop(i).Position.psi + (1-U1) .* (y_psi + rand()*(pop(m).Position.psi - pop(k).Position.psi));
                
                % phi component
                y_phi = (pop(i).Position.phi + pop(k).Position.phi) / 2;
                pop(i).Position.phi = U1 .* pop(i).Position.phi + (1-U1) .* (y_phi + rand()*(pop(m).Position.phi - pop(k).Position.phi));
            end
            
        else  % EXPLOITATION phase
            
            % Cycle parameter Yt
            Yt = 2 * rand() * (1 - t/MaxIt)^(t/MaxIt);
            
            % U2: random direction vector with values -1, 0, or +1
            % NOTE: Fixed from original paper (Python version correction)
            U2 = (rand(VarSize) < 0.5) * 2 - 1;
            S = rand() * U2;
            
            % Sum of fitness for weighting (avoid division by zero)
            allCosts = [pop.Cost];
            sumFitness = sum(allCosts) + eps;
            
            if rand() < Tf  % Strategy 3: Physical Attack
                St = exp(pop(i).Cost / sumFitness);
                S = S .* Yt .* St;
                
                k = randi(nPop);
                m = randi(nPop);
                
                % r component
                pop(i).Position.r = (1-U1).*pop(i).Position.r + U1.*(pop(k).Position.r + St*(pop(m).Position.r - pop(k).Position.r) - S);
                
                % psi component
                pop(i).Position.psi = (1-U1).*pop(i).Position.psi + U1.*(pop(k).Position.psi + St*(pop(m).Position.psi - pop(k).Position.psi) - S);
                
                % phi component
                pop(i).Position.phi = (1-U1).*pop(i).Position.phi + U1.*(pop(k).Position.phi + St*(pop(m).Position.phi - pop(k).Position.phi) - S);
                
            else  % Strategy 4: Lethal Attack
                Mt = exp(pop(i).Cost / sumFitness);
                
                k = randi(nPop);
                vt = pop(i).Position;
                Vtp = pop(k).Position;
                r2_param = rand();
                
                % r component
                Ft_r = rand(VarSize) .* (Mt * (-vt.r + Vtp.r));
                S_r = S .* Yt .* Ft_r;
                pop(i).Position.r = GlobalBest.Position.r + (alpha*(1-r2_param)+r2_param)*(U2.*GlobalBest.Position.r - pop(i).Position.r) - S_r;
                
                % psi component
                Ft_psi = rand(VarSize) .* (Mt * (-vt.psi + Vtp.psi));
                S_psi = S .* Yt .* Ft_psi;
                pop(i).Position.psi = GlobalBest.Position.psi + (alpha*(1-r2_param)+r2_param)*(U2.*GlobalBest.Position.psi - pop(i).Position.psi) - S_psi;
                
                % phi component
                Ft_phi = rand(VarSize) .* (Mt * (-vt.phi + Vtp.phi));
                S_phi = S .* Yt .* Ft_phi;
                pop(i).Position.phi = GlobalBest.Position.phi + (alpha*(1-r2_param)+r2_param)*(U2.*GlobalBest.Position.phi - pop(i).Position.phi) - S_phi;
            end
        end
        
        % ===== Enforce Bounds (SPSO-style: clip to bounds) =====
        % r bounds
        pop(i).Position.r = max(pop(i).Position.r, VarMin.r);
        pop(i).Position.r = min(pop(i).Position.r, VarMax.r);
        % psi bounds
        pop(i).Position.psi = max(pop(i).Position.psi, VarMin.psi);
        pop(i).Position.psi = min(pop(i).Position.psi, VarMax.psi);
        % phi bounds
        pop(i).Position.phi = max(pop(i).Position.phi, VarMin.phi);
        pop(i).Position.phi = min(pop(i).Position.phi, VarMax.phi);
        
        % ===== Evaluation =====
        cartPos = SphericalToCart(pop(i).Position, model);
        % Safety: reject solutions with NaN/Inf or out-of-bounds
        if any(isnan(cartPos.x)) || any(isnan(cartPos.y)) || any(isnan(cartPos.z)) || ...
           any(isinf(cartPos.x)) || any(isinf(cartPos.y)) || any(isinf(cartPos.z))
            newCost = inf;
        else
            try
                newCost = CostFunction(cartPos);
            catch
                newCost = inf;
            end
        end
        
        % ===== Greedy Selection =====
        if pop(i).Cost < newCost
            % Keep old position if new is worse
            pop(i).Position = prev_positions{i};
        else
            % Accept new position
            prev_positions{i} = pop(i).Position;
            pop(i).Cost = newCost;
            
            % Update Global Best
            if pop(i).Cost < GlobalBest.Cost
                GlobalBest.Position = pop(i).Position;
                GlobalBest.Cost = pop(i).Cost;
            end
        end
    end
    
    % Show Iteration Information
    if mod(t, 10) == 0 || t == 1
        disp(['Iteration ' num2str(t) ': Best Cost = ' num2str(BestCost(t))]);
    end
end

%% Plot Results

disp('==================== CPO OPTIMIZATION COMPLETE ====================');
disp(['Final Best Cost (CPO) = ' num2str(GlobalBest.Cost)]);

% Best solution
BestPosition = SphericalToCart(GlobalBest.Position, model);
disp('Best path waypoints:');
disp(BestPosition);

smooth = 0.95;
PlotSolution(BestPosition, model, smooth);
title('CPO: Crested Porcupine Optimizer - UAV Path Planning');

% Convergence curve
figure(2);
plot(BestCost, 'LineWidth', 2, 'Color', [0.2 0.4 0.8]);
xlabel('Iteration');
ylabel('Best Cost');
title('CPO Convergence Curve');
grid on;

% Save results
save('results/CPO_results_v2.mat', 'BestCost', 'GlobalBest', 'BestPosition');
disp('Results saved to results/CPO_results_v2.mat');
