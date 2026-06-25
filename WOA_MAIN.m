%_________________________________________________________________________________%
%  Whale Optimization Algorithm (WOA) for UAV 3D Path Planning                    %
%  Adapted to SPSO framework with spherical vector encoding                       %
%                                                                                 %
%  WOA Reference:                                                                 %
%    Mirjalili, S., Lewis, A.                                                      %
%    "The Whale Optimization Algorithm", Advances in Engineering Software, 2016   %
%_________________________________________________________________________________%

clc;
clear;
close all;

%% Problem Definition (SPSO framework)

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

%% WOA Parameters

MaxIt=200;          % Maximum Number of Iterations
nPop=150;           % Population Size
b=1;                % Spiral shape constant

%% Initialization

empty_whale.Position=[];
empty_whale.Cost=[];

GlobalBest.Cost=inf;
GlobalBest.Position=[];

pop=repmat(empty_whale,nPop,1);

% Initialization with retry
isInit = false;
while (~isInit)
    disp('Initializing WOA population...');
    for i=1:nPop
        pop(i).Position=CreateRandomSolution(VarSize,VarMin,VarMax);
        cartPos = SphericalToCart(pop(i).Position, model);
        if any(isnan(cartPos.x))||any(isnan(cartPos.y))||any(isnan(cartPos.z))
            pop(i).Cost=inf;
        else
            try; pop(i).Cost=CostFunction(cartPos); catch; pop(i).Cost=inf; end
        end
        
        if pop(i).Cost < GlobalBest.Cost
            GlobalBest.Position=pop(i).Position;
            GlobalBest.Cost=pop(i).Cost;
            isInit = true;
        end
    end
end

BestCost=zeros(MaxIt,1);

%% WOA Main Loop

disp('Starting WOA optimization...');

for t=1:MaxIt
    BestCost(t)=GlobalBest.Cost;
    
    % Linearly decrease a from 2 to 0
    a = 2 - t*(2/MaxIt);
    a2 = -1 + t*(-1/MaxIt);  % a2 linearly decreases from -1 to -2
    
    for i=1:nPop
        
        % Random numbers for each component
        r1 = rand(); r2 = rand();
        A = 2*a*r1 - a;
        C = 2*r2;
        l = (a2-1)*rand() + 1;  % l in [-2, 1]
        p = rand();
        
        % === Separate updates for r, psi, phi ===
        
        if p < 0.5
            if abs(A) >= 1
                % Exploration: search from random whale
                k = randi(nPop);
                
                % r component
                D_r = abs(C*pop(k).Position.r - pop(i).Position.r);
                pop(i).Position.r = pop(k).Position.r - A*D_r;
                
                % psi component
                D_psi = abs(C*pop(k).Position.psi - pop(i).Position.psi);
                pop(i).Position.psi = pop(k).Position.psi - A*D_psi;
                
                % phi component
                D_phi = abs(C*pop(k).Position.phi - pop(i).Position.phi);
                pop(i).Position.phi = pop(k).Position.phi - A*D_phi;
                
            else
                % Exploitation: encircling prey
                % r component
                D_r = abs(C*GlobalBest.Position.r - pop(i).Position.r);
                pop(i).Position.r = GlobalBest.Position.r - A*D_r;
                
                % psi component
                D_psi = abs(C*GlobalBest.Position.psi - pop(i).Position.psi);
                pop(i).Position.psi = GlobalBest.Position.psi - A*D_psi;
                
                % phi component
                D_phi = abs(C*GlobalBest.Position.phi - pop(i).Position.phi);
                pop(i).Position.phi = GlobalBest.Position.phi - A*D_phi;
            end
        else
            % Spiral update (bubble-net attack)
            % r component
            D_r = abs(GlobalBest.Position.r - pop(i).Position.r);
            pop(i).Position.r = D_r * exp(b*l) .* cos(2*pi*l) + GlobalBest.Position.r;
            
            % psi component
            D_psi = abs(GlobalBest.Position.psi - pop(i).Position.psi);
            pop(i).Position.psi = D_psi * exp(b*l) .* cos(2*pi*l) + GlobalBest.Position.psi;
            
            % phi component
            D_phi = abs(GlobalBest.Position.phi - pop(i).Position.phi);
            pop(i).Position.phi = D_phi * exp(b*l) .* cos(2*pi*l) + GlobalBest.Position.phi;
        end
        
        % Enforce bounds
        pop(i).Position.r = max(pop(i).Position.r, VarMin.r); pop(i).Position.r = min(pop(i).Position.r, VarMax.r);
        pop(i).Position.psi = max(pop(i).Position.psi, VarMin.psi); pop(i).Position.psi = min(pop(i).Position.psi, VarMax.psi);
        pop(i).Position.phi = max(pop(i).Position.phi, VarMin.phi); pop(i).Position.phi = min(pop(i).Position.phi, VarMax.phi);
        
        % Evaluation
        cartPos = SphericalToCart(pop(i).Position, model);
        if any(isnan(cartPos.x))||any(isnan(cartPos.y))||any(isnan(cartPos.z))
            pop(i).Cost = inf;
        else
            try; pop(i).Cost = CostFunction(cartPos); catch; pop(i).Cost = inf; end
        end
        
        % Update Global Best
        if pop(i).Cost < GlobalBest.Cost
            GlobalBest.Position = pop(i).Position;
            GlobalBest.Cost = pop(i).Cost;
        end
    end
    
    if mod(t,10)==0 || t==1
        disp(['Iteration ' num2str(t) ': Best Cost = ' num2str(BestCost(t))]);
    end
end

%% Plot Results
disp('==================== WOA OPTIMIZATION COMPLETE ====================');
disp(['Final Best Cost (WOA) = ' num2str(GlobalBest.Cost)]);

BestPosition = SphericalToCart(GlobalBest.Position, model);
smooth = 0.95;
PlotSolution(BestPosition, model, smooth);
title('WOA: Whale Optimization Algorithm - UAV Path Planning');

figure(2);
plot(BestCost, 'LineWidth', 2, 'Color', [0.2 0.6 0.6]);
xlabel('Iteration'); ylabel('Best Cost');
title('WOA Convergence Curve'); grid on;

save('results/WOA_results.mat', 'BestCost', 'GlobalBest', 'BestPosition');
disp('Results saved to results/WOA_results.mat');
