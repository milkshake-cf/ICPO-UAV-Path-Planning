%_________________________________________________________________________________%
%  Grey Wolf Optimizer (GWO) for UAV 3D Path Planning                             %
%  Adapted to SPSO framework with spherical vector encoding                       %
%                                                                                 %
%  GWO Reference:                                                                 %
%    Mirjalili, S., Mirjalili, S.M., Lewis, A.                                     %
%    "Grey Wolf Optimizer", Advances in Engineering Software, 2014                %
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

%% GWO Parameters

MaxIt=200;          % Maximum Number of Iterations
nPop=150;           % Population Size (same as ICPO for fair comparison)

%% Initialization

% Create empty wolf structure
empty_wolf.Position=[];
empty_wolf.Cost=[];

% Initialize alpha, beta, delta
Alpha.Cost=inf; Alpha.Position=[];
Beta.Cost=inf;  Beta.Position=[];
Delta.Cost=inf; Delta.Position=[];

% Create population
pack=repmat(empty_wolf,nPop,1);

% Initialization with retry (same as SPSO/CPO)
isInit = false;
while (~isInit)
    disp('Initializing GWO population...');
    for i=1:nPop
        pack(i).Position=CreateRandomSolution(VarSize,VarMin,VarMax);
        cartPos = SphericalToCart(pack(i).Position, model);
        if any(isnan(cartPos.x))||any(isnan(cartPos.y))||any(isnan(cartPos.z))
            pack(i).Cost=inf;
        else
            try; pack(i).Cost=CostFunction(cartPos); catch; pack(i).Cost=inf; end
        end
        
        % Update alpha, beta, delta
        if pack(i).Cost < Alpha.Cost
            Delta = Beta; Beta = Alpha; Alpha.Position = pack(i).Position; Alpha.Cost = pack(i).Cost;
            isInit = true;
        elseif pack(i).Cost < Beta.Cost
            Delta = Beta; Beta.Position = pack(i).Position; Beta.Cost = pack(i).Cost;
        elseif pack(i).Cost < Delta.Cost
            Delta.Position = pack(i).Position; Delta.Cost = pack(i).Cost;
        end
    end
end

% Safety: ensure Beta and Delta are valid (fallback to Alpha)
if isempty(Beta.Position); Beta = Alpha; end
if isempty(Delta.Position); Delta = Alpha; end

BestCost=zeros(MaxIt,1);

%% GWO Main Loop

disp('Starting GWO optimization...');

for t=1:MaxIt
    BestCost(t)=Alpha.Cost;
    
    % Linearly decrease a from 2 to 0
    a = 2 - t*(2/MaxIt);
    
    for i=1:nPop
        % Separate updates for r, psi, phi components
        
        % === r component ===
        [X1_r, X2_r, X3_r] = GWO_update(pack(i).Position.r, Alpha.Position.r, Beta.Position.r, Delta.Position.r, a, VarSize);
        pack(i).Position.r = (X1_r + X2_r + X3_r) / 3;
        
        % === psi component ===
        [X1_psi, X2_psi, X3_psi] = GWO_update(pack(i).Position.psi, Alpha.Position.psi, Beta.Position.psi, Delta.Position.psi, a, VarSize);
        pack(i).Position.psi = (X1_psi + X2_psi + X3_psi) / 3;
        
        % === phi component ===
        [X1_phi, X2_phi, X3_phi] = GWO_update(pack(i).Position.phi, Alpha.Position.phi, Beta.Position.phi, Delta.Position.phi, a, VarSize);
        pack(i).Position.phi = (X1_phi + X2_phi + X3_phi) / 3;
        
        % Enforce bounds
        pack(i).Position.r = max(pack(i).Position.r, VarMin.r); pack(i).Position.r = min(pack(i).Position.r, VarMax.r);
        pack(i).Position.psi = max(pack(i).Position.psi, VarMin.psi); pack(i).Position.psi = min(pack(i).Position.psi, VarMax.psi);
        pack(i).Position.phi = max(pack(i).Position.phi, VarMin.phi); pack(i).Position.phi = min(pack(i).Position.phi, VarMax.phi);
        
        % Evaluation
        cartPos = SphericalToCart(pack(i).Position, model);
        if any(isnan(cartPos.x))||any(isnan(cartPos.y))||any(isnan(cartPos.z))
            pack(i).Cost = inf;
        else
            try; pack(i).Cost = CostFunction(cartPos); catch; pack(i).Cost = inf; end
        end
        
        % Update alpha, beta, delta
        if pack(i).Cost < Alpha.Cost
            Delta = Beta; Beta = Alpha; Alpha.Position = pack(i).Position; Alpha.Cost = pack(i).Cost;
        elseif pack(i).Cost < Beta.Cost
            Delta = Beta; Beta.Position = pack(i).Position; Beta.Cost = pack(i).Cost;
        elseif pack(i).Cost < Delta.Cost
            Delta.Position = pack(i).Position; Delta.Cost = pack(i).Cost;
        end
    end
    
    if mod(t,10)==0 || t==1
        disp(['Iteration ' num2str(t) ': Best Cost = ' num2str(BestCost(t))]);
    end
end

%% Plot Results
disp('==================== GWO OPTIMIZATION COMPLETE ====================');
disp(['Final Best Cost (GWO) = ' num2str(Alpha.Cost)]);

BestPosition = SphericalToCart(Alpha.Position, model);
smooth = 0.95;
PlotSolution(BestPosition, model, smooth);
title('GWO: Grey Wolf Optimizer - UAV Path Planning');

figure(2);
plot(BestCost, 'LineWidth', 2, 'Color', [0.4 0.4 0.4]);
xlabel('Iteration'); ylabel('Best Cost');
title('GWO Convergence Curve'); grid on;

save('results/GWO_results.mat', 'BestCost', 'Alpha', 'BestPosition');
disp('Results saved to results/GWO_results.mat');

%% Helper function: GWO position update
function [X1, X2, X3] = GWO_update(X, A_Pos, B_Pos, D_Pos, a, VarSize)
    % For Alpha
    A1 = 2*a*rand(VarSize) - a;
    C1 = 2*rand(VarSize);
    D_alpha = abs(C1.*A_Pos - X);
    X1 = A_Pos - A1.*D_alpha;
    
    % For Beta
    A2 = 2*a*rand(VarSize) - a;
    C2 = 2*rand(VarSize);
    D_beta = abs(C2.*B_Pos - X);
    X2 = B_Pos - A2.*D_beta;
    
    % For Delta
    A3 = 2*a*rand(VarSize) - a;
    C3 = 2*rand(VarSize);
    D_delta = abs(C3.*D_Pos - X);
    X3 = D_Pos - A3.*D_delta;
end
