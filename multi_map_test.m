%_________________________________________________________________________%
%  Multi-Map Comparison: 5 Algorithms across 3 Terrain Types              %
%  Maps: 1=Mountain (original), 2=Urban, 3=Coastal                       %
%_________________________________________________________________________%

clc; clear; close all;

%% Configuration
N_RUNS = 10;          % Runs per algorithm per map
MAX_IT = 200;         % Max iterations
N_POP  = 150;         % Population for all algorithms
SPSO_POP = 500;        % SPSO uses more particles

MAP_NAMES = {'Mountain', 'Urban', 'Coastal'};
ALG_NAMES = {'SPSO', 'GWO', 'AGWO', 'WOA', 'CPO'};

% Preallocate results: algs × maps × runs
all_results = struct();

%% Run all combinations
for map_id = 1:3
    fprintf('\n========== MAP %d: %s ==========\n', map_id, MAP_NAMES{map_id});

    % Create the model for this map
    model = CreateMap(map_id);

    % Set up problem
    CostFunction = @(x) MyCost(x, model);
    nVar = model.n;
    VarSize = [1 nVar];

    % Bounds
    VarMin.x = model.xmin;  VarMax.x = model.xmax;
    VarMin.y = model.ymin;  VarMax.y = model.ymax;
    VarMin.z = model.zmin;  VarMax.z = model.zmax;
    VarMax.r = 2*norm(model.start - model.end)/nVar;
    VarMin.r = 0;
    AngleRange = pi/4;
    VarMin.psi = -AngleRange;  VarMax.psi = AngleRange;
    dirVector = model.end - model.start;
    phi0 = atan2(dirVector(2), dirVector(1));
    VarMin.phi = phi0 - AngleRange;  VarMax.phi = phi0 + AngleRange;

    % Save model info
    all_results(map_id).map_name = MAP_NAMES{map_id};
    all_results(map_id).model = model;

    for alg_id = 1:5
        alg_name = ALG_NAMES{alg_id};
        fprintf('  %s: ', alg_name);

        best_costs = zeros(1, N_RUNS);
        times = zeros(1, N_RUNS);
        all_convergence = zeros(N_RUNS, MAX_IT);

        for run = 1:N_RUNS
            tic;
            switch alg_name
                case 'SPSO'
                    [best_cost, ~, conv] = RunSPSO(model, SPSO_POP, MAX_IT, ...
                        VarSize, VarMin, VarMax, CostFunction);
                case 'GWO'
                    [best_cost, ~, conv] = RunGWO(model, N_POP, MAX_IT, ...
                        VarSize, VarMin, VarMax, CostFunction);
                case 'AGWO'
                    [best_cost, ~, conv] = RunAGWO(model, N_POP, MAX_IT, ...
                        VarSize, VarMin, VarMax, CostFunction);
                case 'WOA'
                    [best_cost, ~, conv] = RunWOA(model, N_POP, MAX_IT, ...
                        VarSize, VarMin, VarMax, CostFunction);
                case 'CPO'
                    [best_cost, ~, conv] = RunCPO(model, N_POP, MAX_IT, ...
                        VarSize, VarMin, VarMax, CostFunction);
            end
            elapsed = toc;
            best_costs(run) = best_cost;
            times(run) = elapsed;
            all_convergence(run, :) = conv(1:MAX_IT);
            fprintf('.');
        end

        % Store results
        all_results(map_id).algs(alg_id).name = alg_name;
        all_results(map_id).algs(alg_id).best_costs = best_costs;
        all_results(map_id).algs(alg_id).times = times;
        all_results(map_id).algs(alg_id).convergence = all_convergence;
        all_results(map_id).algs(alg_id).mean_cost = mean(best_costs);
        all_results(map_id).algs(alg_id).std_cost = std(best_costs);
        all_results(map_id).algs(alg_id).min_cost = min(best_costs);

        fprintf(' Mean=%.1f ±%.1f (%.1fs)\n', ...
            mean(best_costs), std(best_costs), mean(times));
    end
end

%% Save results
mkdir('results'); mkdir('figures');
save('results/multi_map_results.mat', 'all_results', 'MAP_NAMES', 'ALG_NAMES', 'N_RUNS', 'MAX_IT');

%% Print summary tables
fprintf('\n========== MULTI-MAP SUMMARY ==========\n');
for map_id = 1:3
    fprintf('\n--- %s ---\n', MAP_NAMES{map_id});
    fprintf('%-6s %10s %10s %10s %10s\n', 'Alg', 'Best', 'Mean', 'Std', 'Time');
    fprintf('%s\n', repmat('-', 1, 50));
    for alg_id = 1:5
        a = all_results(map_id).algs(alg_id);
        fprintf('%-6s %10.1f %10.1f %10.1f %10.1f\n', ...
            a.name, min(a.best_costs), a.mean_cost, a.std_cost, mean(a.times));
    end
end

%% Generate multi-map bar chart
figure('Position', [100 100 1200 450]);
for map_id = 1:3
    subplot(1, 3, map_id);
    means = zeros(1,5); stds = zeros(1,5);
    for alg_id = 1:5
        means(alg_id) = all_results(map_id).algs(alg_id).mean_cost;
        stds(alg_id) = all_results(map_id).algs(alg_id).std_cost;
    end
    b = bar(means);
    hold on;
    [ngroups, nbars] = size(means);
    x = nan(nbars, ngroups);
    for i = 1:nbars
        x(i,:) = b(i).XEndPoints;
    end
    errorbar(x', means, stds, 'k', 'LineStyle', 'none', 'LineWidth', 1.5);
    set(gca, 'XTickLabel', ALG_NAMES);
    title(sprintf('%s Terrain', MAP_NAMES{map_id}));
    ylabel('Mean Best Cost');
    grid on;
end
sgtitle('Algorithm Performance Across Different Terrains (10 runs each)');
saveas(gcf, 'figures/multi_map_comparison.png');
fprintf('\nFigures saved: figures/multi_map_comparison.png\n');
fprintf('Results saved: results/multi_map_results.mat\n');

%% ======================= MAP CREATION =======================
function model = CreateMap(map_id)
    switch map_id
        case 1  % Mountain (original)
            H = imread('ChrismasTerrain.tif');
        case 2  % Urban (flat + buildings)
            H = CreateUrbanTerrain();
        case 3  % Coastal (half flat, half hills)
            H = CreateCoastalTerrain();
        otherwise
            error('Unknown map_id: %d', map_id);
    end

    H(H < 0) = 0;
    MAPSIZE_X = size(H,2);
    MAPSIZE_Y = size(H,1);
    [X,Y] = meshgrid(1:MAPSIZE_X, 1:MAPSIZE_Y);

    % Define threats based on map type
    switch map_id
        case 1  % Mountain: 6 scattered cylinders
            threats = [400 500 100 80;
                       600 200 150 70;
                       500 350 150 80;
                       350 200 150 70;
                       700 550 150 70;
                       650 750 150 80];
            start_loc = [200; 100; 150];
            end_loc   = [800; 800; 150];

        case 2  % Urban: 12 buildings (dense, tall, narrow)
            threats = [250 150 150 40;   350 180 160 35;
                       500 120 180 30;   650 160 170 35;
                       180 300 150 40;   300 350 160 30;
                       550 300 180 35;   700 350 170 30;
                       400 500 150 40;   500 550 160 35;
                       250 650 170 30;   600 700 180 40;
                       750 500 150 35;   150 500 160 30;
                       350 750 170 35;   700 250 150 40];
            start_loc = [100; 100; 150];
            end_loc   = [750; 780; 150];

        case 3  % Coastal: fewer but larger threats + water
            threats = [300 300 100 100;
                       650 450 120 90;
                       500 650 110 85;
                       200 600 100 95];
            start_loc = [150; 120; 150];
            end_loc   = [700; 750; 150];
    end

    model.start = start_loc;
    model.end   = end_loc;
    model.n     = 10;
    model.xmin  = 1;
    model.xmax  = MAPSIZE_X;
    model.ymin  = 1;
    model.ymax  = MAPSIZE_Y;
    model.zmin  = 100;
    model.zmax  = 200;
    model.MAPSIZE_X = MAPSIZE_X;
    model.MAPSIZE_Y = MAPSIZE_Y;
    model.X = X;
    model.Y = Y;
    model.H = H;
    model.threats = threats;

    % Show the map briefly (non-blocking)
    figure('Visible', 'off');
    PlotModel(model);
    mkdir('figures');
    saveas(gcf, sprintf('figures/map_%d_%s.png', map_id, ...
        {'mountain','urban','coastal'}{map_id}));
    close;
end

function H = CreateUrbanTerrain()
    % Flat terrain with slight random variation
    H = 50 * ones(800, 800) + 5 * randn(800, 800);
    H(H < 0) = 0;
end

function H = CreateCoastalTerrain()
    % Left half: sea (flat ~0), Right half: hills (sinusoidal)
    H = zeros(800, 800);
    for i = 1:800
        for j = 1:800
            if j < 300
                H(i,j) = max(0, 5 + 2*randn());
            elseif j < 500
                H(i,j) = max(0, (j-300)/2 + 10*sin(i/50)*cos(j/50) + 5*randn());
            else
                x = j - 500;
                H(i,j) = max(0, 100 + 80*sin(i/30 + x/40)*cos(x/60) + 10*randn());
            end
        end
    end
end

%% ======================= ALGORITHM WRAPPERS =======================

function [best_cost, best_sol, conv] = RunSPSO(model, nPop, MaxIt, ...
        VarSize, VarMin, VarMax, CostFunction)
    nVar = model.n;
    CostFunction_wrap = @(x) safeCost(x, CostFunction);

    % PSO parameters
    w = 1; wdamp = 0.98; c1 = 1.5; c2 = 1.5;
    alpha = 0.5;
    VelMax.r = alpha*(VarMax.r - VarMin.r);  VelMin.r = -VelMax.r;
    VelMax.psi = alpha*(VarMax.psi - VarMin.psi);  VelMin.psi = -VelMax.psi;
    VelMax.phi = alpha*(VarMax.phi - VarMin.phi);  VelMin.phi = -VelMax.phi;

    % Empty particle template
    empty_particle.Position = [];
    empty_particle.Velocity = [];
    empty_particle.Cost     = [];
    empty_particle.Best.Position = [];
    empty_particle.Best.Cost     = [];

    % Initialize
    GlobalBest.Cost = inf;
    particle = repmat(empty_particle, nPop, 1);

    isInit = false;
    while ~isInit
        for i = 1:nPop
            particle(i).Position = CreateRandomSolution(VarSize, VarMin, VarMax);
            particle(i).Velocity.r   = zeros(VarSize);
            particle(i).Velocity.psi = zeros(VarSize);
            particle(i).Velocity.phi = zeros(VarSize);
            sol_cart = SphericalToCart(particle(i).Position, model);
            particle(i).Cost = CostFunction_wrap(sol_cart);
            particle(i).Best.Position = particle(i).Position;
            particle(i).Best.Cost = particle(i).Cost;
            if particle(i).Best.Cost < GlobalBest.Cost
                GlobalBest = particle(i).Best;
                isInit = true;
            end
        end
    end

    conv = zeros(1, MaxIt);
    for iter = 1:MaxIt
        for i = 1:nPop
            % r component
            particle(i).Velocity.r = w*particle(i).Velocity.r ...
                + c1*rand(VarSize).*(particle(i).Best.Position.r - particle(i).Position.r) ...
                + c2*rand(VarSize).*(GlobalBest.Position.r - particle(i).Position.r);
            particle(i).Velocity.r = max(min(particle(i).Velocity.r, VelMax.r), VelMin.r);
            particle(i).Position.r = particle(i).Position.r + particle(i).Velocity.r;
            particle(i).Position.r = max(min(particle(i).Position.r, VarMax.r), VarMin.r);

            % psi component
            particle(i).Velocity.psi = w*particle(i).Velocity.psi ...
                + c1*rand(VarSize).*(particle(i).Best.Position.psi - particle(i).Position.psi) ...
                + c2*rand(VarSize).*(GlobalBest.Position.psi - particle(i).Position.psi);
            particle(i).Velocity.psi = max(min(particle(i).Velocity.psi, VelMax.psi), VelMin.psi);
            particle(i).Position.psi = particle(i).Position.psi + particle(i).Velocity.psi;
            particle(i).Position.psi = max(min(particle(i).Position.psi, VarMax.psi), VarMin.psi);

            % phi component
            particle(i).Velocity.phi = w*particle(i).Velocity.phi ...
                + c1*rand(VarSize).*(particle(i).Best.Position.phi - particle(i).Position.phi) ...
                + c2*rand(VarSize).*(GlobalBest.Position.phi - particle(i).Position.phi);
            particle(i).Velocity.phi = max(min(particle(i).Velocity.phi, VelMax.phi), VelMin.phi);
            particle(i).Position.phi = particle(i).Position.phi + particle(i).Velocity.phi;
            particle(i).Position.phi = max(min(particle(i).Position.phi, VarMax.phi), VarMin.phi);

            sol_cart = SphericalToCart(particle(i).Position, model);
            particle(i).Cost = CostFunction_wrap(sol_cart);

            if particle(i).Cost < particle(i).Best.Cost
                particle(i).Best.Position = particle(i).Position;
                particle(i).Best.Cost = particle(i).Cost;
                if particle(i).Best.Cost < GlobalBest.Cost
                    GlobalBest = particle(i).Best;
                end
            end
        end
        w = w * wdamp;
        conv(iter) = GlobalBest.Cost;
    end
    best_cost = GlobalBest.Cost;
    best_sol  = GlobalBest.Position;
end

function [best_cost, best_sol, conv] = RunGWO(model, nPop, MaxIt, ...
        VarSize, VarMin, VarMax, CostFunction)
    CostFunction_wrap = @(x) safeCost(x, CostFunction);

    % Initialize wolves
    empty_wolf.Position = []; empty_wolf.Cost = [];
    pack = repmat(empty_wolf, nPop, 1);

    Alpha.Cost = inf; Alpha.Position = [];
    Beta.Cost  = inf; Beta.Position  = [];
    Delta.Cost = inf; Delta.Position = [];

    isInit = false;
    while ~isInit
        for i = 1:nPop
            pack(i).Position = CreateRandomSolution(VarSize, VarMin, VarMax);
            sol_cart = SphericalToCart(pack(i).Position, model);
            pack(i).Cost = CostFunction_wrap(sol_cart);
            if pack(i).Cost < Alpha.Cost
                Delta = Beta; Beta = Alpha;
                Alpha.Position = pack(i).Position; Alpha.Cost = pack(i).Cost;
                isInit = true;
            elseif pack(i).Cost < Beta.Cost
                Delta = Beta;
                Beta.Position = pack(i).Position; Beta.Cost = pack(i).Cost;
            elseif pack(i).Cost < Delta.Cost
                Delta.Position = pack(i).Position; Delta.Cost = pack(i).Cost;
            end
        end
    end
    % Fallback: if Beta/Delta not set, copy Alpha
    if isempty(Beta.Position),  Beta = Alpha; end
    if isempty(Delta.Position), Delta = Alpha; end

    conv = zeros(1, MaxIt);
    for iter = 1:MaxIt
        a = 2 * (1 - iter/MaxIt);
        for i = 1:nPop
            for comp = {'r', 'psi', 'phi'}
                c = comp{1};
                r1 = rand(VarSize); r2 = rand(VarSize);
                A1 = 2*a*r1 - a; C1 = 2*r2;
                D_alpha = abs(C1 .* Alpha.Position.(c) - pack(i).Position.(c));
                X1 = Alpha.Position.(c) - A1 .* D_alpha;

                r1 = rand(VarSize); r2 = rand(VarSize);
                A2 = 2*a*r1 - a; C2 = 2*r2;
                D_beta = abs(C2 .* Beta.Position.(c) - pack(i).Position.(c));
                X2 = Beta.Position.(c) - A2 .* D_beta;

                r1 = rand(VarSize); r2 = rand(VarSize);
                A3 = 2*a*r1 - a; C3 = 2*r2;
                D_delta = abs(C3 .* Delta.Position.(c) - pack(i).Position.(c));
                X3 = Delta.Position.(c) - A3 .* D_delta;

                pack(i).Position.(c) = (X1 + X2 + X3) / 3;
            end
            % Clip bounds
            pack(i).Position.r   = max(min(pack(i).Position.r,   VarMax.r),   VarMin.r);
            pack(i).Position.psi = max(min(pack(i).Position.psi, VarMax.psi), VarMin.psi);
            pack(i).Position.phi = max(min(pack(i).Position.phi, VarMax.phi), VarMin.phi);

            sol_cart = SphericalToCart(pack(i).Position, model);
            pack(i).Cost = CostFunction_wrap(sol_cart);

            if pack(i).Cost < Alpha.Cost
                Delta = Beta; Beta = Alpha; Alpha.Position = pack(i).Position; Alpha.Cost = pack(i).Cost;
            elseif pack(i).Cost < Beta.Cost
                Delta = Beta; Beta.Position = pack(i).Position; Beta.Cost = pack(i).Cost;
            elseif pack(i).Cost < Delta.Cost
                Delta.Position = pack(i).Position; Delta.Cost = pack(i).Cost;
            end
        end
        conv(iter) = Alpha.Cost;
    end
    best_cost = Alpha.Cost; best_sol = Alpha.Position;
end

function [best_cost, best_sol, conv] = RunAGWO(model, nPop, MaxIt, ...
        VarSize, VarMin, VarMax, CostFunction)
    CostFunction_wrap = @(x) safeCost(x, CostFunction);

    empty_agent.Position = []; empty_agent.Cost = [];
    empty_agent.pBest.Position = []; empty_agent.pBest.Cost = [];

    GlobalBest.Cost = inf; GlobalBest.Position = [];

    pop = repmat(empty_agent, nPop, 1);
    prev_pos = cell(nPop, 1);

    isInit = false;
    while ~isInit
        for i = 1:nPop
            pop(i).Position = CreateRandomSolution(VarSize, VarMin, VarMax);
            sol_cart = SphericalToCart(pop(i).Position, model);
            pop(i).Cost = CostFunction_wrap(sol_cart);
            pop(i).pBest.Position = pop(i).Position;
            pop(i).pBest.Cost = pop(i).Cost;
            prev_pos{i} = pop(i).Position;
            if pop(i).Cost < GlobalBest.Cost
                GlobalBest.Position = pop(i).Position;
                GlobalBest.Cost = pop(i).Cost;
                isInit = true;
            end
        end
    end

    conv = zeros(1, MaxIt);
    for iter = 1:MaxIt
        exploration_rate = 0.5 * exp(-iter/MaxIt);  % Adaptive: 0.5→0.18
        Tf = 0.8;

        for i = 1:nPop
            U1 = rand(VarSize) > rand();
            if rand() < exploration_rate  % Exploration
                if rand() < 0.5  % Strategy 1: Visual deterrence (use pBest)
                    k = randi(nPop);
                    y_r   = (pop(i).pBest.Position.r   + pop(k).Position.r)   / 2;
                    y_psi = (pop(i).pBest.Position.psi + pop(k).Position.psi) / 2;
                    y_phi = (pop(i).pBest.Position.phi + pop(k).Position.phi) / 2;
                    pop(i).Position.r   = pop(i).Position.r   + randn(VarSize) .* abs(2*rand()*GlobalBest.Position.r   - y_r);
                    pop(i).Position.psi = pop(i).Position.psi + randn(VarSize) .* abs(2*rand()*GlobalBest.Position.psi - y_psi);
                    pop(i).Position.phi = pop(i).Position.phi + randn(VarSize) .* abs(2*rand()*GlobalBest.Position.phi - y_phi);
                else  % Strategy 2: Sound deterrence
                    k = randi(nPop); m = randi(nPop);
                    pop(i).Position.r   = U1.*pop(i).Position.r   + (1-U1).*((pop(i).Position.r+pop(k).Position.r)/2 + rand()*(pop(m).Position.r   - pop(k).Position.r));
                    pop(i).Position.psi = U1.*pop(i).Position.psi + (1-U1).*((pop(i).Position.psi+pop(k).Position.psi)/2 + rand()*(pop(m).Position.psi - pop(k).Position.psi));
                    pop(i).Position.phi = U1.*pop(i).Position.phi + (1-U1).*((pop(i).Position.phi+pop(k).Position.phi)/2 + rand()*(pop(m).Position.phi - pop(k).Position.phi));
                end
            else  % Exploitation
                Yt = 2*rand()*(1 - iter/MaxIt)^(iter/MaxIt);
                U2 = (rand(VarSize) < 0.5)*2 - 1;
                S_r = rand()*U2; S_psi = rand()*U2; S_phi = rand()*U2;

                % Fitness-weighted step
                allCosts = zeros(nPop, 1);
                for j = 1:nPop, allCosts(j) = pop(j).pBest.Cost; end
                St = exp(pop(i).pBest.Cost / (sum(allCosts) + 1e-10));

                if rand() < Tf  % Strategy 3: Physical attack
                    k = randi(nPop); m = randi(nPop);
                    S_r = S_r*Yt*St; S_psi = S_psi*Yt*St; S_phi = S_phi*Yt*St;
                    pop(i).Position.r   = (1-U1).*pop(i).Position.r   + U1.*(pop(k).Position.r   + St*(pop(m).Position.r   - pop(k).Position.r)   - S_r);
                    pop(i).Position.psi = (1-U1).*pop(i).Position.psi + U1.*(pop(k).Position.psi + St*(pop(m).Position.psi - pop(k).Position.psi) - S_psi);
                    pop(i).Position.phi = (1-U1).*pop(i).Position.phi + U1.*(pop(k).Position.phi + St*(pop(m).Position.phi - pop(k).Position.phi) - S_phi);
                else  % Strategy 4: Lethal attack
                    k = randi(nPop);
                    alpha = 0.2; r2 = rand();
                    Ft = rand(VarSize) * (St * (-pop(i).pBest.Position.r + pop(k).Position.r));
                    S_r = S_r*Yt*Ft;
                    pop(i).Position.r = GlobalBest.Position.r + (alpha*(1-r2)+r2)*(U2.*GlobalBest.Position.r - pop(i).Position.r) - S_r;
                    Ft = rand(VarSize) * (St * (-pop(i).pBest.Position.psi + pop(k).Position.psi));
                    S_psi = S_psi*Yt*Ft;
                    pop(i).Position.psi = GlobalBest.Position.psi + (alpha*(1-r2)+r2)*(U2.*GlobalBest.Position.psi - pop(i).Position.psi) - S_psi;
                    Ft = rand(VarSize) * (St * (-pop(i).pBest.Position.phi + pop(k).Position.phi));
                    S_phi = S_phi*Yt*Ft;
                    pop(i).Position.phi = GlobalBest.Position.phi + (alpha*(1-r2)+r2)*(U2.*GlobalBest.Position.phi - pop(i).Position.phi) - S_phi;
                end
            end

            % Clip bounds
            pop(i).Position.r   = max(min(pop(i).Position.r,   VarMax.r),   VarMin.r);
            pop(i).Position.psi = max(min(pop(i).Position.psi, VarMax.psi), VarMin.psi);
            pop(i).Position.phi = max(min(pop(i).Position.phi, VarMax.phi), VarMin.phi);

            % Evaluate
            sol_cart = SphericalToCart(pop(i).Position, model);
            newCost = CostFunction_wrap(sol_cart);

            % Greedy selection with pBest
            if newCost <= pop(i).pBest.Cost
                pop(i).pBest.Position = pop(i).Position;
                pop(i).pBest.Cost = newCost;
                pop(i).Cost = newCost;
                prev_pos{i} = pop(i).Position;
                if newCost <= GlobalBest.Cost
                    GlobalBest.Position = pop(i).Position;
                    GlobalBest.Cost = newCost;
                end
            else
                pop(i).Position = prev_pos{i};
                pop(i).Cost = pop(i).pBest.Cost;
            end
        end
        conv(iter) = GlobalBest.Cost;
    end
    best_cost = GlobalBest.Cost; best_sol = GlobalBest.Position;
end

function [best_cost, best_sol, conv] = RunWOA(model, nPop, MaxIt, ...
        VarSize, VarMin, VarMax, CostFunction)
    CostFunction_wrap = @(x) safeCost(x, CostFunction);

    empty_whale.Position = []; empty_whale.Cost = [];
    whales = repmat(empty_whale, nPop, 1);

    GlobalBest.Cost = inf; GlobalBest.Position = [];

    isInit = false;
    while ~isInit
        for i = 1:nPop
            whales(i).Position = CreateRandomSolution(VarSize, VarMin, VarMax);
            sol_cart = SphericalToCart(whales(i).Position, model);
            whales(i).Cost = CostFunction_wrap(sol_cart);
            if whales(i).Cost < GlobalBest.Cost
                GlobalBest.Position = whales(i).Position;
                GlobalBest.Cost = whales(i).Cost;
                isInit = true;
            end
        end
    end

    conv = zeros(1, MaxIt);
    for iter = 1:MaxIt
        a = 2 * (1 - iter/MaxIt);
        for i = 1:nPop
            for comp = {'r', 'psi', 'phi'}
                c = comp{1};
                A = 2*a*rand(VarSize) - a;
                C = 2*rand(VarSize);
                l = rand(VarSize)*2 - 1;
                p = rand();
                b = 1;

                if p < 0.5
                    if abs(mean(A)) < 1
                        D = abs(C.*GlobalBest.Position.(c) - whales(i).Position.(c));
                        whales(i).Position.(c) = GlobalBest.Position.(c) - A.*D;
                    else
                        k = randi(nPop);
                        D = abs(C.*whales(k).Position.(c) - whales(i).Position.(c));
                        whales(i).Position.(c) = whales(k).Position.(c) - A.*D;
                    end
                else
                    D = abs(GlobalBest.Position.(c) - whales(i).Position.(c));
                    whales(i).Position.(c) = D.*exp(b*l).*cos(2*pi*l) + GlobalBest.Position.(c);
                end
            end
            whales(i).Position.r   = max(min(whales(i).Position.r,   VarMax.r),   VarMin.r);
            whales(i).Position.psi = max(min(whales(i).Position.psi, VarMax.psi), VarMin.psi);
            whales(i).Position.phi = max(min(whales(i).Position.phi, VarMax.phi), VarMin.phi);

            sol_cart = SphericalToCart(whales(i).Position, model);
            whales(i).Cost = CostFunction_wrap(sol_cart);
            if whales(i).Cost < GlobalBest.Cost
                GlobalBest.Position = whales(i).Position;
                GlobalBest.Cost = whales(i).Cost;
            end
        end
        conv(iter) = GlobalBest.Cost;
    end
    best_cost = GlobalBest.Cost; best_sol = GlobalBest.Position;
end

function [best_cost, best_sol, conv] = RunCPO(model, nPop, MaxIt, ...
        VarSize, VarMin, VarMax, CostFunction)
    CostFunction_wrap = @(x) safeCost(x, CostFunction);

    empty_agent.Position = []; empty_agent.Cost = [];
    GlobalBest.Cost = inf; GlobalBest.Position = [];

    pop = repmat(empty_agent, nPop, 1);
    prev_pos = cell(nPop, 1);

    isInit = false;
    while ~isInit
        for i = 1:nPop
            pop(i).Position = CreateRandomSolution(VarSize, VarMin, VarMax);
            sol_cart = SphericalToCart(pop(i).Position, model);
            pop(i).Cost = CostFunction_wrap(sol_cart);
            prev_pos{i} = pop(i).Position;
            if pop(i).Cost < GlobalBest.Cost
                GlobalBest.Position = pop(i).Position;
                GlobalBest.Cost = pop(i).Cost;
                isInit = true;
            end
        end
    end

    conv = zeros(1, MaxIt);
    for iter = 1:MaxIt
        Tf = 0.8;
        for i = 1:nPop
            U1 = rand(VarSize) > rand();
            if rand() < rand()  % 50/50 Exploration/Exploitation
                if rand() < 0.5  % Strategy 1: Visual
                    k = randi(nPop);
                    y_r = (pop(i).Position.r + pop(k).Position.r)/2;
                    y_psi = (pop(i).Position.psi + pop(k).Position.psi)/2;
                    y_phi = (pop(i).Position.phi + pop(k).Position.phi)/2;
                    pop(i).Position.r   = pop(i).Position.r   + randn(VarSize).*abs(2*rand()*GlobalBest.Position.r - y_r);
                    pop(i).Position.psi = pop(i).Position.psi + randn(VarSize).*abs(2*rand()*GlobalBest.Position.psi - y_psi);
                    pop(i).Position.phi = pop(i).Position.phi + randn(VarSize).*abs(2*rand()*GlobalBest.Position.phi - y_phi);
                else  % Strategy 2: Sound
                    k = randi(nPop); m = randi(nPop);
                    pop(i).Position.r   = U1.*pop(i).Position.r   + (1-U1).*((pop(i).Position.r+pop(k).Position.r)/2 + rand()*(pop(m).Position.r - pop(k).Position.r));
                    pop(i).Position.psi = U1.*pop(i).Position.psi + (1-U1).*((pop(i).Position.psi+pop(k).Position.psi)/2 + rand()*(pop(m).Position.psi - pop(k).Position.psi));
                    pop(i).Position.phi = U1.*pop(i).Position.phi + (1-U1).*((pop(i).Position.phi+pop(k).Position.phi)/2 + rand()*(pop(m).Position.phi - pop(k).Position.phi));
                end
            else  % Exploitation
                Yt = 2*rand()*(1 - iter/MaxIt)^(iter/MaxIt);
                U2 = (rand(VarSize) < 0.5)*2 - 1;
                S_r = rand()*U2; S_psi = rand()*U2; S_phi = rand()*U2;
                allCosts = zeros(nPop,1);
                for j = 1:nPop, allCosts(j) = pop(j).Cost; end
                St = exp(pop(i).Cost / (sum(allCosts) + 1e-10));

                if rand() < Tf  % Strategy 3
                    k = randi(nPop); m = randi(nPop);
                    S_r = S_r*Yt*St; S_psi = S_psi*Yt*St; S_phi = S_phi*Yt*St;
                    pop(i).Position.r   = (1-U1).*pop(i).Position.r   + U1.*(pop(k).Position.r   + St*(pop(m).Position.r - pop(k).Position.r) - S_r);
                    pop(i).Position.psi = (1-U1).*pop(i).Position.psi + U1.*(pop(k).Position.psi + St*(pop(m).Position.psi - pop(k).Position.psi) - S_psi);
                    pop(i).Position.phi = (1-U1).*pop(i).Position.phi + U1.*(pop(k).Position.phi + St*(pop(m).Position.phi - pop(k).Position.phi) - S_phi);
                else  % Strategy 4
                    k = randi(nPop); alpha=0.2; r2=rand();
                    Ft = rand(VarSize)*(St*(-pop(i).Position.r + pop(k).Position.r));
                    S_r = S_r*Yt*Ft;
                    pop(i).Position.r = GlobalBest.Position.r + (alpha*(1-r2)+r2)*(U2.*GlobalBest.Position.r - pop(i).Position.r) - S_r;
                    Ft = rand(VarSize)*(St*(-pop(i).Position.psi + pop(k).Position.psi));
                    S_psi = S_psi*Yt*Ft;
                    pop(i).Position.psi = GlobalBest.Position.psi + (alpha*(1-r2)+r2)*(U2.*GlobalBest.Position.psi - pop(i).Position.psi) - S_psi;
                    Ft = rand(VarSize)*(St*(-pop(i).Position.phi + pop(k).Position.phi));
                    S_phi = S_phi*Yt*Ft;
                    pop(i).Position.phi = GlobalBest.Position.phi + (alpha*(1-r2)+r2)*(U2.*GlobalBest.Position.phi - pop(i).Position.phi) - S_phi;
                end
            end

            pop(i).Position.r   = max(min(pop(i).Position.r,   VarMax.r),   VarMin.r);
            pop(i).Position.psi = max(min(pop(i).Position.psi, VarMax.psi), VarMin.psi);
            pop(i).Position.phi = max(min(pop(i).Position.phi, VarMax.phi), VarMin.phi);

            sol_cart = SphericalToCart(pop(i).Position, model);
            newCost = CostFunction_wrap(sol_cart);

            if pop(i).Cost < newCost
                pop(i).Position = prev_pos{i};
            else
                prev_pos{i} = pop(i).Position;
                pop(i).Cost = newCost;
                if pop(i).Cost <= GlobalBest.Cost
                    GlobalBest.Position = pop(i).Position;
                    GlobalBest.Cost = pop(i).Cost;
                end
            end
        end
        conv(iter) = GlobalBest.Cost;
    end
    best_cost = GlobalBest.Cost; best_sol = GlobalBest.Position;
end

function cost = safeCost(sol_cart, CostFunction)
    try
        cost = CostFunction(sol_cart);
        if isnan(cost) || isinf(cost)
            cost = 1e10;
        end
    catch
        cost = 1e10;
    end
end
