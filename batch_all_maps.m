%_________________________________________________________________________%
%  Multi-Map Multi-Algorithm Batch Experiment                               %
%  7 maps × 5 algorithms × 10 runs = 350 total runs                       %
%  Maps: Map2(dense) Map3(high-alt) Map4(corridor) Map5(edge)              %
%        Mountain Urban Coastal                                             %
%  Algs: SPSO GWO AGWO AGWOv5 CPO                                          %
%_________________________________________________________________________%

clc; clear; close all;

N_RUNS = 10;
MAX_IT = 200;

%% Map definitions
% Scenario maps (CreateModel_map)
scenario_maps = {2, 3, 4, 5};
scenario_names = {'Map2_Dense', 'Map3_HighAlt', 'Map4_Corridor', 'Map5_Edge'};

% Terrain maps (pre-built .mat models)
terrain_maps = {'model_mountain.mat', 'model_urban.mat', 'model_coastal.mat'};
terrain_names = {'Mountain', 'Urban', 'Coastal'};

all_map_ids = [scenario_maps, {1, 2, 3}];  % numeric IDs: 2,3,4,5,1,2,3
all_map_names = [scenario_names, terrain_names];
all_map_types = [repmat({'scenario'}, 1, 4), repmat({'terrain'}, 1, 3)];
N_MAPS = 7;

ALG_NAMES = {'SPSO', 'GWO', 'AGWO', 'AGWOv5', 'CPO'};
N_ALGS = 5;

%% Results storage
all_results = struct();

for map_idx = 1:N_MAPS
    map_name = all_map_names{map_idx};
    map_type = all_map_types{map_idx};
    fprintf('\n========== MAP %d/7: %s ==========\n', map_idx, map_name);

    % === Load model ===
    if strcmp(map_type, 'scenario')
        map_id = all_map_ids{map_idx};
        model = CreateModel_map(map_id);
    else
        model_file = terrain_maps{all_map_ids{map_idx}};
        loaded = load(model_file);
        model = loaded.model;
    end

    CostFunction = @(x) MyCost(x, model);
    nVar = model.n; VarSize = [1 nVar];

    VarMin.x=model.xmin; VarMax.x=model.xmax;
    VarMin.y=model.ymin; VarMax.y=model.ymax;
    VarMin.z=model.zmin; VarMax.z=model.zmax;
    VarMax.r=2*norm(model.start-model.end)/nVar; VarMin.r=0;
    AngleRange=pi/4; VarMin.psi=-AngleRange; VarMax.psi=AngleRange;
    dirVector=model.end-model.start; phi0=atan2(dirVector(2),dirVector(1));
    VarMin.phi=phi0-AngleRange; VarMax.phi=phi0+AngleRange;

    all_results(map_idx).map_name = map_name;
    all_results(map_idx).model = model;

    for alg_idx = 1:N_ALGS
        alg_name = ALG_NAMES{alg_idx};
        fprintf('  %-8s: ', alg_name);

        best_costs = zeros(1, N_RUNS);
        times = zeros(1, N_RUNS);

        for run = 1:N_RUNS
            tic;
            switch alg_name
                case 'SPSO'
                    best_costs(run) = RunSPSO_batch(model, VarSize, VarMin, VarMax, CostFunction, MAX_IT);
                case 'GWO'
                    best_costs(run) = RunGWO_batch(model, VarSize, VarMin, VarMax, CostFunction, MAX_IT);
                case 'AGWO'
                    best_costs(run) = RunAGWO_batch(model, VarSize, VarMin, VarMax, CostFunction, MAX_IT);
                case 'AGWOv5'
                    best_costs(run) = RunAGWOv5_batch(model, VarSize, VarMin, VarMax, CostFunction, MAX_IT);
                case 'CPO'
                    best_costs(run) = RunCPO_batch(model, VarSize, VarMin, VarMax, CostFunction, MAX_IT);
            end
            times(run) = toc;
            fprintf('.');
        end

        all_results(map_idx).algs(alg_idx).name = alg_name;
        all_results(map_idx).algs(alg_idx).costs = best_costs;
        all_results(map_idx).algs(alg_idx).times = times;
        all_results(map_idx).algs(alg_idx).mean_cost = mean(best_costs);
        all_results(map_idx).algs(alg_idx).std_cost = std(best_costs);
        all_results(map_idx).algs(alg_idx).min_cost = min(best_costs);

        fprintf(' Mean=%.0f+-%.0f (%.1fs)\n', mean(best_costs), std(best_costs), mean(times));
    end
end

%% Save results
mkdir('results'); mkdir('figures');
save('results/all_maps_batch_results.mat', 'all_results', 'all_map_names', 'ALG_NAMES', 'N_RUNS', 'MAX_IT');

%% Print summary tables
fprintf('\n========== SECTION 2: SCENARIO ROBUSTNESS (Maps 2-5) ==========\n');
fprintf('%-10s', 'Map');
for a = 1:N_ALGS; fprintf('%12s', ALG_NAMES{a}); end
fprintf('\n%s\n', repmat('-', 1, 10+12*N_ALGS));

for map_idx = 1:4  % Maps 2-5
    fprintf('%-10s', all_map_names{map_idx});
    for alg_idx = 1:N_ALGS
        a = all_results(map_idx).algs(alg_idx);
        fprintf('%8.0f+-%3.0f', a.mean_cost, a.std_cost);
    end
    fprintf('\n');
end

fprintf('\n========== SECTION 3: TERRAIN GENERALIZATION  ==========\n');
fprintf('%-10s', 'Terrain');
for a = 1:N_ALGS; fprintf('%12s', ALG_NAMES{a}); end
fprintf('\n%s\n', repmat('-', 1, 10+12*N_ALGS));

for map_idx = 5:7  % Mountain/Urban/Coastal
    fprintf('%-10s', all_map_names{map_idx});
    for alg_idx = 1:N_ALGS
        a = all_results(map_idx).algs(alg_idx);
        fprintf('%8.0f+-%3.0f', a.mean_cost, a.std_cost);
    end
    fprintf('\n');
end

%% Determine winner per map
fprintf('\n========== WINNER PER MAP ==========\n');
for map_idx = 1:N_MAPS
    best_mean = inf; best_alg = '';
    for alg_idx = 1:N_ALGS
        if all_results(map_idx).algs(alg_idx).mean_cost < best_mean
            best_mean = all_results(map_idx).algs(alg_idx).mean_cost;
            best_alg = all_results(map_idx).algs(alg_idx).name;
        end
    end
    fprintf('%-15s: %s (%.0f)\n', all_map_names{map_idx}, best_alg, best_mean);
end

fprintf('\nAll results saved to results/all_maps_batch_results.mat\n');

%% ===================== ALGORITHM RUNNERS =====================

function best_cost = RunSPSO_batch(model, VarSize, VarMin, VarMax, CostFunction, MaxIt)
    nVar = model.n; nPop = 500;
    alpha_v = 0.5;
    VelMax.r=alpha_v*(VarMax.r-VarMin.r); VelMin.r=-VelMax.r;
    VelMax.psi=alpha_v*(VarMax.psi-VarMin.psi); VelMin.psi=-VelMax.psi;
    VelMax.phi=alpha_v*(VarMax.phi-VarMin.phi); VelMin.phi=-VelMax.phi;

    w=1; wdamp=0.98; c1=1.5; c2=1.5;
    empty_p.Position=[]; empty_p.Velocity=[]; empty_p.Cost=[]; empty_p.Best.Position=[]; empty_p.Best.Cost=[];
    GBest.Cost=inf; particle=repmat(empty_p,nPop,1);
    isInit=false;
    while ~isInit
        for i=1:nPop
            particle(i).Position=CreateRandomSolution(VarSize,VarMin,VarMax);
            particle(i).Velocity.r=zeros(VarSize); particle(i).Velocity.psi=zeros(VarSize); particle(i).Velocity.phi=zeros(VarSize);
            particle(i).Cost=CostFunction(SphericalToCart(particle(i).Position,model));
            particle(i).Best.Position=particle(i).Position; particle(i).Best.Cost=particle(i).Cost;
            if particle(i).Best.Cost<GBest.Cost; GBest=particle(i).Best; isInit=true; end
        end
    end
    for it=1:MaxIt
        for i=1:nPop
            particle(i).Velocity.r=w*particle(i).Velocity.r+c1*rand(VarSize).*(particle(i).Best.Position.r-particle(i).Position.r)+c2*rand(VarSize).*(GBest.Position.r-particle(i).Position.r);
            particle(i).Velocity.r=max(particle(i).Velocity.r,VelMin.r); particle(i).Velocity.r=min(particle(i).Velocity.r,VelMax.r);
            particle(i).Position.r=particle(i).Position.r+particle(i).Velocity.r;
            Out=(particle(i).Position.r<VarMin.r|particle(i).Position.r>VarMax.r); particle(i).Velocity.r(Out)=-particle(i).Velocity.r(Out);
            particle(i).Position.r=max(particle(i).Position.r,VarMin.r); particle(i).Position.r=min(particle(i).Position.r,VarMax.r);

            particle(i).Velocity.psi=w*particle(i).Velocity.psi+c1*rand(VarSize).*(particle(i).Best.Position.psi-particle(i).Position.psi)+c2*rand(VarSize).*(GBest.Position.psi-particle(i).Position.psi);
            particle(i).Velocity.psi=max(particle(i).Velocity.psi,VelMin.psi); particle(i).Velocity.psi=min(particle(i).Velocity.psi,VelMax.psi);
            particle(i).Position.psi=particle(i).Position.psi+particle(i).Velocity.psi;
            Out=(particle(i).Position.psi<VarMin.psi|particle(i).Position.psi>VarMax.psi); particle(i).Velocity.psi(Out)=-particle(i).Velocity.psi(Out);
            particle(i).Position.psi=max(particle(i).Position.psi,VarMin.psi); particle(i).Position.psi=min(particle(i).Position.psi,VarMax.psi);

            particle(i).Velocity.phi=w*particle(i).Velocity.phi+c1*rand(VarSize).*(particle(i).Best.Position.phi-particle(i).Position.phi)+c2*rand(VarSize).*(GBest.Position.phi-particle(i).Position.phi);
            particle(i).Velocity.phi=max(particle(i).Velocity.phi,VelMin.phi); particle(i).Velocity.phi=min(particle(i).Velocity.phi,VelMax.phi);
            particle(i).Position.phi=particle(i).Position.phi+particle(i).Velocity.phi;
            Out=(particle(i).Position.phi<VarMin.phi|particle(i).Position.phi>VarMax.phi); particle(i).Velocity.phi(Out)=-particle(i).Velocity.phi(Out);
            particle(i).Position.phi=max(particle(i).Position.phi,VarMin.phi); particle(i).Position.phi=min(particle(i).Position.phi,VarMax.phi);

            particle(i).Cost=CostFunction(SphericalToCart(particle(i).Position,model));
            if particle(i).Cost<particle(i).Best.Cost; particle(i).Best.Position=particle(i).Position; particle(i).Best.Cost=particle(i).Cost;
                if particle(i).Best.Cost<GBest.Cost; GBest=particle(i).Best; end
            end
        end
        w=w*wdamp;
    end
    best_cost = GBest.Cost;
end

function best_cost = RunGWO_batch(model, VarSize, VarMin, VarMax, CostFunction, MaxIt)
    nVar = model.n; nPop = 150;
    empty_w.Position=[]; empty_w.Cost=[];
    Alpha.Cost=inf; Alpha.Position=[]; Beta.Cost=inf; Beta.Position=[]; Delta.Cost=inf; Delta.Position=[];
    pack=repmat(empty_w,nPop,1);
    isInit=false;
    while ~isInit
        for i=1:nPop
            pack(i).Position=CreateRandomSolution(VarSize,VarMin,VarMax);
            cp=SphericalToCart(pack(i).Position,model);
            if any(isnan(cp.x))||any(isnan(cp.y))||any(isnan(cp.z)); pack(i).Cost=inf;
            else; try; pack(i).Cost=CostFunction(cp); catch; pack(i).Cost=inf; end; end
            if pack(i).Cost<Alpha.Cost; Delta=Beta; Beta=Alpha; Alpha.Position=pack(i).Position; Alpha.Cost=pack(i).Cost; isInit=true;
            elseif pack(i).Cost<Beta.Cost; Delta=Beta; Beta.Position=pack(i).Position; Beta.Cost=pack(i).Cost;
            elseif pack(i).Cost<Delta.Cost; Delta.Position=pack(i).Position; Delta.Cost=pack(i).Cost; end
        end
    end
    if isempty(Beta.Position); Beta=Alpha; end
    if isempty(Delta.Position); Delta=Alpha; end
    for t=1:MaxIt
        a=2-t*(2/MaxIt);
        for i=1:nPop
            [X1,X2,X3]=GWO_update(pack(i).Position.r,Alpha.Position.r,Beta.Position.r,Delta.Position.r,a,VarSize);
            pack(i).Position.r=(X1+X2+X3)/3;
            [X1,X2,X3]=GWO_update(pack(i).Position.psi,Alpha.Position.psi,Beta.Position.psi,Delta.Position.psi,a,VarSize);
            pack(i).Position.psi=(X1+X2+X3)/3;
            [X1,X2,X3]=GWO_update(pack(i).Position.phi,Alpha.Position.phi,Beta.Position.phi,Delta.Position.phi,a,VarSize);
            pack(i).Position.phi=(X1+X2+X3)/3;
            pack(i).Position.r=max(pack(i).Position.r,VarMin.r); pack(i).Position.r=min(pack(i).Position.r,VarMax.r);
            pack(i).Position.psi=max(pack(i).Position.psi,VarMin.psi); pack(i).Position.psi=min(pack(i).Position.psi,VarMax.psi);
            pack(i).Position.phi=max(pack(i).Position.phi,VarMin.phi); pack(i).Position.phi=min(pack(i).Position.phi,VarMax.phi);
            cp=SphericalToCart(pack(i).Position,model);
            if any(isnan(cp.x))||any(isnan(cp.y))||any(isnan(cp.z)); pack(i).Cost=inf;
            else; try; pack(i).Cost=CostFunction(cp); catch; pack(i).Cost=inf; end; end
            if pack(i).Cost<Alpha.Cost; Delta=Beta; Beta=Alpha; Alpha.Position=pack(i).Position; Alpha.Cost=pack(i).Cost;
            elseif pack(i).Cost<Beta.Cost; Delta=Beta; Beta.Position=pack(i).Position; Beta.Cost=pack(i).Cost;
            elseif pack(i).Cost<Delta.Cost; Delta.Position=pack(i).Position; Delta.Cost=pack(i).Cost; end
        end
    end
    best_cost = Alpha.Cost;
end

function [X1,X2,X3]=GWO_update(X,A_Pos,B_Pos,D_Pos,a,VarSize)
    A1=2*a*rand(VarSize)-a; C1=2*rand(VarSize); X1=A_Pos-A1.*abs(C1.*A_Pos-X);
    A2=2*a*rand(VarSize)-a; C2=2*rand(VarSize); X2=B_Pos-A2.*abs(C2.*B_Pos-X);
    A3=2*a*rand(VarSize)-a; C3=2*rand(VarSize); X3=D_Pos-A3.*abs(C3.*D_Pos-X);
end

function best_cost = RunAGWO_batch(model, VarSize, VarMin, VarMax, CostFunction, MaxIt)
    nVar = model.n; nPop = 150; alpha=0.2; Tf=0.8;
    empty_a.Position=[]; empty_a.Cost=[]; empty_a.pBest.Position=[]; empty_a.pBest.Cost=[];
    GBest.Cost=inf; pop=repmat(empty_a,nPop,1); prev_pos=cell(nPop,1);
    isInit=false;
    while ~isInit
        for i=1:nPop
            pop(i).Position=CreateRandomSolution(VarSize,VarMin,VarMax);
            cp=SphericalToCart(pop(i).Position,model);
            if any(isnan(cp.x))||any(isnan(cp.y))||any(isnan(cp.z)); pop(i).Cost=inf;
            else; try; pop(i).Cost=CostFunction(cp); catch; pop(i).Cost=inf; end; end
            pop(i).pBest.Position=pop(i).Position; pop(i).pBest.Cost=pop(i).Cost; prev_pos{i}=pop(i).Position;
            if pop(i).Cost<GBest.Cost; GBest.Position=pop(i).Position; GBest.Cost=pop(i).Cost; isInit=true; end
        end
    end
    for t=1:MaxIt
        explRatio=0.7*(1-t/MaxIt)^0.5+0.3;
        for i=1:nPop
            U1=rand(VarSize)>rand();
            if rand()<explRatio
                if rand()<rand()
                    k=randi(nPop); m=randi(nPop);
                    gr=(pop(k).pBest.Position.r+pop(m).pBest.Position.r)/2;
                    pop(i).Position.r=pop(i).Position.r+randn(VarSize).*abs(2*rand()*GBest.Position.r-gr);
                    gpsi=(pop(k).pBest.Position.psi+pop(m).pBest.Position.psi)/2;
                    pop(i).Position.psi=pop(i).Position.psi+randn(VarSize).*abs(2*rand()*GBest.Position.psi-gpsi);
                    gphi=(pop(k).pBest.Position.phi+pop(m).pBest.Position.phi)/2;
                    pop(i).Position.phi=pop(i).Position.phi+randn(VarSize).*abs(2*rand()*GBest.Position.phi-gphi);
                else
                    k=randi(nPop); m=randi(nPop);
                    yr=(pop(i).Position.r+pop(k).pBest.Position.r)/2; dr=pop(m).pBest.Position.r-pop(k).pBest.Position.r;
                    pop(i).Position.r=U1.*pop(i).Position.r+(1-U1).*(yr+rand()*dr);
                    ypsi=(pop(i).Position.psi+pop(k).pBest.Position.psi)/2; dpsi=pop(m).pBest.Position.psi-pop(k).pBest.Position.psi;
                    pop(i).Position.psi=U1.*pop(i).Position.psi+(1-U1).*(ypsi+rand()*dpsi);
                    yphi=(pop(i).Position.phi+pop(k).pBest.Position.phi)/2; dphi=pop(m).pBest.Position.phi-pop(k).pBest.Position.phi;
                    pop(i).Position.phi=U1.*pop(i).Position.phi+(1-U1).*(yphi+rand()*dphi);
                end
            else
                Yt=2*rand()*(1-t/MaxIt)^(t/MaxIt); U2=(rand(VarSize)<0.5)*2-1; S=rand()*U2;
                sc=0; for j=1:nPop; sc=sc+pop(j).pBest.Cost; end; sf=sc+eps;
                if rand()<Tf
                    St=exp(pop(i).pBest.Cost/sf); S=S.*Yt.*St; k=randi(nPop); m=randi(nPop);
                    pop(i).Position.r=(1-U1).*pop(i).Position.r+U1.*(pop(k).pBest.Position.r+St*(pop(m).pBest.Position.r-pop(k).pBest.Position.r)-S);
                    pop(i).Position.psi=(1-U1).*pop(i).Position.psi+U1.*(pop(k).pBest.Position.psi+St*(pop(m).pBest.Position.psi-pop(k).pBest.Position.psi)-S);
                    pop(i).Position.phi=(1-U1).*pop(i).Position.phi+U1.*(pop(k).pBest.Position.phi+St*(pop(m).pBest.Position.phi-pop(k).pBest.Position.phi)-S);
                else
                    Mt=exp(pop(i).pBest.Cost/sf); k=randi(nPop); r2_p=rand();
                    Ft_r=rand(VarSize).*(Mt*(-pop(i).Position.r+pop(k).pBest.Position.r)); S_r=S.*Yt.*Ft_r;
                    pop(i).Position.r=GBest.Position.r+(alpha*(1-r2_p)+r2_p)*(U2.*GBest.Position.r-pop(i).Position.r)-S_r;
                    Ft_psi=rand(VarSize).*(Mt*(-pop(i).Position.psi+pop(k).pBest.Position.psi)); S_psi=S.*Yt.*Ft_psi;
                    pop(i).Position.psi=GBest.Position.psi+(alpha*(1-r2_p)+r2_p)*(U2.*GBest.Position.psi-pop(i).Position.psi)-S_psi;
                    Ft_phi=rand(VarSize).*(Mt*(-pop(i).Position.phi+pop(k).pBest.Position.phi)); S_phi=S.*Yt.*Ft_phi;
                    pop(i).Position.phi=GBest.Position.phi+(alpha*(1-r2_p)+r2_p)*(U2.*GBest.Position.phi-pop(i).Position.phi)-S_phi;
                end
            end
            pop(i).Position.r=max(pop(i).Position.r,VarMin.r); pop(i).Position.r=min(pop(i).Position.r,VarMax.r);
            pop(i).Position.psi=max(pop(i).Position.psi,VarMin.psi); pop(i).Position.psi=min(pop(i).Position.psi,VarMax.psi);
            pop(i).Position.phi=max(pop(i).Position.phi,VarMin.phi); pop(i).Position.phi=min(pop(i).Position.phi,VarMax.phi);
            cp=SphericalToCart(pop(i).Position,model);
            if any(isnan(cp.x))||any(isnan(cp.y))||any(isnan(cp.z)); nc=inf;
            else; try; nc=CostFunction(cp); catch; nc=inf; end; end
            if pop(i).Cost<nc; pop(i).Position=prev_pos{i};
            else
                prev_pos{i}=pop(i).Position; pop(i).Cost=nc;
                if nc<pop(i).pBest.Cost; pop(i).pBest.Position=pop(i).Position; pop(i).pBest.Cost=nc; end
                if nc<GBest.Cost; GBest.Position=pop(i).Position; GBest.Cost=nc; end
            end
        end
    end
    best_cost = GBest.Cost;
end

function best_cost = RunAGWOv5_batch(model, VarSize, VarMin, VarMax, CostFunction, MaxIt)
    nVar = model.n; nPop = 150;
    empty_a.Position=[]; empty_a.Cost=[]; empty_a.pBest.Position=[]; empty_a.pBest.Cost=[];
    Alpha.Cost=inf; Alpha.Position=[]; Beta.Cost=inf; Beta.Position=[]; Delta.Cost=inf; Delta.Position=[];
    pop=repmat(empty_a,nPop,1);
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
    for t=1:MaxIt
        a=2*(0.7*(1-t/MaxIt)^0.5+0.3);
        if t<MaxIt*0.3; lp=0.1; elseif t<MaxIt*0.7; lp=0.25; else; lp=0.45; end
        for i=1:nPop
            % GWO targets
            A1=2*a*rand(VarSize)-a; C1=2*rand(VarSize); X1_r=Alpha.Position.r-A1.*abs(C1.*Alpha.Position.r-pop(i).Position.r);
            A2=2*a*rand(VarSize)-a; C2=2*rand(VarSize); X2_r=Beta.Position.r-A2.*abs(C2.*Beta.Position.r-pop(i).Position.r);
            A3=2*a*rand(VarSize)-a; C3=2*rand(VarSize); X3_r=Delta.Position.r-A3.*abs(C3.*Delta.Position.r-pop(i).Position.r);
            wa=1/(Alpha.Cost+eps); wb=1/(Beta.Cost+eps); wd=1/(Delta.Cost+eps); ws=wa+wb+wd;
            tr=(wa*X1_r+wb*X2_r+wd*X3_r)/ws;
            A1=2*a*rand(VarSize)-a; C1=2*rand(VarSize); X1_psi=Alpha.Position.psi-A1.*abs(C1.*Alpha.Position.psi-pop(i).Position.psi);
            A2=2*a*rand(VarSize)-a; C2=2*rand(VarSize); X2_psi=Beta.Position.psi-A2.*abs(C2.*Beta.Position.psi-pop(i).Position.psi);
            A3=2*a*rand(VarSize)-a; C3=2*rand(VarSize); X3_psi=Delta.Position.psi-A3.*abs(C3.*Delta.Position.psi-pop(i).Position.psi);
            tpsi=(wa*X1_psi+wb*X2_psi+wd*X3_psi)/ws;
            A1=2*a*rand(VarSize)-a; C1=2*rand(VarSize); X1_phi=Alpha.Position.phi-A1.*abs(C1.*Alpha.Position.phi-pop(i).Position.phi);
            A2=2*a*rand(VarSize)-a; C2=2*rand(VarSize); X2_phi=Beta.Position.phi-A2.*abs(C2.*Beta.Position.phi-pop(i).Position.phi);
            A3=2*a*rand(VarSize)-a; C3=2*rand(VarSize); X3_phi=Delta.Position.phi-A3.*abs(C3.*Delta.Position.phi-pop(i).Position.phi);
            tphi=(wa*X1_phi+wb*X2_phi+wd*X3_phi)/ws;

            % Local search
            if a<1.0 && rand()<lp
                lr=(VarMax.r-VarMin.r)*0.03*(a/2); tr=Alpha.Position.r+lr*randn(VarSize);
                lpsi=(VarMax.psi-VarMin.psi)*0.03*(a/2); tpsi=Alpha.Position.psi+lpsi*randn(VarSize);
                lphi=(VarMax.phi-VarMin.phi)*0.03*(a/2); tphi=Alpha.Position.phi+lphi*randn(VarSize);
            end

            pop(i).Position.r=tr; pop(i).Position.psi=tpsi; pop(i).Position.phi=tphi;
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
    best_cost = Alpha.Cost;
end

function best_cost = RunCPO_batch(model, VarSize, VarMin, VarMax, CostFunction, MaxIt)
    nVar = model.n; nPop = 150; alpha=0.2; Tf=0.8;
    empty_a.Position=[]; empty_a.Cost=[];
    GBest.Cost=inf; pop=repmat(empty_a,nPop,1); prev_pos=cell(nPop,1);
    isInit=false;
    while ~isInit
        for i=1:nPop
            pop(i).Position=CreateRandomSolution(VarSize,VarMin,VarMax);
            cp=SphericalToCart(pop(i).Position,model);
            if any(isnan(cp.x))||any(isnan(cp.y))||any(isnan(cp.z)); pop(i).Cost=inf;
            else; try; pop(i).Cost=CostFunction(cp); catch; pop(i).Cost=inf; end; end
            prev_pos{i}=pop(i).Position;
            if pop(i).Cost<GBest.Cost; GBest.Position=pop(i).Position; GBest.Cost=pop(i).Cost; isInit=true; end
        end
    end
    for t=1:MaxIt
        for i=1:nPop
            U1=rand(VarSize)>rand();
            if rand()<rand()
                if rand()<rand()
                    k=randi(nPop); y_r=(pop(i).Position.r+pop(k).Position.r)/2;
                    pop(i).Position.r=pop(i).Position.r+randn(VarSize).*abs(2*rand()*GBest.Position.r-y_r);
                    y_psi=(pop(i).Position.psi+pop(k).Position.psi)/2;
                    pop(i).Position.psi=pop(i).Position.psi+randn(VarSize).*abs(2*rand()*GBest.Position.psi-y_psi);
                    y_phi=(pop(i).Position.phi+pop(k).Position.phi)/2;
                    pop(i).Position.phi=pop(i).Position.phi+randn(VarSize).*abs(2*rand()*GBest.Position.phi-y_phi);
                else
                    k=randi(nPop); m=randi(nPop);
                    y_r=(pop(i).Position.r+pop(k).Position.r)/2;
                    pop(i).Position.r=U1.*pop(i).Position.r+(1-U1).*(y_r+rand()*(pop(m).Position.r-pop(k).Position.r));
                    y_psi=(pop(i).Position.psi+pop(k).Position.psi)/2;
                    pop(i).Position.psi=U1.*pop(i).Position.psi+(1-U1).*(y_psi+rand()*(pop(m).Position.psi-pop(k).Position.psi));
                    y_phi=(pop(i).Position.phi+pop(k).Position.phi)/2;
                    pop(i).Position.phi=U1.*pop(i).Position.phi+(1-U1).*(y_phi+rand()*(pop(m).Position.phi-pop(k).Position.phi));
                end
            else
                Yt=2*rand()*(1-t/MaxIt)^(t/MaxIt); U2=(rand(VarSize)<0.5)*2-1; S=rand()*U2;
                sc=0; for j=1:nPop; sc=sc+pop(j).Cost; end; sf=sc+eps;
                if rand()<Tf
                    St=exp(pop(i).Cost/sf); S=S.*Yt.*St; k=randi(nPop); m=randi(nPop);
                    pop(i).Position.r=(1-U1).*pop(i).Position.r+U1.*(pop(k).Position.r+St*(pop(m).Position.r-pop(k).Position.r)-S);
                    pop(i).Position.psi=(1-U1).*pop(i).Position.psi+U1.*(pop(k).Position.psi+St*(pop(m).Position.psi-pop(k).Position.psi)-S);
                    pop(i).Position.phi=(1-U1).*pop(i).Position.phi+U1.*(pop(k).Position.phi+St*(pop(m).Position.phi-pop(k).Position.phi)-S);
                else
                    Mt=exp(pop(i).Cost/sf); k=randi(nPop); r2_p=rand();
                    Ft_r=rand(VarSize).*(Mt*(-pop(i).Position.r+pop(k).Position.r)); S_r=S.*Yt.*Ft_r;
                    pop(i).Position.r=GBest.Position.r+(alpha*(1-r2_p)+r2_p)*(U2.*GBest.Position.r-pop(i).Position.r)-S_r;
                    Ft_psi=rand(VarSize).*(Mt*(-pop(i).Position.psi+pop(k).Position.psi)); S_psi=S.*Yt.*Ft_psi;
                    pop(i).Position.psi=GBest.Position.psi+(alpha*(1-r2_p)+r2_p)*(U2.*GBest.Position.psi-pop(i).Position.psi)-S_psi;
                    Ft_phi=rand(VarSize).*(Mt*(-pop(i).Position.phi+pop(k).Position.phi)); S_phi=S.*Yt.*Ft_phi;
                    pop(i).Position.phi=GBest.Position.phi+(alpha*(1-r2_p)+r2_p)*(U2.*GBest.Position.phi-pop(i).Position.phi)-S_phi;
                end
            end
            pop(i).Position.r=max(pop(i).Position.r,VarMin.r); pop(i).Position.r=min(pop(i).Position.r,VarMax.r);
            pop(i).Position.psi=max(pop(i).Position.psi,VarMin.psi); pop(i).Position.psi=min(pop(i).Position.psi,VarMax.psi);
            pop(i).Position.phi=max(pop(i).Position.phi,VarMin.phi); pop(i).Position.phi=min(pop(i).Position.phi,VarMax.phi);
            cp=SphericalToCart(pop(i).Position,model);
            if any(isnan(cp.x))||any(isnan(cp.y))||any(isnan(cp.z)); nc=inf;
            else; try; nc=CostFunction(cp); catch; nc=inf; end; end
            if pop(i).Cost<nc; pop(i).Position=prev_pos{i};
            else; prev_pos{i}=pop(i).Position; pop(i).Cost=nc;
                if nc<GBest.Cost; GBest.Position=pop(i).Position; GBest.Cost=nc; end
            end
        end
    end
    best_cost = GBest.Cost;
end
