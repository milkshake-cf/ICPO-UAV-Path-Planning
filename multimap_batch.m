%_________________________________________________________________________%
%  Multi-Map Benchmark: SPSO vs GWO vs ICPO across 8 scenarios            %
%  10 runs each, saves all results                                        %
%_________________________________________________________________________%

clc; clear; close all;
N_RUNS = 10; MaxIt = 200; N_MAPS = 8;
algs = {'SPSO', 'GWO', 'ICPO'};
fprintf('Multi-Map: %d algs x %d maps x %d runs = %d total experiments\n', 3, N_MAPS, N_RUNS, 3*N_MAPS*N_RUNS);

allResults = struct();
for a = 1:3
    alg = algs{a};
    for m = 1:N_MAPS
        mk = ['map' num2str(m)];
        allResults.(alg).(mk).finalCosts = zeros(1, N_RUNS);
        allResults.(alg).(mk).times = 0;
    end
end

for m = 1:N_MAPS
    fprintf('\n========== MAP %d/8 ==========\n', m);
    model = CreateModel_map(m);
    CostFunction=@(x) MyCost(x,model);
    nVar=model.n; VarSize=[1 nVar];
    VarMin.x=model.xmin; VarMax.x=model.xmax;
    VarMin.y=model.ymin; VarMax.y=model.ymax;
    VarMin.z=model.zmin; VarMax.z=model.zmax;
    VarMax.r=2*norm(model.start-model.end)/nVar; VarMin.r=0;
    AngleRange=pi/4; VarMin.psi=-AngleRange; VarMax.psi=AngleRange;
    dirVector=model.end-model.start; phi0=atan2(dirVector(2),dirVector(1));
    VarMin.phi=phi0-AngleRange; VarMax.phi=phi0+AngleRange;
    
    for a = 1:3
        alg = algs{a};
        t0 = tic;
        
        for run = 1:N_RUNS
            fc = inf;
            retry = 0;
            while fc == inf && retry < 50
                switch a
                    case 1, fc = run_SPSO_core(model, CostFunction, nVar, VarSize, VarMin, VarMax, MaxIt);
                    case 2, fc = run_GWO_core(model, CostFunction, nVar, VarSize, VarMin, VarMax, MaxIt);
                    case 3, fc = run_ICPO_core(model, CostFunction, nVar, VarSize, VarMin, VarMax, MaxIt);
                end
                retry = retry + 1;
            end
            mk = ['map' num2str(m)];
            allResults.(alg).(mk).finalCosts(run) = fc;
        end
        mk = ['map' num2str(m)];
        allResults.(alg).(mk).times = toc(t0);
        fprintf('  %s: mean=%.1f, min=%.1f, time=%.1fs\n', alg, ...
            mean(allResults.(alg).(mk).finalCosts), ...
            min(allResults.(alg).(mk).finalCosts), ...
            allResults.(alg).(mk).times);
    end
end

%% Summary
fprintf('\n========== MULTI-MAP SUMMARY ==========\n');
fprintf('%-8s', 'Map');
for a = 1:3; fprintf('%14s', algs{a}); end
fprintf('\n');
summary = zeros(N_MAPS, 3);
for m = 1:N_MAPS
    fprintf('Map %-4d', m);
    for a = 1:3
        mk = ['map' num2str(m)];
        summary(m,a) = mean(allResults.(algs{a}).(mk).finalCosts);
        fprintf('%14.1f', summary(m,a));
    end
    fprintf('\n');
end

save('results/multimap_results.mat', 'allResults', 'summary', 'N_RUNS', 'N_MAPS', 'algs');
fprintf('\nSaved to results/multimap_results.mat\n');

%% Helper: SPSO core
function fc = run_SPSO_core(model, CostFunction, nVar, VarSize, VarMin, VarMax, MaxIt)
    nPop=500; w=1; wdamp=0.98; c1=1.5; c2=1.5;
    alpha_v=0.5; VelMax.r=alpha_v*(VarMax.r-VarMin.r); VelMin.r=-VelMax.r;
    VelMax.psi=alpha_v*(VarMax.psi-VarMin.psi); VelMin.psi=-VelMax.psi;
    VelMax.phi=alpha_v*(VarMax.phi-VarMin.phi); VelMin.phi=-VelMax.phi;
    empty_p.Position=[]; empty_p.Velocity=[]; empty_p.Cost=[]; empty_p.Best.Position=[]; empty_p.Best.Cost=[];
    GlobalBest.Cost=inf;
    particle=repmat(empty_p,nPop,1);
    isInit=false;
    while ~isInit
        for i=1:nPop
            particle(i).Position=CreateRandomSolution(VarSize,VarMin,VarMax);
            particle(i).Velocity.r=zeros(VarSize); particle(i).Velocity.psi=zeros(VarSize); particle(i).Velocity.phi=zeros(VarSize);
            particle(i).Cost=CostFunction(SphericalToCart(particle(i).Position,model));
            particle(i).Best.Position=particle(i).Position; particle(i).Best.Cost=particle(i).Cost;
            if particle(i).Best.Cost<GlobalBest.Cost; GlobalBest=particle(i).Best; isInit=true; end
        end
    end
    for it=1:MaxIt
        for i=1:nPop
            particle(i).Velocity.r=w*particle(i).Velocity.r+c1*rand(VarSize).*(particle(i).Best.Position.r-particle(i).Position.r)+c2*rand(VarSize).*(GlobalBest.Position.r-particle(i).Position.r);
            particle(i).Velocity.r=max(particle(i).Velocity.r,VelMin.r); particle(i).Velocity.r=min(particle(i).Velocity.r,VelMax.r);
            particle(i).Position.r=particle(i).Position.r+particle(i).Velocity.r;
            Out=(particle(i).Position.r<VarMin.r|particle(i).Position.r>VarMax.r); particle(i).Velocity.r(Out)=-particle(i).Velocity.r(Out);
            particle(i).Position.r=max(particle(i).Position.r,VarMin.r); particle(i).Position.r=min(particle(i).Position.r,VarMax.r);
            particle(i).Velocity.psi=w*particle(i).Velocity.psi+c1*rand(VarSize).*(particle(i).Best.Position.psi-particle(i).Position.psi)+c2*rand(VarSize).*(GlobalBest.Position.psi-particle(i).Position.psi);
            particle(i).Velocity.psi=max(particle(i).Velocity.psi,VelMin.psi); particle(i).Velocity.psi=min(particle(i).Velocity.psi,VelMax.psi);
            particle(i).Position.psi=particle(i).Position.psi+particle(i).Velocity.psi;
            Out=(particle(i).Position.psi<VarMin.psi|particle(i).Position.psi>VarMax.psi); particle(i).Velocity.psi(Out)=-particle(i).Velocity.psi(Out);
            particle(i).Position.psi=max(particle(i).Position.psi,VarMin.psi); particle(i).Position.psi=min(particle(i).Position.psi,VarMax.psi);
            particle(i).Velocity.phi=w*particle(i).Velocity.phi+c1*rand(VarSize).*(particle(i).Best.Position.phi-particle(i).Position.phi)+c2*rand(VarSize).*(GlobalBest.Position.phi-particle(i).Position.phi);
            particle(i).Velocity.phi=max(particle(i).Velocity.phi,VelMin.phi); particle(i).Velocity.phi=min(particle(i).Velocity.phi,VelMax.phi);
            particle(i).Position.phi=particle(i).Position.phi+particle(i).Velocity.phi;
            Out=(particle(i).Position.phi<VarMin.phi|particle(i).Position.phi>VarMax.phi); particle(i).Velocity.phi(Out)=-particle(i).Velocity.phi(Out);
            particle(i).Position.phi=max(particle(i).Position.phi,VarMin.phi); particle(i).Position.phi=min(particle(i).Position.phi,VarMax.phi);
            particle(i).Cost=CostFunction(SphericalToCart(particle(i).Position,model));
            if particle(i).Cost<particle(i).Best.Cost; particle(i).Best.Position=particle(i).Position; particle(i).Best.Cost=particle(i).Cost;
                if particle(i).Best.Cost<GlobalBest.Cost; GlobalBest=particle(i).Best; end
            end
        end
        w=w*wdamp;
    end
    fc = GlobalBest.Cost;
end

%% Helper: GWO core
function fc = run_GWO_core(model, CostFunction, nVar, VarSize, VarMin, VarMax, MaxIt)
    nPop=150;
    empty_w.Position=[]; empty_w.Cost=[];
    Alpha.Cost=inf; Alpha.Position=[]; Beta.Cost=inf; Beta.Position=[]; Delta.Cost=inf; Delta.Position=[];
    pack=repmat(empty_w,nPop,1);
    isInit=false;
    while ~isInit
        for i=1:nPop
            pack(i).Position=CreateRandomSolution(VarSize,VarMin,VarMax);
            pack(i).Cost=CostFunction(SphericalToCart(pack(i).Position,model));
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
            [X1,X2,X3]=GWOh(Alpha.Position.r,Beta.Position.r,Delta.Position.r,pack(i).Position.r,a,VarSize); pack(i).Position.r=(X1+X2+X3)/3;
            [X1,X2,X3]=GWOh(Alpha.Position.psi,Beta.Position.psi,Delta.Position.psi,pack(i).Position.psi,a,VarSize); pack(i).Position.psi=(X1+X2+X3)/3;
            [X1,X2,X3]=GWOh(Alpha.Position.phi,Beta.Position.phi,Delta.Position.phi,pack(i).Position.phi,a,VarSize); pack(i).Position.phi=(X1+X2+X3)/3;
            pack(i).Position.r=max(pack(i).Position.r,VarMin.r); pack(i).Position.r=min(pack(i).Position.r,VarMax.r);
            pack(i).Position.psi=max(pack(i).Position.psi,VarMin.psi); pack(i).Position.psi=min(pack(i).Position.psi,VarMax.psi);
            pack(i).Position.phi=max(pack(i).Position.phi,VarMin.phi); pack(i).Position.phi=min(pack(i).Position.phi,VarMax.phi);
            pack(i).Cost=CostFunction(SphericalToCart(pack(i).Position,model));
            if pack(i).Cost<Alpha.Cost; Delta=Beta; Beta=Alpha; Alpha.Position=pack(i).Position; Alpha.Cost=pack(i).Cost;
            elseif pack(i).Cost<Beta.Cost; Delta=Beta; Beta.Position=pack(i).Position; Beta.Cost=pack(i).Cost;
            elseif pack(i).Cost<Delta.Cost; Delta.Position=pack(i).Position; Delta.Cost=pack(i).Cost; end
        end
    end
    fc = Alpha.Cost;
end
function [X1,X2,X3]=GWOh(Ap,Bp,Dp,X,a,vs)
    A1=2*a*rand(vs)-a; C1=2*rand(vs); X1=Ap-A1.*abs(C1.*Ap-X);
    A2=2*a*rand(vs)-a; C2=2*rand(vs); X2=Bp-A2.*abs(C2.*Bp-X);
    A3=2*a*rand(vs)-a; C3=2*rand(vs); X3=Dp-A3.*abs(C3.*Dp-X);
end

%% Helper: ICPO core
function fc = run_ICPO_core(model, CostFunction, nVar, VarSize, VarMin, VarMax, MaxIt)
    nPop=150; alpha=0.2; Tf=0.8;
    empty_a.Position=[]; empty_a.Cost=[]; empty_a.pBest.Position=[]; empty_a.pBest.Cost=[];
    GlobalBest.Cost=inf;
    pop=repmat(empty_a,nPop,1); prev_pos=cell(nPop,1);
    isInit=false;
    while ~isInit
        for i=1:nPop
            pop(i).Position=CreateRandomSolution(VarSize,VarMin,VarMax);
            pop(i).Cost=CostFunction(SphericalToCart(pop(i).Position,model));
            pop(i).pBest.Position=pop(i).Position; pop(i).pBest.Cost=pop(i).Cost;
            prev_pos{i}=pop(i).Position;
            if pop(i).Cost<GlobalBest.Cost; GlobalBest.Position=pop(i).Position; GlobalBest.Cost=pop(i).Cost; isInit=true; end
        end
    end
    for t=1:MaxIt
        explRatio=0.7*(1-t/MaxIt)^0.5+0.3;
        for i=1:nPop
            U1=rand(VarSize)>rand();
            if rand()<explRatio
                if rand()<rand()
                    k=randi(nPop); m=randi(nPop);
                    gr=(pop(k).pBest.Position.r+pop(m).pBest.Position.r)/2; pop(i).Position.r=pop(i).Position.r+randn(VarSize).*abs(2*rand()*GlobalBest.Position.r-gr);
                    gpsi=(pop(k).pBest.Position.psi+pop(m).pBest.Position.psi)/2; pop(i).Position.psi=pop(i).Position.psi+randn(VarSize).*abs(2*rand()*GlobalBest.Position.psi-gpsi);
                    gphi=(pop(k).pBest.Position.phi+pop(m).pBest.Position.phi)/2; pop(i).Position.phi=pop(i).Position.phi+randn(VarSize).*abs(2*rand()*GlobalBest.Position.phi-gphi);
                else
                    k=randi(nPop); m=randi(nPop);
                    yr=(pop(i).Position.r+pop(k).pBest.Position.r)/2; dr=pop(m).pBest.Position.r-pop(k).pBest.Position.r; pop(i).Position.r=U1.*pop(i).Position.r+(1-U1).*(yr+rand()*dr);
                    ypsi=(pop(i).Position.psi+pop(k).pBest.Position.psi)/2; dpsi=pop(m).pBest.Position.psi-pop(k).pBest.Position.psi; pop(i).Position.psi=U1.*pop(i).Position.psi+(1-U1).*(ypsi+rand()*dpsi);
                    yphi=(pop(i).Position.phi+pop(k).pBest.Position.phi)/2; dphi=pop(m).pBest.Position.phi-pop(k).pBest.Position.phi; pop(i).Position.phi=U1.*pop(i).Position.phi+(1-U1).*(yphi+rand()*dphi);
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
                    Ft_r=rand(VarSize).*(Mt*(-pop(i).Position.r+pop(k).pBest.Position.r)); S_r=S.*Yt.*Ft_r; pop(i).Position.r=GlobalBest.Position.r+(alpha*(1-r2_p)+r2_p)*(U2.*GlobalBest.Position.r-pop(i).Position.r)-S_r;
                    Ft_psi=rand(VarSize).*(Mt*(-pop(i).Position.psi+pop(k).pBest.Position.psi)); S_psi=S.*Yt.*Ft_psi; pop(i).Position.psi=GlobalBest.Position.psi+(alpha*(1-r2_p)+r2_p)*(U2.*GlobalBest.Position.psi-pop(i).Position.psi)-S_psi;
                    Ft_phi=rand(VarSize).*(Mt*(-pop(i).Position.phi+pop(k).pBest.Position.phi)); S_phi=S.*Yt.*Ft_phi; pop(i).Position.phi=GlobalBest.Position.phi+(alpha*(1-r2_p)+r2_p)*(U2.*GlobalBest.Position.phi-pop(i).Position.phi)-S_phi;
                end
            end
            pop(i).Position.r=max(pop(i).Position.r,VarMin.r); pop(i).Position.r=min(pop(i).Position.r,VarMax.r);
            pop(i).Position.psi=max(pop(i).Position.psi,VarMin.psi); pop(i).Position.psi=min(pop(i).Position.psi,VarMax.psi);
            pop(i).Position.phi=max(pop(i).Position.phi,VarMin.phi); pop(i).Position.phi=min(pop(i).Position.phi,VarMax.phi);
            nc=CostFunction(SphericalToCart(pop(i).Position,model));
            if pop(i).Cost<nc; pop(i).Position=prev_pos{i};
            else; prev_pos{i}=pop(i).Position; pop(i).Cost=nc;
                if nc<pop(i).pBest.Cost; pop(i).pBest.Position=pop(i).Position; pop(i).pBest.Cost=nc; end
                if nc<GlobalBest.Cost; GlobalBest.Position=pop(i).Position; GlobalBest.Cost=nc; end
            end
        end
    end
    fc = GlobalBest.Cost;
end
