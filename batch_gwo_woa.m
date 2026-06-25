%_________________________________________________________________________%
%  Batch: GWO vs WOA for UAV Path Planning                                %
%  Runs each for N_RUNS times, saves results                              %
%_________________________________________________________________________%

clc; clear; close all;
N_RUNS = 20; MaxIt = 200;

results = struct();
results.GWO.bestCosts = zeros(N_RUNS, MaxIt);
results.GWO.finalCosts = zeros(1, N_RUNS);
results.GWO.times = zeros(1, N_RUNS);
results.WOA.bestCosts = zeros(N_RUNS, MaxIt);
results.WOA.finalCosts = zeros(1, N_RUNS);
results.WOA.times = zeros(1, N_RUNS);

%% Run GWO
disp('========== Running GWO ==========');
for run = 1:N_RUNS
    tic;
    model = CreateModel();
    CostFunction=@(x) MyCost(x,model);
    nVar=model.n; VarSize=[1 nVar];
    VarMin.x=model.xmin; VarMax.x=model.xmax;
    VarMin.y=model.ymin; VarMax.y=model.ymax;
    VarMin.z=model.zmin; VarMax.z=model.zmax;
    VarMax.r=2*norm(model.start-model.end)/nVar; VarMin.r=0;
    AngleRange=pi/4; VarMin.psi=-AngleRange; VarMax.psi=AngleRange;
    dirVector=model.end-model.start; phi0=atan2(dirVector(2),dirVector(1));
    VarMin.phi=phi0-AngleRange; VarMax.phi=phi0+AngleRange;
    
    nPop=150;
    empty_wolf.Position=[]; empty_wolf.Cost=[];
    Alpha.Cost=inf; Alpha.Position=[];
    Beta.Cost=inf; Beta.Position=[];
    Delta.Cost=inf; Delta.Position=[];
    pack=repmat(empty_wolf,nPop,1);
    
    isInit=false;
    while ~isInit
        for i=1:nPop
            pack(i).Position=CreateRandomSolution(VarSize,VarMin,VarMax);
            cp=SphericalToCart(pack(i).Position,model);
            if any(isnan(cp.x))||any(isnan(cp.y))||any(isnan(cp.z)); pack(i).Cost=inf;
            else; try; pack(i).Cost=CostFunction(cp); catch; pack(i).Cost=inf; end; end
            if pack(i).Cost<Alpha.Cost
                Delta=Beta; Beta=Alpha; Alpha.Position=pack(i).Position; Alpha.Cost=pack(i).Cost; isInit=true;
            elseif pack(i).Cost<Beta.Cost
                Delta=Beta; Beta.Position=pack(i).Position; Beta.Cost=pack(i).Cost;
            elseif pack(i).Cost<Delta.Cost
                Delta.Position=pack(i).Position; Delta.Cost=pack(i).Cost;
            end
        end
    end
    if isempty(Beta.Position); Beta=Alpha; end
    if isempty(Delta.Position); Delta=Alpha; end
    
    BestCost=zeros(MaxIt,1);
    for t=1:MaxIt
        BestCost(t)=Alpha.Cost;
        a=2-t*(2/MaxIt);
        for i=1:nPop
            [X1_r,X2_r,X3_r]=GWO_update(Alpha.Position.r,Beta.Position.r,Delta.Position.r,pack(i).Position.r,a,VarSize);
            pack(i).Position.r=(X1_r+X2_r+X3_r)/3;
            [X1_p,X2_p,X3_p]=GWO_update(Alpha.Position.psi,Beta.Position.psi,Delta.Position.psi,pack(i).Position.psi,a,VarSize);
            pack(i).Position.psi=(X1_p+X2_p+X3_p)/3;
            [X1_h,X2_h,X3_h]=GWO_update(Alpha.Position.phi,Beta.Position.phi,Delta.Position.phi,pack(i).Position.phi,a,VarSize);
            pack(i).Position.phi=(X1_h+X2_h+X3_h)/3;
            pack(i).Position.r=max(pack(i).Position.r,VarMin.r); pack(i).Position.r=min(pack(i).Position.r,VarMax.r);
            pack(i).Position.psi=max(pack(i).Position.psi,VarMin.psi); pack(i).Position.psi=min(pack(i).Position.psi,VarMax.psi);
            pack(i).Position.phi=max(pack(i).Position.phi,VarMin.phi); pack(i).Position.phi=min(pack(i).Position.phi,VarMax.phi);
            cp=SphericalToCart(pack(i).Position,model);
            if any(isnan(cp.x))||any(isnan(cp.y))||any(isnan(cp.z)); pack(i).Cost=inf;
            else; try; pack(i).Cost=CostFunction(cp); catch; pack(i).Cost=inf; end; end
            if pack(i).Cost<Alpha.Cost
                Delta=Beta; Beta=Alpha; Alpha.Position=pack(i).Position; Alpha.Cost=pack(i).Cost;
            elseif pack(i).Cost<Beta.Cost
                Delta=Beta; Beta.Position=pack(i).Position; Beta.Cost=pack(i).Cost;
            elseif pack(i).Cost<Delta.Cost
                Delta.Position=pack(i).Position; Delta.Cost=pack(i).Cost;
            end
        end
    end
    results.GWO.bestCosts(run,:)=BestCost;
    results.GWO.finalCosts(run)=Alpha.Cost;
    results.GWO.times(run)=toc;
    fprintf('GWO Run %d/%d: Best=%.2f, Time=%.1fs\n',run,N_RUNS,Alpha.Cost,results.GWO.times(run));
end

%% Run WOA
disp('========== Running WOA ==========');
for run = 1:N_RUNS
    tic;
    model = CreateModel();
    CostFunction=@(x) MyCost(x,model);
    nVar=model.n; VarSize=[1 nVar];
    VarMin.x=model.xmin; VarMax.x=model.xmax;
    VarMin.y=model.ymin; VarMax.y=model.ymax;
    VarMin.z=model.zmin; VarMax.z=model.zmax;
    VarMax.r=2*norm(model.start-model.end)/nVar; VarMin.r=0;
    AngleRange=pi/4; VarMin.psi=-AngleRange; VarMax.psi=AngleRange;
    dirVector=model.end-model.start; phi0=atan2(dirVector(2),dirVector(1));
    VarMin.phi=phi0-AngleRange; VarMax.phi=phi0+AngleRange;
    
    nPop=150; b=1;
    empty_whale.Position=[]; empty_whale.Cost=[];
    GlobalBest.Cost=inf;
    pop=repmat(empty_whale,nPop,1);
    isInit=false;
    while ~isInit
        for i=1:nPop
            pop(i).Position=CreateRandomSolution(VarSize,VarMin,VarMax);
            cp=SphericalToCart(pop(i).Position,model);
            if any(isnan(cp.x))||any(isnan(cp.y))||any(isnan(cp.z)); pop(i).Cost=inf;
            else; try; pop(i).Cost=CostFunction(cp); catch; pop(i).Cost=inf; end; end
            if pop(i).Cost<GlobalBest.Cost; GlobalBest.Position=pop(i).Position; GlobalBest.Cost=pop(i).Cost; isInit=true; end
        end
    end
    BestCost=zeros(MaxIt,1);
    for t=1:MaxIt
        BestCost(t)=GlobalBest.Cost;
        a=2-t*(2/MaxIt); a2=-1+t*(-1/MaxIt);
        for i=1:nPop
            r1=rand(); r2=rand(); A=2*a*r1-a; C=2*r2; l=(a2-1)*rand()+1; p=rand();
            if p<0.5
                if abs(A)>=1
                    k=randi(nPop);
                    D_r=abs(C*pop(k).Position.r-pop(i).Position.r); pop(i).Position.r=pop(k).Position.r-A*D_r;
                    D_psi=abs(C*pop(k).Position.psi-pop(i).Position.psi); pop(i).Position.psi=pop(k).Position.psi-A*D_psi;
                    D_phi=abs(C*pop(k).Position.phi-pop(i).Position.phi); pop(i).Position.phi=pop(k).Position.phi-A*D_phi;
                else
                    D_r=abs(C*GlobalBest.Position.r-pop(i).Position.r); pop(i).Position.r=GlobalBest.Position.r-A*D_r;
                    D_psi=abs(C*GlobalBest.Position.psi-pop(i).Position.psi); pop(i).Position.psi=GlobalBest.Position.psi-A*D_psi;
                    D_phi=abs(C*GlobalBest.Position.phi-pop(i).Position.phi); pop(i).Position.phi=GlobalBest.Position.phi-A*D_phi;
                end
            else
                D_r=abs(GlobalBest.Position.r-pop(i).Position.r); pop(i).Position.r=D_r*exp(b*l).*cos(2*pi*l)+GlobalBest.Position.r;
                D_psi=abs(GlobalBest.Position.psi-pop(i).Position.psi); pop(i).Position.psi=D_psi*exp(b*l).*cos(2*pi*l)+GlobalBest.Position.psi;
                D_phi=abs(GlobalBest.Position.phi-pop(i).Position.phi); pop(i).Position.phi=D_phi*exp(b*l).*cos(2*pi*l)+GlobalBest.Position.phi;
            end
            pop(i).Position.r=max(pop(i).Position.r,VarMin.r); pop(i).Position.r=min(pop(i).Position.r,VarMax.r);
            pop(i).Position.psi=max(pop(i).Position.psi,VarMin.psi); pop(i).Position.psi=min(pop(i).Position.psi,VarMax.psi);
            pop(i).Position.phi=max(pop(i).Position.phi,VarMin.phi); pop(i).Position.phi=min(pop(i).Position.phi,VarMax.phi);
            cp=SphericalToCart(pop(i).Position,model);
            if any(isnan(cp.x))||any(isnan(cp.y))||any(isnan(cp.z)); pop(i).Cost=inf;
            else; try; pop(i).Cost=CostFunction(cp); catch; pop(i).Cost=inf; end; end
            if pop(i).Cost<GlobalBest.Cost; GlobalBest.Position=pop(i).Position; GlobalBest.Cost=pop(i).Cost; end
        end
    end
    results.WOA.bestCosts(run,:)=BestCost;
    results.WOA.finalCosts(run)=GlobalBest.Cost;
    results.WOA.times(run)=toc;
    fprintf('WOA Run %d/%d: Best=%.2f, Time=%.1fs\n',run,N_RUNS,GlobalBest.Cost,results.WOA.times(run));
end

%% Stats
fprintf('\n========== FINAL STATISTICS (N_RUNS=%d, MaxIt=%d) ==========\n', N_RUNS, MaxIt);
fprintf('%-10s %12s %12s %12s %12s %12s\n', 'Algorithm', 'Best', 'Worst', 'Mean', 'Std', 'Avg Time(s)');
for a = 1:2
    alg = {'GWO','WOA'}; alg=alg{a};
    fc = results.(alg).finalCosts;
    fprintf('%-10s %12.2f %12.2f %12.2f %12.2f %12.1f\n', alg, min(fc), max(fc), mean(fc), std(fc), mean(results.(alg).times));
end
save('results/batch_gwo_woa_results.mat', 'results', 'N_RUNS', 'MaxIt');
disp('Results saved to results/batch_gwo_woa_results.mat');

%% GWO helper
function [X1,X2,X3]=GWO_update(A_Pos,B_Pos,D_Pos,X,a,VarSize)
    A1=2*a*rand(VarSize)-a; C1=2*rand(VarSize); X1=A_Pos-A1.*abs(C1.*A_Pos-X);
    A2=2*a*rand(VarSize)-a; C2=2*rand(VarSize); X2=B_Pos-A2.*abs(C2.*B_Pos-X);
    A3=2*a*rand(VarSize)-a; C3=2*rand(VarSize); X3=D_Pos-A3.*abs(C3.*D_Pos-X);
end
