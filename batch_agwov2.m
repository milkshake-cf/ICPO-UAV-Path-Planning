%_________________________________________________________________________%
%  AGWO v2 Batch: 20-run comparison vs AGWO v1                             %
%_________________________________________________________________________%

clc; clear; close all;
N_RUNS = 20; MaxIt = 200;
algs = {'AGWOv1', 'AGWOv2'};
results = struct();

for v = 1:2
    alg = algs{v};
    results.(alg).bestCosts = zeros(N_RUNS, MaxIt);
    results.(alg).finalCosts = zeros(1, N_RUNS);
    results.(alg).times = zeros(1, N_RUNS);
    
    disp(['========== Running ' alg ' ==========']);
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
        nPop=150; alpha=0.2; Tf=0.8;
        use_chaos = (v==2);
        use_levy  = (v==2);
        
        empty_a.Position=[]; empty_a.Cost=[]; empty_a.pBest.Position=[]; empty_a.pBest.Cost=[];
        GlobalBest.Cost=inf;
        pop=repmat(empty_a,nPop,1); prev_pos=cell(nPop,1);
        
        ch=rand();
        isInit=false;
        while ~isInit
            for i=1:nPop
                if use_chaos
                    for d=1:nVar; ch=4*ch*(1-ch); pop(i).Position.r(d)=VarMin.r+ch*(VarMax.r-VarMin.r); end
                    for d=1:nVar; ch=4*ch*(1-ch); pop(i).Position.psi(d)=VarMin.psi+ch*(VarMax.psi-VarMin.psi); end
                    for d=1:nVar; ch=4*ch*(1-ch); pop(i).Position.phi(d)=VarMin.phi+ch*(VarMax.phi-VarMin.phi); end
                else
                    pop(i).Position=CreateRandomSolution(VarSize,VarMin,VarMax);
                end
                cp=SphericalToCart(pop(i).Position,model);
                if any(isnan(cp.x))||any(isnan(cp.y))||any(isnan(cp.z)); pop(i).Cost=inf;
                else; try; pop(i).Cost=CostFunction(cp); catch; pop(i).Cost=inf; end; end
                pop(i).pBest.Position=pop(i).Position; pop(i).pBest.Cost=pop(i).Cost;
                prev_pos{i}=pop(i).Position;
                if pop(i).Cost<GlobalBest.Cost; GlobalBest.Position=pop(i).Position; GlobalBest.Cost=pop(i).Cost; isInit=true; end
            end
        end
        BestCost=zeros(MaxIt,1);
        
        for t=1:MaxIt
            BestCost(t)=GlobalBest.Cost;
            explRatio=0.7*(1-t/MaxIt)^0.5+0.3;
            levy_scale=0.01*(1-t/MaxIt)^0.5;
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
                if use_levy && rand()<0.2
                    L=levy(nVar); pop(i).Position.r=pop(i).Position.r+L*levy_scale;
                    pop(i).Position.psi=pop(i).Position.psi+L*0.5*levy_scale;
                    pop(i).Position.phi=pop(i).Position.phi+L*0.5*levy_scale;
                end
                pop(i).Position.r=max(pop(i).Position.r,VarMin.r); pop(i).Position.r=min(pop(i).Position.r,VarMax.r);
                pop(i).Position.psi=max(pop(i).Position.psi,VarMin.psi); pop(i).Position.psi=min(pop(i).Position.psi,VarMax.psi);
                pop(i).Position.phi=max(pop(i).Position.phi,VarMin.phi); pop(i).Position.phi=min(pop(i).Position.phi,VarMax.phi);
                cp=SphericalToCart(pop(i).Position,model);
                if any(isnan(cp.x))||any(isnan(cp.y))||any(isnan(cp.z)); nc=inf;
                else; try; nc=CostFunction(cp); catch; nc=inf; end; end
                if pop(i).Cost<nc; pop(i).Position=prev_pos{i};
                else; prev_pos{i}=pop(i).Position; pop(i).Cost=nc;
                    if nc<pop(i).pBest.Cost; pop(i).pBest.Position=pop(i).Position; pop(i).pBest.Cost=nc; end
                    if nc<GlobalBest.Cost; GlobalBest.Position=pop(i).Position; GlobalBest.Cost=nc; end
                end
            end
        end
        results.(alg).bestCosts(run,:)=BestCost;
        results.(alg).finalCosts(run)=GlobalBest.Cost;
        results.(alg).times(run)=toc;
        fprintf('  %s Run %d/%d: Best=%.1f\n', alg, run, N_RUNS, GlobalBest.Cost);
    end
end

%% Stats
fprintf('\n========== AGWO v1 vs v2 ==========\n');
for v=1:2
    alg=algs{v}; fc=results.(alg).finalCosts;
    fprintf('%-10s Best=%.1f Worst=%.1f Mean=%.1f Std=%.1f\n', alg, min(fc), max(fc), mean(fc), std(fc));
end
[p,h]=ranksum(results.AGWOv1.finalCosts, results.AGWOv2.finalCosts);
fprintf('Wilcoxon p=%.6f (p<0.05 = significant)\n', p);

save('results/AGWOv2_comparison.mat','results','algs','N_RUNS');
disp('Saved to results/AGWOv2_comparison.mat');

function L=levy(d)
    beta=1.5; sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
    L=randn(1,d)*sigma./(abs(randn(1,d)).^(1/beta));
end
