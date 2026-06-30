%_________________________________________________________________________________%
%  ABLATION: AGWO without pBest memory                                           %
%  Same as AGWO but uses current position instead of pBest in all strategies     %
%  This isolates the contribution of the pBest mechanism                         %
%_________________________________________________________________________________%

clc; clear; close all;

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

MaxIt=200; nPop=150; alpha=0.2; Tf=0.8;
empty_agent.Position=[]; empty_agent.Cost=[];
GlobalBest.Cost=inf;
pop=repmat(empty_agent,nPop,1);
prev_pos=cell(nPop,1);

isInit=false;
while ~isInit
    for i=1:nPop
        pop(i).Position=CreateRandomSolution(VarSize,VarMin,VarMax);
        cp=SphericalToCart(pop(i).Position,model);
        if any(isnan(cp.x))||any(isnan(cp.y))||any(isnan(cp.z)); pop(i).Cost=inf;
        else; try; pop(i).Cost=CostFunction(cp); catch; pop(i).Cost=inf; end; end
        prev_pos{i}=pop(i).Position;
        if pop(i).Cost<GlobalBest.Cost; GlobalBest.Position=pop(i).Position; GlobalBest.Cost=pop(i).Cost; isInit=true; end
    end
end
BestCost=zeros(MaxIt,1);

for t=1:MaxIt
    BestCost(t)=GlobalBest.Cost;
    explRatio=0.7*(1-t/MaxIt)^0.5+0.3;  % Adaptive (same as AGWO)
    
    for i=1:nPop
        U1=rand(VarSize)>rand();
        
        if rand()<explRatio
            if rand()<rand()  % Strategy 1: Visual (NO pBest - use current positions)
                k=randi(nPop); m=randi(nPop);
                gr=(pop(k).Position.r+pop(m).Position.r)/2;
                pop(i).Position.r=pop(i).Position.r+randn(VarSize).*abs(2*rand()*GlobalBest.Position.r-gr);
                gpsi=(pop(k).Position.psi+pop(m).Position.psi)/2;
                pop(i).Position.psi=pop(i).Position.psi+randn(VarSize).*abs(2*rand()*GlobalBest.Position.psi-gpsi);
                gphi=(pop(k).Position.phi+pop(m).Position.phi)/2;
                pop(i).Position.phi=pop(i).Position.phi+randn(VarSize).*abs(2*rand()*GlobalBest.Position.phi-gphi);
            else  % Strategy 2: Sound (NO pBest)
                k=randi(nPop); m=randi(nPop);
                yr=(pop(i).Position.r+pop(k).Position.r)/2; dr=pop(m).Position.r-pop(k).Position.r;
                pop(i).Position.r=U1.*pop(i).Position.r+(1-U1).*(yr+rand()*dr);
                ypsi=(pop(i).Position.psi+pop(k).Position.psi)/2; dpsi=pop(m).Position.psi-pop(k).Position.psi;
                pop(i).Position.psi=U1.*pop(i).Position.psi+(1-U1).*(ypsi+rand()*dpsi);
                yphi=(pop(i).Position.phi+pop(k).Position.phi)/2; dphi=pop(m).Position.phi-pop(k).Position.phi;
                pop(i).Position.phi=U1.*pop(i).Position.phi+(1-U1).*(yphi+rand()*dphi);
            end
        else
            Yt=2*rand()*(1-t/MaxIt)^(t/MaxIt); U2=(rand(VarSize)<0.5)*2-1; S=rand()*U2;
            sc=0; for j=1:nPop; sc=sc+pop(j).Cost; end; sf=sc+eps;
            if rand()<Tf  % Strategy 3: Physical (NO pBest)
                St=exp(pop(i).Cost/sf); S=S.*Yt.*St; k=randi(nPop); m=randi(nPop);
                pop(i).Position.r=(1-U1).*pop(i).Position.r+U1.*(pop(k).Position.r+St*(pop(m).Position.r-pop(k).Position.r)-S);
                pop(i).Position.psi=(1-U1).*pop(i).Position.psi+U1.*(pop(k).Position.psi+St*(pop(m).Position.psi-pop(k).Position.psi)-S);
                pop(i).Position.phi=(1-U1).*pop(i).Position.phi+U1.*(pop(k).Position.phi+St*(pop(m).Position.phi-pop(k).Position.phi)-S);
            else  % Strategy 4: Lethal (NO pBest)
                Mt=exp(pop(i).Cost/sf); k=randi(nPop); r2_p=rand();
                Ft_r=rand(VarSize).*(Mt*(-pop(i).Position.r+pop(k).Position.r)); S_r=S.*Yt.*Ft_r;
                pop(i).Position.r=GlobalBest.Position.r+(alpha*(1-r2_p)+r2_p)*(U2.*GlobalBest.Position.r-pop(i).Position.r)-S_r;
                Ft_psi=rand(VarSize).*(Mt*(-pop(i).Position.psi+pop(k).Position.psi)); S_psi=S.*Yt.*Ft_psi;
                pop(i).Position.psi=GlobalBest.Position.psi+(alpha*(1-r2_p)+r2_p)*(U2.*GlobalBest.Position.psi-pop(i).Position.psi)-S_psi;
                Ft_phi=rand(VarSize).*(Mt*(-pop(i).Position.phi+pop(k).Position.phi)); S_phi=S.*Yt.*Ft_phi;
                pop(i).Position.phi=GlobalBest.Position.phi+(alpha*(1-r2_p)+r2_p)*(U2.*GlobalBest.Position.phi-pop(i).Position.phi)-S_phi;
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
            if nc<GlobalBest.Cost; GlobalBest.Position=pop(i).Position; GlobalBest.Cost=nc; end
        end
    end
    
    if mod(t,10)==0||t==1
        disp(['Iteration ' num2str(t) ': Best Cost = ' num2str(BestCost(t))]);
    end
end

disp('=== AGWO-no-pBest COMPLETE ===');
disp(['Final Best Cost = ' num2str(GlobalBest.Cost)]);
BestPosition=SphericalToCart(GlobalBest.Position,model);
PlotSolution(BestPosition,model,0.95); title('AGWO without pBest');
figure(2); plot(BestCost,'LineWidth',2); xlabel('Iteration'); ylabel('Best Cost');
title('AGWO without pBest - Convergence'); grid on;
save('results/AGWO_nopBest_results.mat','BestCost','GlobalBest','BestPosition');
