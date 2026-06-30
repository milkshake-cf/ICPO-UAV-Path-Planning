function cost = runAGWO_mm(model, nPop, MaxIt)
    nVar = model.n; VarSize = [1 nVar];
    VarMin.x = model.xmin; VarMax.x = model.xmax;
    VarMin.y = model.ymin; VarMax.y = model.ymax;
    VarMin.z = model.zmin; VarMax.z = model.zmax;
    VarMax.r = 2*norm(model.start-model.end)/nVar; VarMin.r = 0;
    AngleRange = pi/4; VarMin.psi = -AngleRange; VarMax.psi = AngleRange;
    dirVector = model.end - model.start;
    phi0 = atan2(dirVector(2), dirVector(1));
    VarMin.phi = phi0 - AngleRange; VarMax.phi = phi0 + AngleRange;

    empty_agent.Position=[]; empty_agent.Cost=[];
    empty_agent.pBest.Position=[]; empty_agent.pBest.Cost=[];
    GlobalBest.Cost=inf; GlobalBest.Position=[];
    pop=repmat(empty_agent,nPop,1); prev_pos=cell(nPop,1);

    isInit=false;
    while ~isInit
        for i=1:nPop
            pop(i).Position=CreateRandomSolution(VarSize,VarMin,VarMax);
            pop(i).Cost=MyCost(SphericalToCart(pop(i).Position,model),model);
            pop(i).pBest.Position=pop(i).Position; pop(i).pBest.Cost=pop(i).Cost;
            prev_pos{i}=pop(i).Position;
            if pop(i).Cost < GlobalBest.Cost
                GlobalBest.Position=pop(i).Position; GlobalBest.Cost=pop(i).Cost; isInit=true;
            end
        end
    end

    conv=zeros(1,MaxIt);
    for iter=1:MaxIt
        er=0.5*exp(-iter/MaxIt);
        for i=1:nPop
            U1=rand(VarSize)>rand();
            if rand() < er
                if rand()<0.5
                    k=randi(nPop);
                    pop(i).Position.r=pop(i).Position.r+randn(VarSize).*abs(2*rand()*GlobalBest.Position.r-(pop(i).pBest.Position.r+pop(k).Position.r)/2);
                    pop(i).Position.psi=pop(i).Position.psi+randn(VarSize).*abs(2*rand()*GlobalBest.Position.psi-(pop(i).pBest.Position.psi+pop(k).Position.psi)/2);
                    pop(i).Position.phi=pop(i).Position.phi+randn(VarSize).*abs(2*rand()*GlobalBest.Position.phi-(pop(i).pBest.Position.phi+pop(k).Position.phi)/2);
                else
                    k=randi(nPop); m=randi(nPop);
                    pop(i).Position.r=U1.*pop(i).Position.r+(1-U1).*((pop(i).Position.r+pop(k).Position.r)/2+rand()*(pop(m).Position.r-pop(k).Position.r));
                    pop(i).Position.psi=U1.*pop(i).Position.psi+(1-U1).*((pop(i).Position.psi+pop(k).Position.psi)/2+rand()*(pop(m).Position.psi-pop(k).Position.psi));
                    pop(i).Position.phi=U1.*pop(i).Position.phi+(1-U1).*((pop(i).Position.phi+pop(k).Position.phi)/2+rand()*(pop(m).Position.phi-pop(k).Position.phi));
                end
            else
                Yt=2*rand()*(1-iter/MaxIt)^(iter/MaxIt);
                U2=(rand(VarSize)<0.5)*2-1;
                Sr=rand()*U2; Spsi=rand()*U2; Sphi=rand()*U2;
                allCosts=zeros(nPop,1);
                for j=1:nPop, allCosts(j)=pop(j).pBest.Cost; end
                St=exp(pop(i).pBest.Cost/(sum(allCosts)+1e-10));
                if rand()<0.8
                    k=randi(nPop); m=randi(nPop);
                    Sr=Sr*Yt*St; Spsi=Spsi*Yt*St; Sphi=Sphi*Yt*St;
                    pop(i).Position.r=(1-U1).*pop(i).Position.r+U1.*(pop(k).Position.r+St*(pop(m).Position.r-pop(k).Position.r)-Sr);
                    pop(i).Position.psi=(1-U1).*pop(i).Position.psi+U1.*(pop(k).Position.psi+St*(pop(m).Position.psi-pop(k).Position.psi)-Spsi);
                    pop(i).Position.phi=(1-U1).*pop(i).Position.phi+U1.*(pop(k).Position.phi+St*(pop(m).Position.phi-pop(k).Position.phi)-Sphi);
                else
                    k=randi(nPop); alpha=0.2; r2=rand();
                    Ft=rand(VarSize).*(St*(-pop(i).pBest.Position.r+pop(k).Position.r));
                    Sr=Sr.*Yt.*Ft;
                    pop(i).Position.r=GlobalBest.Position.r+(alpha*(1-r2)+r2)*(U2.*GlobalBest.Position.r-pop(i).Position.r)-Sr;
                    Ft=rand(VarSize).*(St*(-pop(i).pBest.Position.psi+pop(k).Position.psi));
                    Spsi=Spsi.*Yt.*Ft;
                    pop(i).Position.psi=GlobalBest.Position.psi+(alpha*(1-r2)+r2)*(U2.*GlobalBest.Position.psi-pop(i).Position.psi)-Spsi;
                    Ft=rand(VarSize).*(St*(-pop(i).pBest.Position.phi+pop(k).Position.phi));
                    Sphi=Sphi.*Yt.*Ft;
                    pop(i).Position.phi=GlobalBest.Position.phi+(alpha*(1-r2)+r2)*(U2.*GlobalBest.Position.phi-pop(i).Position.phi)-Sphi;
                end
            end
            pop(i).Position.r=max(min(pop(i).Position.r,VarMax.r),VarMin.r);
            pop(i).Position.psi=max(min(pop(i).Position.psi,VarMax.psi),VarMin.psi);
            pop(i).Position.phi=max(min(pop(i).Position.phi,VarMax.phi),VarMin.phi);
            newCost=MyCost(SphericalToCart(pop(i).Position,model),model);
            if newCost<=pop(i).pBest.Cost
                pop(i).pBest.Position=pop(i).Position; pop(i).pBest.Cost=newCost; pop(i).Cost=newCost;
                prev_pos{i}=pop(i).Position;
                if newCost<=GlobalBest.Cost
                    GlobalBest.Position=pop(i).Position; GlobalBest.Cost=newCost;
                end
            else
                pop(i).Position=prev_pos{i}; pop(i).Cost=pop(i).pBest.Cost;
            end
        end
        conv(iter)=GlobalBest.Cost;
    end
    cost = GlobalBest.Cost;
end
