function cost = runWOA_mm(model, nPop, MaxIt)
    nVar = model.n; VarSize = [1 nVar];
    VarMin.x = model.xmin; VarMax.x = model.xmax;
    VarMin.y = model.ymin; VarMax.y = model.ymax;
    VarMin.z = model.zmin; VarMax.z = model.zmax;
    VarMax.r = 2*norm(model.start-model.end)/nVar; VarMin.r = 0;
    AngleRange = pi/4; VarMin.psi = -AngleRange; VarMax.psi = AngleRange;
    dirVector = model.end - model.start;
    phi0 = atan2(dirVector(2), dirVector(1));
    VarMin.phi = phi0 - AngleRange; VarMax.phi = phi0 + AngleRange;

    empty_whale.Position=[]; empty_whale.Cost=[];
    whales=repmat(empty_whale,nPop,1);
    GlobalBest.Cost=inf; GlobalBest.Position=[];
    isInit=false;
    while ~isInit
        for i=1:nPop
            whales(i).Position=CreateRandomSolution(VarSize,VarMin,VarMax);
            whales(i).Cost=MyCost(SphericalToCart(whales(i).Position,model),model);
            if whales(i).Cost<GlobalBest.Cost
                GlobalBest.Position=whales(i).Position; GlobalBest.Cost=whales(i).Cost; isInit=true;
            end
        end
    end

    conv=zeros(1,MaxIt);
    for iter=1:MaxIt
        a=2*(1-iter/MaxIt);
        for i=1:nPop
            for comp={'r','psi','phi'}
                c=comp{1}; A=2*a*rand(VarSize)-a; C=2*rand(VarSize);
                l=rand(VarSize)*2-1; p=rand();
                if p<0.5
                    if abs(mean(A))<1
                        D=abs(C.*GlobalBest.Position.(c)-whales(i).Position.(c));
                        whales(i).Position.(c)=GlobalBest.Position.(c)-A.*D;
                    else
                        k=randi(nPop); D=abs(C.*whales(k).Position.(c)-whales(i).Position.(c));
                        whales(i).Position.(c)=whales(k).Position.(c)-A.*D;
                    end
                else
                    D=abs(GlobalBest.Position.(c)-whales(i).Position.(c));
                    whales(i).Position.(c)=D.*exp(l).*cos(2*pi*l)+GlobalBest.Position.(c);
                end
            end
            whales(i).Position.r=max(min(whales(i).Position.r,VarMax.r),VarMin.r);
            whales(i).Position.psi=max(min(whales(i).Position.psi,VarMax.psi),VarMin.psi);
            whales(i).Position.phi=max(min(whales(i).Position.phi,VarMax.phi),VarMin.phi);
            whales(i).Cost=MyCost(SphericalToCart(whales(i).Position,model),model);
            if whales(i).Cost<GlobalBest.Cost
                GlobalBest.Position=whales(i).Position; GlobalBest.Cost=whales(i).Cost;
            end
        end
        conv(iter)=GlobalBest.Cost;
    end
    cost = GlobalBest.Cost;
end
