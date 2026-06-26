function cost = runGWO_mm(model, nPop, MaxIt)
    nVar = model.n; VarSize = [1 nVar];
    VarMin.x = model.xmin; VarMax.x = model.xmax;
    VarMin.y = model.ymin; VarMax.y = model.ymax;
    VarMin.z = model.zmin; VarMax.z = model.zmax;
    VarMax.r = 2*norm(model.start-model.end)/nVar; VarMin.r = 0;
    AngleRange = pi/4; VarMin.psi = -AngleRange; VarMax.psi = AngleRange;
    dirVector = model.end - model.start;
    phi0 = atan2(dirVector(2), dirVector(1));
    VarMin.phi = phi0 - AngleRange; VarMax.phi = phi0 + AngleRange;

    empty_wolf.Position=[]; empty_wolf.Cost=[];
    pack=repmat(empty_wolf,nPop,1);
    Alpha.Cost=inf; Alpha.Position=[]; Beta.Cost=inf; Beta.Position=[];
    Delta.Cost=inf; Delta.Position=[];
    isInit=false;
    while ~isInit
        for i=1:nPop
            pack(i).Position=CreateRandomSolution(VarSize,VarMin,VarMax);
            pack(i).Cost=MyCost(SphericalToCart(pack(i).Position,model),model);
            if pack(i).Cost < Alpha.Cost
                Delta=Beta; Beta=Alpha; Alpha.Position=pack(i).Position; Alpha.Cost=pack(i).Cost; isInit=true;
            elseif pack(i).Cost < Beta.Cost
                Delta=Beta; Beta.Position=pack(i).Position; Beta.Cost=pack(i).Cost;
            elseif pack(i).Cost < Delta.Cost
                Delta.Position=pack(i).Position; Delta.Cost=pack(i).Cost;
            end
        end
    end
    if isempty(Beta.Position), Beta=Alpha; end
    if isempty(Delta.Position), Delta=Alpha; end

    conv=zeros(1,MaxIt);
    for iter=1:MaxIt
        a=2*(1-iter/MaxIt);
        for i=1:nPop
            for comp={'r','psi','phi'}
                c=comp{1}; r1=rand(VarSize); r2=rand(VarSize);
                A1=2*a*r1-a; C1=2*r2; D_alpha=abs(C1.*Alpha.Position.(c)-pack(i).Position.(c));
                X1=Alpha.Position.(c)-A1.*D_alpha;
                r1=rand(VarSize); r2=rand(VarSize);
                A2=2*a*r1-a; C2=2*r2; D_beta=abs(C2.*Beta.Position.(c)-pack(i).Position.(c));
                X2=Beta.Position.(c)-A2.*D_beta;
                r1=rand(VarSize); r2=rand(VarSize);
                A3=2*a*r1-a; C3=2*r2; D_delta=abs(C3.*Delta.Position.(c)-pack(i).Position.(c));
                X3=Delta.Position.(c)-A3.*D_delta;
                pack(i).Position.(c)=(X1+X2+X3)/3;
            end
            pack(i).Position.r=max(min(pack(i).Position.r,VarMax.r),VarMin.r);
            pack(i).Position.psi=max(min(pack(i).Position.psi,VarMax.psi),VarMin.psi);
            pack(i).Position.phi=max(min(pack(i).Position.phi,VarMax.phi),VarMin.phi);
            pack(i).Cost=MyCost(SphericalToCart(pack(i).Position,model),model);
            if pack(i).Cost<Alpha.Cost
                Delta=Beta; Beta=Alpha; Alpha.Position=pack(i).Position; Alpha.Cost=pack(i).Cost;
            elseif pack(i).Cost<Beta.Cost
                Delta=Beta; Beta.Position=pack(i).Position; Beta.Cost=pack(i).Cost;
            elseif pack(i).Cost<Delta.Cost
                Delta.Position=pack(i).Position; Delta.Cost=pack(i).Cost;
            end
        end
        conv(iter)=Alpha.Cost;
    end
    cost = Alpha.Cost;
end
