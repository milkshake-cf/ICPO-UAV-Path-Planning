function cost = runSPSO_mm(model, nPop, MaxIt)
    nVar = model.n; VarSize = [1 nVar];
    VarMin.x=model.xmin; VarMax.x=model.xmax;
    VarMin.y=model.ymin; VarMax.y=model.ymax;
    VarMin.z=model.zmin; VarMax.z=model.zmax;
    VarMax.r=2*norm(model.start-model.end)/nVar; VarMin.r=0;
    AngleRange=pi/4; VarMin.psi=-AngleRange; VarMax.psi=AngleRange;
    dirVector=model.end-model.start;
    phi0=atan2(dirVector(2),dirVector(1));
    VarMin.phi=phi0-AngleRange; VarMax.phi=phi0+AngleRange;

    w=1; wdamp=0.98; c1=1.5; c2=1.5; alpha=0.5;
    VelMax.r=alpha*(VarMax.r-VarMin.r); VelMin.r=-VelMax.r;
    VelMax.psi=alpha*(VarMax.psi-VarMin.psi); VelMin.psi=-VelMax.psi;
    VelMax.phi=alpha*(VarMax.phi-VarMin.phi); VelMin.phi=-VelMax.phi;

    empty_particle.Position=[]; empty_particle.Velocity=[];
    empty_particle.Cost=[]; empty_particle.Best.Position=[]; empty_particle.Best.Cost=[];
    GlobalBest.Cost=inf; particle=repmat(empty_particle,nPop,1);

    isInit=false;
    while ~isInit
        for i=1:nPop
            particle(i).Position=CreateRandomSolution(VarSize,VarMin,VarMax);
            particle(i).Velocity.r=zeros(VarSize);
            particle(i).Velocity.psi=zeros(VarSize);
            particle(i).Velocity.phi=zeros(VarSize);
            particle(i).Cost = MyCost(SphericalToCart(particle(i).Position,model),model);
            particle(i).Best.Position=particle(i).Position;
            particle(i).Best.Cost=particle(i).Cost;
            if particle(i).Best.Cost < GlobalBest.Cost, GlobalBest=particle(i).Best; isInit=true; end
        end
    end

    for iter=1:MaxIt
        for i=1:nPop
            particle(i).Velocity.r = w*particle(i).Velocity.r + c1*rand(VarSize).*(particle(i).Best.Position.r-particle(i).Position.r) + c2*rand(VarSize).*(GlobalBest.Position.r-particle(i).Position.r);
            particle(i).Velocity.r = max(min(particle(i).Velocity.r,VelMax.r),VelMin.r);
            particle(i).Position.r = max(min(particle(i).Position.r+particle(i).Velocity.r,VarMax.r),VarMin.r);
            particle(i).Velocity.psi = w*particle(i).Velocity.psi + c1*rand(VarSize).*(particle(i).Best.Position.psi-particle(i).Position.psi) + c2*rand(VarSize).*(GlobalBest.Position.psi-particle(i).Position.psi);
            particle(i).Velocity.psi = max(min(particle(i).Velocity.psi,VelMax.psi),VelMin.psi);
            particle(i).Position.psi = max(min(particle(i).Position.psi+particle(i).Velocity.psi,VarMax.psi),VarMin.psi);
            particle(i).Velocity.phi = w*particle(i).Velocity.phi + c1*rand(VarSize).*(particle(i).Best.Position.phi-particle(i).Position.phi) + c2*rand(VarSize).*(GlobalBest.Position.phi-particle(i).Position.phi);
            particle(i).Velocity.phi = max(min(particle(i).Velocity.phi,VelMax.phi),VelMin.phi);
            particle(i).Position.phi = max(min(particle(i).Position.phi+particle(i).Velocity.phi,VarMax.phi),VarMin.phi);
            particle(i).Cost = MyCost(SphericalToCart(particle(i).Position,model),model);
            if particle(i).Cost < particle(i).Best.Cost
                particle(i).Best.Position=particle(i).Position; particle(i).Best.Cost=particle(i).Cost;
                if particle(i).Best.Cost < GlobalBest.Cost, GlobalBest=particle(i).Best; end
            end
        end
        w=w*wdamp;
    end
    cost = GlobalBest.Cost;
end
