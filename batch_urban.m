%_________________________________________________________________________%
%  Cross-Map Validation: Urban Canyon terrain                             %
%  5 algorithms, N runs each                                              %
%_________________________________________________________________________%
clc; clear; close all;
N_RUNS = 5; MaxIt = 200;
algs = {'SPSO','GWO','AGWO','WOA','CPO'};
results = struct();
for a = 1:5; alg=algs{a};
    results.(alg).bestCosts=zeros(N_RUNS,MaxIt);
    results.(alg).finalCosts=zeros(1,N_RUNS);
    results.(alg).times=zeros(1,N_RUNS);
end

for a = 1:5
    alg = algs{a};
    fprintf('=== %s on Urban Map ===\n', alg);
    for run = 1:N_RUNS
        tic;
        model = CreateModel_urban();
        CostFn=@(x) MyCost(x,model);
        nVar=model.n; VarSize=[1 nVar];
        Vmin.x=model.xmin; Vmax.x=model.xmax;
        Vmin.y=model.ymin; Vmax.y=model.ymax;
        Vmin.z=model.zmin; Vmax.z=model.zmax;
        Vmax.r=2*norm(model.start-model.end)/nVar; Vmin.r=0;
        AR=pi/4; Vmin.psi=-AR; Vmax.psi=AR;
        dv=model.end-model.start; p0=atan2(dv(2),dv(1));
        Vmin.phi=p0-AR; Vmax.phi=p0+AR;
        
        switch alg
            case 'SPSO'
                nP=500; w=1; wd=0.98; c1=1.5; c2=1.5;
                av=0.5; VmaxV.r=av*(Vmax.r-Vmin.r); VminV.r=-VmaxV.r;
                VmaxV.psi=av*(Vmax.psi-Vmin.psi); VminV.psi=-VmaxV.psi;
                VmaxV.phi=av*(Vmax.phi-Vmin.phi); VminV.phi=-VmaxV.phi;
                ep.Position=[]; ep.Velocity=[]; ep.Cost=[]; ep.Best.Position=[]; ep.Best.Cost=[];
                GB.Cost=inf; p=repmat(ep,nP,1); ok=false;
                while ~ok; for i=1:nP
                    p(i).Position=CreateRandomSolution(VarSize,Vmin,Vmax);
                    p(i).Velocity.r=zeros(VarSize); p(i).Velocity.psi=zeros(VarSize); p(i).Velocity.phi=zeros(VarSize);
                    cp=SphericalToCart(p(i).Position,model);
                    if any(isnan(cp.x))||any(isnan(cp.y))||any(isnan(cp.z)); p(i).Cost=inf;
                    else; try; p(i).Cost=CostFn(cp); catch; p(i).Cost=inf; end; end
                    p(i).Best.Position=p(i).Position; p(i).Best.Cost=p(i).Cost;
                    if p(i).Best.Cost<GB.Cost; GB=p(i).Best; ok=true; end
                end; end
                BC=zeros(MaxIt,1);
                for t=1:MaxIt; BC(t)=GB.Cost;
                    for i=1:nP
                        p(i).Velocity.r=w*p(i).Velocity.r+c1*rand(VarSize).*(p(i).Best.Position.r-p(i).Position.r)+c2*rand(VarSize).*(GB.Position.r-p(i).Position.r);
                        p(i).Velocity.r=max(p(i).Velocity.r,VminV.r); p(i).Velocity.r=min(p(i).Velocity.r,VmaxV.r);
                        p(i).Position.r=p(i).Position.r+p(i).Velocity.r;
                        OR=(p(i).Position.r<Vmin.r|p(i).Position.r>Vmax.r); p(i).Velocity.r(OR)=-p(i).Velocity.r(OR);
                        p(i).Position.r=max(p(i).Position.r,Vmin.r); p(i).Position.r=min(p(i).Position.r,Vmax.r);
                        p(i).Velocity.psi=w*p(i).Velocity.psi+c1*rand(VarSize).*(p(i).Best.Position.psi-p(i).Position.psi)+c2*rand(VarSize).*(GB.Position.psi-p(i).Position.psi);
                        p(i).Velocity.psi=max(p(i).Velocity.psi,VminV.psi); p(i).Velocity.psi=min(p(i).Velocity.psi,VmaxV.psi);
                        p(i).Position.psi=p(i).Position.psi+p(i).Velocity.psi;
                        OR=(p(i).Position.psi<Vmin.psi|p(i).Position.psi>Vmax.psi); p(i).Velocity.psi(OR)=-p(i).Velocity.psi(OR);
                        p(i).Position.psi=max(p(i).Position.psi,Vmin.psi); p(i).Position.psi=min(p(i).Position.psi,Vmax.psi);
                        p(i).Velocity.phi=w*p(i).Velocity.phi+c1*rand(VarSize).*(p(i).Best.Position.phi-p(i).Position.phi)+c2*rand(VarSize).*(GB.Position.phi-p(i).Position.phi);
                        p(i).Velocity.phi=max(p(i).Velocity.phi,VminV.phi); p(i).Velocity.phi=min(p(i).Velocity.phi,VmaxV.phi);
                        p(i).Position.phi=p(i).Position.phi+p(i).Velocity.phi;
                        OR=(p(i).Position.phi<Vmin.phi|p(i).Position.phi>Vmax.phi); p(i).Velocity.phi(OR)=-p(i).Velocity.phi(OR);
                        p(i).Position.phi=max(p(i).Position.phi,Vmin.phi); p(i).Position.phi=min(p(i).Position.phi,Vmax.phi);
                        cp=SphericalToCart(p(i).Position,model);
                        if any(isnan(cp.x))||any(isnan(cp.y))||any(isnan(cp.z)); p(i).Cost=inf;
                        else; try; p(i).Cost=CostFn(cp); catch; p(i).Cost=inf; end; end
                        if p(i).Cost<p(i).Best.Cost; p(i).Best.Position=p(i).Position; p(i).Best.Cost=p(i).Cost;
                            if p(i).Best.Cost<GB.Cost; GB=p(i).Best; end
                        end
                    end
                    w=w*wd;
                end
            case 'GWO'
                nP=150; ew.Position=[]; ew.Cost=[];
                A.Cost=inf; A.Position=[]; B.Cost=inf; B.Position=[]; D.Cost=inf; D.Position=[];
                pack=repmat(ew,nP,1); ok=false;
                while ~ok; for i=1:nP
                    pack(i).Position=CreateRandomSolution(VarSize,Vmin,Vmax);
                    cp=SphericalToCart(pack(i).Position,model);
                    if any(isnan(cp.x))||any(isnan(cp.y))||any(isnan(cp.z)); pack(i).Cost=inf;
                    else; try; pack(i).Cost=CostFn(cp); catch; pack(i).Cost=inf; end; end
                    if pack(i).Cost<A.Cost; D=B; B=A; A.Position=pack(i).Position; A.Cost=pack(i).Cost; ok=true;
                    elseif pack(i).Cost<B.Cost; D=B; B.Position=pack(i).Position; B.Cost=pack(i).Cost;
                    elseif pack(i).Cost<D.Cost; D.Position=pack(i).Position; D.Cost=pack(i).Cost; end
                end; end
                if isempty(B.Position); B=A; end; if isempty(D.Position); D=A; end
                BC=zeros(MaxIt,1);
                for t=1:MaxIt; BC(t)=A.Cost; a_g=2-t*(2/MaxIt);
                    for i=1:nP
                        [X1,X2,X3]=GWOu(A.Position.r,B.Position.r,D.Position.r,pack(i).Position.r,a_g,VarSize); pack(i).Position.r=(X1+X2+X3)/3;
                        [X1,X2,X3]=GWOu(A.Position.psi,B.Position.psi,D.Position.psi,pack(i).Position.psi,a_g,VarSize); pack(i).Position.psi=(X1+X2+X3)/3;
                        [X1,X2,X3]=GWOu(A.Position.phi,B.Position.phi,D.Position.phi,pack(i).Position.phi,a_g,VarSize); pack(i).Position.phi=(X1+X2+X3)/3;
                        pack(i).Position.r=max(pack(i).Position.r,Vmin.r); pack(i).Position.r=min(pack(i).Position.r,Vmax.r);
                        pack(i).Position.psi=max(pack(i).Position.psi,Vmin.psi); pack(i).Position.psi=min(pack(i).Position.psi,Vmax.psi);
                        pack(i).Position.phi=max(pack(i).Position.phi,Vmin.phi); pack(i).Position.phi=min(pack(i).Position.phi,Vmax.phi);
                        cp=SphericalToCart(pack(i).Position,model);
                        if any(isnan(cp.x))||any(isnan(cp.y))||any(isnan(cp.z)); pack(i).Cost=inf;
                        else; try; pack(i).Cost=CostFn(cp); catch; pack(i).Cost=inf; end; end
                        if pack(i).Cost<A.Cost; D=B; B=A; A.Position=pack(i).Position; A.Cost=pack(i).Cost;
                        elseif pack(i).Cost<B.Cost; D=B; B.Position=pack(i).Position; B.Cost=pack(i).Cost;
                        elseif pack(i).Cost<D.Cost; D.Position=pack(i).Position; D.Cost=pack(i).Cost; end
                    end
                end; GB.Cost=A.Cost;
            case 'AGWO'
                nP=150; al=0.2; Tf=0.8;
                ea.Position=[]; ea.Cost=[]; ea.pBest.Position=[]; ea.pBest.Cost=[];
                GB.Cost=inf; pop=repmat(ea,nP,1); pv=cell(nP,1); ok=false;
                while ~ok; for i=1:nP
                    pop(i).Position=CreateRandomSolution(VarSize,Vmin,Vmax);
                    cp=SphericalToCart(pop(i).Position,model);
                    if any(isnan(cp.x))||any(isnan(cp.y))||any(isnan(cp.z)); pop(i).Cost=inf;
                    else; try; pop(i).Cost=CostFn(cp); catch; pop(i).Cost=inf; end; end
                    pop(i).pBest.Position=pop(i).Position; pop(i).pBest.Cost=pop(i).Cost; pv{i}=pop(i).Position;
                    if pop(i).Cost<GB.Cost; GB.Position=pop(i).Position; GB.Cost=pop(i).Cost; ok=true; end
                end; end
                BC=zeros(MaxIt,1);
                for t=1:MaxIt; BC(t)=GB.Cost; eR=0.7*(1-t/MaxIt)^0.5+0.3;
                    for i=1:nP
                        U1=rand(VarSize)>rand();
                        if rand()<eR
                            if rand()<rand(); k=randi(nP); m=randi(nP);
                                gr=(pop(k).pBest.Position.r+pop(m).pBest.Position.r)/2; pop(i).Position.r=pop(i).Position.r+randn(VarSize).*abs(2*rand()*GB.Position.r-gr);
                                gp=(pop(k).pBest.Position.psi+pop(m).pBest.Position.psi)/2; pop(i).Position.psi=pop(i).Position.psi+randn(VarSize).*abs(2*rand()*GB.Position.psi-gp);
                                gh=(pop(k).pBest.Position.phi+pop(m).pBest.Position.phi)/2; pop(i).Position.phi=pop(i).Position.phi+randn(VarSize).*abs(2*rand()*GB.Position.phi-gh);
                            else; k=randi(nP); m=randi(nP);
                                yr=(pop(i).Position.r+pop(k).pBest.Position.r)/2; dr=pop(m).pBest.Position.r-pop(k).pBest.Position.r; pop(i).Position.r=U1.*pop(i).Position.r+(1-U1).*(yr+rand()*dr);
                                yp=(pop(i).Position.psi+pop(k).pBest.Position.psi)/2; dp=pop(m).pBest.Position.psi-pop(k).pBest.Position.psi; pop(i).Position.psi=U1.*pop(i).Position.psi+(1-U1).*(yp+rand()*dp);
                                yh=(pop(i).Position.phi+pop(k).pBest.Position.phi)/2; dh=pop(m).pBest.Position.phi-pop(k).pBest.Position.phi; pop(i).Position.phi=U1.*pop(i).Position.phi+(1-U1).*(yh+rand()*dh);
                            end
                        else
                            Yt=2*rand()*(1-t/MaxIt)^(t/MaxIt); U2=(rand(VarSize)<0.5)*2-1; S=rand()*U2; sc=0; for j=1:nP; sc=sc+pop(j).pBest.Cost; end; sf=sc+eps;
                            if rand()<Tf; St=exp(pop(i).pBest.Cost/sf); S=S.*Yt.*St; k=randi(nP); m=randi(nP);
                                pop(i).Position.r=(1-U1).*pop(i).Position.r+U1.*(pop(k).pBest.Position.r+St*(pop(m).pBest.Position.r-pop(k).pBest.Position.r)-S);
                                pop(i).Position.psi=(1-U1).*pop(i).Position.psi+U1.*(pop(k).pBest.Position.psi+St*(pop(m).pBest.Position.psi-pop(k).pBest.Position.psi)-S);
                                pop(i).Position.phi=(1-U1).*pop(i).Position.phi+U1.*(pop(k).pBest.Position.phi+St*(pop(m).pBest.Position.phi-pop(k).pBest.Position.phi)-S);
                            else; Mt=exp(pop(i).pBest.Cost/sf); k=randi(nP); r2_p=rand();
                                Fr=rand(VarSize).*(Mt*(-pop(i).Position.r+pop(k).pBest.Position.r)); Sr=S.*Yt.*Fr; pop(i).Position.r=GB.Position.r+(al*(1-r2_p)+r2_p)*(U2.*GB.Position.r-pop(i).Position.r)-Sr;
                                Fp=rand(VarSize).*(Mt*(-pop(i).Position.psi+pop(k).pBest.Position.psi)); Sp=S.*Yt.*Fp; pop(i).Position.psi=GB.Position.psi+(al*(1-r2_p)+r2_p)*(U2.*GB.Position.psi-pop(i).Position.psi)-Sp;
                                Fh=rand(VarSize).*(Mt*(-pop(i).Position.phi+pop(k).pBest.Position.phi)); Sh=S.*Yt.*Fh; pop(i).Position.phi=GB.Position.phi+(al*(1-r2_p)+r2_p)*(U2.*GB.Position.phi-pop(i).Position.phi)-Sh;
                            end
                        end
                        pop(i).Position.r=max(pop(i).Position.r,Vmin.r); pop(i).Position.r=min(pop(i).Position.r,Vmax.r);
                        pop(i).Position.psi=max(pop(i).Position.psi,Vmin.psi); pop(i).Position.psi=min(pop(i).Position.psi,Vmax.psi);
                        pop(i).Position.phi=max(pop(i).Position.phi,Vmin.phi); pop(i).Position.phi=min(pop(i).Position.phi,Vmax.phi);
                        cp=SphericalToCart(pop(i).Position,model);
                        if any(isnan(cp.x))||any(isnan(cp.y))||any(isnan(cp.z)); nc=inf; else; try; nc=CostFn(cp); catch; nc=inf; end; end
                        if pop(i).Cost<nc; pop(i).Position=pv{i};
                        else; pv{i}=pop(i).Position; pop(i).Cost=nc;
                            if nc<pop(i).pBest.Cost; pop(i).pBest.Position=pop(i).Position; pop(i).pBest.Cost=nc; end
                            if nc<GB.Cost; GB.Position=pop(i).Position; GB.Cost=nc; end
                        end
                    end
                end
            case 'WOA'
                nP=150; b_w=1; ew.Position=[]; ew.Cost=[];
                GB.Cost=inf; pop=repmat(ew,nP,1); ok=false;
                while ~ok; for i=1:nP
                    pop(i).Position=CreateRandomSolution(VarSize,Vmin,Vmax);
                    cp=SphericalToCart(pop(i).Position,model);
                    if any(isnan(cp.x))||any(isnan(cp.y))||any(isnan(cp.z)); pop(i).Cost=inf;
                    else; try; pop(i).Cost=CostFn(cp); catch; pop(i).Cost=inf; end; end
                    if pop(i).Cost<GB.Cost; GB.Position=pop(i).Position; GB.Cost=pop(i).Cost; ok=true; end
                end; end
                BC=zeros(MaxIt,1);
                for t=1:MaxIt; BC(t)=GB.Cost; a_=2-t*(2/MaxIt); a2_=-1+t*(-1/MaxIt);
                    for i=1:nP
                        r1=rand(); r2=rand(); A_=2*a_*r1-a_; C_=2*r2; l_=(a2_-1)*rand()+1; pp=rand();
                        if pp<0.5
                            if abs(A_)>=1; k=randi(nP);
                                Dr=abs(C_*pop(k).Position.r-pop(i).Position.r); pop(i).Position.r=pop(k).Position.r-A_*Dr;
                                Dp=abs(C_*pop(k).Position.psi-pop(i).Position.psi); pop(i).Position.psi=pop(k).Position.psi-A_*Dp;
                                Dh=abs(C_*pop(k).Position.phi-pop(i).Position.phi); pop(i).Position.phi=pop(k).Position.phi-A_*Dh;
                            else
                                Dr=abs(C_*GB.Position.r-pop(i).Position.r); pop(i).Position.r=GB.Position.r-A_*Dr;
                                Dp=abs(C_*GB.Position.psi-pop(i).Position.psi); pop(i).Position.psi=GB.Position.psi-A_*Dp;
                                Dh=abs(C_*GB.Position.phi-pop(i).Position.phi); pop(i).Position.phi=GB.Position.phi-A_*Dh;
                            end
                        else
                            Dr=abs(GB.Position.r-pop(i).Position.r); pop(i).Position.r=Dr*exp(b_w*l_).*cos(2*pi*l_)+GB.Position.r;
                            Dp=abs(GB.Position.psi-pop(i).Position.psi); pop(i).Position.psi=Dp*exp(b_w*l_).*cos(2*pi*l_)+GB.Position.psi;
                            Dh=abs(GB.Position.phi-pop(i).Position.phi); pop(i).Position.phi=Dh*exp(b_w*l_).*cos(2*pi*l_)+GB.Position.phi;
                        end
                        pop(i).Position.r=max(pop(i).Position.r,Vmin.r); pop(i).Position.r=min(pop(i).Position.r,Vmax.r);
                        pop(i).Position.psi=max(pop(i).Position.psi,Vmin.psi); pop(i).Position.psi=min(pop(i).Position.psi,Vmax.psi);
                        pop(i).Position.phi=max(pop(i).Position.phi,Vmin.phi); pop(i).Position.phi=min(pop(i).Position.phi,Vmax.phi);
                        cp=SphericalToCart(pop(i).Position,model);
                        if any(isnan(cp.x))||any(isnan(cp.y))||any(isnan(cp.z)); pop(i).Cost=inf;
                        else; try; pop(i).Cost=CostFn(cp); catch; pop(i).Cost=inf; end; end
                        if pop(i).Cost<GB.Cost; GB.Position=pop(i).Position; GB.Cost=pop(i).Cost; end
                    end
                end
            case 'CPO'
                nP=150; al=0.2; Tf=0.8;
                ea.Position=[]; ea.Cost=[]; GB.Cost=inf; pop=repmat(ea,nP,1); pv=cell(nP,1); ok=false;
                while ~ok; for i=1:nP
                    pop(i).Position=CreateRandomSolution(VarSize,Vmin,Vmax);
                    cp=SphericalToCart(pop(i).Position,model);
                    if any(isnan(cp.x))||any(isnan(cp.y))||any(isnan(cp.z)); pop(i).Cost=inf;
                    else; try; pop(i).Cost=CostFn(cp); catch; pop(i).Cost=inf; end; end
                    pv{i}=pop(i).Position;
                    if pop(i).Cost<GB.Cost; GB.Position=pop(i).Position; GB.Cost=pop(i).Cost; ok=true; end
                end; end
                BC=zeros(MaxIt,1);
                for t=1:MaxIt; BC(t)=GB.Cost;
                    for i=1:nP
                        U1=rand(VarSize)>rand();
                        if rand()<rand()
                            if rand()<rand(); k=randi(nP);
                                yr=(pop(i).Position.r+pop(k).Position.r)/2; pop(i).Position.r=pop(i).Position.r+randn(VarSize).*abs(2*rand()*GB.Position.r-yr);
                                yp=(pop(i).Position.psi+pop(k).Position.psi)/2; pop(i).Position.psi=pop(i).Position.psi+randn(VarSize).*abs(2*rand()*GB.Position.psi-yp);
                                yh=(pop(i).Position.phi+pop(k).Position.phi)/2; pop(i).Position.phi=pop(i).Position.phi+randn(VarSize).*abs(2*rand()*GB.Position.phi-yh);
                            else; k=randi(nP); m=randi(nP);
                                yr=(pop(i).Position.r+pop(k).Position.r)/2; pop(i).Position.r=U1.*pop(i).Position.r+(1-U1).*(yr+rand()*(pop(m).Position.r-pop(k).Position.r));
                                yp=(pop(i).Position.psi+pop(k).Position.psi)/2; pop(i).Position.psi=U1.*pop(i).Position.psi+(1-U1).*(yp+rand()*(pop(m).Position.psi-pop(k).Position.psi));
                                yh=(pop(i).Position.phi+pop(k).Position.phi)/2; pop(i).Position.phi=U1.*pop(i).Position.phi+(1-U1).*(yh+rand()*(pop(m).Position.phi-pop(k).Position.phi));
                            end
                        else
                            Yt=2*rand()*(1-t/MaxIt)^(t/MaxIt); U2=(rand(VarSize)<0.5)*2-1; S=rand()*U2; sc=0; for j=1:nP; sc=sc+pop(j).Cost; end; sf=sc+eps;
                            if rand()<Tf; St=exp(pop(i).Cost/sf); S=S.*Yt.*St; k=randi(nP); m=randi(nP);
                                pop(i).Position.r=(1-U1).*pop(i).Position.r+U1.*(pop(k).Position.r+St*(pop(m).Position.r-pop(k).Position.r)-S);
                                pop(i).Position.psi=(1-U1).*pop(i).Position.psi+U1.*(pop(k).Position.psi+St*(pop(m).Position.psi-pop(k).Position.psi)-S);
                                pop(i).Position.phi=(1-U1).*pop(i).Position.phi+U1.*(pop(k).Position.phi+St*(pop(m).Position.phi-pop(k).Position.phi)-S);
                            else; Mt=exp(pop(i).Cost/sf); k=randi(nP); r2_p=rand();
                                Fr=rand(VarSize).*(Mt*(-pop(i).Position.r+pop(k).Position.r)); Sr=S.*Yt.*Fr; pop(i).Position.r=GB.Position.r+(al*(1-r2_p)+r2_p)*(U2.*GB.Position.r-pop(i).Position.r)-Sr;
                                Fp=rand(VarSize).*(Mt*(-pop(i).Position.psi+pop(k).Position.psi)); Sp=S.*Yt.*Fp; pop(i).Position.psi=GB.Position.psi+(al*(1-r2_p)+r2_p)*(U2.*GB.Position.psi-pop(i).Position.psi)-Sp;
                                Fh=rand(VarSize).*(Mt*(-pop(i).Position.phi+pop(k).Position.phi)); Sh=S.*Yt.*Fh; pop(i).Position.phi=GB.Position.phi+(al*(1-r2_p)+r2_p)*(U2.*GB.Position.phi-pop(i).Position.phi)-Sh;
                            end
                        end
                        pop(i).Position.r=max(pop(i).Position.r,Vmin.r); pop(i).Position.r=min(pop(i).Position.r,Vmax.r);
                        pop(i).Position.psi=max(pop(i).Position.psi,Vmin.psi); pop(i).Position.psi=min(pop(i).Position.psi,Vmax.psi);
                        pop(i).Position.phi=max(pop(i).Position.phi,Vmin.phi); pop(i).Position.phi=min(pop(i).Position.phi,Vmax.phi);
                        cp=SphericalToCart(pop(i).Position,model);
                        if any(isnan(cp.x))||any(isnan(cp.y))||any(isnan(cp.z)); nc=inf; else; try; nc=CostFn(cp); catch; nc=inf; end; end
                        if pop(i).Cost<nc; pop(i).Position=pv{i};
                        else; pv{i}=pop(i).Position; pop(i).Cost=nc;
                            if nc<GB.Cost; GB.Position=pop(i).Position; GB.Cost=nc; end
                        end
                    end
                end
        end
        results.(alg).bestCosts(run,:)=BC;
        results.(alg).finalCosts(run)=GB.Cost;
        results.(alg).times(run)=toc;
        fprintf('  %s Run %d/%d: Best=%.2f\n', alg, run, N_RUNS, GB.Cost);
    end
end

%% Stats
fprintf('\n=== URBAN MAP RESULTS ===\n');
fprintf('%-8s %8s %8s %8s %8s\n', 'Alg', 'Best', 'Worst', 'Mean', 'Std');
for a=1:5
    alg=algs{a}; fc=results.(alg).finalCosts;
    fprintf('%-8s %8.0f %8.0f %8.0f %8.0f\n', alg, min(fc), max(fc), mean(fc), std(fc));
end
save('results/urban_map_results.mat','results','algs','N_RUNS','MaxIt');

%% Helper
function [X1,X2,X3]=GWOu(AP,BP,DP,X,a,VS)
    A1=2*a*rand(VS)-a; C1=2*rand(VS); X1=AP-A1.*abs(C1.*AP-X);
    A2=2*a*rand(VS)-a; C2=2*rand(VS); X2=BP-A2.*abs(C2.*BP-X);
    A3=2*a*rand(VS)-a; C3=2*rand(VS); X3=DP-A3.*abs(C3.*DP-X);
end
