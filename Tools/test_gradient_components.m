  %% Get numerical dsdTheta
    
    theta2=Theta2_True + randn(length(Theta2_True),1)/100;
    [expdelta0, flag]=invertmarketshares(ones(length(Model.delta),1),theta2);
    Model.delta=log(expdelta0);
    [S, Shares, M, G_M, dS_ddelta,dS_dTheta2] =sim_model_NFP(expdelta0,theta2);
    
    step=10^-6;
    for i=1:length(theta2)
        
        thetaA=theta2;
        thetaB=theta2;
        
        thetaA(i)=thetaA(i)-step/2;
        thetaB(i)=thetaB(i)+step/2;
        
        [sharesA ]=sim_model_NFP(exp(Model.delta),thetaA);
        [sharesB ]=sim_model_NFP(exp(Model.delta),thetaB);
        
        dSdthetaNUM(:,i)=(sharesB-sharesA)/step;
    end
    
    dif1(1)=max(max(abs(dSdthetaNUM(Model.normID==0,:)-dS_dTheta2(Model.normID==0,:))));
    dif1(2)=max(max(abs((dSdthetaNUM-dS_dTheta2)./dS_dTheta2)));
    
    fprintf('\nMax   difference between numerical dsdTheta and analytical version = \nMax Time = %10.3e  \n',dif1(1))
    fprintf('\nMax percent difference between numerical dsdTheta and analytical version = \nMax Time = %10.3e  \n',dif1(2))
    
    
    %% Get numberical dsdDelta
    step=10^-6;
    dSddeltaNUM=NaN(size(dS_ddelta));
    for i=1:length(Model.delta)
        
        deltaA=Model.delta;
        deltaB=Model.delta;
        
        deltaA(i)=deltaA(i)-step/2;
        deltaB(i)=deltaB(i)+step/2;
        
        [sharesA ]=sim_model_NFP(exp(deltaA),theta2);
        [sharesB ]=sim_model_NFP(exp(deltaB),theta2);
        
        dSddeltaNUM(:,i)=(sharesB-sharesA)/step;
    end
    
    dif2(1)=max(max(abs(dSddeltaNUM-dS_ddelta)))  ;             % error
    dif2(2)=max(max(abs((dSddeltaNUM-dS_ddelta)./dS_ddelta))) ; % percent error
    
    fprintf('\nMax   difference between numerical dsdDelta  and analytical version = \nMax Time = %10.3e  \n',dif2(1))
    fprintf('\nMax percent difference between numerical dsdDelta  and analytical version = \nMax Time = %10.3e  \n',dif2(2))
    
    
    
    %% Get numerical ddelta_dtheta
    step=10^-6;
    ddelta_dthetaNUM=NaN(size(Model.SS,1),12);
    
    for i=1:length(theta2)
        
        thetaA=theta2;
        thetaB=theta2;
        
        thetaA(i)=thetaA(i)-step/2;
        thetaB(i)=thetaB(i)+step/2;
        
        expdeltaA=invertmarketshares(exp(Model.delta),thetaA);
        expdeltaB=invertmarketshares(exp(Model.delta),thetaB);
        
        ddelta_dthetaNUM(:,i)=(log(expdeltaB)-log(expdeltaA))/step;
    end
    
    % Get analytical dDeltadTheta by adjusting the regular formula for the outside option delta changing as well.
    
    dDeltadThetanorm=zeros(size(dS_dTheta));
    for m=1:Model.nMarkets % Parallel work across markets
        for t=1:Model.nTime
            % Setup Index
            index=(Model.SIndex(:,1)==Model.Markets(m) & Model.SIndex(:,2)==Model.Years(t));
            indexNorm=(index==1 & Model.normID==1);
            index(index==1 & Model.normID==1)=0;
            
            TEMP=-inv(dS_ddelta(index,index))*dS_dTheta(index,:); %-dS_ddelta(index,indexNorm)*dS_dTheta(indexNorm,:)./dS_ddelta(indexNorm,indexNorm)
            dDeltadThetanorm(index,:)=TEMP;
        end
    end
    
    dif3(1)=max(max(abs(ddelta_dthetaNUM-dDeltadThetanorm)));
    dif3(2)=max(max(abs(dDeltadThetanorm-ddelta_dthetaNUM)./dDeltadThetanorm));
    
    fprintf('\nMax difference between numerical dDeltadTheta  and analytical version = \nMax Time = %10.3e  \n',dif3(1))
    fprintf('\nMax percent difference between numerical dDeltadTheta  and analytical version = \nMax Time = %10.3e  \n',dif3(2))
    
    
    %% Testing Gradient IV Moments
    step=10^-8;
    
    [expdelta0, flag, norm_maxShares, norm_max,  iter]=invertmarketshares(exp(Model.delta),theta2);
    delta=log(expdelta0);
    Model.delta=delta;
    
    % Analytical Gradient for IV moments
    [Q, Dg, Qm, Qiv, DgM, DgIV ]=gmmObj_NFP(theta2);
    [DgIVb, G, JA1]= gradient_BLP_NFP(theta2);
    
    Jnum1=NaN(size(JA1));
    for i=1:length(theta2)
        
        theta2A=theta2;
        theta2B=theta2;
        
        theta2A(i)=theta2A(i)-step/2;
        theta2B(i)=theta2B(i)+step/2;
        
        expdeltaA=invertmarketshares(exp(Model.delta),theta2A);
        deltaA=log(expdeltaA);
        theta1A=Model.XZ\deltaA;
        residA=deltaA-Model.XZ*theta1A;
        gIVa=(residA'*Model.IV)'-Model.AdjustmentObjIV;
        
        expdeltaB=invertmarketshares(exp(Model.delta),theta2B);
        deltaB=log(expdeltaB);
        theta1B=Model.XZ\deltaB;
        residB=deltaB-Model.XZ*theta1B;
        gIVb=(residB'*Model.IV)'-Model.AdjustmentObjIV;
        
        Jnum1(:,i)=(residB-residA)/step;
        
        Qa=gIVa'*Model.WBLP*gIVa;
        Qb=gIVb'*Model.WBLP*gIVb;
        
        [QA,~,~, QivA]=gmmObj_NFP(theta2A);
        [QB,~,~, QivB]=gmmObj_NFP(theta2B);
        
        DQiv(i,1)=(Qb-Qa)/step;
        DQall(i,1)=(QB-QA)/step;
        
    end
    
    [expdelta0, flag, norm_maxShares, norm_max,  iter]=invertmarketshares(exp(Model.delta),theta2);
    delta=log(expdelta0);
    
    [DgIV, QIV, JacIV]= gradient_BLP_NFP(theta2);
    [DgM, M, QM,JacM] = gradient_MicroM_NFP(expdelta0,theta2, dS_dTheta,dS_ddelta);
    
    dif4(1)=max(abs(DQiv-DgIV))
    dif4(2)=max(abs(DQiv-DgIV)./DgIV)
    
    max(abs(Jnum1-JacIV))
    max(abs(Jnum1-JacIV)./Jnum1)
    
    
    fprintf('\nMax difference between IV numerical gradient and analytical version = \nMax Time = %10.3e  \n',dif4(1))
    fprintf('\nMax percent difference between IV numerical gradient  and analytical version = \nMax Time = %10.3e  \n',dif4(2))
    
    
    
    Jnum2=Model.IPX*ddelta_dthetaNUM;
    Jnum3=Model.IPX*dDeltadThetanorm;
    max(abs(Jnum2-Jnum3))
    max(abs(Jnum2-JacIV)./Jnum2)
    
    
    theta1=Model.XZ\Model.delta;
    residTest=delta-Model.XZ*theta1;
    gIV=(residTest'*Model.IV)'-Model.AdjustmentObjIV;
    DgN=2*(Jnum1)'*Model.IV*Model.WBLP*(gIV);
    DgN-DgIV
    
    DgN=2*(JacIV)'*Model.IV*Model.WBLP*(gIV);
    DgN-DgIV
    
    
    [Q, Dg, Q1, Q2, DgM, DgIV ]=gmmObj_NFP(theta2);
    
    
    dif4(3)=max(abs(DgN-DgIV));
    dif4(4)=max(abs(DgN-DgIV)./DgIV);
    
    fprintf('\nMax difference between IV numerical gradient and analytical version = \nMax Time = %10.3e  \n',dif4(3))
    fprintf('\nMax percent difference between IV numerical gradient  and analytical version = \nMax Time = %10.3e  \n',dif4(4))
    
    
    %% Testing Gradient micromoments
    
    % Try new Theta2 that is further from the optimal.
    theta2=Theta2_True; %zeros(size(Theta2_True));
    
    theta2=Theta2_True + randn(length(Theta2_True),1)/100;
    
    expdelta=invertmarketshares(ones(length(Model.SS),1),theta2);
    Model.delta=log(expdelta);
    [~, Shares, M, ~, dS_ddelta,dS_dTheta] =sim_model_NFP(expdelta,theta2);
    [DgM2, M2, G2,J2] = gradient_MicroM_NFP(expdelta,theta2,dS_dTheta,dS_ddelta);
    g0=M-Model.MM;
    g0((Model.mAllIndex(:,end)==0))=[];
    Q1=g0'*Model.WMM*g0;
    
    expdelta=invertmarketshares(ones(length(Model.SS),1),theta2);
    Model.delta=log(expdelta);
    [Qv, Dgv, Q1v, Q2v, DgMv, DgIVv ]=gmmObj_NFP(theta2);
    
    step=10^-8;
    DGNUM=NaN(length(theta2),1);DM=NaN(length(M),length(theta2));DALLNUM=NaN(length(theta2),1);DAiv=NaN(length(theta2),1);
    for i=1:length(theta2)
        thetaA=theta2;
        thetaB=theta2;
        
        thetaA(i)=thetaA(i)-step/2;
        thetaB(i)=thetaB(i)+step/2;
        
        [expdeltaA flagA]=invertmarketshares(expdelta,thetaA);
        [~, ~, Ma, gA]=sim_model_NFP(expdeltaA,thetaA);
        [~, QivA]= gradient_BLP_NFP(thetaA);
        Qa=gmmObj_NFP(thetaA);
        
        [expdeltaB flagB]=invertmarketshares(expdelta,thetaB);
        [~, ~, Mb, gB]=sim_model_NFP(expdeltaB,thetaB);
        [~, QivB]= gradient_BLP_NFP(thetaB);
        Qb=gmmObj_NFP(thetaB);
        
        DM(:,i)=(Mb-Ma)/step;
        DGNUM(i)=(gB-gA)/step;
        DALLNUM(i,1)=(Qb-Qa)/step;
        
        
        DAiv(i,1)=(QivB-QivA)/step;
        
    end
    
    
    [Q, Dg, Q1, Q2, DgM, DgIV ]=gmmObj_NFP(theta2);
    
    dif5(1)=max(abs(DAiv-DgIV)./DgIV);
    dif5(2)=max(abs(DGNUM-DgM));
    dif5(3)=max(abs(DGNUM-DgM)./DgM);
    
    fprintf('\nMax difference between Micro Moments numerical gradient and analytical version = \nMax Time = %10.3e  \n',dif5(1))
    fprintf('\nMax percent difference between Micro Moments numerical gradient  and analytical version = \nMax Time = %10.3e  \n',dif5(2))
    
    dif6(4)=max(max(abs(DALLNUM-Dg)));
    dif6(5)=max(max(abs((DALLNUM-Dg)./Dg)));
    
    fprintf('\nMax   difference between numerical gradient and analytical version = \nMax Time = %10.3e  \n',dif6(1))
    fprintf('\nMax percent difference between numerical gradient and analytical version = \nMax Time = %10.3e  \n',dif6(2))
    
    
    [DgM2, M2, G2,J2] = gradient_MicroM_NFP(expdelta,theta2,dS_dTheta,dS_ddelta);
    dif5(6)=max(abs(DGNUM-DgM2));