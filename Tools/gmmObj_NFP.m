function [Q, Dg, Qm, Qiv, DgM, DgIV]=gmmObj_NFP(theta2)

global Model 

expdelta=invertmarketshares(exp(Model.delta),theta2);% Invert to get mean utilities

if nargout==1
    [~, ~, ~,Qm] =sim_model_NFP(expdelta,theta2);
else
    [dS_ddelta,dS_dTheta] =sim_model_NFP_DD(expdelta,theta2);
    if Model.AggregateMoments==1
        [DgM, ~, Qm]= gradient_MicroM_NFP_Ag(expdelta,theta2, dS_dTheta,dS_ddelta);
       
    else
        [DgM, ~, Qm]= gradient_MicroM_NFP(expdelta,theta2, dS_dTheta,dS_ddelta);
      
    end

end
      
n=sum(Model.useObs);
delta0=log(expdelta(Model.useObs));
delta0=Model.A_FE*delta0;
%mdliv=fitlm(Model.XZ,delta0,'Weights',Model.wtMeanReg);
%deltaHat = predict(mdliv,Model.XZ);
%resid0=Model.wtMeanReg .* (delta0-deltaHat);
theta1=Model.XZ\delta0;
resid0=delta0 - Model.XXX*theta1;
gIV=(1/n)*(resid0'*Model.IV)' - Model.AdjustmentObjIV;
Qiv=gIV'*Model.WBLP*gIV;

Q=Model.mWm*Qm+Model.mWiv*Qiv;

%Model.delta=log(expdelta); % store current delta
if nargout>=2 
    dDeltadThetanorm=zeros(size(dS_dTheta));
    for m=1:Model.nMarkets % Parallel work across markets
        for t=1:Model.nTime
            % Setup Index
            index=(Model.SIndex(:,1)==Model.Markets(m) & Model.SIndex(:,2)==Model.Years(t));
            index(index==1 & Model.normID==1)=0;
            
            TEMP=-inv(dS_ddelta(index,index))*dS_dTheta(index,:);
            dDeltadThetanorm(index,:)=TEMP;
        end
    end
    JIV=(Model.IPX*Model.A_FE*dDeltadThetanorm(Model.useObs,:))*(1/n);
    DgIV=2*JIV'*Model.IV*Model.WBLP*gIV;
    Dg=Model.mWm*DgM+Model.mWiv*DgIV;
end

