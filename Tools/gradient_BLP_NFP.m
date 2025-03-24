function [DgIV, Qiv, JIV]= gradient_BLP_NFP(theta2)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
 
global Model

expdelta=invertmarketshares(exp(Model.delta),theta2);% Invert to get mean utilities


delta0=log(expdelta);
delta0=Model.A_FE*delta0;
theta1=Model.XZ\delta0;
resid0=delta0-Model.XXX*theta1;
gIV=(resid0'*Model.IV)'-Model.AdjustmentObjIV;
Qiv=gIV'*Model.WBLP*gIV;

if nargout>=2 
    [dS_ddelta,dS_dTheta] =sim_model_NFP_DD(expdelta,theta2);
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
    JIV=(Model.IPX*Model.A_FE*dDeltadThetanorm);
    DgIV=2*JIV'*Model.IV*Model.WBLP*gIV;
   
  
end
 
