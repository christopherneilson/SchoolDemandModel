
function [SE, Vmm, Viv, Wiv, Wmicro,dropParameter,Out_MIndex]=getVarCov(Theta2Star, SModel,students,option)

%% Variance Covariance from IV moments 

[expdelta, flag, norm_maxShares, norm_max,  iter]=invertmarketshares(exp(SModel.delta),Theta2Star);
[eS, eShares, eM, eG_M, dS_ddelta,dS_dTheta]=sim_model_NFP(expdelta,Theta2Star);


% Get analytical dDeltadTheta by adjusting the regular formula for the
% outside option delta changing as well. 
dDeltadThetanorm=zeros(size(dS_dTheta));
for m=1:SModel.nMarkets % Parallel work across markets
    for t=1:SModel.nTime
        % Setup Index
        index=(SModel.SIndex(:,1)==SModel.Markets(m) & SModel.SIndex(:,2)==SModel.Years(t));
        index(index==1 & SModel.normID==1)=0;
        
        TEMP=-inv(dS_ddelta(index,index))*dS_dTheta(index,:); 
        dDeltadThetanorm(index,:)=TEMP;
    end
end


n=length(SModel.S);
delta0=log(expdelta);
SModel.theta1=SModel.XZ\delta0;
SModel.resid0=delta0-SModel.XXX*SModel.theta1;

g=bsxfun(@times,SModel.resid0,SModel.IV);
gstarIV=bsxfun(@minus,g,mean(g));
Viv=gstarIV'*gstarIV/n;
Wiv=inv(Viv);

Jiv=SModel.IV'*[-SModel.XX dDeltadThetanorm];
dropParameter=max(abs(Jiv))<10^-4;
Jiv(:,dropParameter==1)=[];

SEiv=sqrt(diag(inv(Jiv'*SModel.WBLP*Jiv)*Jiv'*SModel.WBLP*Viv*SModel.WBLP*Jiv*inv(Jiv'*SModel.WBLP*Jiv)));

%% Variance from Micro Moments 

% Collect microdata used to calculate moments. 
 %students.Properties.VariableNames;
N1=NaN(length(eM),1);
NN=N1;
V1=zeros(length(eM));
VV=V1;
W1=VV;
Qnorm=SModel.Q(SModel.normIndex);
DropMM=NaN(length(eM),1);
for year=reshape(SModel.Years,1,length(SModel.Years))
    for market=reshape(SModel.Markets,1,length(SModel.Markets))
        
    
    for type=1:6 
    momentIndex=(SModel.mAllIndex(:,1)==market & SModel.mAllIndex(:,2)==year & SModel.mAllIndex(:,3)==type);
    microdataIndex=find(market==students.Market & students.Year==year & students.Type==type);
        
    qnorm=Qnorm(SModel.SIndex(:,1)==market & SModel.SIndex(:,2)==year);
    
    microdataVA=students.VA2_AVE(microdataIndex)-qnorm(1);
    microdataPrice=students.AvePrice(microdataIndex)/1000;
    microdataDist=students.Dist2School(microdataIndex);
    
    % Calculate contribution of each datapoint to the varcov of moments
    
    momentValueVA=eM(momentIndex==1 & SModel.mAllIndex(:,4)==1);
    momentValuePrice=eM(momentIndex==1 & SModel.mAllIndex(:,4)==2);
    momentValueDist=eM(momentIndex==1 & SModel.mAllIndex(:,4)==3);
    
    dev_1=bsxfun(@minus,microdataVA-momentValueVA,nanmean(microdataVA-momentValueVA));
    dev_2=bsxfun(@minus,microdataPrice-momentValuePrice,nanmean(microdataPrice-momentValuePrice));
    dev_3=bsxfun(@minus,microdataDist-momentValueDist,nanmean(microdataDist-momentValueDist));
    
    n1=sum(isnan(dev_1)==0);    
    n2=sum(isnan(dev_2)==0);
    n3=sum(isnan(dev_3)==0);
    
    Simga_i=zeros(3,3);
    Simga_i(1,1)=nanmean(dev_1.*dev_1); d1=(Simga_i(1,1)<10^-8 | n1<10);
    Simga_i(2,2)=nanmean(dev_2.*dev_2); d2=(Simga_i(2,2)<10^-8 | n2<10);   
    Simga_i(3,3)=nanmean(dev_3.*dev_3); d3=(Simga_i(3,3)<10^-8 | n3<10);
    
    
    Simga_i(1,2)=nanmean(dev_1.*dev_2);
    Simga_i(2,1)=Simga_i(1,2);
    
    Simga_i(1,3)=nanmean(dev_1.*dev_3);
    Simga_i(3,1)=Simga_i(1,3);
        
    Simga_i(2,3)=nanmean(dev_2.*dev_3);
    Simga_i(3,2)=Simga_i(2,3);
     
    momentLoc=find(SModel.mAllIndex(:,1)==market & SModel.mAllIndex(:,2)==year & SModel.mAllIndex(:,3)==type);
    
    % Make robust and drop moments
    drop=[d1;d2;d3];
    Sigma=NaN(3,3);
    invSigma=NaN(3,3);
    Stemp=Simga_i(drop==0,drop==0);
    
    if rcond(Stemp)>10^-4
    Sigma(drop==0,drop==0)=Stemp;
    invSigma(drop==0,drop==0)=inv(Stemp);
    W1(momentLoc, momentLoc)=invSigma;
    VV(momentLoc, momentLoc)=Sigma; 
    else
    drop=[1;1;1];    
    end
    
    
    
    
    NN(momentLoc)=[n1;n2;n3];
    DropMM(momentLoc)=drop;
   
    end
    end
end

skipMoment=(SModel.mAllIndex(:,end)==0 | DropMM==1);

% Fix to add off diag 
V1=VV(skipMoment==0,skipMoment==0);
N1=NN(skipMoment==0);
Wmicro=W1(skipMoment==0,skipMoment==0); 
Vmm=V1;

if option==0  
SE=[];
elseif option==1 %% Get Variance-Covariance Matrix (sandwich formula)

% Get numerical J for micromoments 
[Dg, M, G,Jm, dMdEps] = gradient_MicroM_NFP(exp(SModel.delta),Theta2Star, dS_dTheta,dS_ddelta);

Jmicro=[dMdEps*SModel.XX Jm];
Jmicro(:,dropParameter==1)=[];
Jmicro(DropMM(SModel.mAllIndex(:,end)==1)==1,:)=[];
JJ=[Jiv; Jmicro];
WW=blkdiag(SModel.WBLP,SModel.WMM(DropMM(SModel.mAllIndex(:,end)==1)==0,DropMM(SModel.mAllIndex(:,end)==1)==0));
VVV=blkdiag(Viv,Vmm);

Omega=inv(JJ'*WW*JJ)*JJ'*WW*VVV*WW'*JJ*inv(JJ'*WW*JJ);
SE=sqrt(diag(Omega));

Out_MIndex=(skipMoment==0);

% JJ=[Jiv];
% WW=SModel.WBLP;
% VV=Viv;
% 
% 
% Omega=inv(JJ'*WW*JJ)*JJ'*WW*VV*WW'*JJ*inv(JJ'*WW*JJ);
% SE=sqrt(diag(Omega));
% 
% 
% JJ=[Jmicro];
% WW=SModel.WMM;
% VV=Vmm;
% 
% Omega=inv(JJ'*WW*JJ)*JJ'*WW*VV*WW'*JJ*inv(JJ'*WW*JJ);
% SE=sqrt(diag(Omega));


end





