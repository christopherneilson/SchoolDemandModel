function [Dg, eeM, G,J,  DM] = gradient_MicroM_NFP(expdelta0,Theta2, dS_dTheta,dS_ddelta)
% This function calculates the gradient of the GMM objective function of  the random coefficients model
% using only micromoments described in Neilson 2014, estimated via Nested Fixed Point

% Source: Neilson (2013)
% Code Revised: 6/11/2013

global  Model 

% Note: Originally thought to take shares and not recalculate them. Needs to get this update. 

%% Gradient Specific Variables
if nargin==1
    [expdelta]=invertmarketshares(expdelta0,Theta2);
    [dS_ddelta,dS_dTheta] =sim_model_NFP_DD(expdelta,Theta2);
end

nM=size(Model.mAllIndex,1); % Potential Moments
eeM=zeros(nM,1);

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

%% Share equation variables
% Apply Normalization by Market x Year
expdelta=norm_delta(expdelta0,Model.normIndex);
if (expdelta(Model.normID==1)~=1)
    fprintf('Not normed correctly')
end

nTheta2=Model.nTheta2;

BetaEdu=Theta2(Model.BetaEduO);
BetaSEP=Theta2(Model.BetaSEPO);
AlphaEdu=Theta2(Model.AlphaEduO);
AlphaSEP=Theta2(Model.AlphaSEPO);
LambdaEdu=Theta2(Model.LambdaEduO);
LambdaSEP=Theta2(Model.LambdaSEPO);
TypeD = Model.TypeD;
if TypeD==2 
    LambdaSQ = Theta2(Model.LambdaSQ);
else 
    LambdaSQ = 0;    
end
if TypeD==4 || TypeD==6 
    LambdaSplHi = Theta2(Model.LambdaSplHi);
else 
    LambdaSplHi = 0;    
end
if TypeD==5 || TypeD==6 
    LambdaSplLo = Theta2(Model.LambdaSplLo);
else 
    LambdaSplLo = 0;    
end
if TypeD==3
    LambdaSantiago=Theta2(Model.LambdaSantiago);      
else
	LambdaSantiago = 0; 
end
DistanceSplineLo = Model.DistanceSplineLo;
DistanceSplineHi = Model.DistanceSplineHi;
if TypeD==7 || TypeD==8 || TypeD==9
    LambdaEdu12=Theta2(Model.LambdaEduO12);
    LambdaSEP12=Theta2(Model.LambdaSEPO12);    
    LambdaEdu23=Theta2(Model.LambdaEduO23);
    LambdaSEP23=Theta2(Model.LambdaSEPO23);  
    LambdaEdu34=Theta2(Model.LambdaEduO34);
    LambdaSEP34=Theta2(Model.LambdaSEPO34);  
    LambdaEdu4p=Theta2(Model.LambdaEduO4p);
    LambdaSEP4p=Theta2(Model.LambdaSEPO4p);      
    LambdaDiffO = [LambdaEdu12(1)+LambdaSEP12    LambdaEdu23(1)+LambdaSEP23    LambdaEdu34(1)+LambdaSEP34    LambdaEdu4p(1)+LambdaSEP4p;
                   LambdaEdu12(1)                LambdaEdu23(1)                LambdaEdu34(1)                LambdaEdu4p(1)            ;
                   LambdaEdu12(2)+LambdaSEP12    LambdaEdu23(2)+LambdaSEP23    LambdaEdu34(2)+LambdaSEP34    LambdaEdu4p(2)+LambdaSEP4p;    
                   LambdaEdu12(2)                LambdaEdu23(2)                LambdaEdu34(2)                LambdaEdu4p(2)            ;
                   LambdaEdu12(3)                LambdaEdu23(3)                LambdaEdu34(3)                LambdaEdu4p(3)            ; 
                   LambdaEdu12(4)                LambdaEdu23(4)                LambdaEdu34(4)                LambdaEdu4p(4)            ];
else  
    LambdaDiffO = 0;
end

BetaO=[BetaSEP               AlphaEdu(1)+AlphaSEP  LambdaEdu(1)+LambdaSEP   ; % No HS - SEP
    0                     AlphaEdu(1)           LambdaEdu(1)             ; % No HS - Not SEP
    BetaEdu(1)+BetaSEP    AlphaEdu(2)+AlphaSEP  LambdaEdu(2)+LambdaSEP   ; % HS - SEP
    BetaEdu(1)            AlphaEdu(2)           LambdaEdu(2)             ; % HS - Not SEP
    BetaEdu(2)            AlphaEdu(3)           LambdaEdu(3)             ; % Tech Moms
    BetaEdu(3)            AlphaEdu(4)           LambdaEdu(4)             ];% College Moms
SigRC=Theta2(Model.SigRC);

nVi=length(Model.q_Nodes);        % Number of simulated consumers
if Model.lognormalRC==1
    ViNodes=exp(SigRC.'*Model.q_Nodes);
    dViNodes=Model.q_Nodes .* ViNodes;    
else
    ViNodes = SigRC.'*Model.q_Nodes;    
    dViNodes = Model.q_Nodes;     
end
ViWeights=Model.q_Weights;

P=Model.P;
Q=Model.Q;
SEP=Model.SEP;
Q(isnan(Q))=0; % fix this in initial data


% Calc parallel by Market-Year
WNodes=Model.wNodes;
Qnorm=Q(Model.normIndex);
D=Model.D;
locSchoolDist=Model.locSchoolDist;
locNodeDist=Model.locNodeDist;
sList=Model.mapNodeList;


piW=Model.piW;
counterMarket=Model.counterMarket;
nodeMarket = Model.nodeMarket;

nS=size(Model.SIndex,1);

S=zeros(nS,1);
if TypeD==8 || TypeD==9
  eM_temp=zeros(max(counterMarket),24);  
else
  eM_temp=zeros(max(counterMarket),18); 
end
eM=eM_temp;

dMdTheta=cell(length(locNodeDist),1);
dMdEps=cell(length(locNodeDist),1);

%for i=1:length(locNodeDist)
parfor i=1:length(locNodeDist)
    ni=locNodeDist(i);
    slist=sList(ni,:);
    market = nodeMarket(i);   
    
    dist_ni=D.Value(:,ni);
    %dist_ni=Model.distanceMatrix(:,ni);
    wn=WNodes(ni,:);
    weights=wn(:,5:end);
    
    dmdtheta=cell(length(sList),1);
    dmdeps=cell(length(sList),1);
    
  
    
    it=1;
    for s=slist % Parallel work across markets
        pickS=(s==counterMarket);
        
        
        delta_mt=log(expdelta(pickS));
        nFirms=length(delta_mt);
        
        % initialize placeholders
        em=eM_temp;
        
        p=P(pickS);
        q=Q(pickS);
        qq=q-Qnorm(pickS);
        sep=SEP(pickS);
        pi=piW(s,:);
        
        dni=dist_ni(locSchoolDist(pickS));
        if TypeD==4 || TypeD==6 
            dniHi = (dni - DistanceSplineHi) .* (dni>DistanceSplineHi);
        else
            dniHi = 0;
        end
        if TypeD==5 || TypeD==6 
            dniLo = (dni - DistanceSplineLo) .* (dni<DistanceSplineLo);
        else
            dniLo = 0;
        end 
        if TypeD==7 || TypeD==8 || TypeD==9
            dni1p = (dni - 1) .* (dni>=1);
            dni2p = (dni - 2) .* (dni>=2); 
            dni2l = (dni - 2) .* (dni<2);      
            dni3p = (dni - 3) .* (dni>=3);
            sni3p = (dni>=3);            
            dni4p = (dni - 4) .* (dni>=4);
        else
            dni1p = 0; dni2p = 0; dni2l = 0; dni3p = 0; dni4p = 0; sni3p=0;
        end  
        
        if TypeD==8 || TypeD==9
            mType=zeros(6,4);       
            dmoTypes=zeros(6,4*nTheta2);
            dmdepsTypes=zeros(6,4*nFirms);
        else
            mType=zeros(6,3);
            dmoTypes=zeros(6,3*nTheta2);
            dmdepsTypes=zeros(6,3*nFirms);
        end
        
        % Get Aggregate Derivatives for this Market-Time
        ddeltadtheta=dDeltadThetanorm(pickS,:);
        
        
        for type=1:6
            wnode=weights(:,type);
            % Adjust out of pocket price for SEP
            if (type==1 || type==3)
                pp=p;
                pp(sep==1)=0;
            else
                pp=p;
            end
            % Set up interaction between firm Xs and unobservables
            if TypeD==2
                UoType=qq*BetaO(type,1)+pp*BetaO(type,2)+dni*BetaO(type,3) + (dni.^2) * LambdaSQ;
            elseif TypeD==3 && market==312
                UoType=qq*BetaO(type,1)+pp*BetaO(type,2)+dni*BetaO(type,3) + dni*LambdaSantiago;
            elseif TypeD==4 
                UoType=qq*BetaO(type,1)+pp*BetaO(type,2)+dni*BetaO(type,3) + dniLo * LambdaSplHi;                
            elseif TypeD==5 
                UoType=qq*BetaO(type,1)+pp*BetaO(type,2)+dni*BetaO(type,3) + dniLo * LambdaSplLo;                
            elseif TypeD==6 
                UoType=qq*BetaO(type,1)+pp*BetaO(type,2)+dni*BetaO(type,3) + dniLo * LambdaSplLo + dniHi * LambdaSplHi;                                      
            elseif TypeD==7 || TypeD==8 || TypeD==9
                UoType=qq*BetaO(type,1)+pp*BetaO(type,2)+dni*BetaO(type,3);  
                UoType = UoType + LambdaDiffO(type,1) * dni1p;
                UoType = UoType + LambdaDiffO(type,2) * dni2p;
                UoType = UoType + LambdaDiffO(type,3) * dni3p;
                UoType = UoType + LambdaDiffO(type,4) * dni4p;                
            else 
                UoType=qq*BetaO(type,1)+pp*BetaO(type,2)+dni*BetaO(type,3);
            end
            Uv=ViNodes.*qq; % Random Coefficient Part of Utility
            %get maxutility to norm vector of utility to avoid overflow
            maxuij=max(delta_mt)+max(max(max(UoType)))+max(max(max(Uv)));
            
            num=exp(delta_mt+UoType+Uv-maxuij);
            denom=sum(num,1); % Note, no 1 added.
            share_ijv=num./repmat(denom,[nFirms 1 1]);
            
            
            if TypeD==8 || TypeD==9
                dmoVi=zeros(1,4*nTheta2);
                dmdepsVi=zeros(1,4*nFirms);                  
            else
                dmoVi=zeros(1,3*nTheta2);
                dmdepsVi=zeros(1,3*nFirms);
            end
            
            for vik=1:nVi
                
                si=share_ijv(:,vik);
                ds_ddelta_mt=(diag(si)-si*si');
                ds_dtheta2_mt=zeros(nFirms,nTheta2);
                
                if type==1 % No HS: SEP
                    ds_dtheta2_mt(:,4)=si.*(qq-si'*qq);
                    ds_dtheta2_mt(:,5)=si.*(pp-si'*pp);
                    ds_dtheta2_mt(:,9)=si.*(pp-si'*pp);
                    ds_dtheta2_mt(:,10)=si.*(dni-si'*dni);
                    ds_dtheta2_mt(:,14)=si.*(dni-si'*dni);
                    if TypeD==7 || TypeD==8 || TypeD==9
                        ds_dtheta2_mt(:,15)=si.*(dni1p-si'*dni1p);                            
                        ds_dtheta2_mt(:,19)=si.*(dni1p-si'*dni1p);         
                        ds_dtheta2_mt(:,20)=si.*(dni2p-si'*dni2p);                            
                        ds_dtheta2_mt(:,24)=si.*(dni2p-si'*dni2p);       
                        ds_dtheta2_mt(:,25)=si.*(dni3p-si'*dni3p);                            
                        ds_dtheta2_mt(:,29)=si.*(dni3p-si'*dni3p);       
                        ds_dtheta2_mt(:,30)=si.*(dni4p-si'*dni4p);                            
                        ds_dtheta2_mt(:,34)=si.*(dni4p-si'*dni4p);            
                    end
                elseif type==2 % No HS; No SEP
                    ds_dtheta2_mt(:,5)=si.*(pp-si'*pp);
                    ds_dtheta2_mt(:,10)=si.*(dni-si'*dni);
                    if TypeD==7 || TypeD==8 || TypeD==9
                        ds_dtheta2_mt(:,15)=si.*(dni1p-si'*dni1p);                                 
                        ds_dtheta2_mt(:,20)=si.*(dni2p-si'*dni2p);                            
                        ds_dtheta2_mt(:,25)=si.*(dni3p-si'*dni3p);                                  
                        ds_dtheta2_mt(:,30)=si.*(dni4p-si'*dni4p);                                    
                    end
                elseif type==3 % HS; SEP
                    ds_dtheta2_mt(:,1)=si.*(qq-si'*qq);
                    ds_dtheta2_mt(:,4)=si.*(qq-si'*qq);
                    ds_dtheta2_mt(:,6)=si.*(pp-si'*pp);
                    ds_dtheta2_mt(:,9)=si.*(pp-si'*pp);
                    ds_dtheta2_mt(:,11)=si.*(dni-si'*dni);
                    ds_dtheta2_mt(:,14)=si.*(dni-si'*dni);
                    if TypeD==7 || TypeD==8 || TypeD==9
                        ds_dtheta2_mt(:,16)=si.*(dni1p-si'*dni1p);                            
                        ds_dtheta2_mt(:,19)=si.*(dni1p -si'*dni1p);         
                        ds_dtheta2_mt(:,21)=si.*(dni2p-si'*dni2p);                            
                        ds_dtheta2_mt(:,24)=si.*(dni2p-si'*dni2p);       
                        ds_dtheta2_mt(:,26)=si.*(dni3p-si'*dni3p);                            
                        ds_dtheta2_mt(:,29)=si.*(dni3p-si'*dni3p);       
                        ds_dtheta2_mt(:,31)=si.*(dni4p-si'*dni4p);                            
                        ds_dtheta2_mt(:,34)=si.*(dni4p-si'*dni4p);            
                    end
                elseif type==4 % HS; No SEP
                    ds_dtheta2_mt(:,1)= si.*(qq-si'*qq);
                    ds_dtheta2_mt(:,6)= si.*(pp-si'*pp);
                    ds_dtheta2_mt(:,11)= si.*(dni-si'*dni);
                    if TypeD==7 || TypeD==8 || TypeD==9                      
                        ds_dtheta2_mt(:,16)=si.*(dni1p-si'*dni1p);                            
                        ds_dtheta2_mt(:,21)=si.*(dni2p-si'*dni2p);                            
                        ds_dtheta2_mt(:,26)=si.*(dni3p-si'*dni3p);                            
                        ds_dtheta2_mt(:,31)=si.*(dni4p-si'*dni4p);                            
                    end
                elseif type==5 % Tech Moms
                    ds_dtheta2_mt(:,2)=si.*(qq-si'*qq);
                    ds_dtheta2_mt(:,7)=si.*(pp-si'*pp);
                    ds_dtheta2_mt(:,12)=si.*(dni-si'*dni);
                    if TypeD==7 || TypeD==8 || TypeD==9
                        ds_dtheta2_mt(:,17)=si.*(dni1p-si'*dni1p);                            
                        ds_dtheta2_mt(:,22)=si.*(dni2p-si'*dni2p);                            
                        ds_dtheta2_mt(:,27)=si.*(dni3p-si'*dni3p);                            
                        ds_dtheta2_mt(:,32)=si.*(dni4p-si'*dni4p);                            
                    end
                elseif type==6 % College Moms
                    ds_dtheta2_mt(:,3)=si.*(qq-si'*qq);
                    ds_dtheta2_mt(:,8)=si.*(pp-si'*pp);
                    ds_dtheta2_mt(:,13)=si.*(dni-si'*dni);
                    if TypeD==7 || TypeD==8 || TypeD==9
                        ds_dtheta2_mt(:,18)=si.*(dni1p-si'*dni1p);                            
                        ds_dtheta2_mt(:,23)=si.*(dni2p-si'*dni2p);                            
                        ds_dtheta2_mt(:,28)=si.*(dni3p-si'*dni3p);                            
                        ds_dtheta2_mt(:,33)=si.*(dni4p-si'*dni4p);                            
                    end
                end
                if TypeD==2 % nonlinear distance component
                    ds_dtheta2_mt(:,nTheta2-1)=si.*((dni.^2)-si'*(dni.^2));   
                elseif TypeD==3 && market==312
                    ds_dtheta2_mt(:,nTheta2-1)=si.*(dni - si'*dni);          
                elseif TypeD==4 
                    ds_dtheta2_mt(:,nTheta2-1)=si.*(dniHi - si'*dniHi);         
                elseif TypeD==5
                    ds_dtheta2_mt(:,nTheta2-1)=si.*(dniLo - si'*dniLo);   
                elseif TypeD==6
                    ds_dtheta2_mt(:,nTheta2-2)=si.*(dniLo - si'*dniLo);   
                    ds_dtheta2_mt(:,nTheta2-1)=si.*(dniHi - si'*dniHi); 
                end                    
                ds_dtheta2_mt(:,nTheta2)=dViNodes(vik)*(si.*(qq-si'*qq));
                
                %Micro Moments Jacobian--------------------
                ds=ds_ddelta_mt*ddeltadtheta+ds_dtheta2_mt;
                
                if TypeD==8
                    dmoVi=dmoVi+ [[ds'*q]'  [ds'*p]' [ds'*dni2l]' [ds'*dni2p]']*ViWeights(vik)';
                    dmdepsVi=dmdepsVi+[[ds_ddelta_mt'*q]'  [ds_ddelta_mt'*p]'  [ds_ddelta_mt'*dni2l]' [ds_ddelta_mt'*dni2p]']*ViWeights(vik);                    
                elseif TypeD==9
                    dmoVi=dmoVi+ [[ds'*q]'  [ds'*p]' [ds'*dni]' [ds'*sni3p]']*ViWeights(vik)';
                    dmdepsVi=dmdepsVi+[[ds_ddelta_mt'*q]'  [ds_ddelta_mt'*p]'  [ds_ddelta_mt'*dni]' [ds_ddelta_mt'*sni3p]']*ViWeights(vik);                    
                else 
                    dmoVi=dmoVi+ [[ds'*q]'  [ds'*p]' [ds'*dni]']*ViWeights(vik)';
                    dmdepsVi=dmdepsVi+[[ds_ddelta_mt'*q]'  [ds_ddelta_mt'*p]'  [ds_ddelta_mt'*dni]']*ViWeights(vik);                    
                end
            end
            % Store contribution to dm_dtheta of Node ni
            dmoTypes(type,:)=dmoVi*wnode;
            
            % Store contribution to dm_deps of Node ni
            dmdepsTypes(type,:)=dmdepsVi*wnode;
            
            if TypeD==8
                mType(type,:)=ViWeights'*share_ijv'*[q p dni2l dni2p]*wnode;
            elseif TypeD==9   
                mType(type,:)=ViWeights'*share_ijv'*[q p dni sni3p]*wnode;                
            else 
                mType(type,:)=ViWeights'*share_ijv'*[q p dni]*wnode;               
            end
        end
        
        
        if sum(isnan(dmoTypes))>0
            Theta2
            error('NaNs in the moment calculations \n')
        end
        
        
        if TypeD==8 || TypeD==9
            em(s,:)=reshape(mType',24,1)';
            %mType(:,3:4) = mType(:,3:4) ./ shareDistType;
        else
            em(s,:)=reshape(mType',18,1)';
        end
        eM=eM+em;
        
        dmdtheta{it}=dmoTypes;
        dmdeps{it}=dmdepsTypes;
        it=it+1;    
    end
    
    dMdTheta{i,1}=dmdtheta; % dmoTypes below
    dMdEps{i,1}=dmdeps; % dmdepsTypes
end




    
    
    

Jtemp=zeros(nM,Model.nTheta2);
DMtemp=zeros(nM,length(Model.S));

for ni=1:length(locNodeDist)
    slist=Model.mapNodeList(ni,:);
    for y=1:length(Model.Years)
        
        indexM=(Model.mAllIndex(:,1)==WNodes(ni,1) & Model.mAllIndex(:,2)==Model.Years(y));
        pickS=(slist(y) ==counterMarket);
        
        dmoTypes=dMdTheta{ni,1}{y};
        dmdepsTypes=dMdEps{ni}{y};
        
        if TypeD==8 || TypeD==9
            nFirms=size(dmdepsTypes,2)/4;
            dMdTheta_m=NaN(4*6,nTheta2);
            dMdEps_m=NaN(4*6,nFirms);
            for tt=1:6
                dMdTheta_m(1+4*(tt-1),1:nTheta2)=dmoTypes(tt,1:nTheta2);
                dMdTheta_m(2+4*(tt-1),1:nTheta2)=dmoTypes(tt,nTheta2+1:2*nTheta2);
                dMdTheta_m(3+4*(tt-1),1:nTheta2)=dmoTypes(tt,2*nTheta2+1:3*nTheta2);
                dMdTheta_m(4+4*(tt-1),1:nTheta2)=dmoTypes(tt,3*nTheta2+1:end);
                
                dMdEps_m(1+4*(tt-1),1:nFirms)=dmdepsTypes(tt,1:nFirms);
                dMdEps_m(2+4*(tt-1),1:nFirms)=dmdepsTypes(tt,nFirms+1:2*nFirms);
                dMdEps_m(3+4*(tt-1),1:nFirms)=dmdepsTypes(tt,2*nFirms+1:3*nFirms);
                dMdEps_m(4+4*(tt-1),1:nFirms)=dmdepsTypes(tt,3*nFirms+1:4*nFirms);                
            end         
        else
            nFirms=size(dmdepsTypes,2)/3;
            dMdTheta_m=NaN(3*6,nTheta2);
            dMdEps_m=NaN(3*6,nFirms);
            for tt=1:6
                dMdTheta_m(1+3*(tt-1),1:nTheta2)=dmoTypes(tt,1:nTheta2);
                dMdTheta_m(2+3*(tt-1),1:nTheta2)=dmoTypes(tt,nTheta2+1:2*nTheta2);
                dMdTheta_m(3+3*(tt-1),1:nTheta2)=dmoTypes(tt,2*nTheta2+1:end);

                dMdEps_m(1+3*(tt-1),1:nFirms)=dmdepsTypes(tt,1:nFirms);
                dMdEps_m(2+3*(tt-1),1:nFirms)=dmdepsTypes(tt,nFirms+1:2*nFirms);
                dMdEps_m(3+3*(tt-1),1:nFirms)=dmdepsTypes(tt,2*nFirms+1:3*nFirms);
            end
        end
        Jtemp(indexM,:)=Jtemp(indexM,:)+dMdTheta_m;
        DMtemp(indexM,pickS)=DMtemp(indexM,pickS)+dMdEps_m;        
    end
end

%PI=repmat(Model.nPi,1,3);
% Model
% 
% EM=[sum(eM.*PI)./sum(PI)]';
% 
% 
% Model.MM
% 
% for i=1:size(PI,1)
%     m=unique(Model.SIndex(i==Model.counterMarket,1));
%     y=unique(Model.SIndex(i==Model.counterMarket,2));
% 
%     for type=1:6
%     for k=1:3
%     pickM=find(Model.mAllIndex(:,1)==m & Model.mAllIndex(:,2)==y & Model.mAllIndex(:,3)==type & Model.mAllIndex(:,4)==k);
%     end
%     end
%     
% end
% Model.mAllIndex(:,4)==1
% Model.mAllIndex(:,4)==2
% Model.mAllIndex(:,4)==3
% 
% 
% Model.MM(Model.mAllIndex(:,4)==1 )
% 
% 
% reshape(Model.MM,size(PI,1),size(PI,2));

eeM=reshape(eM',size(Model.MM,1),1);
g0=eeM-Model.MM-Model.AdjustmentObjMM;
g0=g0(Model.mAllIndex(:,end)==1);
G=g0'*Model.WMM*g0;


DM=DMtemp(Model.mAllIndex(:,end)==1,:);

J=Jtemp(Model.mAllIndex(:,end)==1,:);
Dg=2*J'*Model.WMM*g0;

end


