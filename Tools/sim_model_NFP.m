function [S, Shares, eM, G_M, dS_ddelta,dS_dTheta] =sim_model_NFP(expdelta,Theta2)

global  Model D 

% This function takes as inputs Theta2 and mean utilities and calculates
% shares, giving as output both aggregate (vector) and individual shares
% (structure)
% Other inputs are derivatives of shares with repect to theta2 and mean
% utilities.

% Source: Neilson (2013)
% Code Revised: 02/06/2021

% Determine flow 
if nargout==1 % request only for S, Shares, empirical moments 
    skipp=0;
elseif nargout>1 && nargout<=4  % request for S, Shares, empirical moments 
    skipp=1;
elseif nargout>4 
    skipp=2;        
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

Shares_Temp=zeros(nS,6);
Shares=Shares_Temp;

dS_dTh=cell(length(locNodeDist),1);
dS_dd=cell(length(locNodeDist),1);

parfor i=1:length(locNodeDist)
    ni=locNodeDist(i);
    slist_node=sList(ni,:);
    market = nodeMarket(i);
    
    dist_ni=D.Value(:,ni); 
    wn=WNodes(ni,:);
    weights=wn(:,5:end);
    
    ds_dths=cell(length(sList),1);
    ds_dds=cell(length(sList),1);
    
    it=1;
    for s=slist_node % Parallel work across markets
        
    
        pickS=(s==counterMarket);
        stemp=zeros(nS,1);
         
        delta_mt=log(expdelta(pickS));
        nFirms=length(delta_mt);
        
        % initialize placeholders
        em=eM_temp;
        sharestemp=Shares_Temp;
        
        ds_dthetaTypes=zeros(nFirms,nTheta2);
        ds_ddeltaTypes=zeros(nFirms,nFirms);
   
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
        else
            mType=zeros(6,3);
        end
        shares=NaN(length(delta_mt),6);
        
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
            share_ij=share_ijv*ViWeights;

            if TypeD==8 % moments are dist*(dist<2) AND dist*(dist>2)
                mType(type,:)=ViWeights'*share_ijv'*[q p dni2l dni2p]*wnode;
            elseif TypeD==9 % moments are distance and sh(dist>3)
                mType(type,:)=ViWeights'*share_ijv'*[q p dni sni3p]*wnode;                
            else
                mType(type,:)=ViWeights'*share_ijv'*[q p dni]*wnode;               
            end
            
            if skipp==2
                ds_dthetaVi=zeros(nFirms,nTheta2);
                ds_ddeltaVi=zeros(nFirms,nFirms);

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
                    
                    % Store contribution to ds_dtheta of Vi vik
                    ds_dthetaVi=ds_dthetaVi+ds_dtheta2_mt*ViWeights(vik);

                    % Store contribution to ds_ddelta of Vi vik
                    ds_ddeltaVi=ds_ddeltaVi+ds_ddelta_mt*ViWeights(vik);
                end
                
                % Store contribution to ds_dtheta of Node ni
                ds_dthetaTypes=ds_dthetaTypes+ds_dthetaVi*pi(type)*wnode;
                ds_ddeltaTypes=ds_ddeltaTypes+ds_ddeltaVi*pi(type)*wnode;
            end
            
            shares(:,type)= share_ij;
        end
        
        if skipp==2
            ds_dths{it}=ds_dthetaTypes;
            ds_dds{it}=ds_ddeltaTypes;
            it=it+1;    
        end
        
        sharestemp(pickS,:)=shares.*weights;
        Shares=Shares+sharestemp;
        if TypeD==8 || TypeD==9
            em(s,:)=reshape(mType',24,1)';            
        else
            em(s,:)=reshape(mType',18,1)';
        end

        eM=eM+em;
        
        stemp(pickS,1)=(shares.*weights)*pi';
        S=S+stemp;
    end
    
    if skipp==2
        dS_dTh{i,1}=ds_dths;
        dS_dd{i,1}=ds_dds;
    end

        
end

if nargout>2
    
    if Model.AggregateMoments==1
        microMoments=reshape(eM',size(Model.MM,1),1);
        eM=collapseMoments(microMoments,Model);
        pick=(Model.mAllIndex_Ag(:,end)==1);
        WWA=Model.WMM_Ag;
        g0=eM-Model.MM_Ag-Model.AdjustmentObjMM;
        gMM=g0(pick);
        G_M=gMM'*WWA*gMM;
    else
        eM=reshape(eM',size(Model.MM,1),1);
        g0=eM-Model.MM-Model.AdjustmentObjMM;
        g0=g0(Model.mAllIndex(:,end)==1);
        G_M=g0'*Model.WMM*g0;
    end
end


if skipp==2
    dS_dTheta=zeros(length(expdelta),Model.nTheta2);
    dS_ddelta=zeros(length(expdelta),length(expdelta));
    
    for ni=1:length(locNodeDist)
        slist_node=Model.mapNodeList(ni,:);
        for y=1:length(Model.Years)
            
            pickS=(slist_node(y) ==counterMarket);
            
            dS_dTheta(pickS,:) = dS_dTheta(pickS,:) + dS_dTh{ni,1}{y} ;
            dS_ddelta(pickS,pickS) = dS_ddelta(pickS,pickS)+ dS_dd{ni}{y};
        end
    end
end





