function [expdelta2, flag, norm_maxShares, norm_max,  iter]=invertmarketshares(expdelta0,Theta2)

global Model



%% Unpack model data and parameters
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
if Model.lognormalRC==1
    ViNodes=exp(SigRC.'*Model.q_Nodes);
else
    ViNodes = SigRC.'*Model.q_Nodes;    
end
ViWeights=Model.q_Weights;

P=Model.P;
Q=Model.Q;
SEP=Model.SEP;

% Calc parallel by Market-Year
WNodes=Model.wNodes;
Qnorm=Q(Model.normIndex);
Model.D=parallel.pool.Constant(Model.distanceMatrix);
D = Model.D;
locSchoolDist=Model.locSchoolDist;
locNodeDist=Model.locNodeDist;

piW=Model.piW;
counterMarket=Model.counterMarket;
nodeMarket = Model.nodeMarket;
nS=size(Model.SIndex,1);
sList=Model.mapNodeList;
yList=Model.yList;  



%% Iterate inner loop        
        
iter=1;
expdelta0=norm_delta(expdelta0,Model.normIndex);
norm_maxSharesOld=1;
norm_max=1; norm_maxShares=1;showit=1;normaxold=10;

count=0;countM=0;countS=0;flag=0;

while norm_max > Model.tol_inner && norm_maxShares>Model.tol_inner  && iter < 1000
 
    sharesE =   sim_Slim(expdelta0);
    expdelta1 = expdelta0.*(Model.S./sharesE);
    expdelta2 = norm_delta(expdelta1,Model.normIndex); % Apply Normalization
    
    sharesE2 =sim_Slim(expdelta2);
    expdelta3 = expdelta2.*(Model.S./sharesE2);
    expdelta4 = norm_delta(expdelta3,Model.normIndex); % Apply Normalization
    
    v=log(expdelta4.*expdelta0./(expdelta2.*expdelta2));
    r=log(expdelta2./expdelta0);
    alph=(v'*r)./(v'*v);
    
    delta2=log(expdelta0)-2*alph*r+alph^2*v;
    
    [norm_max, ~]= max(abs(exp(delta2) -expdelta0 ));
    [norm_maxShares, ~]= max(abs(Model.S-sharesE2));
    
    if iter/100>=showit && iter/1001<showit
        fprintf('Dif of shares of %17.12f and diference of expdelta of %17.12f  at iteration %10.0f   \n ', [norm_maxShares norm_max  iter] )
        showit=showit+1;
    end
    
    % START ROBUSTNESS SECTION
    % If iteration gets stuck, get out after some time.
    if (norm_max==normaxold) || (norm_maxShares==norm_maxSharesOld)
        count=count+1;
        %  fprintf('Not updating it seems %4.0f \n',count)
        if count>=11
            flag=1;
            disp('jumping ship - flag 1')
            break
        end
    end
    % If iteration goes wrong way, get out after some time.
    if norm_max>normaxold
        countS=countS+1;
        if countM>=11
            flag=2;
            disp('jumping ship - flag 2')
            break
        end
    end
    
    if max(isnan(expdelta0))>=1  || max(isnan(delta2))>=1
        expdelta2 = expdelta0 ;
        flag=3;
        disp('jumping ship - flag 3')
        break
    end
    
    
    % ----------------------------------------------------------------------------
    % END ROBUSTNESS SECTION
    
    expdelta0 = exp(delta2);
    normaxold=norm_max;
    norm_maxSharesOld=norm_maxShares;
    iter = iter + 1;
end




    function S =sim_Slim(expdelta)


        S=zeros(nS,1);
        parfor i=1:length(locNodeDist)
            ni=locNodeDist(i);
            slist=sList(ni,:);
            market = nodeMarket(i);
     
            dist_ni=D.Value(:,ni);
            wn=WNodes(ni,:);
            weights=wn(:,5:end);
            
            for s=slist % Parallel work across markets
                
                pickS=(s==counterMarket);
                stemp=zeros(nS,1);
                
                delta_mt=log(expdelta(pickS));
                nFirms=length(delta_mt);
                
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
                    dni3p = (dni - 3) .* (dni>=3);
                    dni4p = (dni - 4) .* (dni>=4);
                else
                    dni1p = 0; dni2p = 0; dni3p = 0; dni4p = 0;            
                end            
        
                shares=NaN(length(delta_mt),6);
                for type=1:6
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
                    shares(:,type)=share_ijv*ViWeights;
                end
                stemp(pickS,1)=(shares.*weights)*pi';
                S=S+stemp;
            end
        end
        
        
        
        
    end




end