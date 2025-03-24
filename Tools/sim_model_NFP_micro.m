function [S, Shares,SharesAll] =sim_model_NFP_micro(expdelta,Theta2)

global Markets Model 

% This function takes as inputs Theta2 and mean utilities and calculates
% shares, giving as output both aggregate (vector) and individual shares
% (structure)
% Other inputs are derivatives of shares with repect to theta2 and mean
% utilities.

% Source: Neilson (2013)
% Code Revised: 5/11/2019

%%
%{

%}


% Apply Normalization by Market x Year
expdelta=norm_delta(expdelta,Model.normIndex);
if (expdelta(Model.normID==1)~=1)
    fprintf('Not normed correctly')
end

BetaEdu=Theta2(Model.BetaEduO);
BetaSEP=Theta2(Model.BetaSEPO);
AlphaEdu=Theta2(Model.AlphaEduO);
AlphaSEP=Theta2(Model.AlphaSEPO);
LambdaEdu=Theta2(Model.LambdaEduO);
LambdaSEP=Theta2(Model.LambdaSEPO);

BetaO=[BetaSEP               AlphaEdu(1)+AlphaSEP  LambdaEdu(1)+LambdaSEP   ; % No HS - SEP
       0                     AlphaEdu(1)           LambdaEdu(1)             ; % No HS - Not SEP
       BetaEdu(1)+BetaSEP    AlphaEdu(2)+AlphaSEP  LambdaEdu(2)+LambdaSEP   ; % HS - SEP
       BetaEdu(1)            AlphaEdu(2)           LambdaEdu(2)             ; % HS - Not SEP
       BetaEdu(2)            AlphaEdu(3)           LambdaEdu(3)             ; % Tech Moms
       BetaEdu(3)            AlphaEdu(4)           LambdaEdu(4)             ];% College Moms 
SigRC=Theta2(Model.SigRC);
      
% Number of simulated consumers
if Model.lognormalRC==1
    ViNodes=exp(SigRC.'*Model.q_Nodes);  
else
    ViNodes = SigRC.'*Model.q_Nodes;    
end
ViWeights=Model.q_Weights;


P=Model.P;
Q=Model.Q;
SEP=Model.SEP;


% By Market, by year
S=NaN(length(expdelta),1);
Shares=NaN(length(expdelta),6);
SharesAll={};
for m=1:length(Model.Markets) % Parallel work across markets
    
    sIndex=Markets(m).sIndex;
    dIndex=Markets(m).dIndex;
    dist=Markets(m).dist;
    nNodes=Markets(m).nNodes;
    wNodes=Model.wNodes(Model.wNodes(:,1)==Model.Markets(m),5:end);
    pi=reshape(Model.wPi(Model.wPi(:,1)==Model.Markets(m),4),[6 length(Model.Years)]);
    
    for t=1:length(Model.Years)     
        % Setup Index
        index=sIndex(:,t);
        index(index==0)=[];
        p=P(index);
        q=Q(index);
        qq=q-Q(Model.normIndex(index));
        sep=SEP(index);
        d=dist(:,dIndex(:,t)==1)';
        
        expdelta_mt=expdelta(index);
        delta_mt=log(expdelta_mt);
        
        nFirms=length(index);
        for type=1:6
        wnodes=wNodes(:,type);        
            % Adjust out of pocket price for SEP 
            if (type==1 || type==3 ) 
               pp=p;
               pp(sep==1)=0;
            else
               pp=p;
            end            
            % Set up interaction between firm Xs and unobservables
            UoType=qq*BetaO(type,1)+pp*BetaO(type,2)+d*BetaO(type,3);
            Uv=ViNodes.*qq; % Random Coefficient Part of Utility
            %get maxutility to norm vector of utility to avoid overflow
            maxuij=max(delta_mt)+max(max(max(UoType)))+max(max(max(Uv)));
            shares_m=zeros(nNodes,nFirms);
            for ni=1:nNodes
                num=exp(delta_mt+UoType(:,ni)+Uv-maxuij);
                denom=sum(num,1); % Note, no 1 added.
                share_ijv=num./repmat(denom,[nFirms 1 1]);
                share_ij=share_ijv*ViWeights;
                shares_m(ni,:)=share_ij; 
            end
            SharesAll{m,t,type}=shares_m;
            Shares(index,type)= shares_m'*wnodes;
        end   
    S(index)=Shares(index,:)*pi(:,t);
    end
end

end






