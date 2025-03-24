Theta2Guess=Theta2_True;
schools.ShareMarketMicro=Model.S;
Model.delta=ones(size(Model.S,1),1);
tic
[expdelta1, flag, norm_maxShares, norm_max,  iter]=invertmarketshares(ones(length(Model.delta),1),Theta2Guess);
toc
Model.delta=log(expdelta1);
tic
S0 =sim_model_NFP(expdelta1,Theta2Guess);
toc

tic
[S2, Shares, M, G_M, dS_ddelta,dS_dTheta] =sim_model_NFP(expdelta1,Theta2Guess);
toc
if max(abs(schools.ShareMarketMicro-S0))>10^-8   || max(abs(schools.ShareMarketMicro-S2))>10^-8
    error('Check share equations')
else
    fprintf('Share equation checks out ')
end


