%% Testing Inversion Function
% Test invertmarket shares function. Starting at a random delta, return exp(deltaTrueNormed)
[expdelta0, flag, norm_maxShares, norm_max,  iter]=invertmarketshares(ones(size(Model.delta,1),1),Theta2_True);
% case with normed delta 
% simple case
adjSigma=1;
Eps=randn(length(Model.S),1)*2;
beta=[ones(size(Eps,1),1) Model.IV]\Eps;
Eps=Eps-[ones(size(Eps,1),1) Model.IV]*beta; % Make these uncorrelated to IV.
Eps=(Eps-mean(Eps))./(adjSigma*std(Eps));

XX=randn(size(Model.XX,1),length(Theta1_True));
XX=(XX-mean(XX))./(adjSigma*std(XX));
XX(Model.normID==1,:)=0;
XXX=[Model.Q XX(:,2:end)];
deltaTrue=XXX*Theta1_True+Eps;
deltaTrueNormed=deltaTrue-deltaTrue(Model.normIndex);
expdeltaTrue=exp(deltaTrueNormed);
temp=expdelta0-exp(deltaTrueNormed);
fprintf('\nWe test the inversion function by starting at a random vector of mean utilities.\n')
fprintf('\nWe calculate the max difference between true and estimated mean utility from inversion.\nMax dif = %3.2e  \n',max(max(abs(temp./expdelta0))))
