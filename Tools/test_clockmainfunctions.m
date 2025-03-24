











[expdelta0, flag, norm_maxShares, norm_max,  iter]=invertmarketshares(expdelta0,Theta2_True);
Model.delta=log(expdelta0);

for k=1:100
tic
S =sim_model_NFP(expdelta0,Theta2Guess);
time1(k)=toc;
end
time(1)=median(time1);

tic
[S3, Shares, M, G_M, dS_ddelta,dS_dTheta] =sim_model_NFP(expdelta0,Theta2_True);
time(2)=toc;
tic
[Dg0, M2, G2,J] = gradient_MicroM_NFP(expdelta0,Theta2_True,dS_dTheta,dS_ddelta);
time(3)=toc;
tic
[expdelta0 flag]=invertmarketshares(ones(length(Model.delta),1),Theta2_True);
time(4)=toc;

% Evaluating the Objective Function
tic
Q0=gmmObj_NFP(Theta2_True);
time(5)=toc;
tic
[Q, Dg, Qm, Qblp, DgM, DgIV]=gmmObj_NFP(Theta2_True);
time(6)=toc;

fprintf('\nWe time evaluating the main functions used in the program.\n\n')
fprintf('\nSimulating the model aggregate shares = \nMax Time = %10.8f  \n',time(1))
fprintf('\nSimulating the model micro shares, and moments = \nMax Time = %10.8f  \n',time(2))
fprintf('\nCalculating the gradient of the micromoments = \nMax Time = %10.8f  \n',time(3))
fprintf('\nTiming share inversion = \nMax Time = %10.8f  \n',time(4))
fprintf('\nEvaluating the objective function = \nMax Time = %10.8f  \n',time(5))
fprintf('\nEvaluating the gradient of the objective function = \nMax Time = %10.8f  \n',time(6))