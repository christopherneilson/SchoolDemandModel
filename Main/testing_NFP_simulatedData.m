%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TITLE: Testing Code Using Simulated Data
%
% DESCRIPTION:
% This script loads (or creates) simulated data to estimate the model
% described in Neilson (2025). It defines model parameters, loads data, 
% performs sanity checks, and runs Monte Carlo simulations for evaluation.
%
% CREATED BY: Christopher Neilson
% CREATION DATE: 2021-04-13
%
% LAST MODIFIED BY: Ignacio Lepe
% LAST MODIFIED DATE: 2025-02-25
%
% VERSION: 1.0
%
% NOTES:
% - This script requires access to real data for generating simulated data.
% - Ensure all necessary toolboxes and dependencies are installed.
% - Adjust the test parameters as needed for different experimental setups.
%
% REQUIREMENTS:
% - MATLAB R202X or later
% - Global variables: Model, pathData, Theta1_True, Theta2_True
% - Dependencies: `create_fake_data`, `test_checkshareEquation`, `test_inversion`
%                 `test_clockmainfunctions`, `test_gradient_components`
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Setup

    clear
    global Model pathData pathFakeData Theta1_True Theta2_True % Define global variables used across functions
    warning('on') % Enable warnings

%% Paths

% 1. Data
    pathData = '/Users/ignaciolepe/ConsiliumBots Dropbox/Ignacio Lepe/NFP-FakeData/'; % Path to the real data
    pathFakeData = '/Users/ignaciolepe/Documents/GitHub/NFP-FakeData/Data/'; % Path to the simulated data

% 2. Codes
    % This code automatically locates the root directory
    rootPath = fileparts(fileparts(matlab.desktop.editor.getActiveFilename())); % Go one level up
    pathMain = fullfile(rootPath, 'Main'); % Path to the Main directory
    pathTools = fullfile(rootPath, 'Tools'); % Path to the Tools
    addpath(genpath(pathTools)); % Add Tools directory to MATLAB path

%% Options

    % Generating new simulated data or using last created
    Test.simFakedata=1; % if equal to 1, create_fake_data will create a new simulated database. NOTE: You need access to the full read-data to run this code.

%% Set Model Parameters

    % Define true structural parameters for the model
    Theta1_True=[
         1.5
        -0.1
        -0.3
        -0.6
        -0.3
         0.9
        -0.2
         1.0
         0.0
        -1.4
        -0.4
        -0.5
         0.5
        -0.5
        ];
    Theta2_True=[
        1.20
        1.60
        2.10
       -0.40
       -1.10
       -0.30
       -0.10
       -0.50
       -0.80
       -1.40
       -0.90
       -0.90
       -0.70
       -0.20
        1.00
    ];

%% Load Data

    if Test.simFakedata==0
       load([pathFakeData 'TestDataNFP_2025-02-25_1825.mat']) % Last output of create_fake_data. 
    else
        create_fake_data % Based on real data moments, this code creates only the necesary data to simulate the model. NOTE: You need access to the full read-data to run this code.
    end

%% Sanity Cheks

    % Choose what to test 
    Test.ShareEquation=1;
    Test.Inversion=1;
    Test.ClockFun=1;
    Test.Gradient=1;

    % Check Share Equation
    if Test.ShareEquation==1
        test_checkshareEquation
    end

    % Check Inversion
    if Test.Inversion==1
        test_inversion
    end

    % Clock Main Functions
    if Test.ClockFun==1
        test_clockmainfunctions
    end

    % Evaluating Gradients
    if Test.Gradient==1
       test_gradient_components
    end


%% Run Monte Carlo (MC) Simulation for Full Objective Function

%1. Set optimization options for different solvers
    opts1 = optimset('Display','iter','FinDiffType','central','MaxFunEvals',1000,...
              'TolFun',1E-4,'TolX',1E-4); 

    opts2 = optimset('Algorithm', 'interior-point', 'Display','iter', ...
            'GradObj','on','GradConstr','off', ...
            'MaxIter',1000, ...
            'TolX', 1e-15, 'TolFun', 1e-8, 'TolCon', 1e-8, 'UseParallel', true);

%2. Define parameter limits
    % Set lower bound for parameters
    x_L=-15*ones(Model.nTheta2,1);
    positive=[1 2 Model.nTheta2]; % Indices of parameters that must be positive (e.g., standard deviations)
    negative=[ 4 5 6 8 9 10]; % Indices of parameters that must be negative (e.g., price coefficients)
    x_L(positive)=0; % std dev must be positive

    % Set upper bound for parameters
    x_U=5*ones(Model.nTheta2,1); 
    x_U(negative)=0; % price coefficient must be negative 

    % Store bounds in the model structure
    Model.x_L=x_L;
    Model.x_U=x_U;

%3. Generate random initial parameter values
    RandTheta2=NaN(Model.nTheta2,nsim); % Initialize random parameter matrix
    noise=rand(Model.nTheta2,nsim);  % Generate uniform random noise
    
    % Scale noise to fit within parameter bounds
    for i=1:nsim
        for k=1:Model.nTheta2
            RandTheta2(k,i)=noise(k,i)*(x_U(k)-x_L(k))+x_L(k);
        end
    end
    
    % Include true parameter values as the first column
    RandTheta2=[Theta2_True RandTheta2];

%4. Initialize variables for storing results

    Model.Parallel = false; % Disable parallel computing
    timeA = NaN(nsim,1); % Store computation time for first optimization
    timeB = NaN(nsim,1); % Store computation time for second optimization
    Sol0 = NaN(size(RandTheta2,1), nsim); % Store initial solutions
    Sol1 = NaN(size(RandTheta2,1), nsim); % Store final solutions
    G0 = timeA; % Store GMM objective function values for first optimization
    G1 = timeA; % Store GMM objective function values for second optimization
    Qo = NaN(nsim,1); % Store initial GMM objective function evaluations

    % Compute initial delta values for market shares inversion
    [expdelta0, flag, norm_maxShares, norm_max,  iter]=invertmarketshares(exp(Model.deltaTrue),RandTheta2(:,1));
    Model.delta=log(expdelta0);

    
%5. Run Monte Carlo Simulations

    for i=1:size(RandTheta2,2)
        fval0 = NaN; % Placeholder for first objective value
        fval1 = NaN; % Placeholder for second objective value

        randTheta=RandTheta2(:,i); % Select current set of random parameters
        [expdelta0, flag, norm_maxShares, norm_max,  iter]=invertmarketshares(exp(Model.deltaTrue),randTheta);
        Model.delta=log(expdelta0);
        
        % Initialize weighting matrices
        Model.WBLP = eye(size(Model.IV,2)); %
        Model.WMM = eye(sum(Model.mAllIndex(:,end))); %
        
        % Evaluate initial GMM objective function
        x=randTheta;
        tic
        [Qo(i), Dgo, Q1o, Q2o, DgMo, DgIVop]=gmmObj_NFP(x);
        tic
        
        % Check if Knitro solver is available, otherwise use fmincon
        if exist('knitromatlab.m')>0
            % First optimization step
            tic
            [Sol0(:,i),G0(i),exitflag0(i)]=knitromatlab(@gmmObj_NFP,x,[],[],[],[],Model.x_L,Model.x_U,[],[],opts2);
            timeA(i)=toc/60; % Store computation time in minutes
            
            % Update weighting matrices
            [SE, Vmm, Vblp, Wblp, Wmicro]=getVarCov(Sol0(:,i), Model,students,0);
            Model.WBLP=Wblp;
            Model.WMM=Wmicro;
            
            % Second optimization step using updated weights
            tic
            [Sol1(:,i),G1(i),exitflag(i)]=knitromatlab(@gmmObj_NFP,Sol0(:,i),[],[],[],[],Model.x_L,Model.x_U,[],opts2);
            timeB(i)=toc/60;
        else
            % Use fmincon if Knitro is not available
            tic
            [Sol0(:,i),G0(i),exitflag0(i)]=fmincon(@gmmObj_NFP,x,[],[],[],[],Model.x_L,Model.x_U,[],opts2);
            timeA(i)=toc/60;
            
            % Update weighting matrices
            [SE, Vmm, Vblp, Wblp, Wmicro]=getVarCov(Sol0(:,i), Model,students,0);
            Model.WBLP=Wblp;
            Model.WMM=Wmicro;
            
            % Second optimization step using updated weights
            tic
            [Sol1(:,i),G1(i),exitflag(i)]=fmincon(@gmmObj_NFP,Sol0(:,i),[],[],[],[],Model.x_L,Model.x_U,[],opts2);
            timeB(i)=toc/60;
        end

    end

%6. Identify the best solution
    
    gtemp=G1;
    gtemp(1)=Inf; % Ignore the first solution when finding the minimum
    [Gmin0, loc0]=min(gtemp); % Identify the best parameter set
    xstar0=Sol1(:,loc0); % Extract optimal parameter estimates

    % Compute final delta values based on the best parameter estimates
    [expdelta0, flag, norm_maxShares, norm_max,  iter]=invertmarketshares(exp(Model.delta),xstar0);
    Model.delta=log(expdelta0); 
    delta0=log(expdelta0);

    % Compute OLS estimates for Theta1
    Model.theta1=inv(Model.XZ'*Model.XZ)*Model.XZ'*delta0;
    Theta1Star=Model.theta1(1:length(Theta1_True)); % Extract Theta1 estimates
    Theta2Star=xstar0; % Extract Theta2 estimates

%7. Evaluate the GMM Objective Function
    
    QSol=gmmObj_NFP(Theta2Star);
    QTrue=gmmObj_NFP(Theta2_True);

%8. Display results

    info.rnames=strvcat('Coefficients 1',strvcat(Pack.Theta1_Names(1:length(Theta1_True))));
    mprint([Theta1Star Theta1_True abs(Theta1_True-Theta1Star)./Theta1_True   ] ,info)
    norm(Theta1_True-Theta1Star)

    info.rnames=strvcat('Coefficients 2',strvcat(Pack.Theta2_Names));
    mprint([Theta2Star Theta2_True abs(Theta2_True-Theta2Star)./Theta2_True],info)
    norm(Theta2_True-Theta2Star)

%9. Compute Gradients and Model Fit

    expdelta=invertmarketshares(ones(length(Model.S),1),Theta2_True);
    Model.delta=log(expdelta);
    [~, Shares, M, G_M, dS_ddelta,dS_dTheta] =sim_model_NFP(expdelta,Theta2_True);
    [DgM2, M2, G2,J2] = gradient_MicroM_NFP(expdelta,Theta2_True,dS_dTheta,dS_ddelta);
    [QSol,  DgSol,  Q1Sol,  Q2Sol,  DgMSol,  DgIVSol ]=gmmObj_NFP(Theta2Star);
    [QTrue, DgTrue, Q1True, Q2True, DgMTrue, DgIVTrue ]=gmmObj_NFP(Theta2_True);

    % Store final estimated parameters
    Model.Theta1Star=Theta1Star;
    Model.Theta2Star=Theta2Star;


%% Evaluating Moments

%1. Compute market shares using the estimated parameters (Theta2Star)
    [expdelta, flag, norm_maxShares, norm_max, iter] = invertmarketshares(exp(Model.delta), Theta2Star);
    [S, Shares, M] = sim_model_NFP(expdelta, Theta2Star); % Simulate model with estimated parameters
    printMoments(Model, M); % Print the moments computed with estimated parameters

%2. Compute market shares using the true parameters (Theta2_True)
    [expdelta, flag, norm_maxShares, norm_max, iter] = invertmarketshares(exp(Model.delta), Theta2_True);
    [S, Shares, MTrue] = sim_model_NFP(expdelta, Theta2_True); % Simulate model with true parameters
    printMoments(Model, MTrue); % Print the moments computed with true parameters

%3. Compute Variance from Simulation

    % Compute standard errors and variance-covariance matrices based on simulated data
    [SE, Vmm, Vblp, Wblp, Wmicro] = getVarCov(Theta2_True, Model, students, 1);

%4. Print Estimated and True Parameters with Standard Errors

    % Create a structured table with parameter names
    infoT2.rnames = strvcat('Starting Point Coefficients Theta2', ...
                            strvcat(strvcat(Pack.Theta1_Names(1:length(Theta1_True))), strvcat(Pack.Theta2_Names)));

    % Define column headers for the output table
    infoT2.cnames = strvcat('True Coeff', 'Estimated Coeff', 'StdError');

    % Print a comparison table of true and estimated coefficients, including standard errors
    mprint([[Theta1_True(end-length(Theta1_True)+1:end); Theta2_True] ...  % True parameters
            [Model.Theta1Star(end-length(Theta1_True)+1:end); Model.Theta2Star] ... % Estimated parameters
            SE], ... % Standard Errors
            infoT2);

%5. Format Estimates for LaTeX Table Output

    % Extract the last five Theta1 estimates and all Theta2 estimates
    ThetaStar = [Model.Theta1Star(end-4:end); Model.Theta2Star];

    % Initialize a cell array to store formatted LaTeX table rows
    tableTestEstimates = {};

    % Loop through each estimated parameter and format it for LaTeX output
    for i = 1:length(ThetaStar)
        tableTestEstimates{i,1} = infoT2.rnames(i+1,:); % Parameter name
        tableTestEstimates{i,2} = ' & '; % Column separator
        tableTestEstimates{i,3} = num2str(ThetaStar(i), '%5.4f'); % Estimated coefficient
        tableTestEstimates{i,4} = ' & '; % Column separator
        tableTestEstimates{i,5} = num2str(SE(i), '%5.3f'); % Standard error
        tableTestEstimates{i,6} = ' \\ '; % LaTeX row ending
    end
