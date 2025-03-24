%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TITLE: Main Simulation Script for Market and School Choice Model
%
% DESCRIPTION:
% This script sets up, configures, and runs a simulation of the market 
% and school choice model. It loads necessary data, defines model parameters, 
% generates instrumental variables (IVs), normalizes school characteristics, 
% and simulates microdata. The script also computes estimated parameters 
% and moments for comparison.
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
% - This script relies on simulated and real data moments for estimation.
% - Ensure all necessary toolboxes and dependencies are installed.
% - Adjust market selection, years, and test parameters as needed.
% - Simulated microdata is generated for estimation validation.
%
% REQUIREMENTS:
% - MATLAB R202X or later
% - Global variables: Model, Markets, pathData
% - Dependencies: `setup_modeldata`, `sim_model_NFP`, `gradient_BLP_NFP`, 
%                 `gradient_MicroM_NFP`, `gmmObj_NFP`, `invertmarketshares`
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Setup

    clear
    global Markets Model pathData pathFakeData Theta1_True Theta2_True % Define global variables used across functions
    warning('on') % Enable warnings

%% Options

% 1. Select dates to load/save datasets
    Pack.dateMarkets = '2021_04_13'; % Specifies the date of the real dataset to upload. Ensure the dataset was generated using v2-codes.

% 2. Select global parameters
    Years = [2007 2011]; % Defines the range of years considered in the simulation.
    listMarkets = [1:100]; % Specifies the range of markets included in the simulation.

% 3. Define model parameters
    Pack.pathData = pathData; % Assign the data path to load the real data
    Pack.pathFakeData = pathFakeData; % Assign the data path to save the simulated data
    Pack.Model = 'Model 0'; % Select the model type to use in the simulation. Options: 'Model 0', 'Model 1', 'Model 2', 'Model 1b', 'Model 2b'
    Pack.Distance = 'Drive Distance'; % Define the distance metric for computations. Options: 'Linear', 'Drive Distance', 'Drive Time'
    Pack.StartingValues = 0; % Flag to indicate whether to use starting values (0 = no)
    Pack.DistanceRoof = 10; % Maximum allowable distance for school assignment
    Pack.TypeQ = 0; % Type of quantity parameter (model-specific)
    Pack.TypeD = 1; % Type of distance parameter (model-specific)
    Pack.DistanceSplineLo = 0; % Lower bound for distance splines (if used)
    Pack.DistanceSplineHi = 0; % Upper bound for distance splines (if used)

%% Set Route and Model Parameters

% 1. Load and Initialize Model Data
    % Load and preprocess model data based on the specified years, markets, and parameters.
    [Model, Pack, nodes, schools, students, moments, pi, iv, xx, Dist] = setup_modeldata(Years, listMarkets, Pack);
    
    % Assign the distance data to the Markets variable for further use.
    Markets = Dist; 

%2. Define Model Configuration and Simulation Parameters

    % Model-Specific Settings
    % Configure model-level settings such as aggregate moments and random coefficients.
    Model.AggregateMoments = 0; % Disable aggregate moments in the model.
    Model.lognormalRC = 0; % Disable log-normal random coefficients.

    % Numerical and Optimization Parameters
    % Set tolerance levels and weight parameters for estimation.
    Model.mWiv = 0; % Weighting parameter for instrumental variables.
    Model.mWm = 1; % Weighting parameter for moments in GMM estimation.
    Model.tol_inner = 10^-14; % Tolerance for inner optimization loop.
    Model.tol_outer = 10^-8;  % Tolerance for outer optimization loop.

    % Simulation Configuration
    % Define the number of simulations (sets of random starting points).
    nsim = 2; % Number of Monte Carlo simulations.

    % Enable or disable gradient testing.
    TestGradient = 0; % Set to 1 to enable gradient testing, 0 to disable.

%3. Set Up Randomization and Quadrature Integration

    % Set Random Seed for Reproducibility
    rng(1199, 'twister'); % Sets the random seed.

    % Initialize Delta (Unobserved School Characteristics)
    Model.delta = zeros(length(Model.S), 1); % Set all initial values to zero before running the simulation.

    % Define Quadrature Parameters
    % Set the number of quadrature points for Gaussian-Hermite quadrature.
    Model.nVi = 7; % Number of points used in numerical integration.

    % Generate quadrature nodes and weights for numerical integration.
    [q_Nodes, q_Weights] = GHQuadInit(1, Model.nVi);
    Model.q_Nodes = q_Nodes; % Store quadrature nodes in the model structure.
    Model.q_Weights = q_Weights; % Store quadrature weights in the model structure.

    % Define a colormap for any visualization needs during analysis.
    C = colormap(parula(5)); % Use 'parula' colormap with 5 levels.
    
%% Parameter Setup

% 1. Define and Display Model Parameters

    % Clear previous parameter structures
    clear infoT1 infoT2 infoT3

    % Print the true coefficients for Theta1 (first set of model parameters)
    infoT1.rnames = strvcat('True Coefficients Theta1', strvcat(Pack.Theta1_Names));
    mprint([Theta1_True], infoT1)

    % Print the true coefficients for Theta2 (second set of model parameters)
    infoT2.rnames = strvcat('True Coefficients Theta2', strvcat(Pack.Theta2_Names));
    mprint([Theta2_True], infoT2)

% 2. Generate Instrumental Variables (IVs) and Model Inputs

    % Generate and normalize explanatory variables (XX)
    adjSigma=1;
    XX=randn(size(Model.XX,1),length(Theta1_True));
    XX=(XX-mean(XX))./(adjSigma*std(XX));
    XX(Model.normID==1,:)=0; % Ensure specific normalization constraints

    % Generate IVs (Instrumental Variables) given market size
    IV=randn(size(Model.XX,1),500);      % Generating IVs
    IV=(IV-mean(IV))./(adjSigma*std(IV));  % Normalize IVs

    % Define the relationship between IVs and endogenous variables (Price & Quality)
    GammaTrue=rand(size(IV,2),2); 

    % Generate unobserved school characteristics that are endogenous to P and Q
    Eps=randn(length(Model.S),1)*2; % Random error term
    beta=[ones(size(Eps,1),1) IV]\Eps; % Regression for orthogonalization
    Eps=Eps-[ones(size(Eps,1),1) IV]*beta; % Make these uncorrelated to IV.
    Eps=(Eps-mean(Eps))./(adjSigma*std(Eps)); % Normalize Eps


    % Generate endogenous variables: Quality (Q) and Price (P)
    QP=IV*GammaTrue + XX*ones(size(XX,2),1)+ Eps*[-2 2];
    QP=(QP-mean(QP))./(adjSigma*std(QP));

    % Verify meaningful correlations
    corr(Eps,IV)
    corr(Eps,QP)
    corr(QP,IV)

    % Assign computed variables to the model structure
    Model.IV=IV;
    Model.Q=QP(:,1);
    Model.P=QP(:,2);



% 3. Estimate Parameters Using OLS, IV, and GMM with Fixed Effects

    % Compute the initial delta (utility) for students
    XXX=[Model.Q XX(:,2:end)];
    deltaTrue=XXX*Theta1_True+Eps;
    
    % Ordinary Least Squares (OLS) estimation (biased)
    theta1_ols=XXX\deltaTrue;

    % Instrumental Variables (IV) estimation (unbiased)
    W=eye(size(Model.IV,2)); % Identity weighting matrix
    theta1_iv=(XXX'*IV * W * IV'*XXX)\(XXX'*IV * W * IV'*deltaTrue);

    % Normalize delta using a reference school in each market
    deltaTrueNormed=deltaTrue-deltaTrue(Model.normIndex);
    QQ=Model.Q - Model.Q(Model.normIndex); % Adjust quality variable
    Model.XX=[QQ XX(:,2:end)]; % Update model structure
    Model.WMM=diag(ones(length(Model.MM),1)*100);
    
    XZ=[Model.IV ];
    EtaHat=XZ\Model.XX;
    Xhat=XZ*EtaHat;
    XZ=[Xhat ];
    theta1_iv=(Model.XX'*Model.IV * W * Model.IV'*Model.XX)\(Model.XX'*Model.IV * W * Model.IV'*deltaTrueNormed);

    % Estimate fixed effects for market-year combinations
    FE = [];
    FE_g = [];
    g = 1;
    for i = 1:length(Model.Markets)
        for j = 1:length(Model.Years)
            pick = (Model.SIndex(:,1) == Model.Markets(i) & Model.SIndex(:,2) == Model.Years(j));
            FE = [FE pick];
            FE_g = [FE_g g * pick];
            g = g + 1;
        end
    end
    FE_g = sum(FE_g, 2); % Aggregate fixed effects

    % Compute IV regression including fixed effects
    IVFE=[Model.IV ];
    Model.XXX=[Model.XX FE(:,2:end)];
    EtaHat=IVFE\[Model.XX FE(:,2:end)];
    Xhat=IVFE*EtaHat;
    
    Model.XZ=Xhat;
    theta1all=Model.XZ\deltaTrueNormed;
    theta1all(1:length(Theta1_True))

    % Compute unbiased IV estimates with fixed effects
    theta1_iv2=(Model.XXX'*Model.IV * W * Model.IV'*Model.XXX)\(Model.XXX'*Model.IV * W * Model.IV'*deltaTrueNormed);

    % Print final unbiased estimates
    theta1_iv2(1:length(Theta1_True)); % Correctly estimated


%% Concentrating Out Fixed Effects and Computing Model Estimates

% 1. Iteratively Remove Fixed Effects (FE) from Delta and Covariates

    % Initialize variables for iteration
    tempdelta = deltaTrueNormed; % Temporary storage for the normalized delta
    deltaFE = deltaTrueNormed; % Initial delta value before demeaning
    tempXX = Model.XX; % Copy of explanatory variables (XX) for iteration
    tempIV = Model.IV; % Copy of instrumental variables (IV)
    dif = 1; count = 1; % Convergence criteria and iteration counter

    while dif > 10^-15  % Continue until convergence threshold is met
        % Compute mean XX per market-year group and demean explanatory variables
        meanXX = grpstats(tempXX, FE_g); 
        XXc = tempXX - meanXX(FE_g, :); 
        XXc(Model.normID == 1, :) = 0; % Normalize specific schools

        % Compute mean delta per market-year group and demean delta values
        meanDelta = grpstats(tempdelta, FE_g);
        deltaFE = tempdelta - meanDelta(FE_g);
        deltaFE(Model.normID == 1, :) = 0; % Normalize delta for specific schools

        % Compute the maximum absolute difference to check convergence
        dif = max(max(abs(tempXX - XXc))) + max(max(abs(tempdelta - deltaFE)));

        % Update temp variables for the next iteration
        tempdelta = deltaFE;
        tempXX = XXc; 
        count = count + 1;
    end

    % Compute the final demeaned delta values using market-year groups
    meanD = grpstats(deltaTrueNormed, FE_g);
    deltaFE = deltaTrueNormed - meanD(FE_g);

% 2. Estimate Model Parameters Using IV and Fixed Effects

    % Demean explanatory variables (XX) per market-year
    meanXX = grpstats(Model.XX, FE_g);
    XXc = Model.XX - meanXX(FE_g, :);

    % Define instrumental variables for estimation
    IVc = Model.IV; % Instrumental Variables remain unchanged

    % Prepare regression variables
    IVFE = [IVc]; % Store instrumental variables with fixed effects included
    XXXc = [XXc]; % Store demeaned explanatory variables

    % Compute projection of explanatory variables onto IVs
    EtaHat = IVFE \ XXXc; 
    Xhat = IVFE * EtaHat;

    % Store projected variables for estimation
    XZc = Xhat;

    % Compute final IV-based estimates for model parameters
    theta1all = XZc \ deltaFE; 
    theta1_iv2 = (XXXc' * IVc * W * IVc' * XXXc) \ (XXXc' * IVc * W * IVc' * deltaFE);

% 3. Compute Simulated Market Shares and Validate Model Predictions

    % Define weighting matrices for GMM estimation
    Model.WBLP = eye(size(Model.IV,2)); % Weighting matrix for IV-based estimation
    Model.WMM = eye(sum(Model.mAllIndex(:, end))); % Weighting matrix for moments

    % Compute projection matrix for XZ
    Model.IPX = eye(size(Model.XZ,1)) - Model.XZ * inv(Model.XZ' * Model.XZ) * Model.XZ';

    % Store estimated parameters and residuals
    Model.theta1all = Model.XZ \ deltaTrueNormed;
    Model.residTrue = deltaTrueNormed - Model.XXX * Model.theta1all;

    % Compute simulated market shares based on estimated parameters
    [STrue, SharesTrue, MTrue] = sim_model_NFP(exp(deltaTrueNormed), Theta2_True);

    % Plot histogram of simulated market shares
    figure(1)
    hist(STrue)
    h = findobj(gca, 'Type', 'patch');
    h.FaceColor = Pack.Colors(3, :);
    h.EdgeColor = 'w';
    box on;

    % Update model with simulated shares
    Model.SS = STrue;
    Model.S = Model.SS;
    schools.ShareMarketMicro = Model.S;
    
    % Update model moments based on simulated data
    Model.MM = MTrue;
    printMoments(Model, MTrue);

% 4. Compute Objective Function Adjustments and Gradients

    % Print estimated coefficients for Theta1 and Theta2
    infoT1.rnames = strvcat('True Coefficients Theta1', strvcat(Pack.Theta1_Names(1:length(Theta1_True))));
    mprint([Theta1_True], infoT1);

    infoT2.rnames = strvcat('True Coefficients Theta2', strvcat(Pack.Theta2_Names));
    mprint([Theta2_True], infoT2);

    % Compute adjustments for moment-based estimation
    Model.AdjustmentObjMM = Model.MM - MTrue;
    Model.AdjustmentObjIV = Model.IV' * Model.residTrue; % Should be close to zero

    % Compute final delta estimates by inverting market shares
    [expdelta, flag, norm_maxShares, norm_max, iter] = invertmarketshares(exp(deltaTrueNormed), Theta2_True);
    [S, Shares, M, G_M, dS_ddelta, dS_dTheta] = sim_model_NFP(expdelta, Theta2_True);
    Model.delta = log(expdelta);
    Model.deltaTrue = Model.delta;

    % Compute gradients for GMM estimation
    [DgIV, QIV, J] = gradient_BLP_NFP(Theta2_True);
    [S, Shares, M, G_M, dS_ddelta, dS_dTheta] = sim_model_NFP(expdelta, Theta2_True);
    [DgMM, M2, QM, J] = gradient_MicroM_NFP(expdelta, Theta2_True, dS_dTheta, dS_ddelta);
    [Q, Dg, Qm, Qiv, DgM, DgIV] = gmmObj_NFP(Theta2_True);

%% Simulate Microdata and Compute Moments

% 1. Generate Simulated Student-Level Data

    % Extract unique market identifiers
    marketList = unique(Model.SIndex(:,1));

    % Compute predicted shares using the model
    [~, ~, SharesAll] = sim_model_NFP_micro(exp(deltaTrueNormed), Theta2_True);

    % Define the number of observations per market
    nObs = 100000; 

    % Initialize student data storage
    clear simStudents;
    iObs = 1;

    % Loop through each market to generate synthetic student data
    for m = 1:length(Markets)
        market = marketList(m);
        nNodes = Markets(m).nNodes; % Number of nodes in the market
        wNodes = Model.wNodes(Model.wNodes(:,1) == Model.Markets(m), 5:end); % Node weights
        pi = reshape(Model.wPi(Model.wPi(:,1) == Model.Markets(m), 4), [6, length(Model.Years)]); % Type proportions
        sIndex = Markets(m).sIndex; % School index
        dIndex = Markets(m).dIndex; % Distance index
        dist = Markets(m).dist; % Distance matrix

        for t = 1:length(Years)   
            index = sIndex(:,t);
            index(index == 0) = []; % Remove zero entries
            p = Model.P(index); % Prices
            q = Model.Q(index); % Quality measures
            d = dist(:, dIndex(:,t) == 1)'; % Distances
            rbd = Model.SIndex(index,3); % School identifiers

            for type = 1:6  % Iterate over student types
                for ni = 1:nNodes % Iterate over nodes within the market
                    choiceProb = SharesAll{m,t,type}(ni,:)'; % Extract choice probabilities
                    nBlock = round(nObs * wNodes(ni) * pi(type, t)); % Compute the number of students to simulate

                    for ijt = 1:nBlock
                        % Simulate school choice using a multinomial draw
                        choice = mnrnd(1, choiceProb)';
                        
                        % Assign simulated student attributes
                        simStudents(iObs).Market = market;
                        simStudents(iObs).Year = Years(t);
                        simStudents(iObs).Type = type;
                        simStudents(iObs).School_RBD = rbd(choice == 1);
                        simStudents(iObs).VA2_AVE = q(choice == 1);
                        simStudents(iObs).AvePrice = p(choice == 1) * 1000; % Convert price to consistent scale
                        simStudents(iObs).Dist2School = d(choice == 1, ni);
                        iObs = iObs + 1; % Move to the next observation
                    end
                end
            end
        end
    end

    % Convert structured data to table format
    students = struct2table(simStudents);

    % Sort data by Year, Market, and School
    students = sortrows(students, {'Year', 'Market', 'School_RBD'});

% 2. Validate Simulated Data by Comparing Market Shares

    iter = 1; % Initialize iteration counter
    for i = 1:length(Model.Markets)
        market = Model.Markets(i);
        for nYear = 1:length(Model.Years)
            y = Model.Years(nYear);
            for type = 1:6
                % Extract predicted and simulated proportions for student types
                piModel=Model.wPi(Model.wPi(:,1)==market & Model.wPi(:,2)==y & Model.wPi(:,3)==type,end);
                piSim=sum(students.Type==type & students.Market==market & students.Year==y)/sum(students.Market==market & students.Year==y); 
                
                % Compute the percentage deviation from the predicted values
                difPi(iter) = (piSim - piModel) / piModel;
                iter = iter + 1;
            end
        end
    end

% 3. Compute Aggregate Moments for Validation

    % Initialize moment arrays
    MomentList = []; 
    eMM = NaN(length(Model.mAllIndex), 1); % Expected moments
    eNMM = NaN(length(Model.mAllIndex), 1); % Number of observations per moment

    for i = 1:length(Model.Markets)
        market = Model.Markets(i);
        for nYear = 1:length(Model.Years)
            y = Model.Years(nYear);
            for type = 1:6
                % Identify relevant students in the dataset
                r = find(students.Type == type & students.Market == market & students.Year == y);

                % Compute average statistics for the given type and market
                ave_q = mean(students.VA2_AVE(r));
                ave_p = mean(students.AvePrice(r) / 1000); % Convert price back to original scale
                ave_d = mean(students.Dist2School(r));

                % Store computed moments in the expected moments matrix
                eMM(Model.mAllIndex(:,1) == market & Model.mAllIndex(:,2) == y & Model.mAllIndex(:,3) == type & Model.mAllIndex(:,4) == 1) = ave_q;
                eMM(Model.mAllIndex(:,1) == market & Model.mAllIndex(:,2) == y & Model.mAllIndex(:,3) == type & Model.mAllIndex(:,4) == 2) = ave_p;
                eMM(Model.mAllIndex(:,1) == market & Model.mAllIndex(:,2) == y & Model.mAllIndex(:,3) == type & Model.mAllIndex(:,4) == 3) = ave_d;

                % Store observation counts for normalization
                eNMM(Model.mAllIndex(:,1) == market & Model.mAllIndex(:,2) == y & Model.mAllIndex(:,3) == type & Model.mAllIndex(:,4) == 1) = length(r);
                eNMM(Model.mAllIndex(:,1) == market & Model.mAllIndex(:,2) == y & Model.mAllIndex(:,3) == type & Model.mAllIndex(:,4) == 2) = length(r);
                eNMM(Model.mAllIndex(:,1) == market & Model.mAllIndex(:,2) == y & Model.mAllIndex(:,3) == type & Model.mAllIndex(:,4) == 3) = length(r);
            end
        end
    end

% 4. Compute and Print Deviations from Expected Moments

    % Compute absolute deviations from the model's predicted moments
    dif1 = abs(Model.MM - eMM);
    dif2 = abs((Model.MM - eMM) ./ Model.MM);

    % Print computed moments and deviations for validation
    printMoments(Model, eMM);
    printMoments(Model, [dif1]);
    printMoments(Model, [dif2]);

% 5. Save the Model with Updated Moments

    % Update moment weighting matrix
    Model.mAllIndex(:, end) = 1;
    Model.WMM = eye(sum(Model.mAllIndex(:, end)));

    % Store the number of observations per moment in the model structure
    Model.NMM = eNMM;

    % Save the updated model data to a MAT file
    save([Pack.pathFakeData 'TestDataNFP_' datestr(now, 'yyyy-mm-dd_HHMM') '.mat'], 'Model');
