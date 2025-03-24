function [Model,Pack,nodes,schools,students,moments,pi,iv,xx,Dist]=setup_modeldata(Years,listMarkets,Pack)
% =========================================================================
% FUNCTION: setup_modeldata
%
% DESCRIPTION:
% This function loads and processes the necessary data for the market and 
% school choice model. It filters datasets based on selected years and markets 
% and initializes relevant model variables.
%
% INPUTS:
% - Years: Vector specifying the years to include in the analysis.
% - listMarkets: Vector specifying the markets to include in the analysis.
% - Pack: Structure containing metadata such as paths, model settings, etc.
%
% OUTPUTS:
% - Model: Structure storing processed model variables.
% - Pack: Updated Pack structure with any modifications during setup.
% - nodes: Table containing information about geographic nodes.
% - schools: Table containing school-level data.
% - students: Table containing student-level data.
% - moments: Table storing aggregated market moments.
% - pi: Table with proportions of different student types per market/year.
% - iv: Table of instrumental variables (IVs).
% - xx: Feature matrix (to be populated later in the script).
% - Dist: Structure containing distance matrices for each market.
% =========================================================================

%% Load Data
    % Load preprocessed data files that contain information about nodes, schools, 
    % students, market moments, student type proportions (Pi), instrumental variables (IVs), 
    % and the list of available markets.

    load([Pack.pathData 'Data_' Pack.dateMarkets '.mat'],'Nodes','SchoolsAll','Students','MarketMoments','Pi','IVs','finalMarketList')
            
%% Validate Input Arguments

    % If no specific markets are provided, use all available markets from the dataset.
    if isempty(listMarkets) || nargin == 0
        listMarkets = unique(finalMarketList);
    end

    % If no specific years are provided, use all available years from the Schools dataset.
    if isempty(Years) || nargin == 0
        Years = unique(SchoolsAll.Year)'; 
    end

%% Select Subset of Years and Markets to Work On

    % Ensure only markets that exist in the final dataset are considered
    listMarkets=listMarkets(ismember(listMarkets,finalMarketList));


    % 1. Process Nodes Data
        % Extract relevant nodes for the selected markets
        nodes=Nodes(ismember(Nodes.Market,listMarkets),:);

        % Ensure the variable name for node identifiers is consistent
        if sum(strcmp(nodes.Properties.VariableNames,'NodeID'))==0 || sum(strcmp(nodes.Properties.VariableNames,'NodeNumber'))==1
            nodes.Properties.VariableNames{'NodeNumber'} = 'NodeID';
        end

        % Clear unused variable to free up memory
        clear Nodes

    % 2. Process Schools Data
        Schools=SchoolsAll;

        % Extract schools for the selected years
        tempSchools=Schools(ismember(Schools.Year,Years),:);

        % Further filter schools to include only those from selected markets
        schools=tempSchools(ismember(tempSchools.Market,listMarkets),:);

        % Remove schools with missing quality measure (Q2)
        dropSchool=isnan(schools.Q2);
        schools(dropSchool,:)=[];

        % Sort schools by Market, Year, and School ID (RBD) for consistency
        schools=sortrows(schools,{'Market','Year','School_RBD'}); 

    % 3. Process Students, Market Moments, and Pi Data
        % Filter student data to include only selected years and markets
        tempStudents=Students(ismember(Students.Year,Years),:);
        students=tempStudents(ismember(tempStudents.Market,listMarkets),:);

        % Filter market moments data for selected years and markets
        tempMoments=MarketMoments(ismember(MarketMoments.Year,Years),:);
        moments=tempMoments(ismember(tempMoments.Market,listMarkets),:);

        % Filter student type proportions (Pi) data for selected years and markets        
        tempPi=Pi(ismember(Pi.Year,Years),:);
        pi=tempPi(ismember(tempPi.Market,listMarkets),:);


    % 4. Process Instrumental Variables (IVs)
        tempiv=IVs(ismember(IVs.Year,Years),:);
        iv=tempiv(ismember(tempiv.Market,listMarkets),:);

        % Ensure IV data aligns with the filtered school dataset
        iv(dropSchool,:)=[];
        iv=sortrows(iv,{'Market','Year','School_RBD'}); 

        % Verify that IV and school data match correctly
        if max(abs(iv.School_RBD-schools.School_RBD))~=0
        error('IV and Schools do not match')
        end


    % 5. Assign Period Labels to IV Data
        % Define periods based on year ranges:
        % - Period 1: Before 2008
        % - Period 2: 2008-2009
        % - Period 3: 2010-2012
        % - Period 4: 2013 and beyond
        iv.Period=NaN(size(iv,1),1);
        iv.Period=(iv.Year<2008)+(iv.Year==2008 | iv.Year==2009)*2 +(iv.Year<2013 & iv.Year>2009)*3 +(iv.Year>2012)*4 ;

%% Compute Distance Matrix for Schools and Nodes

    % Load precomputed distance matrices and node information
    load([Pack.pathData 'Node_Schools_TravelTimes_Imputed.mat'], ...
        'DurationMatrix', 'DistanceMatrix', 'LineDistanceMatrix', 'Nodes', 'ListAllSchools');

    % Store Nodes information separately to avoid overwriting later
    NodesDistance=Nodes;
    clear Nodes % Free memory by removing the original variable

%1. Iterate over each selected market to compute distances
    for i=1:length(listMarkets)
        fprintf('On Market %3.0f \n',i)
        market=listMarkets(i);
        
        % Get unique school identifiers for the market
        listSchools=unique(schools.School_RBD(schools.Market==market));
        
        % Extract geographical locations (latitude and longitude) of schools
        SchoolLocations=NaN(length(listSchools),2);
        for j=1:length(listSchools)
            r=find(listSchools(j)==schools.School_RBD);
            SchoolLocations(j,1:2)=[schools.Geo_X(r(1)) schools.Geo_Y(r(1))];
        end
        
        % Extract nodes corresponding to the market
        node_m=nodes(nodes.Market==market,:) ;
        
        % Map school IDs to their positions in the ListAllSchools dataset
        sLoc=NaN(size(listSchools,1),1);
        for j=1:length(listSchools)
            sLoc(j)=find(listSchools(j)==ListAllSchools(:,1));
        end
        
        % Map node IDs to their positions in the NodesDistance dataset
        nLoc=NaN(size(node_m,1),1);
        for nNode=1:size(node_m,1)
            nLoc(nNode)=find(ismember(NodesDistance.Market,market) & ismember(NodesDistance.NodeNumber,node_m.NodeID(nNode)));
        end
        
        % Compute distance matrix based on selected distance metric
        switch Pack.Distance
            case 'Linear'
                dist_ij=LineDistanceMatrix(sLoc,nLoc)';
                dist_ij(dist_ij>Pack.DistanceRoof)=Pack.DistanceRoof;
                dist_ij(dist_ij==0)=Pack.DistanceRoof;

           case 'Drive Distance'
                dist_ij=DistanceMatrix(sLoc,nLoc)';
                dist_ij(dist_ij>Pack.DistanceRoof)=Pack.DistanceRoof;
                dist_ij(dist_ij==0)=Pack.DistanceRoof;

            case 'Drive Time'
                dist_ij=DurationMatrix(sLoc,nLoc)';
                dist_ij(dist_ij>Pack.DistanceRoof)=Pack.DistanceRoof;
                dist_ij(dist_ij==0)=Pack.DistanceRoof;
        end
        
        % Cap distances at the specified maximum distance threshold (DistanceRoof)
        dist_ij(dist_ij > Pack.DistanceRoof) = Pack.DistanceRoof;
        dist_ij(dist_ij == 0) = Pack.DistanceRoof;
    
        % Create school-to-distance index mapping for different years
        years=unique(schools.Year);
        sIndex=zeros(length(listSchools),length(years));
        dIndex=zeros(length(listSchools),length(years));
        
        for nYear=1:length(years)
            y=years(nYear);
            for j=1:length(listSchools)
                r=find(listSchools(j)==schools.School_RBD & schools.Year==y);
                if isempty(r)
                else
                     dIndex(j,nYear)=1;
                     sIndex(j,nYear)=r;
                end
            end
        end
        
        % Store computed distance matrices and indices in the Dist structure
        Dist(i).sIndex=sIndex;
        Dist(i).dIndex=dIndex;
        Dist(i).dist=dist_ij;
        Dist(i).nNodes=size(node_m,1);
        Dist(i).nSchools=length(listSchools);
    end


%2. Alternative Sparse Representation of Distance Data

    % Identify schools that exist in the full list of schools
    inlistAllSchools=ismember(ListAllSchools(:,1),schools.School_RBD);
    listAllSchools=ListAllSchools(inlistAllSchools,1);

    % Ensure the node identifiers use a consistent variable name
    NodesDistance.Properties.VariableNames{'NodeNumber'} = 'NodeID';

    % Identify nodes that belong to selected markets
    inlistAllNodeMarkets=ismember(NodesDistance(:,{'Market','NodeID'}),nodes(:,{'Market','NodeID'}));
    allnodes=NodesDistance(inlistAllNodeMarkets,:);
    
    % Map each school to its corresponding index in the distance matrix
    locDistance=NaN(size(schools,1),1);
    for i=1:size(schools,1)
        locDistance(i,1)=find(schools(i,3).School_RBD==listAllSchools(:,1));
    end
    Model.locSchoolDist=locDistance;    
    
    % Map each node to its corresponding index in the distance matrix
    locNodes=NaN(size(nodes,1),1);
    for i=1:size(locNodes,1)
        locNodes(i,1)=find(nodes.Market(i)==allnodes.Market & nodes.NodeID(i)==allnodes.NodeID);
    end
    Model.locNodeDist=locNodes;  

%3. Store the Final Distance Matrix Based on the Selected Distance Type
    switch Pack.Distance
        case 'Linear'  
            % Use the precomputed linear distance matrix
            distanceMatrix = sparse(LineDistanceMatrix(inlistAllSchools,inlistAllNodeMarkets));
            students.Dist2School=students.Dist2SchoolMeasures(:,1);
            students.Dist2School(students.Dist2School>Pack.DistanceRoof)=Pack.DistanceRoof; 
        case 'Drive Distance'
            % Use the precomputed driving distance matrix
            distanceMatrix = sparse(DistanceMatrix(inlistAllSchools,inlistAllNodeMarkets));
            students.Dist2School=students.Dist2SchoolMeasures(:,2);
            students.Dist2School(students.Dist2School>Pack.DistanceRoof)=Pack.DistanceRoof;
        case 'Drive Time'
            % Use the precomputed driving time matrix
            distanceMatrix = sparse(DurationMatrix(inlistAllSchools,inlistAllNodeMarkets));        
            students.Dist2School=students.Dist2SchoolMeasures(:,3);
            students.Dist2School(students.Dist2School>Pack.DistanceRoof)=Pack.DistanceRoof;
    end
    
    % Apply the maximum distance cap (DistanceRoof)
    distanceMatrix(distanceMatrix>Pack.DistanceRoof) = Pack.DistanceRoof;
    
    % Store the final distance matrix in the model structure
    Model.distanceMatrix=distanceMatrix;
    Model.D=parallel.pool.Constant(distanceMatrix);

%% Add Shrunk Quality Measures to Schools Data

    % Load shrunk quality dataset, which includes value-added (VA) estimates and their adjusted (shrunk) versions.
    Qshrunk=readtable([Pack.pathData 'ShrunkQuality_09_03_2021.csv']);

    % Extract and sort the relevant columns from the dataset, ensuring no duplicates
    qshrunk=unique(Qshrunk(:,{'Year','Period','School_RBD','VAshrunk_j','VA_j','ShrinkFactor_j','VAshrunk_jt','VA_Hat_jt'}));
    qshrunk=sortrows(qshrunk,{'Year','Period','School_RBD','VAshrunk_j','VA_j','ShrinkFactor_j','VAshrunk_jt','VA_Hat_jt'});

    % Initialize a temporary matrix to store matched quality values for each school
    temp=NaN(size(schools,1),6);

    % Loop through each school to match its quality data from the shrunk dataset
    for i=1:size(schools,1)
        rbd=schools.School_RBD(i); % Get the school's unique ID
        y=schools.Year(i);         % Get the corresponding year
        
        % Find the matching row in the shrunk quality dataset
        loc=find(qshrunk.School_RBD==rbd & qshrunk.Year==y);
        
        % If a match is found, extract the relevant values
        if isempty(loc)==0
            if length(loc)>1
                % If multiple matches exist, take the mean values    
                temp(i,:)=mean([qshrunk.VAshrunk_j(loc) qshrunk.VAshrunk_jt(loc) qshrunk.VA_j(loc)  qshrunk.ShrinkFactor_j(loc) qshrunk.VA_Hat_jt(loc) qshrunk.Period(loc)]);
            else
                % If a single match exists, directly assign values    
                temp(i,:)=[qshrunk.VAshrunk_j(loc) qshrunk.VAshrunk_jt(loc) qshrunk.VA_j(loc)  qshrunk.ShrinkFactor_j(loc) qshrunk.VA_Hat_jt(loc) qshrunk.Period(loc)];
            end
        end
    end
    
    % Assign extracted shrunk quality measures to the schools dataset
    schools.Qshrunk_jt = temp(:, 2);   % Shrunk VA for school-year type
    schools.Qshrunk_p  = temp(:, 1);   % Shrunk VA for period
    schools.Qjt        = temp(:, 3);   % Raw VA measure
    schools.Qshrinkfactor = temp(:, 4); % Shrinking factor applied to VA
    schools.Q_Hat_jt   = temp(:, 5);   % Predicted VA estimate
    schools.Period     = temp(:, 6);   % Time period classification

%% Store Key Model Variables  

    % Assign school market shares to the model
    Model.S=schools.ShareMarketMicro;
    
%1. Normalize school quality (Q) based on different methods specified by Pack.TypeQ
    if Pack.TypeQ==1
        % Normalize shrunk value-added measure (Qshrunk_jt) using mean and std from 2007
        muVA=nanmean(Schools.Q2(Schools.Year==2007));
        sigma2=nanstd(Schools.Q2(Schools.Year==2007));
        Model.Q=(schools.Qshrunk_jt-muVA)/sigma2; %VA2_AVE (Standardized value-added estimate)
    elseif Pack.TypeQ==2
        % Similar to TypeQ=1 but imputes missing values using Q_Hat_jt or Q2
        muVA=nanmean(Schools.Q2(Schools.Year==2007));
        sigma2=nanstd(Schools.Q2(Schools.Year==2007));
        temp = schools.Qshrunk_p;
        temp(isnan(temp)) = schools.Q_Hat_jt(isnan(temp)); % Fill missing with predicted VA
        temp(isnan(temp)) = schools.Q2(isnan(temp)); % Fill remaining missing values with raw VA
        Model.Q=(temp-muVA)/sigma2; %VA2_AVE;
    elseif Pack.TypeQ==3
        % Uses an alternative quality measure Q3, replacing missing values with Q2
        muVA=nanmean(Schools.Q3(Schools.Year==2007));
        sigma2=nanstd(Schools.Q3(Schools.Year==2007));
        temp=schools.Q3;
        temp(isnan(temp))=schools.Q2(isnan(temp));  % Fill missing values with Q2
        Model.Q=(temp-muVA)/sigma2; 
    elseif Pack.TypeQ==4 
        % Uses AVE (average scores) instead of VA for quality normalization
        muVA=nanmean(Schools.AVE(Schools.Year==2007));
        sigma2=nanstd(Schools.AVE(Schools.Year==2007));
        temp=schools.AVE;
        
        % Impute missing values by computing the mean for each market-year-type combination
        for mm = listMarkets 
           for yy = Years 
               for dep = 1:3 % Loop over different dependency categories (e.g., education levels)
                   idx = find(schools.Market==mm & schools.Year==yy & schools.DEP==dep);
                   if size(idx,1)<=1
                      mAVE = nanmean(schools.AVE(schools.Market==mm & schools.Year==yy));  
                   else 
                      mAVE = nanmean(schools.AVE(idx));                     
                   end

                   temp(isnan(temp) & schools.Market==mm & schools.Year==yy & schools.DEP==dep) = mAVE;
               end
           end
        end
        Model.Q=(temp-muVA)/sigma2;    
        
    elseif Pack.TypeQ==0
        % Standardizes Q2 directly using mean and standard deviation from 2007
        muVA=nanmean(Schools.Q2(Schools.Year==2007));
        sigma2=nanstd(Schools.Q2(Schools.Year==2007));
        Model.Q=(schools.Q2 - muVA)/sigma2; %VA2_AVE

    elseif Pack.TypeQ==992
        % Directly assigns Q2 as quality measure without standardization
        Model.Q=schools.Q2; %VA2_AVE
        
    elseif Pack.TypeQ==993
        % Uses Q3 but fills missing values with Q2
        temp=schools.Q3;
        temp(isnan(temp))=schools.Q2(isnan(temp));
        Model.Q=temp; %VA2_AVE    
    else
        % Default to using Q2 as the quality measure
        Model.Q=schools.Q2; %VA2_AVE    
    end
    
    % Check for missing quality values and raise an error if any exist
    if sum(isnan(Model.Q))>0
       error('Missing Quality') 
    end

    % Assign socioeconomic priority (SEP) index from schools to the model
    Model.SEP=schools.SEP;

%2. Assign School Prices to the Model Based on Selected Model Type
    switch Pack.Model
        case {'Model 0','Model 1', 'Model 2'}   
        % Use average price per school, converting from local currency (assuming in thousands)    
        PriceMoments= moments.AVEP;
        Model.P=schools.AvePrice/1000;
        Model.P(isnan(schools.AvePrice))=0;
        
        case {'Model 1b', 'Model 2b'}
        % Uses trimmed average prices instead of raw prices
        PriceMoments= moments.AVEPt;
        Model.P=schools.AvePriceTrim/1000;
        Model.P(isnan(schools.AvePrice))=0;
        
        otherwise
        % Default: Use standard average price    
        PriceMoments= moments.AVEP;
        Model.P=schools.AvePrice/1000;
        Model.P(isnan(schools.AvePrice))=0; 
    end

%% Make Index and Fixed Effects (FE)
    % This section defines various fixed effects (FEs) and indexes used in the model.
    
    % Extract unique market IDs and years from the schools dataset
    Model.Markets=unique(schools.Market);
    Model.Years=unique(schools.Year);
    
    % Create a composite index combining Market, Year, and School RBD (unique school identifier)
    Model.SIndex=[schools.Market schools.Year schools.School_RBD];
    
    % Create an index for firms (schools) by assigning unique numerical values to school RBDs
    [RBD,~,indx_j]=unique(Model.SIndex(:,3));
    Model.FirmIndex=indx_j;
    Model.FirmID=RBD;       

%1. Create Market-Year Fixed Effects    
    Rfe=[];
    Cfe=[];
    MarketID=unique(Model.SIndex(:,1)); % List of unique markets
    TimeID=unique(Model.SIndex(:,2)); % List of unique years
    counter=1;xMarketTimeNames=[];
    
    % Loop through all market-year combinations to construct fixed effects
    for j=1:length(MarketID)
        for t=1:length(TimeID)
        idm=MarketID(j);
        idt=TimeID(t);
        xMarketTimeNames=strvcat(strvcat(xMarketTimeNames),['Market ' num2str(idm) ' - Year ' num2str(idt)]);
        rj=find(idm==Model.SIndex(:,1) & idt==Model.SIndex(:,2));
        Rfe=[Rfe;rj];
        Cfe=[Cfe; counter*ones(length(rj),1)];

        counter=counter+1;
        end
    end
    
    % Create a sparse matrix for Market-Year Fixed Effects
    MarketTimeFE=sparse(Rfe,Cfe,ones(length(Rfe),1));
    Model.MarketTimeFE=full(MarketTimeFE);
    
%2. Create Time Fixed Effects
    Rfe=[];
    Cfe=[];
    TimeID=unique(Model.SIndex(:,2));
    counter=1;
    T=NaN(length(Model.SIndex(:,2)),1);
    
    % Loop through all unique years to construct time fixed effects    
    for t=1:length(TimeID)
        idt=TimeID(t);
        rj=find(idt==Model.SIndex(:,2));
        T(rj)=t;
        Rfe=[Rfe;rj];
        Cfe=[Cfe; counter*ones(length(rj),1)];
        counter=counter+1;
    end
    
    % Create a sparse matrix for Time Fixed Effects
    TimeFE=sparse(Rfe,Cfe,ones(length(Rfe),1));

%3. Construct Municipality (Comuna) Fixed Effects    
    Dcomuna=[];
    ComunaName={};
    killcomuna=unique(schools.Geo_ComunaCode(schools.normID==1)); % Identify excluded municipalities
    listComuna=unique(schools.Geo_ComunaCode);
    dropcomuna=ismember(listComuna,killcomuna);
    listComuna(dropcomuna)=[];
    
    % Generate a fixed effect matrix for municipalities (comunas)
    for y=listComuna(1:end)'
       Dcomuna=[Dcomuna schools.Geo_ComunaCode==y]; 
       ComunaName=strvcat(ComunaName,['Comuna ' num2str(y)]); 
    end
    Model.Dcomuna=Dcomuna;

%4. Construct School Fixed Effects    
    killschool=unique(schools.School_RBD(schools.normID==1)); % Identify excluded schools
    listSchool=unique(schools.School_RBD);
    dropschool=ismember(listSchool,killschool);
    listSchool(dropschool)=[];
    FirmFE=zeros(size(schools,1),length(listSchool)); % School Fixed Effect Matrix
    FirmFE_Pre=FirmFE; % Before 2008
    FirmFE_Post=FirmFE; % After 2008
    
    % Loop through schools and assign fixed effects
    for i=1:length(listSchool)
        firm=listSchool(i);
        FirmFE(Model.SIndex(:,3)==firm,i)=1;
        FirmFE_Pre(Model.SIndex(:,3)==firm & Model.SIndex(:,2)<2008,i)=1;  
        FirmFE_Post(Model.SIndex(:,3)==firm  & Model.SIndex(:,2)>=2008,i)=1;  
    end
    Model.FirmFE=FirmFE;

%5. Construct Normalization Index and Fixed Effect    
    Model.normID=schools.normID;it=1;
    Model.normIndex=zeros(size(schools,1),1);
    normFE=zeros(size(schools,1),length(Model.Markets)*length(Model.Years));
    
    % Create a normalization index for each market-year combination
    for i=1:length(Model.Markets)
        market=Model.Markets(i);
        for nYear=1:length(Model.Years)
            y=Model.Years(nYear);

            pickSchools=find(schools.Market==market & schools.Year==y);
            r=find(schools.normID==1 & schools.Market==market & schools.Year==y);
            if length(r)~=1
                error('Norm Wrong')
            end
            normFE(r,it)=1;it=it+1;
            Model.normIndex(pickSchools)=r;
        end
    end
    Model.normFE=normFE;
    Model.wNodes=table2array(nodes);

%6. Construct Market-Probability Matrix
    tempPi=table2array(pi);
    finalPi=[];
    
    % Compute transition probabilities for each market-year-type
    for i=1:length(Model.Markets)
        market=Model.Markets(i);
        for nYear=1:length(Model.Years)
            y=Model.Years(nYear);
            for type=1:6
                r=find(tempPi(:,1)==market & tempPi(:,2)==y & tempPi(:,3)==type);
                if  isempty(r)==0
                    pi_ij=tempPi(r,4);
                else
                    pi_ij=0;
                end
                finalPi=[finalPi;  Model.Markets(i) y  type pi_ij];
            end
        r=find(finalPi(:,1)==market & finalPi(:,2)==y);
        finalPi(r,4)=finalPi(r,4)./sum(finalPi(r,4)); % Normalize probabilities        
        end

    end
    Model.wPi=finalPi;

%7. Specify Fixed Effects for Different Models
    Model.useObs=true(size(Model.S,1),1); 
    % IV Regression Elements
    switch Pack.Model
        case {'Model 0'}
        % Admin FE for public school districts (comunas)
        % X : Voucher and Private school characteristics: Religion, Traditional, High School 
        spec_Model0

        case {'Model 1','Model 1b'}
        % Brand FE for elite non private schools
        % Admin FE for public school districts (comunas)
        spec_Model1                

        case {'Model 2','Model 2b'}
        % Model with school FE and quality / price changing over time     
        spec_Model2
        %spec_Model2

        case 'Model 3'    
        %spec_Model1
         spec_Model3

        case 'Model 4'    
        %spec_Model1
         spec_Model4    
    end
    Model.nSchool=schools.nSchool;

%8. Update Student Data with School Quality and Price Imputations    

    % Initialize columns in the students table to store quality (Q2) and price (P2) values
    students.Q2=NaN(size(students.Market,1),1); % School quality measure
    students.P2=NaN(size(students.Market,1),1); % School price measure
    
    % Assign school quality and price data to students based on their enrolled school
    for i=1:length(Model.Markets)
        market=Model.Markets(i);
        for nYear=1:length(Model.Years)
            y=Model.Years(nYear);    
            
            % Identify schools operating in the current market and year
            pickSchools=(schools.Market(:,1)==market & schools.Year==y);
            listRBD=unique(schools.School_RBD(pickSchools));
            
            % Identify students in the same market and year
            pickStudents=(students.Market(:,1)==market & students.Year==y);
            TotalStudents=sum(pickStudents);     
            
            % Assign school quality and price data to each student
            for j=1:length(listRBD)
                pickStudents_j=(students.Market(:,1)==market & students.Year==y & students.School_RBD==listRBD(j));
                pickSchool_j=(pickSchools==1 & schools.School_RBD==listRBD(j));
                
                % Validate that the school's market share is correctly assigned
                if abs(sum(pickStudents_j)/TotalStudents-schools.ShareMarketMicro(pickSchool_j))>10^8
                    error('MicroShares Wrong')
                end 
                
                % Find the corresponding school index in the schools dataset
                locQ=find(schools.School_RBD==listRBD(j) & schools.Market(:,1)==market & schools.Year==y);    
                
                % Assign the school's quality and price to the students enrolled there
                students.Q2(pickStudents_j)=Model.Q(locQ);
                students.P2(pickStudents_j)=Model.P(locQ);
            end
        end
    end
    
%9. Compute Student Distance Weights (DW)
    
    % Initialize distance weights (DW) for students
    students.DW=zeros(size(students.Market,1),1);
    
    % Extract node share data (weights for different student types per node)
    nShares=[nodes.w1 nodes.w2 nodes.w3 nodes.w4 nodes.w5 nodes.w6];
    
    
    % Compute weighted distances for each student based on their node location
    for i=1:length(Model.Markets) % Loop through all markets
        market=Model.Markets(i);
        for nYear=1:length(Model.Years)  % Loop through all years
            y=Model.Years(nYear);
            
            % Extract node information for the current market
            mNodeID=nodes.NodeID(nodes.Market==market,:);
            mNodeShares=nShares(nodes.Market==market,:);
            
            % Compute distance weights separately for each student type (1 to 6)
            for type=1:6
                pickStudents=(students.Market(:,1)==market & students.Year==y & students.Type==type);
                TotalStudentsType=sum(pickStudents);

                if sum(TotalStudentsType)>0 % Only proceed if there are students of this type
                        
                    % Extract distances for selected students
                    dist2school=students.Dist2School(pickStudents);
                    hasDist=~isnan(dist2school);   % Identify students with valid distance data
                    
                    % Get student node locations
                    sNodes=students.NodeID(pickStudents==1,:);
                    
                    % Tabulate occurrences of each node
                    tab=tabulate(sNodes);
                    
                    % Identify relevant nodes in the current market
                    rep=ismember(mNodeID,tab(:,1));
                    Rep=sum(mNodeShares(rep,type)); % Compute total node share for the student type
                    
                    % Initialize weights
                    weights=zeros(length(sNodes),1);                
                    empirical=zeros(length(mNodeID),1);
                    
                    % Compute weights for each node
                    for ni=1:length(mNodeID)
                        
                        pickStudentNode=(sNodes==mNodeID(ni));

                        % Calculate the empirical share of students at this node
                        empirical(ni)=sum(pickStudentNode)/sum(hasDist);
                        
                        % Compute the weight as the node share relative to the empirical share
                        weights(pickStudentNode,1)=(mNodeShares(ni,type)/Rep)/empirical(ni);

                    end
                    
                    % Assign computed weights to students
                    students.DW(pickStudents)=weights;

                end
            end
        end
    end


%% Get Index for Parallel Market Calculation
    
    % Initialize counters and index variables for market and school data
    counter=1;
    it_s=1; % Index for schools
    it_w=1; % Index for nodes
    
    % Initialize matrices to store indexing information
    indexMarket=NaN(size(Model.SIndex,1),1); % Index for markets
    counterMarket=NaN(size(Model.SIndex,1),1); % Counter for market iterations
    cutNodes=[]; % Stores node index ranges for each market-year
    piW=[]; % Stores adjusted probability distributions for markets
    
%1. Loop through all markets to assign indices
    for m=1:length(Model.Markets)
        
       % Get the list of nodes associated with the current market
       pickNodes=find(Model.wNodes(:,1)==Model.Markets(m));
    
    % Iterate over all available years
    for t=1:length(Model.Years)
        
       % Select schools belonging to the current market and year 
       pickSchools=find(Model.SIndex(:,1)==Model.Markets(m) & Model.Years(t)==Model.SIndex(:,2));
       indexMarket(it_s:it_s+length(pickSchools)-1,1)=pickSchools;
       counterMarket(it_s:it_s+length(pickSchools)-1,1)=counter; % Store market counter 
          
       % Store the node index range for the current market-year combination
       cutNodes=[cutNodes; [it_w, it_w+length(pickNodes)-1]];
        
       % Retrieve weight probabilities (pi) for the current market-year
       pickPi=find(Model.wPi(:,1)==Model.Markets(m) & Model.Years(t)==Model.wPi(:,2));
       piW=[piW; Model.wPi(pickPi,4)'];
        
       % Increment counters
       it_s=it_s+length(pickSchools); 
       counter=counter+1;
    end
    it_w=it_w+length(pickNodes);  
    end
    
    Model.SIndex(counterMarket==2,:)
    
    % Normalize the weight probability matrix piW to ensure it sums to 1
    for s=1:max(counterMarket)

        pis= piW(s,:);
        pis(pis==0)=10^-4; % Avoid zero probabilities
        pis=round(pis,6)./sum(round(pis,6)); % Normalize probabilities
        
        % Ensure probabilities sum exactly to 1
        if sum(pis)<1
        pis=add2one(pis);
        end
        piW(s,:)=pis;
        
        % Update weight probabilities in the Model structure
        marketLine=unique(Model.SIndex(counterMarket==s,1:2))';
        pickPi=find(Model.wPi(:,1)==marketLine(1) & marketLine(2)==Model.wPi(:,2));
        Model.wPi(pickPi,4)=piW(s,:)';
        
        % Normalize node weights within the market
        w=Model.wNodes(cutNodes(s,1):cutNodes(s,2),5:end);
        w(w==0)=10^-4;
        w=round(w,6)./sum(round(w,6),1);

        for k=1:6
        if sum(w(:,k))<1   
            w(:,k)=add2one(w(:,k));
            %sum(w(:,k))<1;
        end
        end
        Model.wNodes(cutNodes(s,1):cutNodes(s,2),5:end)=w;
    end
    
    % Store final calculated values in the model structure
    Model.piW=piW;
    Model.counterMarket=counterMarket;
    Model.nodeMarket = nodes.Market;
    Model.cutNodes=cutNodes;
    Model.GeoSchools=[schools.Geo_X schools.Geo_Y];

%2. Map Nodes to Corresponding Market-Year Entries    
    % Initialize a matrix to store mapping information
    sList=NaN(length(Model.locNodeDist),length(Years));
    
    % Loop through all nodes and assign them to market-year index
    for i=1:length(Model.locNodeDist)
        ni=Model.locNodeDist(i);
        listi=find(cutNodes(:,1)<=ni & cutNodes(:,2)>=ni);
        if ~isempty(listi)
        sList(i,:)=listi';
        end
    end
    Model.mapNodeList=sList;

%3. Extract Market and Year Indices for Each Counter Market Entry
    
    % Initialize matrices to store market and year indices
    mList=NaN(max(Model.counterMarket),1);
    yList=NaN(max(Model.counterMarket),1);
    
    % Loop through all market counters to extract corresponding market-year pairs
    for i=1:max(Model.counterMarket)
        yList(i)=unique(Model.SIndex(i==Model.counterMarket,2));

        mList(i)=unique(Model.SIndex(i==Model.counterMarket,1));

    end
    Model.yList=yList;
    Model.mList=mList;
    
%4. Compute the Number of Students for Each Type in Each Market-Year    
    
    % Initialize matrix to store student counts by type
    nPi=NaN(max(Model.counterMarket),6);
    
    % Loop through each market-year and count the number of students in each type
    for i=1:max(Model.counterMarket)
        pickS=find(i==Model.counterMarket);

        npi(1)=sum(schools.nType_1(pickS));
        npi(2)=sum(schools.nType_2(pickS));
        npi(3)=sum(schools.nType_3(pickS));
        npi(4)=sum(schools.nType_4(pickS));
        npi(5)=sum(schools.nType_5(pickS));
        npi(6)=sum(schools.nType_6(pickS));

        nPi(i,:)=npi;
    end
    Model.nPi=nPi;        

%5. Construct List of Market Moments
    
    % Initialize an empty matrix to store market moments
    MomentList=[];
    
    % Iterate through each market in the dataset
    for i=1:length(Model.Markets)
        market=Model.Markets(i);
        for nYear=1:length(Model.Years)
            y=Model.Years(nYear);

            % Iterate through all student types
            for type=1:6
                r=find(moments.Market(:,1)==market & moments.Year==y & moments.Type==type);
                if isempty(r)==1
                    MomentList=[MomentList ; market y type 1 0 NaN NaN];
                    MomentList=[MomentList ; market y type 2 0 NaN NaN];
                    MomentList=[MomentList ; market y type 3 0 NaN NaN];
                    if Pack.TypeD==8 || Pack.TypeD==9
                        MomentList=[MomentList ; market y type 4 0 NaN NaN];                   
                    end
                else
                    pickStudents=(students.Market(:,1)==market & students.Year==y & students.Type==type);
                    if moments.nObsVA2(r)>30 && sum(pickStudents)>0
                        MomentList=[MomentList ; market y type 1 1 nanmean(students.Q2(pickStudents)) sum(pickStudents)];
                    else
                        MomentList=[MomentList ; market y type 1 0 NaN NaN];
                    end
                    if moments.nObsPrice(r)>30 && sum(pickStudents)>0

                        MomentList=[MomentList ; market y type 2 1 nanmean(students.P2(pickStudents)) sum(pickStudents)];
                    else
                        MomentList=[MomentList ; market y type 2 0 NaN NaN];
                    end
                    if moments.nObsDist(r)>30 && sum(students.DW(pickStudents))>0
                        if Pack.TypeD==8
                            dist = students.Dist2School(pickStudents);
                            Dmean_g2=nanmean(dist(dist>=2));
                            Dmean_l2=nanmean(dist(dist<2));
                            MomentList=[MomentList ; market y type 3 1 Dmean_l2 moments.nObsDist(r)];
                            MomentList=[MomentList ; market y type 4 1 Dmean_g2 moments.nObsDist(r)];
                        elseif Pack.TypeD==9
                            dist = students.Dist2School(pickStudents);
                            Dmean=nanmean(dist);
                            DShmean=nanmean(dist>=3);
                            MomentList=[MomentList ; market y type 3 1 Dmean moments.nObsDist(r)];
                            MomentList=[MomentList ; market y type 4 1 DShmean moments.nObsDist(r)];            
                        else 
                            Dmean=nanmean(students.Dist2School(pickStudents));
                            MomentList=[MomentList ; market y type 3 1 Dmean moments.nObsDist(r)];
                        end
                    else
                        MomentList=[MomentList ; market y type 3 0 NaN NaN];
                        if Pack.TypeD==8 || Pack.TypeD==9
                            MomentList=[MomentList ; market y type 4 0 NaN NaN];
                        end
                    end
                end
            end
        end
    end
    
    % Store computed moment data in the model structure
    Model.NMM=MomentList(:,end);
    Model.MM=MomentList(:,end-1);
    Model.mAllIndex=MomentList(:,1:end-2);
    
    % Compute weighting matrix WMM using identity matrix (diagonal weights)
    V=diag(ones(sum(Model.mAllIndex(:,end)==1),1));
    WMM=inv(V);
    Model.WMM=WMM;

 %6. Aggregate Market Moments by Period
    
    % Determine the number of moment types to aggregate
    numMomentType = 3;
    if Pack.TypeD==8 || Pack.TypeD==9
        numMomentType = 4;
    end
    
    % Initialize aggregated moment storage
    NMM=Model.NMM;
    NMM(isnan(NMM) | Model.mAllIndex(:,end)==0)=0;
    NMMTOT=[];
    mAllIndex_Ag=[];
    MM_Ag=[];
    
    % Compute aggregate moments for periods before 2008
    for m=1:numMomentType
        for t=1:6
            period=1; % Pre-2008 period
            
            % Select observations before 2008
            nObs=NMM(Model.mAllIndex(:,end)==1 & Model.mAllIndex(:,3)==t & Model.mAllIndex(:,4)==m & Model.mAllIndex(:,2)<2008);
            mM=Model.MM(Model.mAllIndex(:,end)==1 & Model.mAllIndex(:,3)==t & Model.mAllIndex(:,4)==m & Model.mAllIndex(:,2)<2008);
            
            % Compute weighted means for moments
            MM_Ag=[MM_Ag; sum(mM.*nObs)/sum(nObs)];
            mAllIndex_Ag=[mAllIndex_Ag ;[0 period t m sum(nObs)>0]];    
            NMMTOT=[NMMTOT ;sum(nObs)];        
        end
    end
    
    % Compute aggregate moments for periods after 2008
    for m=1:numMomentType
        for t=1:6
            period=2; % Post-2008 period
            
            % Select observations after 2008
            nObs=NMM(Model.mAllIndex(:,end)==1 & Model.mAllIndex(:,3)==t & Model.mAllIndex(:,4)==m & Model.mAllIndex(:,2)>2008);
            mM=Model.MM(Model.mAllIndex(:,end)==1 & Model.mAllIndex(:,3)==t & Model.mAllIndex(:,4)==m & Model.mAllIndex(:,2)>2008);
            
             % Compute weighted means for moments
            MM_Ag=[MM_Ag; sum(mM.*nObs)/sum(nObs)];
            mAllIndex_Ag=[mAllIndex_Ag ;[0 period t m sum(nObs)>0]];   
            NMMTOT=[NMMTOT ;sum(nObs)];         
        end
    end
    
    % Store aggregated moment data in the model structure
    Model.MM_Ag=MM_Ag;
    Model.NMM_Ag=NMMTOT;
    Model.mAllIndex_Ag=mAllIndex_Ag;
    Model.WMM_Ag=diag(Model.NMM_Ag)/sum(Model.NMM_Ag);

%7. Compute Market-Level Moment Proportions
    
    % Initialize market-level probability matrix
    piMarkets=zeros(max(Model.counterMarket),6*3);
    
    % Compute moment probabilities for each market-year combination
    for i=1:max(Model.counterMarket)

        year=yList(i);
        market=mList(i);
        
        % Assign period based on the year
        if year<2008
            period=1;
        elseif year>2008
            period=2;
        else
            period=0;
        end
        
        % Iterate over moment types and student types
        temp=[];
        for mm=1:3
        for type=1:6
        pick=(Model.mAllIndex(:,1)==market & Model.mAllIndex(:,2)==year & Model.mAllIndex(:,3)==type & Model.mAllIndex(:,4)==mm);
        pickAll=(Model.mAllIndex_Ag(:,2)==period & Model.mAllIndex_Ag(:,3)==type & Model.mAllIndex_Ag(:,4)==mm);

        num=nansum(Model.NMM(pick));
        dem=nansum(Model.NMM_Ag(pickAll));
        
        % Compute moment probability
        if isnan(num) | isnan(dem) | dem==0
            pitemp=0;
        else
            pitemp=nansum(Model.NMM(pick))./nansum(Model.NMM_Ag(pickAll));
        end
        temp=[temp pitemp];
        end

        end
        piMarkets(i,:)=temp;

    end
    
    % Store market-level moment probabilities
    Model.piTypes_Ag=piMarkets;
    Model.WBLP = inv(Model.IV'*Model.IV)/size(Model.IV,1);


%% Adjust Market Shares to add to one exactly
    for i=1:length(Model.Markets)
        market=Model.Markets(i);
        for nYear=1:length(Model.Years)
            y=Model.Years(nYear);
            index=find(schools.Year==y & schools.Market==market);
            s=schools.ShareMarketMicro(index);
            sum(s);
            s1=s./sum(s);
            if sum(s1)<1   
            s1=add2one(s);
            end
            schools.ShareMarketMicro(index)=s1;
        end
    end
    Model.S=schools.ShareMarketMicro;

    % Adjust Objective Function 
    Model.AdjustmentObjMM=0;
    Model.AdjustmentObjIV=0;

%% Set up Names and starting values for optimization routines.

    Model.nMarkets=length(Model.Markets);
    Model.nTypes=6;
    Model.nTime= length(Years);
    Model.nTheta1= size(Model.XX,2);
    Model.nVi=5;

    Xij_Names    ={
        'Quality x High School Mother ',...
        'Quality x 2 year Technical Degree Mother ',...
        'Quality x 4 year College Degree Mother',...
        'Quality x Poor',...
        'Price x Non High School Mother ',...
        'Price x High School Mother ',...
        'Price x 2 year Technical Degree Mother ',...
        'Price x 4 year College Degree Mother',...
        'Price x Poor'};

    if Pack.TypeD==1
        Model.nTheta2= 15;
        Dist_Names={
            'Distance x Non High School Mother',...
            'Distance x High School Mother ',...
            'Distance x 2 year Technical Degree Mother ',...
            'Distance x 4 year College Degree Mother ',...
            'Distance x Poor',...
            };
    elseif Pack.TypeD==2
        Model.nTheta2= 16;
        Dist_Names={
            'Distance x Non High School Mother',...
            'Distance x High School Mother ',...
            'Distance x 2 year Technical Degree Mother ',...
            'Distance x 4 year College Degree Mother ',...
            'Distance x Poor',...
            'Distance Squared'
            };    
    elseif Pack.TypeD==3
        Model.nTheta2= 16;
        Dist_Names={
            'Distance x Non High School Mother',...
            'Distance x High School Mother ',...
            'Distance x 2 year Technical Degree Mother ',...
            'Distance x 4 year College Degree Mother ',...
            'Distance x Poor',...
            'Distance - Santiago'
            };    
    elseif Pack.TypeD==4
        Model.nTheta2= 16;
        Dist_Names={
            'Distance x Non High School Mother',...
            'Distance x High School Mother ',...
            'Distance x 2 year Technical Degree Mother ',...
            'Distance x 4 year College Degree Mother ',...
            'Distance x Poor',...
            strcat('Distance > ',Pack.DistanceSplineHi),...
            };    
    elseif Pack.TypeD==5
        Model.nTheta2= 16;
        Dist_Names={
            'Distance x Non High School Mother',...
            'Distance x High School Mother ',...
            'Distance x 2 year Technical Degree Mother ',...
            'Distance x 4 year College Degree Mother ',...
            'Distance x Poor',...
            strcat('Distance < ',Pack.DistanceSplineLo),...
            };    
    elseif Pack.TypeD==6
        Model.nTheta2= 17;
        Dist_Names={
            'Distance x Non High School Mother',...
            'Distance x High School Mother ',...
            'Distance x 2 year Technical Degree Mother ',...
            'Distance x 4 year College Degree Mother ',...
            'Distance x Poor',...
            strcat('Distance > ',Pack.DistanceSplineHi),...
            strcat('Distance < ',Pack.DistanceSplineLo),...
            };    
    elseif Pack.TypeD==7 || Pack.TypeD==8 || Pack.TypeD==9
        Model.nTheta2= 35;
        Dist_Names={
            'Distance x Non High School Mother',...
            'Distance x High School Mother ',...
            'Distance x 2 year Technical Degree Mother ',...
            'Distance x 4 year College Degree Mother ',...
            'Distance x Poor',...
            'Distance x Distance >=1 x Non High School Mother',...
            'Distance x Distance >=1 x High School Mother ',...
            'Distance x Distance >=1 x 2 year Technical Degree Mother ',...
            'Distance x Distance >=1 x 4 year College Degree Mother ',...
            'Distance x Distance >=1 x Poor',...
            'Distance x Distance >=2 x Non High School Mother',...
            'Distance x Distance >=2 x High School Mother ',...
            'Distance x Distance >=2 x 2 year Technical Degree Mother ',...
            'Distance x Distance >=2 x 4 year College Degree Mother ',...
            'Distance x Distance >=2 x Poor',...
            'Distance x Distance >=3 x Non High School Mother',...
            'Distance x Distance >=3 x High School Mother ',...
            'Distance x Distance >=3 x 2 year Technical Degree Mother ',...
            'Distance x Distance >=3 x 4 year College Degree Mother ',...
            'Distance x Distance >=3 x Poor',...
            'Distance x Distance >= 4 x Non High School Mother',...
            'Distance x Distance >= 4 x High School Mother ',...
            'Distance x Distance >= 4 x 2 year Technical Degree Mother ',...
            'Distance x Distance >= 4 x 4 year College Degree Mother ',...
            'Distance x Distance >= 4 x Poor',...
            };    
    end

    % Random Coefficients
    %RC_Names={'Sigma Preference - Quality','Correlation Preference x Talent'};
    RC_Names={'Sigma Preference - Quality'};
    % SigGuess=1;

    Pack.Theta1_Names=Xj_Names;
    Pack.Theta2_Names={Xij_Names{:},Dist_Names{:},RC_Names{:}};
    %Pack.Theta3_Names=Test_Names;

    Pack.IV_Names=strvcat(strvcat(ivVarnames),strvcat(Pack.Theta1_Names(2:end,:)));

    % Index for finding parameters
    if Pack.TypeD==1
        Model.BetaO          = logical([ 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 ]);   % 9 parameters
        Model.Lambda         = logical([ 0 0 0 0 0 0 0 0 0 1 1 1 1 1 0 ]);      % 4 parameters
        Model.SigRC          = logical([ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 ]);     % 3 parameters

        Model.BetaEduO       = logical([ 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0]);   % 2 parameters
        Model.BetaSEPO       = logical([ 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0]);   % 1 parameters
        Model.AlphaEduO      = logical([ 0 0 0 0 1 1 1 1 0 0 0 0 0 0 0]);   % 2 parameters
        Model.AlphaSEPO      = logical([ 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0]);   % 1 parameters
        Model.LambdaEduO     = logical([ 0 0 0 0 0 0 0 0 0 1 1 1 1 0 0]);   % 2 parameters
        Model.LambdaSEPO     = logical([ 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0]);   % 1 parameters
        Model.LambdaSQ       = logical([ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]);   % 1 parameters     
        Model.LambdaSplLo    = logical([ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]);   % 1 parameters 
        Model.LambdaSplHi    = logical([ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]);   % 1 parameters       
        Model.LambdaSantiago = logical([ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]);   % 0 parameters 
        Model.SigRC          = logical([ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1]);   % 3 parameters
    elseif Pack.TypeD==2
        Model.BetaO       = logical([ 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 ]);   % 9 parameters
        Model.Lambda      = logical([ 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 0 ]);      % 5 parameters
        Model.SigRC       = logical([ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 ]);     % 3 parameters

        Model.BetaEduO       = logical([ 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0]);   % 2 parameters
        Model.BetaSEPO       = logical([ 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0]);   % 1 parameters
        Model.AlphaEduO      = logical([ 0 0 0 0 1 1 1 1 0 0 0 0 0 0 0 0]);   % 2 parameters
        Model.AlphaSEPO      = logical([ 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0]);   % 1 parameters
        Model.LambdaEduO     = logical([ 0 0 0 0 0 0 0 0 0 1 1 1 1 0 0 0]);   % 2 parameters
        Model.LambdaSEPO     = logical([ 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0]);   % 1 parameters
        Model.LambdaEduO12   = logical([ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]);   % 2 parameters
        Model.LambdaSEPO12   = logical([ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]);   % 1 parameters    
        Model.LambdaEduO23   = logical([ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]);   % 2 parameters
        Model.LambdaSEPO23   = logical([ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]);   % 1 parameters    
        Model.LambdaEduO34   = logical([ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]);   % 2 parameters
        Model.LambdaSEPO34   = logical([ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]);   % 1 parameters    
        Model.LambdaEduO4p   = logical([ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]);   % 2 parameters
        Model.LambdaSEPO4p   = logical([ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]);   % 1 parameters       
        Model.LambdaSQ       = logical([ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0]);   % 1 parameters         
        Model.LambdaSplLo    = logical([ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]);   % 1 parameters 
        Model.LambdaSplHi    = logical([ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]);   % 1 parameters    
        Model.LambdaSantiago = logical([ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]);   % 0 parameters     
        Model.SigRC          = logical([ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1]);   % 3 parameters    
    elseif Pack.TypeD==3
        Model.BetaO          = logical([ 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 ]);   % 9 parameters
        Model.Lambda         = logical([ 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 0 ]);      % 5 parameters
        Model.SigRC          = logical([ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 ]);     % 3 parameters

        Model.BetaEduO       = logical([ 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0]);   % 2 parameters
        Model.BetaSEPO       = logical([ 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0]);   % 1 parameters
        Model.AlphaEduO      = logical([ 0 0 0 0 1 1 1 1 0 0 0 0 0 0 0 0]);   % 2 parameters
        Model.AlphaSEPO      = logical([ 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0]);   % 1 parameters
        Model.LambdaEduO     = logical([ 0 0 0 0 0 0 0 0 0 1 1 1 1 0 0 0]);   % 2 parameters
        Model.LambdaSEPO     = logical([ 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0]);   % 1 parameters
        Model.LambdaEduO12   = logical([ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]);   % 2 parameters
        Model.LambdaSEPO12   = logical([ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]);   % 1 parameters    
        Model.LambdaEduO23   = logical([ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]);   % 2 parameters
        Model.LambdaSEPO23   = logical([ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]);   % 1 parameters    
        Model.LambdaEduO34   = logical([ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]);   % 2 parameters
        Model.LambdaSEPO34   = logical([ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]);   % 1 parameters    
        Model.LambdaEduO4p   = logical([ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]);   % 2 parameters
        Model.LambdaSEPO4p   = logical([ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]);   % 1 parameters       
        Model.LambdaSQ       = logical([ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]);   % 0 parameters             
        Model.LambdaSplLo    = logical([ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0]);   % 1 parameters 
        Model.LambdaSplHi    = logical([ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]);   % 1 parameters    
        Model.LambdaSantiago = logical([ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0]);   % 0 parameters     
        Model.SigRC          = logical([ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1]);   % 3 parameters    
    elseif Pack.TypeD==4
        Model.BetaO       = logical([ 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 ]);   % 9 parameters
        Model.Lambda      = logical([ 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 0 ]);      % 5 parameters
        Model.SigRC       = logical([ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 ]);     % 3 parameters

        Model.BetaEduO       = logical([ 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0]);   % 2 parameters
        Model.BetaSEPO       = logical([ 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0]);   % 1 parameters
        Model.AlphaEduO      = logical([ 0 0 0 0 1 1 1 1 0 0 0 0 0 0 0 0]);   % 2 parameters
        Model.AlphaSEPO      = logical([ 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0]);   % 1 parameters
        Model.LambdaEduO     = logical([ 0 0 0 0 0 0 0 0 0 1 1 1 1 0 0 0]);   % 2 parameters
        Model.LambdaSEPO     = logical([ 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0]);   % 1 parameters
        Model.LambdaEduO12   = logical([ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]);   % 2 parameters
        Model.LambdaSEPO12   = logical([ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]);   % 1 parameters    
        Model.LambdaEduO23   = logical([ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]);   % 2 parameters
        Model.LambdaSEPO23   = logical([ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]);   % 1 parameters    
        Model.LambdaEduO34   = logical([ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]);   % 2 parameters
        Model.LambdaSEPO34   = logical([ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]);   % 1 parameters    
        Model.LambdaEduO4p   = logical([ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]);   % 2 parameters
        Model.LambdaSEPO4p   = logical([ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]);   % 1 parameters       
        Model.LambdaSQ       = logical([ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]);   % 0 parameters            
        Model.LambdaSplLo    = logical([ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]);   % 1 parameters 
        Model.LambdaSplHi    = logical([ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0]);   % 1 parameters     
        Model.LambdaSantiago = logical([ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]);   % 0 parameters    
        Model.SigRC          = logical([ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1]);   % 3 parameters    
    elseif Pack.TypeD==5
        Model.BetaO       = logical([ 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 ]);   % 9 parameters
        Model.Lambda      = logical([ 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 0 ]);      % 5 parameters
        Model.SigRC       = logical([ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 ]);     % 3 parameters

        Model.BetaEduO       = logical([ 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0]);   % 2 parameters
        Model.BetaSEPO       = logical([ 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0]);   % 1 parameters
        Model.AlphaEduO      = logical([ 0 0 0 0 1 1 1 1 0 0 0 0 0 0 0 0]);   % 2 parameters
        Model.AlphaSEPO      = logical([ 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0]);   % 1 parameters
        Model.LambdaEduO     = logical([ 0 0 0 0 0 0 0 0 0 1 1 1 1 0 0 0]);   % 2 parameters
        Model.LambdaSEPO     = logical([ 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0]);   % 1 parameters
        Model.LambdaEduO12   = logical([ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ]);   % 2 parameters
        Model.LambdaSEPO12   = logical([ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ]);   % 1 parameters    
        Model.LambdaEduO23   = logical([ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ]);   % 2 parameters
        Model.LambdaSEPO23   = logical([ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ]);   % 1 parameters    
        Model.LambdaEduO34   = logical([ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ]);   % 2 parameters
        Model.LambdaSEPO34   = logical([ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ]);   % 1 parameters    
        Model.LambdaEduO4p   = logical([ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ]);   % 2 parameters
        Model.LambdaSEPO4p   = logical([ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ]);   % 1 parameters        
        Model.LambdaSQ       = logical([ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]);   % 0 parameters            
        Model.LambdaSplLo    = logical([ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0]);   % 1 parameters 
        Model.LambdaSplHi    = logical([ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]);   % 1 parameters     
        Model.LambdaSantiago = logical([ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]);   % 0 parameters    
        Model.SigRC          = logical([ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1]);   % 3 parameters  
    elseif Pack.TypeD==6
        Model.BetaO          = logical([ 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 ]);   % 9 parameters
        Model.Lambda         = logical([ 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 0 ]);      % 5 parameters
        Model.SigRC          = logical([ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 ]);     % 3 parameters

        Model.BetaEduO       = logical([ 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ]);   % 2 parameters
        Model.BetaSEPO       = logical([ 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 ]);   % 1 parameters
        Model.AlphaEduO      = logical([ 0 0 0 0 1 1 1 1 0 0 0 0 0 0 0 0 0 ]);   % 2 parameters
        Model.AlphaSEPO      = logical([ 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 ]);   % 1 parameters
        Model.LambdaEduO     = logical([ 0 0 0 0 0 0 0 0 0 1 1 1 1 0 0 0 0 ]);   % 2 parameters
        Model.LambdaSEPO     = logical([ 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 ]);   % 1 parameters
        Model.LambdaEduO12   = logical([ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ]);   % 2 parameters
        Model.LambdaSEPO12   = logical([ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ]);   % 1 parameters    
        Model.LambdaEduO23   = logical([ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ]);   % 2 parameters
        Model.LambdaSEPO23   = logical([ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ]);   % 1 parameters    
        Model.LambdaEduO34   = logical([ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ]);   % 2 parameters
        Model.LambdaSEPO34   = logical([ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ]);   % 1 parameters    
        Model.LambdaEduO4p   = logical([ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ]);   % 2 parameters
        Model.LambdaSEPO4p   = logical([ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ]);   % 1 parameters      
        Model.LambdaSQ       = logical([ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ]);   % 0 parameters            
        Model.LambdaSplLo    = logical([ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 ]);   % 1 parameters 
        Model.LambdaSplHi    = logical([ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 ]);   % 1 parameters     
        Model.LambdaSantiago = logical([ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ]);   % 0 parameters    
        Model.SigRC          = logical([ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 ]);   % 3 parameters    
    elseif Pack.TypeD==7 || Pack.TypeD==8 || Pack.TypeD==9
        Model.BetaO          = logical([ 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ]);   % 9 parameters
        Model.Lambda         = logical([ 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 ]);      % 5 parameters
        Model.SigRC          = logical([ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 ]);     % 3 parameters

        Model.BetaEduO       = logical([ 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ]);   % 2 parameters
        Model.BetaSEPO       = logical([ 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ]);   % 1 parameters
        Model.AlphaEduO      = logical([ 0 0 0 0 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ]);   % 2 parameters
        Model.AlphaSEPO      = logical([ 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ]);   % 1 parameters
        Model.LambdaEduO     = logical([ 0 0 0 0 0 0 0 0 0 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ]);   % 2 parameters
        Model.LambdaSEPO     = logical([ 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ]);   % 1 parameters
        Model.LambdaEduO12   = logical([ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ]);   % 2 parameters
        Model.LambdaSEPO12   = logical([ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ]);   % 1 parameters    
        Model.LambdaEduO23   = logical([ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 ]);   % 2 parameters
        Model.LambdaSEPO23   = logical([ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 ]);   % 1 parameters    
        Model.LambdaEduO34   = logical([ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 0 0 0 0 0 0 0 ]);   % 2 parameters
        Model.LambdaSEPO34   = logical([ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 ]);   % 1 parameters    
        Model.LambdaEduO4p   = logical([ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 0 0 ]);   % 2 parameters
        Model.LambdaSEPO4p   = logical([ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 ]);   % 1 parameters    
        Model.LambdaSQ       = logical([ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ]);   % 0 parameters            
        Model.LambdaSplLo    = logical([ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ]);   % 1 parameters 
        Model.LambdaSplHi    = logical([ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ]);   % 1 parameters     
        Model.LambdaSantiago = logical([ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ]);   % 0 parameters    
        Model.SigRC          = logical([ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 ]);   % 3 parameters    
    end
    Model.TypeD = Pack.TypeD;
    Model.DistanceSplineLo = Pack.DistanceSplineLo;
    Model.DistanceSplineHi = Pack.DistanceSplineHi;
    Pack.Colors=colormap(parula(6));



    switch Pack.Model
        case {'Model none','Model 0','Model 1','Model 1b'}   
            Model.Theta1_Guess=[1.3919,...
                                -0.0875,...
                                -0.2721,...
                                -0.61,...
                                -0.2808,...
                                0.883,...
                                -0.1573,...
                                1.0608,...
                                0.8924,...
                                -1.4070,...
                                -0.4678,...
                                -0.5285];


            Model.Theta2_Guess=[ 1.2301  ,...
                                1.6159,...
                                2.1341,...
                                -0.3870,...
                                -1.1062,...
                                -0.3123,...
                                -0.0975,...
                                -0.5,...
                                -0.7883,...
                                -1.4268,...
                                -0.9387,...
                                -0.8800,...
                                -0.7157,...
                                -0.2420,...
                                1.0012];


        case {'Model 2','Model 2b'}    
                    Model.Theta1_Guess=1;


            Model.Theta2_Guess=[    
                2.7294,...
        3.4956,...
        3.6738,...
       -0.9463,...
       -1.0474,...
       -0.4670,...
       -0.2194,...
       -0.0138,...
       -0.8383,...
       -2.0599,...
       -1.4411,...
       -1.2975,...
       -1.2478,...
       -0.2583,...
        1.0000];


    end



end
