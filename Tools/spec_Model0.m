  % Simple model with minimal Xs 
                xVarnames={ 
                'Voucher School                      ',...
                'For Profit                 x Voucher',... 
                'Religious - Non Catholic   x Voucher',...   
                'Religious - Catholic       x Voucher',...   
                'Has High School            x Voucher',...       
                'Old                        x Voucher',...   
                'Brand New                  x Voucher',...
                'Private Non Voucher School          ',...    
                'Religious                  x Private',...    
                'Religious - Catholic       x Private',...  
                'Has High School            x Private',...     
                'Old                        x Private',...
                'Brand New                  x Private',...
                }';     
                NewSchool=(isnan(schools.ShareSEP(:,1))==1 & isnan(schools.ShareSEP(:,2))==0);
                ClosedSchool=(isnan(schools.ShareSEP(:,1))==0 & isnan(schools.ShareSEP(:,2))==1);
                % sum([NewSchool ClosedSchool])
               
                xx=[schools.DEP==2 schools.ForProfit.*schools.DEP==2  schools.Religious.*schools.DEP==2 schools.Religion_Catholic.*schools.DEP==2 schools.Continuation.*schools.DEP==2 schools.History(:,2) NewSchool.*schools.DEP==2,...
                    schools.DEP==3                                    schools.Religious.*schools.DEP==3 schools.Religion_Catholic.*schools.DEP==3 schools.Continuation.*schools.DEP==3 schools.History(:,3) NewSchool.*schools.DEP==3,...
                    ];
                if rank(xx)<size(xx,2)
                        error('Rank of xx')
                end
                
                Model.xxVarnames=xVarnames;
                
                
                %% Fixed Effects to concentrate out 
                % Model Spec 0 - Market and Time FE 
                
                % Setup matrix for FWL and concentrate out FEs
                FE=full(Model.MarketTimeFE);
                fe=[FE   ];%DComunaPublic
                A_FE=eye(size(fe,1))-fe*inv(fe'*fe)*fe';
                Model.A_FE=sparse(A_FE);   %       spy(Model.A_FE)
                
                
              %% Additional IV Instruments 
                
                % Setup IVs
                ivVarnames='';ivTransfer=[];
                tempLabels={ 'Change Mg Transfers (fixed shares)','Change Mg Transfers (fixed shares) x High Exposure'};
                tempNames={'Public Schools ','Voucher Schools '};
               
                for k=1:length(tempNames) 
                for l=1:length(tempLabels)
                ivVarnames=strvcat(ivVarnames,  [tempLabels{l} ' X ' tempNames{k}]);
                end
                end
                
                
                ivTransfer=[ (NewSchool==0).*[ schools.dMgTransfer schools.dMgTransfer.*schools.HighFixedExposure(:,1) ].*(schools.DEP==1)    (NewSchool==0).*[ schools.dMgTransfer schools.dMgTransfer.*schools.HighFixedExposure(:,1)].*(schools.DEP==2) ];  
                ivTransfer(isnan(ivTransfer))=0;                
                %ivVarnames=strvcat(ivVarnames,'Not in Both');
%                 ivVarnames=strvcat(ivVarnames,'PriorShare');
%                 ivTransfer=[ivTransfer schools.SharesepPre];
                
                % Policy Exposure to Poor Students  
                ivPolicy=table2array([schools(:,strcmp(schools.Properties.VariableNames','PrePolicy Share SEP (1km)')),...
                                      schools(:,strcmp(schools.Properties.VariableNames','Post Policy x PrePolicy Share SEP (1km)'))]);
                tempLabels={'Share SEP (1km)','Post Policy x PrePolicy Share SEP (1km)'};
                tempNames={'Public Schools ','Voucher Schools '};
                for k=1:2
                ivVarnames=strvcat(ivVarnames,  [tempLabels{k} ' X ' tempNames{1}]);
                ivVarnames=strvcat(ivVarnames,  [tempLabels{k} ' X ' tempNames{2}]);
                end
                ivPolicy=schools.HighFixedExposure;
                %ivPolicy=ivPolicy.*(ivPolicy>prctile(ivPolicy(schools.Year==2007,1),80) )>0; 
                % prctile(ivPolicy(schools.Year==2007,1),80) 
                % prctile(table2array(schools(schools.Year==2007,{'PrePolicy Share SEP (1km)'})),80)
                PolicyIV=[ivPolicy(:,1).*[schools.DEP==1 schools.DEP==2 ] ivPolicy(:,2).*[schools.DEP==1 schools.DEP==2 ]];
 
                %sum(schools.HighExposure)
%                 tempLabels={'Dist to closest SEP'};   
                 tempLabels={'N Public SEP (0.5km)','N Voucher SEP (0.5km)','N Voucher Non SEP (0.5km)'};       %'Transfers Private (1km)','Transfers to Public (1km)','Market Growth','Market Share / N schools (0.5km)',...
                 tempNames={'Public Schools  ','Voucher Schools '};
                 tempLabelsV2={'Post Policy'};    
                 OtherXsIV=[];
%                 for k=1:length(tempNames)
%                     for l=1:length(tempLabels)
%                         
%                     strLabel=[tempLabels{l} repmat(' ',1,40-length(tempLabels{l}))   ' x ' tempNames{k}];    
%                     ivVarnames=strvcat(ivVarnames,  strLabel);
%                     d=table2array(iv(:,strcmp(iv.Properties.VariableNames',tempLabels{l})));
%                     d=(d>0);
%                     var=d.*(schools.DEP==k);
%                     OtherXsIV=[OtherXsIV var];
%                     
%                     ivVarnames=strvcat(ivVarnames,  [strLabel ' - ' tempLabelsV2{1}]);
%                     OtherXsIV=[OtherXsIV var.*(schools.Year>=2008)];
%                     end
%                 end

%                %Utilities 
                 ivMgUtilitiesCosts=[];UtilitiesIV=[];
%                 ivMgUtilitiesCosts(:,1)=table2array([schools(:,strcmp(schools.Properties.VariableNames','CostAgua'))]);                          
%                 tempLabels={'Utilities - Water'};
%                 tempNames={'Public Schools ','Voucher Schools '};
%                 for k=1:2
%                 ivVarnames=strvcat(ivVarnames,  [tempLabels{1} '_X_' tempNames{k}]);
%                 %ivVarnames=strvcat(ivVarnames,  [tempLabels{2} '_X_' tempNames{k}]);
%                 end
%                 UtilitiesIV=[ivMgUtilitiesCosts.*(schools.DEP==1) ivMgUtilitiesCosts.*(schools.DEP==2)];

                % Labor Markets 
                 listR={'FE_OtherIndustries',...                                     
                        %'WageBill_OtherIndustries',...% 'WageBill_Education',...  % 'xb_OtherIndustries',...};'FE_Education',...  'xb_OtherIndustries',...  FE_OtherIndustries       
%                        'Revenue_OtherIndustries',... % 'Revenue_Education' 'xb_OtherIndustries',...
%                        'MedianWages_Other',...       % 'MedianWages_Teachers',...
                       };               
                ivLaborCosts=[];
                for jk=1:length(listR)
                ivLaborCosts(:,jk)=table2array([schools(:,strcmp(schools.Properties.VariableNames',listR{jk}))]);           
                ivLaborCosts((isnan(ivLaborCosts(:,jk)) |  ivLaborCosts(:,jk)==0),jk)=nanmean(ivLaborCosts(:,jk));            
                end
                ivLaborCosts=(ivLaborCosts-nanmean(ivLaborCosts))./std(ivLaborCosts);

                tempVars={'Locality Earnings Deviation - nonEducation'}; 
                tempNames={ 'Voucher Schools','Private Schools'};
                for k=1:length(tempNames)
                    for jk=1:length(listR)
                        ivVarnames=strvcat(ivVarnames,  [tempVars{jk} ' X ' tempNames{k}]);
                    end
                end
                LaborCostsIV=[ivLaborCosts.*(schools.DEP==2) ivLaborCosts.*(schools.DEP==3) ];
                
                tempIV=[ivTransfer,...
                        PolicyIV,...  
                        OtherXsIV,...
                        UtilitiesIV,...
                        LaborCostsIV,...
                                  ]; %    
                              
                Model.ivVarnames=strvcat(ivVarnames, strvcat(Model.xxVarnames));
                info.rnames=strvcat(' ',ivVarnames);
                mprint(sum(isnan(tempIV))',info)
                
                
                % Normalize VA 
                QQ=Model.Q - Model.Q(Model.normIndex);
                IVV=tempIV;
                IVV(isnan(IVV) | isinf(IVV))=0;
                dropIV=(sum(IVV)==0);
                ivVarnames(dropIV,:)=[];
                IVV(:,dropIV)=[];
          
                Model.ivVarnames=strvcat(strvcat(ivVarnames),strvcat(xVarnames));           
                XIV=[IVV xx schools.Year>2007];
                listVars=strvcat(strvcat(Model.ivVarnames),'Post Policy');
               
%                 md1=fitlm(XIV,Model.Q,'PredictorVars',listVars )
%                 md1s=fitlm(XIV,schools.Qshrunk_jt,'PredictorVars',listVars )
%                 
                XIVx=XIV;
                XIVx=XIVx-XIVx(Model.normIndex,:);
        
                Model.XIV=XIV;
                Model.XIVnames=listVars;
               
             

                % Normalize VA 
                XX=[Model.Q xx]; % rank(xx)
                Xj_Names=format_vNames(strvcat('Quality',strvcat(xVarnames(1:size(xx,2),:))));
                Model.XX=XX-XX(Model.normIndex,:); % normalized for OO school
                Model.XXX=Model.A_FE*Model.XX; % 
                IVVx=IVV-IVV(Model.normIndex,:);
                ivAll=Model.A_FE*[IVVx Model.XX(:,2:end) ]; 
                Model.IV=ivAll;
            
                rank(Model.A_FE*Model.XX(:,2:end))
                if rank(ivAll)<size(ivAll,2)
                error('Rank Issues with IVs')
                end

                

                %Model.IV=(Model.IV-mean(Model.IV))./std(Model.IV);
                fprintf('\n Value Added on XX and IV in Levels \n')
                mdlivLevels=fitlm([IVV Model.XX(:,2:end) ],Model.Q,'PredictorVars',Model.ivVarnames)
                
                fprintf('\n Value Added on XX and IV with Market Time FE concentrated out \n')
                mdlivFWL=fitlm([Model.IV ],Model.A_FE*QQ,'PredictorVars',Model.ivVarnames,'RobustOpts','on')
                % ypred = predict(mdlivFWL,Model.IV);
                % corr(ypred,Model.XXX(:,1))

                EtaHat=[Model.IV]\Model.XXX; Model.XXX(:,1);
                Xhat=[Model.IV]*EtaHat;   
                Model.XZ=Xhat ;
                Model.IPX=eye(size(Model.XZ,1))-Model.XXX*inv(Model.XZ'*Model.XZ)*Model.XZ';
