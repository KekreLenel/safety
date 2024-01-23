%% estimateIRFs.m
% Estimates effects of safety shocks using recursive VAR

runTime = tic;
rng('default');

%% Import and preliminaries

warning('off','MATLAB:catenate:DimensionMismatch');

inputNames = {'output/monthly_var.csv','output/quarterly_var.csv'};
fileNames = {'output/fig_9','output/fig_10'};

safetyCell = {'cip_govt_3m_2020','cip_govt_3m_2020'};

otherVarsCell = {{'log_q_eom','diffy_3m_2020','excessequity','log_ip_oecd_USD','diff_log_ip'},...
    {'log_q_eom','diffy_3m_2020','excessequity3mo','log_ip_oecd_USD','diff_log_ip'}};

varNamesToPlotCell = {{'3-mo swapped G10-Tbill','Real exchange rate','3-mo G10-Tbill','1-mo real excess MSCI ACWI','U.S. IP','G10-U.S. IP'},...
    {'3-mo swapped G10-Tbill','Real exchange rate','3-mo G10-Tbill','3-mo real excess MSCI ACWI','U.S. IP','G10-U.S. IP'}};
varYears = [1995,2019;1995,2019];
varMonths = [1,12;3,12];

orderOtherCell = {[2 3 4 5 6],[2 3 4 5 6]}; % where other vars should be ordered in VAR
% (precise ordering after safety doesn't matter)
lagsVec = [4,1]; % number of lags of variables to include in VAR
numBoots = 10000;
TVec = [48,12]; % number of periods to be plotted in IRFs
lCI = 5;
uCI = 95;

sizeSafety = 0.25; % sign of innovation to safety to plot
scaleToPlot = 100; % plots increase in safety with all vars in basis points 
rowsToPlot = [2,2];
colsToPlot = [3,3];
plotIRFs = [1,0];

%% Loop over specs

for spec = 1:length(lagsVec)  
    inputName = inputNames{spec};
    fileName = fileNames{spec};
    
    delete([fileName,'.txt']); % clear any old diary
    diary([fileName,'.txt']);
    
    disp(['year range: [',num2str(varYears(spec,1)),',',num2str(varYears(spec,2)),']']);
    disp(['month range: [',num2str(varMonths(spec,1)),',',num2str(varMonths(spec,2)),']']);
    disp(['lags: ',num2str(lagsVec(spec))]);    
    
    safety = safetyCell{spec};
    lags = lagsVec(spec);
    otherVars = otherVarsCell{spec};
    varNamesToPlot = varNamesToPlotCell{spec};
    orderOther = orderOtherCell{spec};
    T = TVec(spec);
    
    % Define data and variable names
    datafile = importdata(inputName);
    data = datafile.data;
    names = datafile.colheaders; 
    
    % Extract relevant subsets of data
    varData = data(find((varYears(spec,1) == data(:,strcmp('year',names))) & ...
        (varMonths(spec,1) == data(:,strcmp('month',names)))):...
        find((varYears(spec,2) == data(:,strcmp('year',names))) & ...
        (varMonths(spec,2) == data(:,strcmp('month',names)))),:);
    numVars2 = size(otherVars,2);
    numVars = numVars2+1;
              
    %% Run recursive VAR with safety ordered first
    
    % Generate data for recursive VAR  
    varInds = nan(1,numVars);
    varInds(1) = find(strcmp(safety,names),1);
    for j=1:numVars2
        varInds(orderOther(j)) = find(strcmp(otherVars{j},names),1);
    end    
    var2Inds = varInds(2:end);
    
    rvar1Dep = varData(lags+1:size(varData,1),find(strcmp(safety,names),1));
    rvar1Indep = [lagmatrix(varData(:,find(strcmp(safety,names),1)),1:lags),...
        lagmatrix(varData(:,var2Inds),1:lags)];
    rvar1Indep = rvar1Indep(lags+1:size(rvar1Indep,1),:);
    rvar2Dep = varData(lags+1:size(varData,1),var2Inds);
    rvar2Indep = [lagmatrix(varData(:,find(strcmp(safety,names),1)),0:lags),...
        lagmatrix(varData(:,var2Inds),1:lags)];
    rvar2Indep = rvar2Indep(lags+1:size(rvar2Indep,1),:);
    
    % Estimate recursive VAR 
    rvar1Beta = [rvar1Indep,ones(size(rvar1Indep,1),1)]\rvar1Dep;
    rvar1Resid = rvar1Dep - [rvar1Indep,ones(size(rvar1Indep,1),1)]*rvar1Beta;
    rvar2Beta = [rvar2Indep,ones(size(rvar2Indep,1),1)]\rvar2Dep;
    rvar2Resid = rvar2Dep - [rvar2Indep,ones(size(rvar2Indep,1),1)]*rvar2Beta;    
    
    % Stack VAR into a first-order system
    % Let rA1, rA2 be matrices of coefficients on non-constants
    rA1 = nan(lags+1,lags*numVars);
    rA1(1,:) = rvar1Beta(1:(size(rvar1Beta,1)-1))';
    rA1(2:(lags+1),:) = [eye(lags,lags),zeros(lags,lags*(numVars-1))];
    rA2 = nan(numVars2*lags,(lags+1)+lags*numVars2);
    rA2(1:numVars2,:) = rvar2Beta(1:(size(rvar2Beta,1)-1),:)';
    rA2(numVars2+1:numVars2*lags,:) = [zeros(numVars2*(lags-1),(lags+1)),eye(numVars2*(lags-1),numVars2*(lags-1)),zeros(numVars2*(lags-1),numVars2)];
    
    % Generate impulse responses
    rirf = nan(10000,numVars);
    rendo1 = zeros(10000,lags+1);
    rendo2 = zeros(10000,numVars2*lags);
    rendo1(1,:) = [sizeSafety,zeros(1,lags)];
    rendo2(1,:) = (rA2*[rendo1(1,:)';zeros(lags*numVars2,1)])';
    rirf(1,:) = [rendo1(1,1),rendo2(1,1:numVars2)];
    for t=2:10000
        rendo1(t,:) = (rA1*[rendo1(t-1,1:lags)';rendo2(t-1,:)'])';
        rendo2(t,:) = (rA2*[rendo1(t,:)';rendo2(t-1,:)'])';
        rirf(t,:) = [rendo1(t,1),rendo2(t,1:numVars2)];
    end    
        
    %% Run bootstrap following Christiano-Eichenbaum-Evans (1999) for IRF standard errors 

    rirfBoots = nan(size(rirf,1),size(rirf,2),numBoots);
    residToPick = floor(rand(size(rvar1Resid,1),numBoots)*size(rvar1Resid,1))+1;
    for nBoot=1:numBoots
        % Generate data for recursive VAR
        rendo1Boot = zeros(size(rvar1Dep,1)+1,lags+1);
        rendo2Boot = zeros(size(rvar1Dep,1)+1,numVars2*lags);
        rendo1Boot(1,1:lags) = rvar2Indep(1,2:(lags+1));
        rendo2Boot(1,1:numVars2*lags) = rvar2Indep(1,(lags+1)+1:size(rvar2Indep,2));
        ru1 = rvar1Resid(residToPick(:,nBoot),:);
        ru2 = rvar2Resid(residToPick(:,nBoot),:); 
        for t=2:size(rendo1Boot,1)
            rendo1Boot(t,:) = (rA1*[rendo1Boot(t-1,1:lags)';rendo2Boot(t-1,:)'])' + ...
                [rvar1Beta(end),zeros(1,lags)] + ...
                [ru1(t-1,:),zeros(1,lags)];       
            rendo2Boot(t,:) = (rA2*[rendo1Boot(t,:)';rendo2Boot(t-1,:)'])' + ...
                [rvar2Beta(end,:),zeros(1,numVars2*(lags-1))] + ...
                [ru2(t-1,:),zeros(1,numVars2*(lags-1))];
        end
        rvar1DataBootAftInit = rendo1Boot(2:size(rendo1Boot,1),1);
        rvar1DataBoot = [varData(1:lags,find(strcmp(safety,names),1));rvar1DataBootAftInit];
        rvar2DataBootAftInit = rendo2Boot(2:size(rendo2Boot,1),1:numVars2);
        rvar2DataBoot = [varData(1:lags,var2Inds);rvar2DataBootAftInit];
        
        % Run recursive VAR on simulated data
        rvar1DepBoot = rvar1DataBoot(lags+1:size(rvar1DataBoot,1));
        rvar1IndepBoot = [lagmatrix(rvar1DataBoot,1:lags),...
            lagmatrix(rvar2DataBoot,1:lags)];
        rvar1IndepBoot = rvar1IndepBoot(lags+1:size(rvar1IndepBoot,1),:);
        rvar2DepBoot = rvar2DataBoot(lags+1:size(rvar2DataBoot,1),:);
        rvar2IndepBoot = [lagmatrix(rvar1DataBoot,0:lags),...
            lagmatrix(rvar2DataBoot,1:lags)];
        rvar2IndepBoot = rvar2IndepBoot(lags+1:size(rvar2IndepBoot,1),:);

        rvar1BetaBoot = [rvar1IndepBoot,ones(size(rvar1IndepBoot,1),1)]\rvar1DepBoot;
        rvar1ResidBoot = rvar1DepBoot - [rvar1IndepBoot,ones(size(rvar1IndepBoot,1),1)]*rvar1BetaBoot;
        rvar2BetaBoot = [rvar2IndepBoot,ones(size(rvar2IndepBoot,1),1)]\rvar2DepBoot;        
        
        % Generate impulse responses and save
        rA1Boot = nan(lags+1,lags*numVars);
        rA1Boot(1,:) = rvar1BetaBoot(1:(size(rvar1BetaBoot,1)-1))';
        rA1Boot(2:(lags+1),:) = [eye(lags,lags),zeros(lags,lags*(numVars-1))];
        rA2Boot = nan(numVars2*lags,(lags+1)+lags*numVars2);
        rA2Boot(1:numVars2,:) = rvar2BetaBoot(1:(size(rvar2Beta,1)-1),:)';
        rA2Boot(numVars2+1:numVars2*lags,:) = [zeros(numVars2*(lags-1),lags+1),eye(numVars2*(lags-1),numVars2*(lags-1)),zeros(numVars2*(lags-1),numVars2)];   
        
        rirfBoot = nan(size(rirf,1),numVars);
        rendo1Boot = zeros(size(rirf,1),lags+1);
        rendo2Boot = zeros(size(rirf,1),numVars2*lags);
        rendo1Boot(1,:) = [sizeSafety,zeros(1,lags)];
        rendo2Boot(1,:) = (rA2Boot*[rendo1Boot(1,:)';zeros(lags*numVars2,1)])';
        rirfBoot(1,:) = [rendo1Boot(1,1),rendo2Boot(1,1:numVars2)];
        for t=2:size(rirf,1)
            rendo1Boot(t,:) = (rA1Boot*[rendo1Boot(t-1,1:lags)';rendo2Boot(t-1,:)'])';
            rendo2Boot(t,:) = (rA2Boot*[rendo1Boot(t,:)';rendo2Boot(t-1,:)'])';
            rirfBoot(t,:) = [rendo1Boot(t,1),rendo2Boot(t,1:numVars2)];
        end            
        rirfBoots(:,:,nBoot) = rirfBoot;
        
        if mod(nBoot,100) == 0
            disp(num2str(nBoot));
        end
    end
    rirflCI = nan(size(rirf));
    rirfuCI = nan(size(rirf));
    for t=1:size(rirf,1)
        rirflCI(t,:) = prctile(squeeze(rirfBoots(t,:,:))',lCI);
        rirfuCI(t,:) = prctile(squeeze(rirfBoots(t,:,:))',uCI);
    end    
    fprintf('\n');

    %% Plot

    Interpreter = 'latex';
    fontsize = 14;
    % Define colors
    blue        = [0       0.4470   0.7410];
    lightblue   = [107,174,214]./255;
    darkblue    = [8,69,148]./255;
    lightred    = [0.9500  0.4250   0.2];
    darkred     = [0.6350  0.0780   0.1840];
    red         = [0.8500  0.3250   0.0980];
    yellow      = [0.70    0.640    0.1250];
    green       = [0.1     0.75     0.2];
    grey        = [0.4     0.4      0.4];    
    
    mapIRFtoPlot = [1,orderOther];
           
    if plotIRFs(spec)
        rvarfig = figure;
        for i=1:numVars
            subplot(rowsToPlot(spec),colsToPlot(spec),i);
            yplot = scaleToPlot*[rirflCI(1:T,mapIRFtoPlot(i))', flipud(rirfuCI(1:T,mapIRFtoPlot(i)))'];
            plot(0:(T-1),zeros(1,T),'-k');
            hold on;
            fill([0:(T-1),(T-1):-1:0],yplot,1,'facecolor', lightblue, 'edgecolor', 'none','facealpha',0.5);
            plot(0:(T-1),scaleToPlot*rirf(1:T,mapIRFtoPlot(i))','LineWidth',2,'Color',blue);    
            hold off;
            xlim([0,(T-1)]);
            ylabel('bp','Fontsize',fontsize,'Interpreter',Interpreter);
            title(varNamesToPlot{i},'Fontsize',fontsize,'Interpreter',Interpreter);
            ax2 = gca;
            ax2.FontSize = fontsize;
            ax2.TickLabelInterpreter = Interpreter;
        end
        set(gcf,'PaperType','usletter','PaperOrientation','portrait','PaperPosition',[0.0 0.0 8.5 2.5*rowsToPlot(spec)],'PaperSize',[8.5 5.0]);
        print(rvarfig, fileName,'-dpdf', '-r600');
    end
        
    %% Export IRFs to excel
    toexcel = scaleToPlot*[rirf(1:T,mapIRFtoPlot),rirflCI(1:T,mapIRFtoPlot),rirfuCI(1:T,mapIRFtoPlot)];
    tabletoexcel = array2table(toexcel);
    tabletoexcel.Properties.VariableNames = ...
        [[safetyCell{spec},otherVarsCell{spec}],...
        strcat([safetyCell{spec},otherVarsCell{spec}],'_lb'),...
        strcat([safetyCell{spec},otherVarsCell{spec}],'_ub')];
    writetable(tabletoexcel,[fileName,'.xlsx'],'Sheet',1,'Range','A1');
end

diary off;

disp(['Total runtime: ',num2str(toc(runTime)),'s']);