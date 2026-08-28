clear;
clc;
close all;

%% =========================================================================
%  EXPERIMENTAL LEVELS
% =========================================================================

T_levels  = [17 23 30 45 50];
RH_levels = [0 10 25 40 50 60 70];

nT   = numel(T_levels);
nRH  = numel(RH_levels);
nRep = 3;

%% =========================================================================
%  RAW REPLICATE DATA
%
%  Dimensions:
%  Row    = temperature
%  Column = RH condition
%  Page   = replicate
% =========================================================================

WU_raw = zeros(nT,nRH,nRep);

% -------------------------------------------------------------------------
% Temperature = 17 C
% -------------------------------------------------------------------------

WU_raw(1,:,1) = [0 0.124 0.240 0.352 0.411 0.492 0.521];
WU_raw(1,:,2) = [0 0.123 0.243 0.353 0.412 0.494 0.523];
WU_raw(1,:,3) = [0 0.127 0.244 0.355 0.418 0.495 0.524];

% -------------------------------------------------------------------------
% Temperature = 23 C
% -------------------------------------------------------------------------

WU_raw(2,:,1) = [0 0.151 0.332 0.451 0.904 0.854 1.242];
WU_raw(2,:,2) = [0 0.150 0.332 0.452 0.907 0.857 1.244];
WU_raw(2,:,3) = [0 0.152 0.336 0.454 0.908 0.857 1.247];

% -------------------------------------------------------------------------
% Temperature = 30 C
% -------------------------------------------------------------------------

WU_raw(3,:,1) = [0 0.208 0.384 0.602 1.013 1.153 1.332];
WU_raw(3,:,2) = [0 0.208 0.385 0.604 1.014 1.155 1.335];
WU_raw(3,:,3) = [0 0.210 0.387 0.605 1.017 1.156 1.337];

% -------------------------------------------------------------------------
% Temperature = 45 C
% -------------------------------------------------------------------------

WU_raw(4,:,1) = [0 0.102 0.261 0.322 0.597 0.734 0.780];
WU_raw(4,:,2) = [0 0.101 0.263 0.323 0.599 0.737 0.782];
WU_raw(4,:,3) = [0 0.103 0.263 0.325 0.600 0.737 0.783];

% -------------------------------------------------------------------------
% Temperature = 50 C
% -------------------------------------------------------------------------

WU_raw(5,:,1) = [0 0.053 0.189 0.260 0.313 0.551 0.572];
WU_raw(5,:,2) = [0 0.052 0.190 0.261 0.318 0.552 0.575];
WU_raw(5,:,3) = [0 0.054 0.192 0.263 0.319 0.554 0.576];

%% =========================================================================
%  CONVERT DATA TO LONG FORMAT
% =========================================================================

nConditions = nT*nRH;
n = nConditions*nRep;

T   = zeros(n,1);
RH  = zeros(n,1);
Rep = zeros(n,1);
WU  = zeros(n,1);

row = 0;

for i = 1:nT

    for j = 1:nRH

        for k = 1:nRep

            row = row + 1;

            T(row)   = T_levels(i);
            RH(row)  = RH_levels(j);
            Rep(row) = k;
            WU(row)  = WU_raw(i,j,k);

        end

    end

end

Aw = RH./100;

%% Display experimental design

DesignTable = table(T,RH,Aw,Rep,WU, ...
    'VariableNames',{'Temperature_C','RH_percent', ...
    'WaterActivity','Replicate','Observed_WU'});

fprintf('\n============================================================\n');
fprintf('EXPERIMENTAL DESIGN\n');
fprintf('============================================================\n');
fprintf('Temperature levels       = %d\n',nT);
fprintf('RH levels                = %d\n',nRH);
fprintf('Independent conditions   = %d\n',nConditions);
fprintf('Replicates per condition = %d\n',nRep);
fprintf('Total observations       = %d\n',n);

%% =========================================================================
%  DESIGN MATRICES
% =========================================================================

X0 = ones(n,1);

Xlinear = [ ...
    ones(n,1), ...
    T, ...
    Aw];

Xquadratic = [ ...
    ones(n,1), ...
    T, ...
    Aw, ...
    T.^2, ...
    Aw.^2, ...
    T.*Aw];

Xcubic = [ ...
    ones(n,1), ...
    T, ...
    Aw, ...
    T.^2, ...
    Aw.^2, ...
    T.*Aw, ...
    T.^3, ...
    Aw.^3, ...
    T.^2.*Aw, ...
    T.*Aw.^2];

termNames = { ...
    'Intercept'
    'T'
    'Aw'
    'T^2'
    'Aw^2'
    'T*Aw'
    'T^3'
    'Aw^3'
    'T^2*Aw'
    'T*Aw^2'};

%% =========================================================================
%  FIT CUBIC MODEL USING ALL 105 OBSERVATIONS
% =========================================================================

b = Xcubic \ WU;

WU_pred = Xcubic*b;
residuals = WU - WU_pred;

p = size(Xcubic,2);

%% =========================================================================
%  GOODNESS-OF-FIT AND OVERALL REGRESSION ANOVA
% =========================================================================

SSEcubic = sum(residuals.^2);
SStotal  = sum((WU-mean(WU)).^2);
SSmodel  = SStotal-SSEcubic;

dfModel    = p-1;
dfResidual = n-p;
dfTotal    = n-1;

MSmodel    = SSmodel/dfModel;
MSresidual = SSEcubic/dfResidual;

Fmodel = MSmodel/MSresidual;

% Upper-tail F-distribution probability
xModel = dfResidual/(dfResidual + dfModel*Fmodel);

Pmodel = betainc( ...
    xModel, ...
    dfResidual/2, ...
    dfModel/2);

R2 = 1-SSEcubic/SStotal;

AdjR2 = 1- ...
    (SSEcubic/dfResidual)/(SStotal/dfTotal);

% Residual standard error
RMSE = sqrt(MSresidual);

% Direct prediction-error RMSE
RMSE_prediction = sqrt(mean(residuals.^2));

%% =========================================================================
%  PRINT COEFFICIENTS
% =========================================================================

fprintf('\n============================================================\n');
fprintf('CUBIC RESPONSE-SURFACE COEFFICIENTS\n');
fprintf('============================================================\n');

for i = 1:numel(b)

    fprintf('%-12s = %+0.10f\n',termNames{i},b(i));

end

%% =========================================================================
%  PRINT OVERALL REGRESSION ANOVA
% =========================================================================

fprintf('\n============================================================\n');
fprintf('OVERALL REGRESSION ANOVA: CUBIC RESPONSE-SURFACE MODEL\n');
fprintf('============================================================\n');

fprintf('%-15s %14s %8s %14s %14s %14s\n', ...
    'Source','SS','DF','MS','F','P value');

fprintf('%-15s %14.6f %8d %14.6f %14.4f %14.4e\n', ...
    'Cubic model', ...
    SSmodel, ...
    dfModel, ...
    MSmodel, ...
    Fmodel, ...
    Pmodel);

fprintf('%-15s %14.6f %8d %14.6f\n', ...
    'Residual', ...
    SSEcubic, ...
    dfResidual, ...
    MSresidual);

fprintf('%-15s %14.6f %8d\n', ...
    'Total', ...
    SStotal, ...
    dfTotal);

fprintf('\nR-squared                  = %.6f\n',R2);
fprintf('Adjusted R-squared         = %.6f\n',AdjR2);
fprintf('Residual RMSE              = %.6f g/g\n',RMSE);
fprintf('Direct prediction RMSE     = %.6f g/g\n',RMSE_prediction);

%% =========================================================================
%  HIERARCHICAL MODEL CONTRIBUTIONS
% =========================================================================

b0 = X0\WU;
bLinear = Xlinear\WU;
bQuadratic = Xquadratic\WU;

pred0 = X0*b0;
predLinear = Xlinear*bLinear;
predQuadratic = Xquadratic*bQuadratic;

SSE0 = sum((WU-pred0).^2);
SSElinear = sum((WU-predLinear).^2);
SSEquadratic = sum((WU-predQuadratic).^2);

% Linear terms added to intercept model
SSlinearAdded = SSE0-SSElinear;
dfLinearAdded = size(Xlinear,2)-size(X0,2);
MSlinearAdded = SSlinearAdded/dfLinearAdded;
FlinearAdded = MSlinearAdded/MSresidual;

xLinear = dfResidual/(dfResidual + ...
    dfLinearAdded*FlinearAdded);

PlinearAdded = betainc( ...
    xLinear, ...
    dfResidual/2, ...
    dfLinearAdded/2);

% Quadratic terms added to linear model
SSquadraticAdded = SSElinear-SSEquadratic;
dfQuadraticAdded = size(Xquadratic,2)-size(Xlinear,2);
MSquadraticAdded = SSquadraticAdded/dfQuadraticAdded;
FquadraticAdded = MSquadraticAdded/MSresidual;

xQuadratic = dfResidual/(dfResidual + ...
    dfQuadraticAdded*FquadraticAdded);

PquadraticAdded = betainc( ...
    xQuadratic, ...
    dfResidual/2, ...
    dfQuadraticAdded/2);

% Cubic terms added to quadratic model
SScubicAdded = SSEquadratic-SSEcubic;
dfCubicAdded = size(Xcubic,2)-size(Xquadratic,2);
MScubicAdded = SScubicAdded/dfCubicAdded;
FcubicAdded = MScubicAdded/MSresidual;

xCubic = dfResidual/(dfResidual + ...
    dfCubicAdded*FcubicAdded);

PcubicAdded = betainc( ...
    xCubic, ...
    dfResidual/2, ...
    dfCubicAdded/2);

fprintf('\n============================================================\n');
fprintf('HIERARCHICAL CONTRIBUTION OF POLYNOMIAL TERMS\n');
fprintf('============================================================\n');

fprintf('%-15s %14s %8s %14s %14s %14s\n', ...
    'Component','SS','DF','MS','F','P value');

fprintf('%-15s %14.6f %8d %14.6f %14.4f %14.4e\n', ...
    'Linear', ...
    SSlinearAdded, ...
    dfLinearAdded, ...
    MSlinearAdded, ...
    FlinearAdded, ...
    PlinearAdded);

fprintf('%-15s %14.6f %8d %14.6f %14.4f %14.4e\n', ...
    'Quadratic', ...
    SSquadraticAdded, ...
    dfQuadraticAdded, ...
    MSquadraticAdded, ...
    FquadraticAdded, ...
    PquadraticAdded);

fprintf('%-15s %14.6f %8d %14.6f %14.4f %14.4e\n', ...
    'Cubic', ...
    SScubicAdded, ...
    dfCubicAdded, ...
    MScubicAdded, ...
    FcubicAdded, ...
    PcubicAdded);

fprintf('%-15s %14.6f %8d %14.6f\n', ...
    'Residual', ...
    SSEcubic, ...
    dfResidual, ...
    MSresidual);

%% =========================================================================
%  COEFFICIENT STATISTICS AND 95% CONFIDENCE INTERVALS
% =========================================================================

CovB = MSresidual*pinv(Xcubic'*Xcubic);

SE = sqrt(diag(CovB));

tStatistic = b./SE;

% Two-sided p-values for t statistics
Pcoeff = betainc( ...
    dfResidual./(dfResidual+tStatistic.^2), ...
    dfResidual/2, ...
    0.5);

alpha = 0.05;

% Two-sided 95% t critical value without tinv
zCritical = betaincinv(alpha,dfResidual/2,0.5);

tCritical = sqrt( ...
    dfResidual*(1-zCritical)/zCritical);

CI_lower = b-tCritical.*SE;
CI_upper = b+tCritical.*SE;

CoefficientTable = table( ...
    termNames, ...
    b, ...
    SE, ...
    tStatistic, ...
    Pcoeff, ...
    CI_lower, ...
    CI_upper, ...
    'VariableNames',{ ...
    'Term', ...
    'Estimate', ...
    'StandardError', ...
    'tStatistic', ...
    'PValue', ...
    'CI95_Lower', ...
    'CI95_Upper'});

fprintf('\n============================================================\n');
fprintf('COEFFICIENT STATISTICS\n');
fprintf('============================================================\n');

disp(CoefficientTable);

%% =========================================================================
%  PURE ERROR AND LACK-OF-FIT ANALYSIS
% =========================================================================

conditionMean = zeros(n,1);

for i = 1:nT

    for j = 1:nRH

        index = T==T_levels(i) & RH==RH_levels(j);

        conditionMean(index) = mean(WU(index));

    end

end

% Pure-error sum of squares
SSpureError = sum((WU-conditionMean).^2);

% Lack-of-fit sum of squares
SSlackOfFit = SSEcubic-SSpureError;

dfPureError = n-nConditions;
dfLackOfFit = nConditions-p;

MSPureError = SSpureError/dfPureError;
MSLackOfFit = SSlackOfFit/dfLackOfFit;

FlackOfFit = MSLackOfFit/MSPureError;

xLOF = dfPureError/(dfPureError + ...
    dfLackOfFit*FlackOfFit);

PlackOfFit = betainc( ...
    xLOF, ...
    dfPureError/2, ...
    dfLackOfFit/2);

fprintf('\n============================================================\n');
fprintf('LACK-OF-FIT ANALYSIS\n');
fprintf('============================================================\n');

fprintf('%-15s %14s %8s %14s %14s %14s\n', ...
    'Source','SS','DF','MS','F','P value');

fprintf('%-15s %14.6f %8d %14.6f %14.4f %14.4e\n', ...
    'Lack of fit', ...
    SSlackOfFit, ...
    dfLackOfFit, ...
    MSLackOfFit, ...
    FlackOfFit, ...
    PlackOfFit);

fprintf('%-15s %14.6f %8d %14.8f\n', ...
    'Pure error', ...
    SSpureError, ...
    dfPureError, ...
    MSPureError);

fprintf('%-15s %14.6f %8d %14.6f\n', ...
    'Residual', ...
    SSEcubic, ...
    dfResidual, ...
    MSresidual);

%% =========================================================================
%  ADD PREDICTIONS AND RESIDUALS TO DESIGN TABLE
% =========================================================================

DesignTable.Predicted_WU = WU_pred;
DesignTable.Residual = residuals;
DesignTable.ConditionMean_WU = conditionMean;

fprintf('\nFirst rows of the complete design matrix:\n');
disp(DesignTable(1:min(15,height(DesignTable)),:));

% Optional export
writetable(DesignTable,'Cubic_RSM_complete_design_matrix.xlsx');
writetable(CoefficientTable,'Cubic_RSM_coefficient_statistics.xlsx');

%% =========================================================================
%  PREDICTION FUNCTION
% =========================================================================

WU_cubic = @(Temp,RHpercent) ...
    b(1) + ...
    b(2).*Temp + ...
    b(3).*(RHpercent./100) + ...
    b(4).*Temp.^2 + ...
    b(5).*(RHpercent./100).^2 + ...
    b(6).*Temp.*(RHpercent./100) + ...
    b(7).*Temp.^3 + ...
    b(8).*(RHpercent./100).^3 + ...
    b(9).*Temp.^2.*(RHpercent./100) + ...
    b(10).*Temp.*(RHpercent./100).^2;

fprintf('\nTest prediction: WU at 30 C and 50%% RH = %.4f g/g\n', ...
    WU_cubic(30,50));

%% =========================================================================
%  CONDITION MEANS FOR SURFACE-PLOT MARKERS
% =========================================================================

[Tcondition,RHcondition] = ndgrid(T_levels,RH_levels);

Tcondition = Tcondition(:);
RHcondition = RHcondition(:);

WUconditionMean = zeros(nConditions,1);
WUconditionSD = zeros(nConditions,1);

row = 0;

for i = 1:nT

    for j = 1:nRH

        row = row+1;

        values = squeeze(WU_raw(i,j,:));

        WUconditionMean(row) = mean(values);
        WUconditionSD(row) = std(values,0);

    end

end

ConditionSummary = table( ...
    Tcondition, ...
    RHcondition, ...
    WUconditionMean, ...
    WUconditionSD, ...
    'VariableNames',{ ...
    'Temperature_C', ...
    'RH_percent', ...
    'Mean_WU', ...
    'SD_WU'});

writetable(ConditionSummary,'Cubic_RSM_condition_summary.xlsx');

%% =========================================================================
%  RESPONSE-SURFACE PLOT
% =========================================================================

[Tg,RHg] = meshgrid( ...
    linspace(min(T_levels),max(T_levels),80), ...
    linspace(min(RH_levels),max(RH_levels),80));

Awg = RHg./100;

WU_grid = ...
    b(1) + ...
    b(2).*Tg + ...
    b(3).*Awg + ...
    b(4).*Tg.^2 + ...
    b(5).*Awg.^2 + ...
    b(6).*Tg.*Awg + ...
    b(7).*Tg.^3 + ...
    b(8).*Awg.^3 + ...
    b(9).*Tg.^2.*Awg + ...
    b(10).*Tg.*Awg.^2;


%% =========================================================================
%  MAXIMUM WATER UPTAKE WITHIN EXPERIMENTAL DOMAIN
% =========================================================================

% Highest measured condition based on the mean of the three replicates
[maxMeasuredWU,idxMeasured] = max(WUconditionMean);

maxMeasuredT  = Tcondition(idxMeasured);
maxMeasuredRH = RHcondition(idxMeasured);

fprintf('\n============================================================\n');
fprintf('MAXIMUM MEASURED WATER UPTAKE\n');
fprintf('============================================================\n');
fprintf('Temperature = %.2f C\n',maxMeasuredT);
fprintf('RH          = %.2f %%\n',maxMeasuredRH);
fprintf('Mean uptake = %.4f g/g\n',maxMeasuredWU);

% Locate the highest point on the RSM grid, then refine the temperature
% continuously at that RH using the fitted cubic model.
[maxGridWU,idxGridMax] = max(WU_grid(:));

RHmaxModelGrid = RHg(idxGridMax);

objectiveT = @(Temp) -WU_cubic(Temp,RHmaxModelGrid);

[TmaxModel,negativeMaxWU] = fminbnd( ...
    objectiveT, ...
    min(T_levels), ...
    max(T_levels));

RHmaxModel = RHmaxModelGrid;
maxModelWU = -negativeMaxWU;

fprintf('\n============================================================\n');
fprintf('MODEL-PREDICTED MAXIMUM WATER UPTAKE\n');
fprintf('============================================================\n');
fprintf('Temperature = %.4f C\n',TmaxModel);
fprintf('RH          = %.2f %%\n',RHmaxModel);
fprintf('Predicted uptake = %.4f g/g\n',maxModelWU);

figure('Color','w');

surf(Tg,RHg,WU_grid,'EdgeColor','none');

hold on;

scatter3( ...
    Tcondition, ...
    RHcondition, ...
    WUconditionMean, ...
    45, ...
    'k', ...
    'filled');

xlabel('Temperature (°C)','FontWeight','bold');
ylabel('Relative Humidity (%)','FontWeight','bold');
zlabel('Water Uptake (g/g)','FontWeight','bold');

colorbar;
view(120,30);

grid on;
box on;

set(gca,'FontWeight','bold','FontSize',14);
