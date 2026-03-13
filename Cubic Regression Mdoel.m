% Cubic Response Surface Model
clear; clc; close all;

% Raw data
T = [17 17 17 17 17 17 17 23 23 23 23 23 23 23 30 30 30 30 30 30 30 45 45 45 45 45 45 45 50 50 50 50 50 50 50]';   % Temperature (°C)
RH = [0 25 50 70 10 40 60 0 25 50 70 10 40 60 0 25 50 70 10 40 60 0 25 50 70 10 40 60 0 25 50 70 10 40 60]';             % Relative Humidity (% RH)
WU = [0 0.24 0.41 0.52 0.12 0.35 0.49 0 0.33 0.9 1.24 0.15 0.45 0.85 0 0.38 1.01 1.33 0.2 0.6 1.15 0 0.26 0.597 0.78 0.1 0.32 0.73 0 0.19 0.31 0.57 0.05 0.26 0.55]';   % Water Uptake (g/g)

% RH to water activity
Aw = RH / 100;  

% Design matrix
X = [ones(size(T)) , T , Aw , T.^2 , Aw.^2 , T.*Aw , T.^3 , Aw.^3 , T.^2.*Aw , T.*Aw.^2];
termNames = {'Intercept','T','aw','T^2','aw^2','T·aw','T^3','aw^3','T^2·aw','T·aw^2'};

% Solve OLS
b = X \ WU;
y_pred = X * b;
resid = WU - y_pred;

% Goodness of fit
SSres = sum(resid.^2);
SStot = sum((WU - mean(WU)).^2);
R2 = 1 - SSres/SStot;
n = numel(WU); p = numel(b);
AdjR2 = 1 - (1-R2)*(n-1)/(n-p-1);
RMSE = sqrt(mean(resid.^2));

% Print results
fprintf('\nCubic RSM Fit for Water Uptake Modelling Parameters\n');
for i = 1:numel(b)
    fprintf('%-8s = %+0.8f\n', termNames{i}, b(i));
end
fprintf('\nR² = %.4f | Adj R² = %.4f | RMSE = %.4f g/g\n', R2, AdjR2, RMSE);

% Plotting
[Tg, RHg] = meshgrid(linspace(min(T),max(T),70), linspace(min(RH),max(RH),70));
awg = RHg/100;

WU_pred = b(1) + b(2).*Tg + b(3).*awg + b(4).*Tg.^2 + b(5).*awg.^2 + ...
          b(6).*Tg.*awg + b(7).*Tg.^3 + b(8).*awg.^3 + ...
          b(9).*Tg.^2.*awg + b(10).*Tg.*awg.^2;

figure('Color','w');
surf(Tg, RHg, WU_pred, 'EdgeColor','none');
xlabel('Temperature (°C)','FontWeight','bold');
ylabel('Relative Humidity (%)','FontWeight','bold');
zlabel('Water Uptake (g/g)','FontWeight','bold');
colormap(jet); colorbar;
view(120,30);
hold on;
%scatter3(T, RH, WU, 70, 'k', 'filled', 'MarkerEdgeColor','w', 'LineWidth',1.2);
grid on; box on;
set(gca,'FontWeight','bold','FontSize',14);

% Test
fprintf('\nTest prediction: WU(30°C, 50%% RH) = %.3f g/g\n', WU_cubic(30,50));
