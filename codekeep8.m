%% SVAR with FFR monetary policy shock and Jordanian macro variables
clear; clc;

%% Settings
dataFile      = 'data.xlsx';   
sheetName     = 1;              
maxLagFFR     = 8;              % maximum AR lag for FFR
maxLagVAR     = 8;              % maximum VAR lag considerid
signHorizon   = 4;              % horizons for sign restrictions
irfHorizon    = 20;             % horizons for IRF plots
fevdHorizons  = [1 4 8 20];     % horizons for FEVD

%% Load data from Excel
tbl   = readtable(dataFile,'Sheet',sheetName);
T     = size(tbl,1);

quarterRaw = tbl{:,1};
ffr        = tbl{:,2};
cbj        = tbl{:,3};
claims     = tbl{:,4};
infJ       = tbl{:,5};
gdp        = tbl{:,6};
cons       = tbl{:,7};

%% Extracting FFR monetary policy shock from AR model
% AR(p) for FFR, lag length chosen by AIC, residuals are FFR shocks

ffr = double(ffr);
bestAIC      = Inf;
bestP        = 1;
bestResid    = [];
bestValidIdx = [];

for p = 1:maxLagFFR
    X = lagmatrix(ffr,1:p);
    y = ffr;
    validIdx = all(~isnan(X),2);
    yEst = y(validIdx);
    XEst = [ones(sum(validIdx),1), X(validIdx,:)];
    
    % OLS
    beta  = XEst \ yEst;
    resid = yEst - XEst*beta;
    
    % Gaussian log-likelihood and AIC
    sigma2 = var(resid,1); % ML variance
    k      = size(XEst,2); % number of estimated parameters
    logL   = -0.5*sum(log(2*pi*sigma2) + (resid.^2)/sigma2);
    aic    = -2*logL + 2*k;
    
    if aic < bestAIC
        bestAIC      = aic;
        bestP        = p;
        bestResid    = resid;
        bestValidIdx = validIdx;
    end
end

fprintf('\nSelected AR lag for FFR (by AIC): p = %d\n',bestP); % check

% Full-length FFR shock series
ffrShockFull = NaN(T,1);
ffrShockFull(bestValidIdx) = bestResid;

%% Stationarity tests, transformations for Jordanian variables
% If positive series, log difference, otherwise just first difference

jordanNames = { ...
    'CBJ policy rate', ...
    'Inflation', ...
    'GDP', ...
    'Total consumption', ...
    'Claims (other depository corps)'};

rawJordan   = [cbj, infJ, gdp, cons, claims];   % T x 5
transJordan = NaN(T,5);
transformType = strings(5,1);

for i = 1:5
    y = double(rawJordan(:,i));
    
    % ADF on levels
    [hLevel,pLevel] = adftest(y);
    fprintf('\nADF on %s in levels: h = %d, p = %.4f\n', ...
        jordanNames{i}, hLevel, pLevel);
    
    %
    yTrans = NaN(T,1);
    switch i
        case 2  % Inflation
            yTrans(2:end) = diff(y);
            transformType(i) = "diff";
        case {3,4,5}  % GDP, Total consumption, Claims
            if any(y <= 0)
                error('Series %s has non-positive values; cannot take log-difference.', ...
                    jordanNames{i});
            end
            yTrans(2:end) = diff(log(y));
            transformType(i) = "log-diff";
        otherwise  % i = 1, CBJ policy rate: original rule
            if hLevel == 0
                if all(y > 0)
                    yTrans(2:end) = diff(log(y));
                    transformType(i) = "log-diff";
                else
                    yTrans(2:end) = diff(y);
                    transformType(i) = "diff";
                end
            else
                yTrans = y;
                transformType(i) = "level";
            end
    end
    
    % Robustness check
    yTransNoNaN = yTrans(~isnan(yTrans));
    [hTrans,pTrans] = adftest(yTransNoNaN);
    fprintf('ADF on %s after %s: h = %d, p = %.4f\n', ...
        jordanNames{i}, transformType(i), hTrans, pTrans);
    
    transJordan(:,i) = yTrans;
end

%% VAR final sample

firstShockIdx = find(~isnan(ffrShockFull),1,'first');
if isempty(firstShockIdx)
    error('FFR shock series is empty. Check AR specification or data.');
end

tStart = max(2, firstShockIdx);
idx    = (tStart:T)';
nObs   = numel(idx);

Y = [ ...
    ffrShockFull(idx), ...
    transJordan(idx,1), ...
    transJordan(idx,2), ...
    transJordan(idx,3), ...
    transJordan(idx,4), ...
    transJordan(idx,5)];

varLabels = { ...
    'FFR shock', ...
    'CBJ policy rate', ...
    'Inflation', ...
    'GDP', ...
    'Total consumption', ...
    'Claims'};

nVars = size(Y,2);

%% VAR lag-order selection, estimation
% NOTE!!! VAR in Y_t = (FFR shock, CBJ, Inflation, GDP, Cons, Claims)'

maxLagVAR = min(maxLagVAR, floor(nObs/3));
if maxLagVAR < 1
    error('Not enough observations for VAR. Reduce lag length or adjust sample.');
end

aic = Inf(maxLagVAR,1);
bic = Inf(maxLagVAR,1);
validLag = false(maxLagVAR,1);

for p = 1:maxLagVAR
    Mdl = varm(nVars,p);
    try
        [EstMdlTmp,~,logL] = estimate(Mdl,Y,'Display','off');
        numParams = nVars^2 * p + nVars; % AR coefficients + intercept
        [aic(p),bic(p)] = aicbic(logL,numParams,nObs);
        validLag(p) = true;
    catch ME
        warning('VAR(%d) estimation failed (%s). Skipping this lag.',p,ME.message);
        aic(p) = Inf;
        bic(p) = Inf;
    end
end

if ~any(validLag)
    error('All candidate VAR lag lengths failed. Consider lowering maxLagVAR or checking the data.');
end

[~,pAIC] = min(aic);
[~,pBIC] = min(bic);
pVAR     = pBIC;   % when i chose AIC VAR was not stable????? used BIC instead

fprintf('\nSelected VAR lag length: pVAR = %d (AIC suggests %d)\n',pVAR,pAIC);

Mdl = varm(nVars,pVAR);
[EstMdl,~,~] = estimate(Mdl,Y,'Display','off');

% companion matrix eigenvalues stability check
n    = EstMdl.NumSeries;
pEst = EstMdl.P;
ARcoeff = EstMdl.AR;

if pEst > 0
    companionDim = n * pEst;
    F = zeros(companionDim);
    
    % fill first block row with [A1, A2, ... Ap]
    for k = 1:pEst
        if iscell(ARcoeff)
            Ak = ARcoeff{k};
        else
            Ak = ARcoeff(:,:,k);
        end
        F(1:n,(k-1)*n+1:k*n) = Ak;
    end
    
    % fill sub-identity for lag shifting
    if pEst > 1
        F(n+1:end,1:n*(pEst-1)) = eye(n*(pEst-1));
    end
    
    eigF = eig(F);
    if any(abs(eigF) >= 1)
        warning('Estimated VAR is not stable (some eigenvalues >= 1 in modulus).');
    end
end

% innovations and covariance
U       = infer(EstMdl,Y);
SigmaU  = cov(U,1);           

%% Structural identification with sign restrictions
P          = chol(SigmaU,'lower');
shockIndex = 1;      % FFR monetary policy shock, 1 s.d.
maxTries   = 20000;
accepted   = false;

% Fix random seed (reproducible), kept changing without...
rng(12345);

% Loop
for draw = 1:maxTries
    % Random orthonormal Q (Haar-distributed)
    [Q,~] = qr(randn(nVars));
    if det(Q) < 0
        Q(:,1) = -Q(:,1);
    end
    
    Bcand = P * Q;
    
    % irfs
    IRFtemp = irf_var(EstMdl,Bcand,signHorizon);
    
    % enforce sign restrictions on impact only (h = 0), weird result
    % otherwise...
    rCBJ  = IRFtemp(2,shockIndex,1); % +
    rINF  = IRFtemp(3,shockIndex,1); % -
    rGDP  = IRFtemp(4,shockIndex,1); % -
    rCONS = IRFtemp(5,shockIndex,1); % -
    
    % sign restrictions: strict inequalities on impact
    if (rCBJ <= 0) || (rINF >= 0) || (rGDP >= 0) || (rCONS >= 0)
        continue;
    end
    
    % test if satisfied restrictions
    B        = Bcand;
    accepted = true;
    fprintf('\nAccepted structural rotation after %d draws.\n',draw);
    break;
end

if ~accepted
    error('No rotation satisfied the sign restrictions. Increase maxTries or review restrictions.');
end

%% Impulse response functions (IRFs)

IRF = irf_var(EstMdl,B,irfHorizon);   % nVars x nVars x (irfHorizon+1)
h   = 0:irfHorizon;

% 65% CI
ciLevel   = 0.65;
lowerPct  = 100*(1-ciLevel)/2;
upperPct  = 100*(1+ciLevel)/2;
nCIdraws  = 500;
IRF_draws = zeros(nVars,nVars,irfHorizon+1,nCIdraws);

% RNG seed for CI draws
rng(54321);
nFound   = 0;
nTriesCI = 0;
maxTriesCI = 200000;

while (nFound < nCIdraws) && (nTriesCI < maxTriesCI)
    nTriesCI = nTriesCI + 1;
    [Qci,~] = qr(randn(nVars));
    if det(Qci) < 0
        Qci(:,1) = -Qci(:,1);
    end
    
    Bci    = P * Qci;
    IRFtmp = irf_var(EstMdl,Bci,signHorizon);
    
    % sign restrictions uniform impact
    rCBJ  = IRFtmp(2,shockIndex,1); % +
    rINF  = IRFtmp(3,shockIndex,1); % -
    rGDP  = IRFtmp(4,shockIndex,1); % -
    rCONS = IRFtmp(5,shockIndex,1); % -
    
    if (rCBJ <= 0) || (rINF >= 0) || (rGDP >= 0) || (rCONS >= 0)
        continue;
    end
    
    % when draw accepted: store full IRFs up to irfHorizon, average irf
    % over all accepted draws
    IRFfull = irf_var(EstMdl,Bci,irfHorizon);
    nFound  = nFound + 1;
    IRF_draws(:,:,:,nFound) = IRFfull;
end

% INDIVIDUAL IRFS

if nFound < nCIdraws
    warning('Only %d accepted draws out of %d attempts for IRF confidence intervals.',nFound,nTriesCI);
    IRF_draws = IRF_draws(:,:,:,1:nFound);
end

IRF_low  = zeros(nVars,irfHorizon+1);
IRF_med  = zeros(nVars,irfHorizon+1);
IRF_high = zeros(nVars,irfHorizon+1);

for iVar = 1:nVars
    for t = 1:irfHorizon+1
        draws = squeeze(IRF_draws(iVar,shockIndex,t,:));
        if isempty(draws)
            IRF_low(iVar,t)  = NaN;
            IRF_med(iVar,t)  = NaN;
            IRF_high(iVar,t) = NaN;
        else
            pct = prctile(draws,[lowerPct 50 upperPct]);
            IRF_low(iVar,t)  = pct(1);
            IRF_med(iVar,t)  = pct(2);
            IRF_high(iVar,t) = pct(3);
        end
    end
end

for iVar = 1:nVars
    figure;
    % Grey 65% CI band (if available)
    if ~all(isnan(IRF_low(iVar,:)))
        xBand = [h, fliplr(h)];
        yBand = [IRF_low(iVar,:), fliplr(IRF_high(iVar,:))];
        fill(xBand,yBand,[0.4 0.4 0.4],'EdgeColor','none','FaceAlpha',0.4);
        hold on;
    else
        hold on;
    end
    % Point estimate (median IRF across accepted draws)
    plot(h, IRF_med(iVar,:),'LineWidth',1.5,'Color','b');
    yline(0,'k-');
    hold off;
    xlabel('Horizon (quarters)');
    ylabel('Response');
    title(sprintf('IRF of %s to FFR monetary policy shock',varLabels{iVar}));
end

% combined IRF figur
figure;
rows = 3; cols = 2;
for iVar = 1:nVars
    subplot(rows,cols,iVar);
    if ~all(isnan(IRF_low(iVar,:)))
        xBand = [h, fliplr(h)];
        yBand = [IRF_low(iVar,:), fliplr(IRF_high(iVar,:))];
        fill(xBand,yBand,[0.4 0.4 0.4],'EdgeColor','none','FaceAlpha',0.4);
        hold on;
    else
        hold on;
    end
    plot(h, IRF_med(iVar,:),'LineWidth',1.5,'Color','b');
    yline(0,'k-');
    hold off;
    xlabel('Horizon (quarters)');
    ylabel('Response');
    title(varLabels{iVar});
end

%% Forecast error variance decomposition
fevdHorizons = fevdHorizons(:)';
if any(fevdHorizons > irfHorizon)
    error('All FEVD horizons must be <= irfHorizon.');
end

nH   = numel(fevdHorizons);
FEVD = zeros(nVars,nVars,nH);   % (variable, shock, horizonIndex)

for hIdx = 1:nH
    H = fevdHorizons(hIdx);
    
    var_i_H = zeros(nVars,1);
    contrib = zeros(nVars,nVars);
    
    % k = 0... H-1 corresponds to IRF(:,:,k+1), i.e. 1..H
    for k = 1:H
        Tk = IRF(:,:,k);      % Theta_{k-1}
        var_i_H = var_i_H + sum(Tk.^2,2);
        contrib = contrib + Tk.^2;
    end
    
    for iVar = 1:nVars
        FEVD(iVar,:,hIdx) = contrib(iVar,:) ./ var_i_H(iVar);
    end
end

% print FEVD results
for hIdx = 1:nH
    H = fevdHorizons(hIdx);
    fprintf('\nForecast error variance decomposition at horizon %d quarters:\n',H);
    for iVar = 1:nVars
        fprintf('%s: ',varLabels{iVar});
        for j = 1:nVars
            fprintf('%s=%.3f ',varLabels{j}, FEVD(iVar,j,hIdx));
        end
        fprintf('\n');
    end
end

%% Local function for IRFs (useful excercise for me)
function IRF = irf_var(EstMdl,B,horizon)
% compute IRFs for a VAR model with structural impact matrix
    n = EstMdl.NumSeries;
    p = EstMdl.P;
    ARcoeff = EstMdl.AR;

    IRF = zeros(n,n,horizon+1);
    IRF(:,:,1) = B;    % impact response at horizon 0
    
    if p == 0
        return;
    end
    
    % companion form
    companionDim = n*p;
    F = zeros(companionDim);
    for j = 1:p
        if iscell(ARcoeff)
            Aj = ARcoeff{j};
        else
            Aj = ARcoeff(:,:,j);
        end
        F(1:n,(j-1)*n+1:j*n) = Aj;
    end
    if p > 1
        F(n+1:end,1:n*(p-1)) = eye(n*(p-1));
    end
    J = [eye(n), zeros(n,n*(p-1))];
    
    Fpow = eye(companionDim);
    for k = 2:horizon+1
        Fpow = Fpow * F;
        IRF(:,:,k) = J * Fpow * J' * B;
    end
end
