function kparam = f_variogram_fit(x, y, lb, ub, fittype, SCFtype)

% function kparam = f_variogram_fit(x, y, lb, ub, fittype, SCFtype)
% 
% This function file takes in a data set and creates the kriging metamodels
% for each response.  The resulting information is stored in a structure
% array named kparam.
% 
% It uses the theoretical variogram to fit the experimental variogram not
% the likelihood method.
% 
% -------
% INPUTS:
% -------
% 
% x - [n X d] matrix of sample point input values
% y - [n X m] matrix of sample point output values
% lb - optional [1 x d] vector of lower bounds on the design variables
%      (default=min(x))
% ub - optional [1 x d] vector of upper bounds on the design variables
%      (default=max(x))
% fittype - 
% SCFtype - 
% 
% --------
% OUTPUTS:
% --------
% kparam - the structure array containing the following kriging parameters:
%   inputs - [n x d] matrix of sample point input values 
%   outputs - [n x m] matrix of sample point output values
%   theta - the [m x d] matrix containing the theta values for prediction
%   p - the [m x d] matrix of p values (for the general exponential SCF)
%   nugget - the [m x 1] vector for the nugget effect
%   lb - the [m x d] matrix of lower bounds on the design variables
%   ub - the [m x d] matrix of upper bounds on the design variables
%   beta - the [(?) x m] matrix of global model terms
%   sigma2 - the [1 x m] vector of variance terms used for prediction
%   K - the [n+? x n+? x m] array of expanded correlation matrices
%   mse - the least square error between the. and exp. variogram
%   fittype - how many parameter to be fiting
%   SCFtype - SCF type selected from f_SCF for fitting kriging model
% 
% NCKU, System Optimization Lab.,
% Coded by Y.-C. Huang
% last updated 2009/08/26

% ----------------------
% some property varilble
% ----------------------
n = size(x,1); % the number of input sampled points
m = size(y,2); % the number of output dimensions

%% model fitting
% ---------------------------------------
% fit variogram for each output dimension
% ---------------------------------------
for met = 1:m
    % -------------------------
    % set the fitting parameter
    % -------------------------
    if met == 1, mse = []; sigma2 = []; theta = []; p = []; nugget = []; end
    
    % ---------------------------------
    % check exp_vario_rt has NaN or not
    % ---------------------------------
    clear exp_vario_hst exp_vario_rt exp_vario_hs exp_vario_r
    [exp_data, exp_vario_hst, exp_vario_rt] = f_variogram_exp(x,y(:,met),lb,ub);
    iter = 1;
    for i = 1:size(exp_vario_rt,1)
        if ~isnan(exp_vario_rt(i,1))
            exp_vario_hs(iter,1) = exp_vario_hst(i,1);
            exp_vario_r(iter,1) = exp_vario_rt(i,1);
            iter = iter + 1;
        end
    end

    % ---------------------------------------
    % beginning fitting kriging with UMDIRECT
    % ---------------------------------------
    if nargin <= 4
        fittype = 1; % select how many parameter to be fiting
        SCFtype = 1; % SCF type for fitting
    end
    
    switch fittype
        case 1
            % -------------
            % sigma2, theta
            % -------------
            vlb = [0.01, 0.01];
            vub = [max(exp_vario_r)*2, 1e3];
            
        case 2
            % ----------------
            % sigma2, theta, p
            % ----------------
            vlb = [0.01, 0.01, 0.01];
            vub = [max(exp_vario_r)*2, 1e3, 1.99];
            
        case 3
            % ------------------------
            % sigma2, theta, p, nugget
            % ------------------------
            vlb = [0.01, 0.01, 1.5, 0.01];
            vub = [max(exp_vario_r)*2, 1e3, 1.99, 0.99];
            
        otherwise
    end
    
    % ---------------------------------------------------------------------
    fileInfo.fName = 'f_variogram_mse';
    fileInfo.fParams = {exp_vario_hs, exp_vario_r, SCFtype};
    UMDoptions.display = 1;
    UMDoptions.maxIter = 100;
    UMDoptions.maxfCount = 2e2*size(vlb,2);
    UMDoptions.termType = 3;
    UMDoptions.saveFile = 'noSave';
    UMDoptions.localSearch = 0; 
    outstruct = UMDIRECT(fileInfo, vlb, vub, UMDoptions);

    mse = [mse; outstruct.fBest];
    sigma2 = [sigma2; outstruct.xBest(1,1)];
    theta = [theta; outstruct.xBest(1,2)];
    
    switch fittype
        case 1
            p = [p; 1.99];
            nugget = [nugget; 1e-12];
        case 2
            p = [p; outstruct.xBest(1,3)];
            nugget = [nugget; 1e-12];
        case 3
            p = [p; outstruct.xBest(1,3)];
            nugget = [nugget; outstruct.xBest(1,4)];
        otherwise
    end
end

%% kriging K-matrix
% ---------------------------
% rescale input data to [0,1]
% ---------------------------
xs = (x-ones(n,1)*lb)./(ones(n,1)*(ub-lb));

for met = 1:m
    R = ones(n); % correlation matrix, R
    for i = 1:n-1
        for j = i+1:n
            R(i,j) = ...
                f_SCF(norm(xs(i,:)-xs(j,:)),theta(met,:),p(met,:),SCFtype);
            R(j,i) = R(i,j);
        end
    end
    R = (1-nugget(met))*R + nugget(met)*eye(n);
    F = ones(n,1); % functional form f(x)
    K(:,:,met) = [[R,F]; [F',0]]; % expanded correlation matrix
    beta(:,met) = ((F'/R)*y(:,met))/((F'/R)*F); % y = fx*beta + z
end

kparam.inputs = x;
kparam.outputs = y;
kparam.theta = theta;
kparam.p = p;
kparam.nugget = nugget;
kparam.lb = lb;
kparam.ub = ub;
kparam.beta = beta;
kparam.sigma2 = sigma2;
kparam.K = K;
kparam.mse = mse;
kparam.fittype = fittype;
kparam.SCFtype = SCFtype;