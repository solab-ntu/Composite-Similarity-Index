function [ynew, yvar, lambda] = f_predictkrige(xnew, kparam)

% function [ynew, yvar, lambda] = f_predictkrige(xnew, kparam)
%
% This function file uses the kparam values from fitkrige.m to predict the value
% of the response at several locations using a kriging model.
%
% -------
% INPUTS:
% -------
% xnew - [nn x d] matrix containing the input values to predict
%
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
% --------
% OUTPUTS:
% --------
% ynew - the [nn x m] matrix of predicted responses where the different
%        responses are separated by columns.
%
% yvar - the [nn x m] matrix of variances in the predicted responses.
%
%
% NCKU, System Optimization Lab.,
% Coded by Y.-C. Huang
% last updated 2009/08/26

% --------------------------------------
% pull parameters out of structure array
% --------------------------------------
x = kparam.inputs;
y = kparam.outputs;
theta = kparam.theta;
p = kparam.p;
nugget = kparam.nugget;
lb = kparam.lb;
ub = kparam.ub;
beta = kparam.beta;
sigma2 = kparam.sigma2;
K = kparam.K;
fittype = kparam.fittype;
SCFtype = kparam.SCFtype;

% ------------------------------
% determine the size of data set
% ------------------------------
[n,d] = size(x); 
[ny,m] = size(y);

% -------------------------------
% compute some stuff needed later
% -------------------------------
o = ones(n,1);
a = o*lb; b=o*ub;
xs = (x-a)./(b-a); % rescale sample data to [0,1]
nn = size(xnew,1);
xnews = (xnew-ones(nn,1)*lb)./ ...
    (ones(nn,1)*ub-ones(nn,1)*lb); % rescale xnew to [0,1] also

% --------------------------
% initial ynew, yvar, lambda
% --------------------------
ynew = ones(nn,m);
if nargout > 1, yvar = ones(nn,m); end
if nargout > 2, lambda = ones(nn,n,m); end

% -------------------------------------------------------------
% begin looping through all response models and untried samples
% -------------------------------------------------------------
for met = 1:m % number of response models
%     Kinv = inv(K(:,:,met)); % kriging weights K = [[R,F]; [F',0]]
    xnewsmat = shiftdim(xnews',-1); % nn,n,1 -> 1,n,nn
    hsmat = sqrt(sum((repmat(xs,[1,1,nn])-repmat(xnewsmat,[n,1])).^2,2));
    rxmat_tmp = (1-nugget(met))*f_SCF(hsmat,theta(met,:),p(met,:),SCFtype);
    if nn == 1, rxmat = rxmat_tmp;
    else rxmat = shiftdim(rxmat_tmp,2)'; end % n,1,nn -> n,nn,1
%     lamKT = [rxmat', ones(size(rxmat,2),1)]*Kinv;
    lamKT = [rxmat', ones(size(rxmat,2),1)]/K(:,:,met);
    ynew(:,met) = lamKT(:,1:n)*y(:,met); % mean estimate
    if (nargout > 1)
        temp = 1 - sum(lamKT.*[rxmat; ones(1,size(rxmat,2))]',2);
        yvar(:,met) = max(eps, sigma2(met)*temp); % variance estimate
        if (nargout > 2), lambda = lamKT; end
    end
end