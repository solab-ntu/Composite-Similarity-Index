function r = f_SCF(hs, theta, p, SCFtype)
% Spatial Correlation Coefficient (w/o sigma2)

theta = repmat(theta,[size(hs,1),1,size(hs,3)]);
p = repmat(p,[size(hs,1),1,size(hs,3)]);

switch SCFtype
    case 1
        % -----------------------------------------------------
        % General exponential model (Gaussian model when p = 2)
        % -----------------------------------------------------
        r(:,1,:) = exp(-theta.*hs.^p);
        
    case 2
        % -------------------
        % Matern linear model
        % -------------------
        r(:,1,:) = (1+theta.*hs).*exp(-theta.*hs);
        
    case 3
        % ------------------
        % Matern cubic model
        % ------------------
        r(:,1,:) = (1+theta.*hs+(theta.^2).*(hs.^2)).*exp(-(theta.*hs));
        
    otherwise
end