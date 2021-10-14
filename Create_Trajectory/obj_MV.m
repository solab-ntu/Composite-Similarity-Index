%% Objective function: Maximum Variance
function [sig] = obj_MV(x)
global kparam
    [y, sig] = f_predictkrige(x, kparam); % kriging
 sig = -sig;
end