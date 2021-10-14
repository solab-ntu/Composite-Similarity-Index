function obj = f_get_kriging_variance(x, kparam)
[~, yvar] = f_predictkrige(x, kparam);
obj = - yvar;