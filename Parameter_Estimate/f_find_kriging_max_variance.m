function [x_next, max_variance] = f_find_kriging_max_variance(kparam)

options.conTol = 1e-3;
options.display = 1;
fileInfo.fName = 'f_get_kriging_variance';
fileInfo.fParams = {kparam};

result = UMDIRECT(fileInfo, kparam.lb, kparam.ub, options);

x_next = result.xBest;
max_variance = result.fBest;
