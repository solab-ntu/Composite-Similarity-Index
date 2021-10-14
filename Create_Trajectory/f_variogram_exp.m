function [exp_data, exp_vario_hs, exp_vario_r] = f_variogram_exp(x, y, lb, ub)
% caution: only deal with one dimension in y

% ------------------------------
% rescale to input data to [0,1]
% ------------------------------
[n, d] = size(x);
xs = (x-ones(n,1)*lb)./(ones(n,1)*(ub-lb));

% ----------------------
% variogram data produce
% ----------------------
data = [];
for i = 1:n-1
    shs = 0; % square hs
    for j = 1:d, shs = shs + (xs(i,j)-xs(i+1:size(xs,1),j)).^2; end
    % hs, sample 1 number, sample 2 number, r(rh)
    data = [data; [sqrt(shs), ones(size(xs,1)-i,1)*i, ((i+1):size(xs,1))', 0.5*(y(i,1) - y(i+1:size(xs,1),1)).^2] ];
end
cmin = min(data(:,1));
cmax = max(data(:,1))/2;
div_reg = 20;

% -------
% set lag
% -------
for i = 1:div_reg+1
    region_index(i,1) = cmin + (i-1)*(cmax - cmin)/div_reg;
end

% --------
% lag zone
% --------
region_zone(1,1:div_reg) = region_index(1:div_reg,1);
region_zone(2,1:div_reg) = region_index(2:div_reg+1,1);

region_value = zeros(div_reg,2);
for j = 1:div_reg
    mask = double(data(:,1) >= region_zone(1,j)) & (data(:,1) <= region_zone(2,j));
    region_value(j,1) = sum(mask.*data(:,4)); % sum of r(h)
    region_value(j,2) = sum(mask); % sum of sample number
end

exp_data = data;
exp_vario_r(1:div_reg,1) = region_value(:,1)./region_value(:,2);
exp_vario_hs(1:div_reg,1) = 0.5*sum(region_zone,1);