function mse = f_variogram_mse(x, exp_vario_hs, exp_vario_r, SCFtype)
% x is [sigma2, theta, p, nugget], not input

sigma2 = x(1,1);
theta = x(1,2);
if size(x,2) > 2, p = x(1,3); else p = 1.99; end
if size(x,2) > 3, nugget = x(1,4); else nugget = 1e-12; end

sigma2mat = repmat(sigma2,[size(exp_vario_r,1),1]);
SCFmat = f_SCF(exp_vario_hs(:,1), theta, p, SCFtype);
the_vario_r(:,1) = sigma2mat.*(1-(1-nugget)*SCFmat);

se(:,1) = (exp_vario_r - the_vario_r).^2;
mse = mean(se);