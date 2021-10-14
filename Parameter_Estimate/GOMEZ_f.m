function [f] = GOMEZ_f(x)
% Gomez #3
% P_def.Problem_eq = @GOMEZ;
% P_def.lb = [-1,-1];
% P_def.ub = [1,1];
% xopt = [0.10925714458181  -0.62344776471809]
% fopt = -0.9711 (g is active)

if size(x,2)~=2, x=x'; end
if size(x,2)~=2, error('not a 2-D input vector!'), end

f = (4-2.1*(x(:,1).^2)+(x(:,1).^4)./3).*(x(:,1).^2) + ...
   	x(:,1).*x(:,2) + (-4+4*(x(:,2).^2)).*(x(:,2).^2);