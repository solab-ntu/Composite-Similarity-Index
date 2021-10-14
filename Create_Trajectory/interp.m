function y_interp = interp(x,y,n_samples)
x_interp = linspace(x(1),x(end),n_samples);
y_interp = interp1(x,y,x_interp,'spline');
end