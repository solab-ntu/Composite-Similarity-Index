function [degs] = degree360(equi_x, equi_y)
difx = diff(equi_x);
dify = diff(equi_y);

degs = zeros(1,length(difx));
for i = 1:length(difx)
    if (difx(i) >= 0) && (dify(i) > 0)
        degs(i) = atan(dify(i)/difx(i));
    elseif difx(i) < 0
        degs(i) = atan(dify(i)/difx(i)) + pi;
    elseif (difx(i) >= 0) && (dify(i) < 0)
        degs(i) = atan(dify(i)/difx(i)) + 2*pi;
    else
        degs(i) = 0;
    end
end
degs = rad2deg(degs');

end