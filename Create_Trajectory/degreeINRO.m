function [degs] = degreeINRO(degs)
for i = 1:length(degs)-1
    if degs(i+1) - degs(i) < -180
        degs(i+1) = degs(i+1) + 360;
    elseif degs(i+1) - degs(i) > 180
        degs(i+1) = degs(i+1) - 360;
    else
        degs(i+1) = degs(i+1);
    end
end

end