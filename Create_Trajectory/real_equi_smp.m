function [real_eq] = real_equi_smp(real_x, real_y, num)
% num:the number of equi sample
% savgol filter
% sg_x =  sgolayfilt(real_x,3,11);% 3,5
% sg_y =  sgolayfilt(real_y,3,11);% 3,5

% equi_smp
% real_eq = interparc(num,sg_x, sg_y,'linear');
real_eq = interparc(num,real_x,real_y,'linear');
end