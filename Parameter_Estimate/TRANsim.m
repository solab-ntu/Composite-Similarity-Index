function [sim_result] = TRANsim(x,ideal_eq, real_eq_x, real_eq_y,th_TRAN)
angle = deg2rad(x(4));
temp_ideal_eq_x = x(3)*(cos(angle)*ideal_eq(:,1) - sin(angle)*ideal_eq(:,2)) + x(1);
temp_ideal_eq_y = x(3)*(sin(angle)*ideal_eq(:,1) + cos(angle)*ideal_eq(:,2)) + x(2);

score = 0;
deltas = sqrt((temp_ideal_eq_x - real_eq_x).^2+(temp_ideal_eq_y - real_eq_y).^2);
score = 1./(1+((deltas./th_TRAN).^(2*3)));
sim_result = (sum(score)/length(score))*100;

end