function [sim_result] = SHAPEsim(x,ideal_degINRO, real_degINRO,th_SHAPE)
score = 0;
temp_ideal_degINRO = ideal_degINRO + x;
deltas = abs(temp_ideal_degINRO - real_degINRO);
score = 1./(1+((deltas./th_SHAPE).^(2*3)));
sim_result = (sum(score)/length(score))*100;

end