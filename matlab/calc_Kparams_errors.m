function [err_alpha err_beta err_u err_v err_gamma] = calc_Kparams_errors(K0, K)
err_u = abs(K(1,3)-K0(1,3));
err_v = abs(K(2,3)-K0(2,3));
err_alpha = abs(K(1,1)-K0(1,1))/K0(1,1)*100; % relative error for alpha and beta
err_beta = abs(K(2,2)-K0(2,2))/K0(2,2)*100;
err_gamma = abs(K(1,2)-K0(1,2));