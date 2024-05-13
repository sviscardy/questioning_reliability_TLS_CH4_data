% This subscript calculates the CH4 vmr in the Herriott cell (eta_H) and in
% the Martian atmosphere (eta) using the weighted average values (Wefg).
% 
% First, run:
%   SS_MSL_full_data_Webster_2015
% or
%   SS_MSL_full_data_Webster_2021
% 
% Data taken from:
% 
% Webster et al. (2015), Mars methane detection and variability at Gale
% crater, Science 347, 415-417, doi:10.1126/science.1261713
% 
% Webster et al. (2021), Day-night differences in Mars methane suggest
% nighttime containment at Gale crater, A&A 650, A166,
% doi:10.1051/0004-6361/202040030
%
% Author: sebastien.viscardy@aeronomie.be
%
%% Full-cell runs [ppbv]
F_mu  = mean(F_Wefg_CH4);
F_sig = std(F_Wefg_CH4);

%% Empty-cell runs [ppbv]
E_mu  = mean(E_Wefg_CH4);
E_sig = std(E_Wefg_CH4);

%% eta_H: CH4 vmr in the Herriott cell [ppbv]
eta_H = F_mu-E_mu;
sig_H = sqrt( (F_sig^2)/length(F_Wefg_CH4) + ...
    (E_sig^2)/length(E_Wefg_CH4) );

%% CH4 vmr in Martian atmosphere [ppbv]
eta   = eta_H/enr_fct;
sig   = sig_H/enr_fct;