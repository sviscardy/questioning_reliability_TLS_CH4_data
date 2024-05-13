%% Load full TLS data
if ( sol_index <  2442 )
    SS_MSL_full_data_Webster_2015
else
    SS_MSL_full_data_Webster_2021
end

%% Statistical analysis of TLS data (3 lines + mean)
SS_stat_TLS_data

%% Determine CH4 abundance in HCell (TLS_eta0) and atmosphere (TLS_sig0) + errors

% 1) CH4 abundance in HCell [ppbv] + standard error [ppbv]
TLS_eta0 = F_lines_mean(4) - E_lines_mean(4);
TLS_sig0 = sqrt( (F_lines_std(4)^2/nFpts + E_lines_std(4)^2/nEpts) );

% 2) CH4 abundance in the Martian atmosphere [ppbv] + standard error [ppbv]
TLS_eta  = TLS_eta0/ER_mu;
TLS_sig  = TLS_eta0/ER_mu*(TLS_sig0/TLS_eta0 + ER_sig/ER_mu);
