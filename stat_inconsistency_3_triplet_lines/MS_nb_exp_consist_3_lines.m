% This script calculates the number of experiments where the 3 CH4
% abundances obtained from the 3 triplet lines are consistent at the
% 1-sigma and 2-sigma levels
% 
%   Elapsed time: 0.1 second
% 
% Author: sebastien.viscardy@aeronomie.be
%
%% Link to functions
addpath('../functions/');

%%
clearvars

tic

%% Main parameters
E_sol_list  = [2442 2446 2615 2627 2644];
t_exp       = 'E';
enr_fct     = 25;
wl          = [1 1 2]; % weight of each triplet line ('g' line twice stronger)
lab_consist = {'no consistency' 'consistency'};

%% Selected experiments
sol_list   = E_sol_list;
nsol       = length(sol_list);

for isig = 1:2 % # sigma
    nsig = isig;
    
    disp(' ')
    disp('=======')
    disp([num2str(nsig),'-sigma'])
    disp('=======')
    disp(' ')
    
    %% isconsist = 1: 3 lines consistent / = 0: not consistent
    isconsist = zeros(1,nsol);
    
    %% Loop over run sols
    for isol = 1:nsol
        sol_index       = sol_list(isol);
        isconsist(isol) = 0;
        
        %% Load full data
        if ( sol_index <  2442 )
            SS_MSL_full_data_Webster_2015
        else
            SS_MSL_full_data_Webster_2021
        end
        
        %% means and uncertainties
        nFpts     = length(F_e_line_CH4);
        nEpts     = length(E_e_line_CH4);
        
        F_mean(1) = mean(F_e_line_CH4);
        F_mean(2) = mean(F_f_line_CH4);
        F_mean(3) = mean(F_g_line_CH4);
        
        F_std(1)  = std(F_e_line_CH4);
        F_std(2)  = std(F_f_line_CH4);
        F_std(3)  = std(F_g_line_CH4);
        
        E_mean(1) = mean(E_e_line_CH4);
        E_mean(2) = mean(E_f_line_CH4);
        E_mean(3) = mean(E_g_line_CH4);
        
        E_std(1)  = std(E_e_line_CH4);
        E_std(2)  = std(E_f_line_CH4);
        E_std(3)  = std(E_g_line_CH4);
        
        eta_line  = zeros(1,3);
        sig_line  = zeros(1,3);
        
        for j_line = 1:3
            eta_line(j_line) = F_mean(j_line)-E_mean(j_line);
            sig_line(j_line) = sqrt( F_std(j_line)^2/nFpts + E_std(j_line)^2/nEpts );
        end
        
        eta_line = eta_line/enr_fct;
        sig_line = sig_line/enr_fct;
        
        %% Wefg
        F_Wefg_CH4  = (wl(1)*F_e_line_CH4 + wl(2)*F_f_line_CH4 + wl(3)*F_g_line_CH4)/sum(wl);
        E_Wefg_CH4  = (wl(1)*E_e_line_CH4 + wl(2)*E_f_line_CH4 + wl(3)*E_g_line_CH4)/sum(wl);
        
        F_Wefg_mean = mean(F_Wefg_CH4);
        E_Wefg_mean = mean(E_Wefg_CH4);
        
        F_Wefg_std  = std(F_Wefg_CH4);
        E_Wefg_std  = std(E_Wefg_CH4);
        
        eta         = (F_Wefg_mean - E_Wefg_mean)/enr_fct;
        sig         = (sqrt((F_Wefg_std)^2/nFpts + (E_Wefg_std)^2/nEpts))/enr_fct;
        
        disp(['eta +/- sig = ',num2str(eta,'%2.2f'),' +/- ',num2str(nsig*sig,'%2.2f'),' ppbv'])
        
        %% Check consistency
        maxval   = max(eta_line - nsig*sig_line);
        minval   = min(eta_line + nsig*sig_line);
        
        if (maxval < minval), isconsist(isol) = 1; end
    end
    
    %% print exp / consistency: (0: not consistent / 1: consistent)
    disp(' ')
    for isol = 1:nsol
        disp(['Sol ',num2str(sol_list(isol)),' : ',lab_consist{isconsist(isol)+1}])
    end
    
    fract_consist = sum(isconsist)/nsol*100;
    disp(['fraction of consistent cases = ',num2str(fract_consist),'%'])
    k = sum(isconsist); % number of 'consistent' set of CH4 abundances
    
    disp(['k = ',num2str(k)])
    disp(['N = ',num2str(nsol)])
    
end

toc