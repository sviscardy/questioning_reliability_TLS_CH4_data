% Probability 'Prob' of obtaining 'k successes' (i.e. 3
% consistent CH4 abudances) in 'N' (= 5, i.e., # TLS experiments)
% independent trials, where 'p' is the probability of success on a single
% trial:
% 
%   Prob(k,N,p) = C(N,k) p^k (1-p)^(N-k)
%
% Running 'MS_nb_exp_consist_3_lines.m' gives:
% 
%   1-sigma: k = 0;
%   2-sigma: k = 3;
% 
% This script generates 'nSIM' = 100,000 instances. For each instance, 3
% sets of random data points (one set for each line) for the empty-cell
% runs and 3 other sets for the full-cell runs. The probability 'p' is
% estimated by analyzing the full sample of nSIM instances. For further
% details, see the manuscript and supplementary information.
%
%   Elapsed time: 2 seconds
% 
% Author: sebastien.viscardy@aeronomie.be
% 
%% Link to functions
addpath('../functions/');

%%
clearvars

tic

disp(' ')
disp('Simulation of instances')
disp('=======================')
disp(' ')

%% Main parameters
nsig_list = [1 2];             % # sigma
nnsig     = length(nsig_list); 
wl        = [1 1 2];           % weight of each triplet line ('g' line twice stronger)
nSIM      = 1e5;               % number of instances

disp(['number of instances = ',num2str(nSIM)])
disp(' ')

%% Type of experiment
t_exp      = 'E';
E_sol_list = [2442 2446 2615 2627 2644];
sol_list   = E_sol_list; nsol = length(sol_list);

%% Loop (# sigma)
for insig = 1:nnsig
    nsig = nsig_list(insig);
    disp(' ')
    disp('-------')
    disp([num2str(nsig),'-sigma'])
    disp('-------')
    disp(' ')
    
    %% isconsist (one value per Sol): 0 = no consistency / 1 = consistency
    isconsist = zeros(1,nsol);
    
    %% p: probability of 'succeeding' (i.e., 3 consistent values from each triplet line)
    p_list    = zeros(1,nsol);
    
    %% P(k,N): Probability of finding k 'good' instances among N instances
    Prob_list = zeros(1,nsol);
    
    %% Loop (experiments)
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
        F_3_lines   = [F_e_line_CH4 F_f_line_CH4 F_g_line_CH4]; % full-cell runs (3 triplet lines)
        E_3_lines   = [E_e_line_CH4 E_f_line_CH4 E_g_line_CH4]; % empty-cell runs (3 triplet lines)
        
        F_Wefg_CH4  = sum(wl.*F_3_lines,2)/sum(wl);             % full-cell runs (Wefg)
        E_Wefg_CH4  = sum(wl.*E_3_lines,2)/sum(wl);             % empty-cell runs (Wefg)
        
        F_Wefg_mean = mean(F_Wefg_CH4);
        E_Wefg_mean = mean(E_Wefg_CH4);
        
        F_Wefg_std  = std(F_Wefg_CH4);
        E_Wefg_std  = std(E_Wefg_CH4);
        
        eta         = (F_Wefg_mean - E_Wefg_mean)/enr_fct;
        sig         = (sqrt((F_Wefg_std)^2/nFpts + (E_Wefg_std)^2/nEpts))/enr_fct;
        
        disp(['Sol ',num2str(sol_index),': eta +/- sig = ', ...
            num2str(eta,'%2.2f'),' +/- ',num2str(nsig*sig,'%2.2f'),' ppbv'])
        
        eta_H       = F_Wefg_mean - E_Wefg_mean; % CH4 vmr in Herriott cell
        eta_bar_E   = E_Wefg_mean;               % mean CH4 vmr (empty-cell runs)
        
        %% Generating data points randomly: Normal distribution N(mu, sigma^2)
        % 'SIM_' stands for 'simulated'
        % ex.: CH4 vmr ('e' line, full-cell runs):
        %       SIM_eta_star_e_F(k=1, ..., 26) ~ N(eta_bar_E+eta_H, F_std(1)^2)
        SIM_eta_star_e_F = eta_bar_E + eta_H + F_std(1)*randn(nSIM,nFpts);
        SIM_eta_star_f_F = eta_bar_E + eta_H + F_std(2)*randn(nSIM,nFpts);
        SIM_eta_star_g_F = eta_bar_E + eta_H + F_std(3)*randn(nSIM,nFpts);
        
        % ex.: CH4 vmr ('e' line, empty-cell runs):
        %       SIM_eta_star_e_E(k=1, ..., 26) ~ N(eta_bar_E, E_std(1)^2)
        SIM_eta_star_e_E = eta_bar_E         + E_std(1)*randn(nSIM,nFpts);
        SIM_eta_star_f_E = eta_bar_E         + E_std(2)*randn(nSIM,nFpts);
        SIM_eta_star_g_E = eta_bar_E         + E_std(3)*randn(nSIM,nFpts);
        
        %% Estimate means and standard deviations
        % means (full-cell runs)
        SIM_eta_bar_e_F  = mean(SIM_eta_star_e_F,2);
        SIM_eta_bar_f_F  = mean(SIM_eta_star_f_F,2);
        SIM_eta_bar_g_F  = mean(SIM_eta_star_g_F,2);
        
        % means (empty-cell runs)
        SIM_eta_bar_e_E  = mean(SIM_eta_star_e_E,2);
        SIM_eta_bar_f_E  = mean(SIM_eta_star_f_E,2);
        SIM_eta_bar_g_E  = mean(SIM_eta_star_g_E,2);
        
        % standard deviations (full-cell runs)
        SIM_sig_bar_e_F  = std(SIM_eta_star_e_F,0,2);
        SIM_sig_bar_f_F  = std(SIM_eta_star_f_F,0,2);
        SIM_sig_bar_g_F  = std(SIM_eta_star_g_F,0,2);
        
        % standard deviations (empty-cell runs)
        SIM_sig_bar_e_E  = std(SIM_eta_star_e_E,0,2);
        SIM_sig_bar_f_E  = std(SIM_eta_star_f_E,0,2);
        SIM_sig_bar_g_E  = std(SIM_eta_star_g_E,0,2);
        
        %
        SIM_sig_bar_lines_F = [SIM_sig_bar_e_F SIM_sig_bar_f_F SIM_sig_bar_g_F];
        SIM_sig_bar_lines_E = [SIM_sig_bar_e_E SIM_sig_bar_f_E SIM_sig_bar_g_E];
        
        % (mean) CH4 vmr in Herriott cell inferred from the 3 triplet lines
        SIM_eta_line = [SIM_eta_bar_e_F SIM_eta_bar_f_F SIM_eta_bar_g_F] - ...
            [SIM_eta_bar_e_E SIM_eta_bar_f_E SIM_eta_bar_g_E];
        
        % associated errors
        SIM_sig_line = zeros(nSIM,3);
        for j_line = 1:3
            SIM_sig_line(:,j_line) = sqrt( ...
                SIM_sig_bar_lines_F(:,j_line).^2/nFpts + ...
                SIM_sig_bar_lines_E(:,j_line).^2/nEpts ...
                );
        end
        
        %% maxval and minval
        % Determine 'maxval' and 'minval' used to check consistency for
        % each set of 3 simulated CH4 abundances
        SIM_maxval   = max(SIM_eta_line - nsig*SIM_sig_line,[],2);
        SIM_minval   = min(SIM_eta_line + nsig*SIM_sig_line,[],2);
        
        %% SIM_isconsistent = 1: 3 lines consistent / = 0: inconsistent
        SIM_isconsistent                = zeros(nSIM,1);
        SIM_isconsistent(SIM_maxval<SIM_minval) = 1;
        
        nb_consistent = sum(SIM_isconsistent);
        p             = nb_consistent/nSIM;
        p_list(isol)  = p;
        disp(['        p = ',num2str(p*100),'%'])
    end
    
    %% mean probability 'Prob', where pmean = mean(p_list)
    if (nsig == 1), k = 0; elseif (nsig == 2), k = 3; end
    N     = nsol; % 5 experiments reported in Webster et al. (2021)
    pmean = mean(p_list);
    Prob  = factorial(N)/(factorial(k)*factorial(N-k)) * ...
        pmean^k * (1-pmean)^(N-k);
    
    disp(' ')
    disp([num2str(nsig),'-sigma: Prob(k=',num2str(k),', N = ', ...
        num2str(N),', p = ',num2str(pmean,'%2.2f'),') = ',num2str(Prob,'%2.0e')])
    disp(' ')
end

toc