clearvars

tic

disp(' ')
disp('Simulation of experiments')
disp('=========================')
disp(' ')

%% Main parameters
nsig = 2;
wl   = [1 1 2]; % weight of each triplet line ('g' line twice stronger)
nSIM = 1e5;     % number of simulated instances

%% Type of experiment
t_exp      = 'E';
E_sol_list = [2442 2446 2615 2627 2644];
% E_sol_list = 2627;
sol_list   = E_sol_list; nsol = length(sol_list);

%% p: probability of finding 3 consistent results (from each of the 3 triplet lines)
p_list = zeros(1,nsol);

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
    
    %     F_Wefg_CH4  = (wl(1)*F_e_line_CH4 + wl(2)*F_f_line_CH4 + wl(3)*F_g_line_CH4)/sum(wl);
    %     E_Wefg_CH4  = (wl(1)*E_e_line_CH4 + wl(2)*E_f_line_CH4 + wl(3)*E_g_line_CH4)/sum(wl);
    
    F_Wefg_mean = mean(F_Wefg_CH4);
    E_Wefg_mean = mean(E_Wefg_CH4);
    
    F_Wefg_std  = std(F_Wefg_CH4);
    E_Wefg_std  = std(E_Wefg_CH4);
    
    eta         = (F_Wefg_mean - E_Wefg_mean)/enr_fct;
    sig         = (sqrt((F_Wefg_std)^2/nFpts + (E_Wefg_std)^2/nEpts))/enr_fct;
    
    disp(['eta +/- sig = ',num2str(eta,'%2.2f'),' +/- ',num2str(nsig*sig,'%2.2f'),' ppbv'])
    
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
    %       eta_star_e_E(k=1, ..., 26) ~ N(eta_bar_E, E_std(1)^2)
    SIM_eta_star_e_E = eta_bar_E         + E_std(1)*randn(nSIM,nFpts);
    SIM_eta_star_f_E = eta_bar_E         + E_std(2)*randn(nSIM,nFpts);
    SIM_eta_star_g_E = eta_bar_E         + E_std(3)*randn(nSIM,nFpts);
    
    SIM_eta_bar_e_F  = mean(SIM_eta_star_e_F,2);
    SIM_eta_bar_f_F  = mean(SIM_eta_star_f_F,2);
    SIM_eta_bar_g_F  = mean(SIM_eta_star_g_F,2);
    
    SIM_eta_bar_e_E  = mean(SIM_eta_star_e_E,2);
    SIM_eta_bar_f_E  = mean(SIM_eta_star_f_E,2);
    SIM_eta_bar_g_E  = mean(SIM_eta_star_g_E,2);
    
    SIM_sig_bar_e_F  = std(SIM_eta_star_e_F,0,2);
    SIM_sig_bar_f_F  = std(SIM_eta_star_f_F,0,2);
    SIM_sig_bar_g_F  = std(SIM_eta_star_g_F,0,2);
    
    SIM_sig_bar_e_E  = std(SIM_eta_star_e_E,0,2);
    SIM_sig_bar_f_E  = std(SIM_eta_star_f_E,0,2);
    SIM_sig_bar_g_E  = std(SIM_eta_star_g_E,0,2);
    
    SIM_sig_bar_lines_F = [SIM_sig_bar_e_F SIM_sig_bar_f_F SIM_sig_bar_g_F];
    SIM_sig_bar_lines_E = [SIM_sig_bar_e_E SIM_sig_bar_f_E SIM_sig_bar_g_E];
    
    SIM_eta_line = [SIM_eta_bar_e_F SIM_eta_bar_f_F SIM_eta_bar_g_F] - ...
        [SIM_eta_bar_e_E SIM_eta_bar_f_E SIM_eta_bar_g_E];
    
    
    SIM_sig_line = zeros(nSIM,3);
    for j_line = 1:3
        SIM_sig_line(:,j_line) = sqrt( ...
            SIM_sig_bar_lines_F(:,j_line).^2/nFpts + ...
            SIM_sig_bar_lines_E(:,j_line).^2/nEpts ...
            );
    end
    
    SIM_maxval   = max(SIM_eta_line - nsig*SIM_sig_line,[],2);
    SIM_minval   = min(SIM_eta_line + nsig*SIM_sig_line,[],2);
    
    
    %% SIM_isconsistent = 1: 3 lines consistent / = 0: inconsistent
    SIM_isconsistent                = zeros(nSIM,1);
    SIM_isconsistent(SIM_maxval<SIM_minval) = 1;
    
    nb_consistent = sum(SIM_isconsistent);
    
    p             = nb_consistent/nSIM;
    p_list(isol)  = p;
    
    disp(['  p = ',num2str(p*100),'%'])
    
    
    %% Estimate of probability of finding k consistency cases among N experiments
    
    if (nsig == 1)
        k = 0;
    elseif (nsig == 2)
        k = 3;
    end
    
    N = 5; % Webster et al. (2021)
    
    Prob = factorial(N)/(factorial(k)*factorial(N-k)) * ...
        p^k * (1-p)^(N-k);
    disp(['Prob(k=',num2str(k),') = ',num2str(Prob)])
    
    Prob_list(isol) = Prob;
    
    
end

meanProb = mean(Prob_list);
disp(' ')
disp(['mean Prob = ',num2str(meanProb)])
disp(' ')

%% mean probability pbar = mean(p_list)
pbar = mean(p_list);

if (nsig == 1)
    k = 0;
elseif (nsig == 2)
    k = 3;
end

N = nsol; % 5 experiments reported in Webster et al. (2021)

Probbar = factorial(N)/(factorial(k)*factorial(N-k)) * ...
    pbar^k * (1-pbar)^(N-k);

disp(['Probbar(k=',num2str(k),') = ',num2str(Probbar)])


