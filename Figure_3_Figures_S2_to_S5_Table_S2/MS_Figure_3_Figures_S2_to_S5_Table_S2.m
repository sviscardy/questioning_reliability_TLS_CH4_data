% ==============================
% Figures 3 and S2-S5 + Table S2
% ==============================
%
% This script makes a figure showing:
% - (left column) the time series of data points (CH4 vmr) from each
% triplet line ('e', 'f', and 'g' lines) and their weighted average (Wefg);
% - (right column) the histogram of data points + normal distribution
%
% The CH4 abundances (inferred to be in the Martian atmosphere) above each
% panel (right column, first three rows) are obtained considering the
% triplet lines separately.
%
% The abundances derived from the three triplet lines are reported
% exclusively in Webster et al. (2021). Therefore, only the data from five
% experiments are considered here, corresponding to Sols 2442, 2446, 2615,
% 2627, and 2644.
%
% Webster et al. (2021), Day-night differences in Mars methane suggest
% nighttime containment at Gale crater, A&A 650, A166,
% doi:10.1051/0004-6361/202040030
%
%   Elapsed time: 55 seconds
% 
% Author: sebastien.viscardy@aeronomie.be
% 
%% Link to functions
addpath('../functions/');

%%
clearvars

tic

%% Type of experiment ('D': direct-ingest' / 'E': enrichment)
t_exp_list = {'E'}; n_exp = length(t_exp_list);
E_sol_list = [2442 2446 2615 2627 2644];
lab_figs   = {'S02' 'S03' 'S04' '03'  'S05'};

%% Sols selected
t_exp    = 'E';
sol_list = E_sol_list;
nsol     = length(sol_list);

%% Infos on figure
figid     = 3;
sfigtype  = 'png';
ftsz      = 14;
ftnm      = 'times';
savefig   = 1; % = 1: save Figure
savedata  = 1; % = 1: save Data (Table)
collistF  = [0.8 0 0.8];
collistE  = [0.2 0.8 0.6];
mksz      = 6;
lab_data  = {'e line' 'f line' 'g line' 'Wefg  '};
lab_line  = {'e' 'f' 'g'};
lab_panel = {'(a)' '(b)' '(c)' '(d)' '(e)' '(f)' '(g)' '(h)'};

%% Define data to be stored (==> Table S2)
sol_idx = zeros(1,nsol);
eta_e_F_list = zeros(1,nsol);
sig_e_F_list = zeros(1,nsol);
eta_e_E_list = zeros(1,nsol);
sig_e_E_list = zeros(1,nsol);
eta_f_F_list = zeros(1,nsol);
sig_f_F_list = zeros(1,nsol);
eta_f_E_list = zeros(1,nsol);
sig_f_E_list = zeros(1,nsol);
eta_g_F_list = zeros(1,nsol);
sig_g_F_list = zeros(1,nsol);
eta_g_E_list = zeros(1,nsol);
sig_g_E_list = zeros(1,nsol);
eta_H_list   = zeros(1,nsol);
sig_H_list   = zeros(1,nsol);

%% Loop (experiments)
for isol = 1:nsol
    sol_index   = sol_list(isol);
    
    %% Load full data
    if ( sol_index <  2442 )
        SS_MSL_full_data_Webster_2015
    else
        SS_MSL_full_data_Webster_2021
    end
    
    %% Statistical analysis of data (3 lines + mean)
    SS_stat_data
    
    %% Calculation of eta and sigma
    SS_TLS_CH4_eta_sig
    
    disp(['Sol ',num2str(sol_index),' : ',num2str(eta,'%2.2f'),' +/- ', ...
        num2str(sig,'%2.2f'),' ppbv'])
    
    %% Xlim and Ylim for time series (left) and histograms (right)
    SS_axes_ranges
    
    %% Make figure
    figure(figid)
    delete(get(gcf,'children'))
    set(gcf,'paperpositionmode','auto')
    set(gcf,'position',[200 -600 1200 1400])
    set(gca,'fontsize',ftsz,'fontname',ftnm)
    nrow = 4; ncol = 2;
    
    %% 1) time series
    for i_ts = 1:4 % (3 triplet lines + Wefg)
        subplot(nrow,ncol,2*i_ts-1)
        ileg = 1; clear legname
        
        yline(0,'--k','handlevisibility','off'); hold on
        
        t_run = 'full';
        SS_plot_ts_data
        
        t_run = 'empty';
        SS_plot_ts_data
        
        SS_mk_label_panel(lab_panel{2*i_ts-1},xlim,'lin',ylim,'lin',ftsz,ftnm)
    end
    
    %% 2) histograms
    maxcount = NaN;
    for i_hist = 1:4 % (3 triplet lines + Wefg)
        subplot(nrow,ncol,2*i_hist)
        ileg = 1; clear legname
        
        t_run = 'full';
        SS_plot_histogram
        SS_plot_normal_distr
        
        t_run = 'empty';
        SS_plot_histogram
        SS_plot_normal_distr
        
        hl = legend(legname,'location','northwest', ...
            'fontsize',ftsz,'fontname',ftnm,'interpreter','latex');
        
        eta_H_line    = F_lines_mean(i_hist) - E_lines_mean(i_hist);
        sig_H_line    = sqrt( (E_lines_std(i_hist).^2)/nEpts + ...
            (F_lines_std(i_hist).^2)/nFpts );
        disp([lab_data{i_hist},':   eta_H  +/- sig_H = ',num2str(eta_H_line,'%2.2f'), ...
            ' +/- ',num2str(sig_H_line,'%2.2f'),' ppbv'])
        
        % eta, sig: CH4 vmr in Martian atmosphere
        eta = eta_H_line/enr_fct;
        sig  = sqrt( (E_lines_std(i_hist).^2)/nEpts + ...
            (F_lines_std(i_hist).^2)/nFpts )/enr_fct;
        
        if (i_hist <= 3) % triplet line
            var_title = ['\eta_',lab_line{i_hist},' \pm \sigma_', ...
                lab_line{i_hist},' = '];
        elseif (i_hist == 4) % Wefg
            var_title = '\eta \pm \sigma = ';
        end
        
        texttitle = [var_title,num2str(eta,'%2.2f'),' \pm ', ...
            num2str(sig,'%2.2f'),' ppbv'];
        title(texttitle,'fontname',ftnm,'fontsize',ftsz)
    end
    
    % Labels
    for i_hist = 1:4 % (3 triplet lines + Wefg)
        subplot(nrow,ncol,2*i_hist)
        xlim(xl_distr)
        ylim([0 maxcount+4])
        if (sol_index == 2442), ylim([0 9]), end
        SS_mk_label_panel(lab_panel{2*i_hist},xlim,'lin',ylim,'lin',ftsz,ftnm)
    end
    
    %% Save figure
    if (savefig == 1)
        sfpath     = '../Figures/';
        sffilename = ['Figure_',lab_figs{isol},'.png'];
        sfflnm     = fullfile(sfpath,sffilename);
        exportgraphics(gcf,sfflnm,'Resolution',400)
    end
    
    %% Store data ==> Table S2
    sol_idx(isol) = sol_index;
    
    % 'e' line
    eta_e_F_list(isol) = F_lines_mean(1); sig_e_F_list(isol) = F_lines_std(1); % full-cell runs
    eta_e_E_list(isol) = E_lines_mean(1); sig_e_E_list(isol) = E_lines_std(1);
    
    % 'f' line
    eta_f_F_list(isol) = F_lines_mean(2); sig_f_F_list(isol) = F_lines_std(2); % full-cell runs
    eta_f_E_list(isol) = E_lines_mean(2); sig_f_E_list(isol) = E_lines_std(2);
    
    % 'g' line
    eta_g_F_list(isol) = F_lines_mean(3); sig_g_F_list(isol) = F_lines_std(3); % full-cell runs
    eta_g_E_list(isol) = E_lines_mean(3); sig_g_E_list(isol) = E_lines_std(3);
    
    % Wefg
    eta_H_list(isol)   = eta_H_line;      sig_H_list(isol)   = sig_H_line; % CH4 vmr (Herriott cell)
end

%% Save data
if (savedata == 1)
    T = table(sol_idx', ...
        eta_e_F_list',sig_e_F_list', ... % full-cell runs: 'e' line
        eta_f_F_list',sig_f_F_list', ... % full-cell runs: 'f' line
        eta_g_F_list',sig_g_F_list', ... % full-cell runs: 'g' line
        eta_e_E_list',sig_e_E_list', ... % empty-cell runs: 'e' line
        eta_f_E_list',sig_f_E_list', ... % empty-cell runs: 'f' line
        eta_g_E_list',sig_g_E_list', ... % empty-cell runs: 'g' line
        eta_H_list',sig_H_list');        % Wefg (Herriott cell)
    sdpath     = '../output_tables/';
    sdfilename = 'table_S2.xlsx';
    sdflnm     = fullfile(sdpath,sdfilename);
    writetable(T,sdflnm,'Sheet',1,'Range','A1')
end

toc