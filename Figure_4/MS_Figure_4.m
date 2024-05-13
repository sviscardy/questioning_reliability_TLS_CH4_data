% =========
% Figure 4
% =========
% 
% This script computes the CH4 vmr in the Martian atmosphere as inferred
% from each of the three triplet lines (e, f, g) and compared to that from
% the weighted mean (Wefg), which corresponds to the abundance reported in
% the literature (Webster et al., 2021).
% 
% The abundances derived from the three triplet lines are reported
% exclusively in Webster et al. (2021). Therefore, only the data from five
% experiments are considered here, corresponding to Sols 2442, 2446, 2615,
% 2627, and 2644.
% 
% Data taken from:
% 
% Webster et al. (2021), Day-night differences in Mars methane suggest
% nighttime containment at Gale crater, A&A 650, A166,
% doi:10.1051/0004-6361/202040030
% 
%   Elapsed time: 6 seconds
% 
% Author: sebastien.viscardy@aeronomie.be
% 
%% Link to functions
addpath('../functions/');

%%
clearvars

tic

%% Type of experiment
t_exp_list = {'E'}; n_exp = length(t_exp_list);
E_sol_list = [2442 2446 2615 2627 2644];

%% Parameters
n_E       = 25; % Enrichment factor

%% Infos on figure
figid     = 4;
sfigtype  = 'png';
ftsz      = 16;
ftnm      = 'times';
savefig   = 1; % = 1: save Figure
savedata  = 1; % = 1: save Data (Table)
collistF  = [0.8 0 0.8];
collistE  = [0.2 0.8 0.6];
mksz      = 6;

lab_data  = {'e line' 'f line' 'g line' 'Wefg'};
lab_line  = {'e' 'f' 'g'};
lab_panel = {'(a)' '(b)' '(c)' '(d)' '(e)' '(f)' '(g)' '(h)'};

%% X-Labels
xlab = cell(1,4);
for iline = 1:3
    xlab{iline} = ['\eta_{\rm',lab_line{iline},'} \pm \sigma_{\rm',lab_line{iline},'}'];
end
xlab{4} = '\eta \pm \sigma';

%% Make figure
figure(figid)
delete(get(gcf,'children'))
set(gcf,'paperpositionmode','auto')
set(gcf,'position',[200 100 1600 800])
set(gca,'fontsize',ftsz,'fontname',ftnm)
nrow = 2; ncol = 3; isub = 1;

%% Loop over the five experiments
for iexp = 1:n_exp
    t_exp = t_exp_list{iexp};
    switch t_exp
        case 'D', sol_list = D_sol_list; enr_fct = 1;
        case 'E', sol_list = E_sol_list; enr_fct = 25;
    end
    nsol = length(sol_list);
    
    %% Loop over run sols
    for isol = 1:nsol
        sol_index   = sol_list(isol);
        
        %% Load full data
        if ( sol_index <  2442 )
            SS_MSL_full_data_Webster_2015
        else
            SS_MSL_full_data_Webster_2021
        end
        
        %% Calculation of eta and sigma
        SS_TLS_CH4_eta_sig
        disp(['Sol ',num2str(sol_index),' : ',num2str(eta,'%2.2f'),' +/- ', ...
            num2str(sig,'%2.2f'),' ppbv'])
        
        %% Number of (full- and empty-cell) runs
        nFpts     = length(F_e_line_CH4); % Number of full-cell runs
        nEpts     = length(E_e_line_CH4); % Number of empty-cell runs
        
        %% Full-cell runs
        mean_F(1) = mean(F_e_line_CH4); % mean CH4 vmr from e line     [ppbv]
        std_F(1)  = std(F_e_line_CH4);  % error on CH4 vmr from e line [ppbv]
        mean_F(2) = mean(F_f_line_CH4); % mean CH4 vmr from f line     [ppbv]
        std_F(2)  = std(F_f_line_CH4);  % error on CH4 vmr from f line [ppbv]
        mean_F(3) = mean(F_g_line_CH4); % mean CH4 vmr from g line     [ppbv]
        std_F(3)  = std(F_g_line_CH4);  % error on CH4 vmr from g line [ppbv]
        mean_F_W  = mean(F_Wefg_CH4);   % mean CH4 vmr (Wefg)          [ppbv]
        std_F_W   = std(F_Wefg_CH4);    % error on CH4 vmr (Wefg)      [ppbv]
        
        %% Empty-cell runs
        mean_E(1) = mean(E_e_line_CH4); % mean CH4 vmr from e line     [ppbv]
        std_E(1)  = std(E_e_line_CH4);  % error on CH4 vmr from e line [ppbv]
        mean_E(2) = mean(E_f_line_CH4); % mean CH4 vmr from f line     [ppbv]
        std_E(2)  = std(E_f_line_CH4);  % error on CH4 vmr from f line [ppbv]
        mean_E(3) = mean(E_g_line_CH4); % mean CH4 vmr from g line     [ppbv]
        std_E(3)  = std(E_g_line_CH4);  % error on CH4 vmr from g line [ppbv]
        mean_E_W  = mean(E_Wefg_CH4);   % mean CH4 vmr (Wefg)          [ppbv]
        std_E_W   = std(E_Wefg_CH4);    % error on CH4 vmr (Wefg)      [ppbv]
        
        %% CH4 abundance in Martian atmosphere (from each of the 3 lines)
        eta_line = zeros(1,3);
        sig_line = zeros(1,3);
        for iline = 1:3 % Loop over the 3 lines
            eta_line(iline) = (mean_F(iline) - mean_E(iline))/n_E;
            sig_line(iline) = (sqrt( std_F(iline)^2/nFpts + std_E(iline)^2/nEpts ))/n_E;
        end
        
        %% CH4 abundance in Herriott cell (Wefg)
        eta_H = mean_F_W - mean_E_W;
        sig_H = sqrt( std_F_W^2/nFpts + std_E_W^2/nEpts );
        
        %% CH4 abundance in Martian atmosphere (eta +/- sig) [ppbv]
        eta   = eta_H/n_E;
        sig   = sig_H/n_E;
        
        %% Error bar
        subplot(nrow,ncol,isub); isub = isub + 1;
        set(gca,'fontsize',ftsz,'fontname',ftnm)
        collist = hsv(3);
        
        if ne(sol_index,2442)
            yline(0,'--k','handlevisibility','off'); hold on
        end
        
        % Error bars (3 lines)
        for iline = 1:3 % Loop over the 3 lines
            he = errorbar(iline,eta_line(iline),sig_line(iline),'r'); hold on
            set(he,'Marker','o', ...
                'MarkerSize',6, ...
                'Color',collist(iline,:), ...
                'MarkerFaceColor',collist(iline,:), ...
                'MarkerEdgeColor',collist(iline,:), ...
                'LineStyle','none','LineWidth',1)
        end
        
        % Error bar (Wefg)
        collist = 'k';
        he = errorbar(4,eta,sig,'r'); hold on
        set(he,'Marker','o', ...
            'MarkerSize',6, ...
            'Color',collist, ...
            'MarkerFaceColor',collist, ...
            'MarkerEdgeColor',collist, ...
            'LineStyle','none','LineWidth',1)
        
        set(gca,'fontsize',ftsz,'fontname',ftnm)
        
        xlim([0 5])
        yl = ylim;
        
        if (sol_index > 2442)
        if yl(1)>0, yl(1) = yl(1)*0.9; else, yl(1) = yl(1)*1.5; end
        yl(2) = yl(2)*1.1;
        ylim(yl)
        end
        
        set(gca,'XTick',1:4,'XTickLabel',lab_data);
        ylabel('CH_4 abundance [ppbv]','fontsize',ftsz,'fontname',ftnm)
        
        title(['Sol ',num2str(sol_index)],'fontsize',ftsz,'fontname',ftnm)
        
        SS_mk_label_panel(lab_panel{isol},xlim,'lin',ylim,'lin',ftsz,ftnm);
        
        set(gca,'YMinorTick','on')
        set(gca,'fontsize',ftsz,'fontname',ftnm,'box','on')
    end
end

%% Legend
collist = [hsv(3); [0 0 0]];
for iline = 1:4
    ht1 = text(9.25,0.7 - 0.5*(iline-1),0,xlab{iline}, ...
        'color',collist(iline,:), ...
        'fontsize',ftsz+2,'fontname',ftnm,'HorizontalAlignment','center');
end

%% Save figure
if (savefig == 1)
    sfpath     = '../Figures/';
    sffilename = 'Figure_04.png';
    sfflnm     = fullfile(sfpath,sffilename);
    exportgraphics(gcf,sfflnm,'Resolution',400)
end

toc