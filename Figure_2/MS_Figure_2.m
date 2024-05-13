% =========
% Figure 2
% =========
% 
% This script makes Figure 2 showing the pressure measurements during full-
% and empty-cell runs on Sols 306, 526, 684, and 2644.
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
%   Elapsed time: 10 seconds
% 
% Author: sebastien.viscardy@aeronomie.be
% 
%% Link to functions
addpath('../functions/');

%%
clearvars

tic

%% Type of experiment ('D': direct-ingest' / 'E': enrichment)
t_exp_list = {'D' 'E'}; n_exp = length(t_exp_list);

%% list of selected experiments
D_sol_list = [306 526];
E_sol_list = [684 2644];

%% Xlim and Ylim
SS_axes_ranges

%% Infos on figure
figid     = 2;
sfigtype  = 'png';
ftsz      = 15;
ftnm      = 'times';
savefig   = 1; % = 1: save Figure
savedata  = 1; % = 1: save Data (Table)
collistFF = [0.8 0 0.8];
collistFH = [1 0 0];
collistE  = [0.2 0.8 0.6];
alpha_lev = 0.2;
mksz      = 6;
xsc       = 'lin';
ysc       = 'lin';
lab_panel = {'(a)' '(b)' '(c)' '(d)' '(e)' '(f)' '(g)' '(h)'};
leg_loc   = {'east' 'east' 'southeast' 'southeast'};

%% Make figure
figure(figid)
delete(get(gcf,'children'))
set(gcf,'paperpositionmode','auto')
set(gcf,'position',[200 -500 1200 1400])
set(gca,'fontsize',ftsz,'fontname',ftnm)
nrow = 4; ncol = 2; isub = 1;

jsol = 0;

%% Loop over types of experiments
for iexp = 1:n_exp
    t_exp = t_exp_list{iexp};
    switch t_exp
        case 'D', sol_list = D_sol_list;
        case 'E', sol_list = E_sol_list;
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
        
        jsol     = jsol + 1;
        
        %% X- and Y-ranges
        xlim_sel = xlim_list(jsol,:);
        ylim_sel = ylim_list(jsol,:);
        
        %% 1) FO pressure
        subplot(nrow,ncol,isub); isub = isub + 1;
        clear legname
        legname = cell(1,2);
        ileg = 1;
        
        % Full runs
        prs_shift  = min(min(F_FO_prs),0);
        F_FO_prs = F_FO_prs - prs_shift;
        E_FO_prs = E_FO_prs - prs_shift;
        x      = F_elapsed_t/3600; % elapsed time [h]
        P      = F_FO_prs;         % FO pressure  [Pa]
        y      = P;                % FO pressure  [Pa]
        F_R_FO = corrcoef(x,y);    % Correlation coefficient
        
        % Linear regression
        X      = [ones(length(x),1) x];
        b      = X\y;
        F_y_lr = b(2)*x+b(1);
        
        F_delta = F_y_lr(end)-F_y_lr(1);
        
        % plot TLS data
        hp1 = plot(x,y); hold on
        set(hp1,'LineStyle','none','Marker','o','MarkerFaceColor',collistFF, ...
            'MarkerSize',mksz,'MarkerEdgeColor',collistFF,'HandleVisibility','off')
        
        % plot linear fit
        plot(x,F_y_lr,'color',collistFF/1.5,'linewidth',4); hold on
        b2 = b(2);
        if ( b(1) >= 0 )
            b1 = ['+ ',num2str(b(1),'%2.1f')];
        elseif (b(1) < 0 )
            bcst = num2str(abs(b(1)),'%2.1f');
            b1 = ['- ',bcst(1:end)];
        end
        legname{ileg} = ['F: {\itp}_{FO}({\itt}) = ', ...
            num2str(b2,'%2.1f'),' {\itt} ',b1]; ileg = ileg + 1;
        
        dp_text = ['\Delta{\itp}_{FO} = ',num2str(F_delta,'%2.1f'),' Pa'];
        x1 = xlim_sel(1); x2 = xlim_sel(2);
        y1 = ylim_sel(1); y2 = ylim_sel(2);
        ht1 = text(x1+(x2-x1)/4,y2-(y2-y1)/10,dp_text,'Color',collistFF, ...
            'FontWeight','bold','fontsize',ftsz-3,'fontname',ftnm, ...
            'HorizontalAlignment','center');
        
        % Empty runs
        % 1) FO
        x       = E_elapsed_t/3600; % elapsed time [h]
        P       = E_FO_prs;         % FO pressure  [Pa]
        y       = P;                % FO pressure  [Pa]
        E_R_FO  = corrcoef(x,y);    % Correlation coefficient
        
        % Linear regression
        X       = [ones(length(x),1) x];
        b       = X\y;
        E_y_lr  = b(2)*x+b(1);
        E_delta = E_y_lr(end)-E_y_lr(1);
        
        % plot TLS data
        hp1 = plot(x,y); hold on
        set(hp1,'LineStyle','none','Marker','o','MarkerFaceColor',collistE, ...
            'MarkerSize',mksz,'MarkerEdgeColor',collistE,'HandleVisibility','off')
        
        % plot linear fit
        plot(x,E_y_lr,'color',collistE/1.5,'linewidth',4); hold on
        b2 = b(2);
        if ( b(1) >= 0 )
            b1 = ['+ ',num2str(b(1),'%2.1f')];
        elseif (b(1) < 0 )
            bcst = num2str(abs(b(1)),'%2.1f');
            b1 = ['- ',bcst(1:end)];
        end
        legname{ileg} = ['E: {\itp}_{FO}({\itt}) = ',num2str(b2,'%2.1f'), ...
            ' {\itt} ',b1]; ileg = ileg + 1;
        
        dp_text = ['\Delta{\itp}_{FO} = ',num2str(E_delta,'%2.1f'),' Pa'];
        x1  = xlim_sel(1); x2 = xlim_sel(2);
        y1  = ylim_sel(1); y2 = ylim_sel(2);
        ht2 = text(x1+(x2-x1)*3/4,y2-(y2-y1)/10,dp_text,'Color',collistE, ...
            'FontWeight','bold','fontsize',ftsz-3,'fontname',ftnm, ...
            'HorizontalAlignment','center');
        
        % Legend
        hl1 = legend(legname,'location','south','fontsize',ftsz-3,'fontname',ftnm);
        
        % Axes, labels
        ylim([ylim_sel(1) ylim_sel(2)])
        
        set(gca,'XMinorTick','on','YMinorTick','on')
        
        if (isub >=nrow*ncol)
            xlabel('elapsed time {\itt} [h]','fontsize',ftsz,'fontname',ftnm)
        end
        ylabel('pressure {\itp}_{FO} [Pa]','fontsize',ftsz,'fontname',ftnm)
        
        % Label of panel
        xlim(xlim_sel)
        xl = xlim;
        yl = ylim;
        SS_mk_label_panel(lab_panel{isub-1},xl,xsc,yl,ysc,ftsz,ftnm);
        
        set(gca,'XMinorTick','on','YMinorTick','on','YTick',yl(1):20:yl(2))
        
        if (iexp == 1 && isol == 1)
            htvert = text((x1+x2)/2,y2+(y2-y1)/2.5,'FO chamber', ...
                'Rotation',0, ...
                'HorizontalAlignment','center','VerticalAlignment','top', ...
                'fontsize',ftsz+2,'fontname',ftnm,'FontWeight','bold');
        end
        set(gca,'fontsize',ftsz,'fontname',ftnm)
        
        texttitle = ['Sol ',num2str(sol_index),': ', ...
            num2str(eta,'%2.2f'),' \pm ',num2str(sig,'%2.2f'),' ppbv'];
        title(texttitle,'fontsize',ftsz-1,'fontname',ftnm)
        
        %% 2) HCell
        subplot(nrow,ncol,isub); isub = isub + 1;
        clear legname
        
        x       = F_elapsed_t/3600; % elapsed time   [h]
        P       = F_HC_prs;         % HCell pressure [Pa]
        y       = P;                % HCell pressure [Pa]
        F_R_HC  = corrcoef(x,y);    % Correlation coefficient
        
        % Linear regression
        X       = [ones(length(x),1) x];
        b       = X\y;
        F_y_lr  = b(2)*x+b(1);
        F_delta = F_y_lr(end)-F_y_lr(1);
        
        % plot TLS data
        hp1 = plot(x,y); hold on
        set(hp1,'LineStyle','none','Marker','v','MarkerFaceColor',collistFF, ...
            'MarkerSize',mksz,'MarkerEdgeColor',collistFF,'HandleVisibility','off')
        
        % plot linear fit
        plot(x,F_y_lr,'color',collistFF/1.5,'linewidth',4); hold on
        b2 = b(2);
        if ( b(1) >= 0 )
            b1 = ['+ ',num2str(b(1),'%2.1f')];
        elseif (b(1) < 0 )
            bcst = num2str(abs(b(1)),'%2.1f');
            b1 = ['- ',bcst(1:end)];
        end
        legname = ['F: {\itp}_{H}({\itt}) = ',num2str(b2,'%2.1f'), ...
            ' {\itt} ',b1];
        
        dp_text = ['\Delta{\itp}_H = ',num2str(F_delta,'%2.1f'),' Pa'];
        x1  = xlim_sel(1); x2 = xlim_sel(2);
        y1  = ylim_sel(3); y2 = ylim_sel(4);
        ht3 = text(x1+(x2-x1)/4,y2-(y2-y1)/10,dp_text,'Color',collistFF, ...
            'FontWeight','bold','fontsize',ftsz-3,'fontname',ftnm, ...
            'HorizontalAlignment','center');
        
        hl3 = legend(legname,'location','south','fontsize',ftsz-3,'fontname',ftnm);
        
        ylim([ylim_sel(3) ylim_sel(4)])
        yl = ylim;
        
        % Axes, labels
        xlim(xl)
        
        set(gca,'XMinorTick','on','YMinorTick','on','YTick',yl(1):20:yl(2))
        
        if (isub >=nrow*ncol)
            xlabel('elapsed time {\itt} [h]','fontsize',ftsz,'fontname',ftnm)
        end
        ylabel('pressure {\itp}_H [Pa]','fontsize',ftsz,'fontname',ftnm)
        
        % Label of panel
        SS_mk_label_panel(lab_panel{isub-1},xlim,xsc,ylim,ysc,ftsz,ftnm);
        
        if (iexp == 1 && isol == 1)
            htvert = text((x1+x2)/2,y2+(y2-y1)/2.5,'Herriott cell', ...
                'Rotation',0, ...
                'HorizontalAlignment','center','VerticalAlignment','top', ...
                'fontsize',ftsz+2,'fontname',ftnm,'FontWeight','bold');
        end
        set(gca,'fontsize',ftsz,'fontname',ftnm)
        
        title(texttitle,'fontsize',ftsz-1,'fontname',ftnm)
    end
end

%% Save figure
if (savefig == 1)
    sfpath     = '../Figures/';
    sffilename = 'Figure_02.png';
    sfflnm     = fullfile(sfpath,sffilename);
    exportgraphics(gcf,sfflnm,'Resolution',400)
end

toc