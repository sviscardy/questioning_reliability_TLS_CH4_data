% ===================
% Figure 1 + Table S1
% ===================
% 
% This script calculates the CH4 amounts/abundances in the foreoptics (FO)
% chamber and Herriott cell and investigates the leak hypothesis from the
% FO chamber and sample cell.
% 
% This script makes Figure 1 in the manuscript and Table S1 in the SI.
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
%   Elapsed time: 15 seconds
% 
% Author: sebastien.viscardy@aeronomie.be (August 01, 2023)
%
%% Link to functions
addpath('../functions/');

%% 
clearvars
tic

%% Physical constants
kboltz     = 1.380649e-23; % Boltzmann constant [m2 kg s-2 K-1]
NA         = 6.022137e23;  % Avogadro's number  [mol-1]

%% Type of experiment
D_sol_list = [79 81 106 292 306 313 466 474 504 526 684];
E_sol_list = [573 684 2442 2446 2615 2627 2644];
t_exp_list = {'D' 'D' 'D' 'D' 'D' 'D' 'D' 'D' 'D' 'D' 'E' 'E' 'D' 'E' 'E' 'E' 'E' 'E'};
sol_list   = [79 81 106 292 306 313 466 474 504 526 573 684 684 2442 2446 2615 2627 2644];

%% Infos on figure
sfigtype  = 'png';
figid     = 1;
ftsz      = 14;
ftnm      = 'times';
savefig   = 1; % = 1: save Figure
savedata  = 1; % = 1: save Data (Table)
mksz      = 6;
lab_panel = {'(a)' '(b)' '(c)' '(d)' '(e)' '(f)' '(g)' '(h)'};
jsol      = 0;

%% Initialize the arrays
ntotsol        = length(D_sol_list) + length(E_sol_list);
sol_idx        = zeros(1,ntotsol);
eta_list       = zeros(1,ntotsol);
sig_list       = zeros(1,ntotsol);
F_pval_list    = zeros(1,ntotsol);
E_pval_list    = zeros(1,ntotsol);
NCH4_FO        = zeros(1,ntotsol);
NCH4_HC        = zeros(1,ntotsol);
etaF_list      = zeros(1,ntotsol);
etaH_list      = zeros(1,ntotsol);
f_list         = zeros(1,ntotsol);
H_rel_diff_p   = zeros(1,ntotsol);
H_dp_list      = zeros(1,ntotsol);
F_rel_diff_p   = zeros(1,ntotsol);
F_dp_list      = zeros(1,ntotsol);
type_exp_list  = cell(1,ntotsol);
F_prs          = zeros(1,ntotsol);
H_prs          = zeros(1,ntotsol);

nsol           = length(sol_list);

%% Loop over run sols
for isol = 1:nsol
    t_exp       = t_exp_list{isol};
    sol_index   = sol_list(isol);
    
    %% Load full data
    if ( sol_index < 2442 )
        SS_MSL_full_data_Webster_2015
    else
        SS_MSL_full_data_Webster_2021
    end
    
    %% Calculation of eta and sigma
    SS_TLS_CH4_eta_sig
    
    jsol      = jsol + 1;
    
    %% Save eta and sigma
    sol_idx(jsol)       = sol_index; % [/]
    type_exp_list{jsol} = t_exp;     % [/]
    eta_list(jsol)      = eta;       % [ppbv]
    sig_list(jsol)      = sig;       % [ppbv]
    
    %% Hypothesis: Leak from FO chamber into Herriott cell
    SS_leak_hypothesis
    
    %% Key variables (FO: Foreoptics / HC: Herriott cell)
    NCH4_FO(jsol)        = N_FO*1e9/NA;             % [nmol]
    NCH4_HC(jsol)        = NH(2)*1e9/NA;            % [nmol]
    etaF_list(jsol)      = etaF(1)*1e9;             % [ppbv]
    etaH_list(jsol)      = eta_H;                   % [ppbv]
    f_list(jsol)         = f*100;                   % [%]
    H_dp_list(jsol)      = pH(2)-pH(1);             % [Pa]
    H_rel_diff_p(jsol)   = (pH(2)-pH(1))/pH(1)*100; % [%]
    F_dp_list(jsol)      = pF(2)-pF(1);             % [Pa]
    F_rel_diff_p(jsol)   = (pF(2)-pF(1))/pF(1)*100; % [%]
    F_prs(jsol)          = pF(1);                   % [Pa]
    H_prs(jsol)          = pH(1);                   % [Pa]
end

%% Save data
if (savedata == 1)
    % Save leak output
    T = table(sol_idx',eta_list',sig_list',NCH4_FO',NCH4_HC',etaF_list',etaH_list', ...
        f_list',H_dp_list',H_rel_diff_p',F_dp_list',F_rel_diff_p');
    
    sdpath     = '../output_tables/';
    sdfilename = 'table_S1.xlsx';
    sdflnm     = fullfile(sdpath,sdfilename);
    writetable(T,sdflnm,'Sheet',1,'Range','A1')
end

%% Make figure
figure(figid)
delete(get(gcf,'children'))
set(gcf,'paperpositionmode','auto')
set(gcf,'position',[200 100 1600 1400])
set(gca,'fontsize',ftsz,'fontname',ftnm)
nrow = 3; ncol = 2; isub = 1;

% Colors
collistHC = [178 34 34]/178;  % Herriott cell
collistFO = [30 144 255]/255; % FO chamber
collistfL = [0.8 0 0.8];      % Ratios and leaking fraction f_L

% Experiments (direct-ingest 'D' and enrichment 'E')
xtlab = cell(1,ntotsol);
for ipts = 1:ntotsol
    xtlab{ipts} = [char(type_exp_list{ipts}),' ',num2str(sol_idx(ipts))];
end

pts_idx = 1:ntotsol;

%% Number of CH4 molecules
subplot(nrow,ncol,isub); isub = isub + 1;
ileg = 1;

hp = plot(pts_idx,NCH4_FO,'or'); hold on
set(hp,'marker','o','MarkerSize',mksz,'MarkerFaceColor',collistFO, ...
    'MarkerEdgeColor',collistFO,'LineStyle','none')
legname{ileg} = 'FO chamber: {\itN}_{FO}'; ileg = ileg + 1;

hp = plot(pts_idx,NCH4_HC,'or'); hold on
set(hp,'marker','o','MarkerSize',mksz,'MarkerFaceColor',collistHC, ...
    'MarkerEdgeColor',collistHC,'LineStyle','none')
legname{ileg} = 'Herriott cell: {\itN}_{H}';

hleg = legend(legname,'fontsize',ftsz-2,'fontname',ftnm);
set(hleg,'position',[0.22642 0.93361 0.10615 0.040625])

xlim([0 ntotsol+1])
ylim([0.99999e-5 1e1])

set(gca,'YScale','log','YMinorTick','on','YTick',10.^(-5:1))

set(gca,'XTick',pts_idx,'XTickLabel',xtlab)
xticklabel_rotate([],45,[],'Fontsize',ftsz,'fontname',ftnm);

grid

yll = ylabel('CH_4 amount {\itN} [nmole]','fontsize',ftsz,'fontname',ftnm);
set(yll,'Position',[-0.13 0.5 0],'VerticalAlignment','middle')

SS_mk_label_panel(lab_panel{isub-1},xlim,'lin',ylim,'log',ftsz,ftnm);

set(gca,'fontsize',ftsz,'fontname',ftnm)

%% Ratio N_FO/N_HC
subplot(nrow,ncol,isub); isub = isub + 1;

hl = plot(pts_idx,NCH4_FO./NCH4_HC,'color',collistfL); hold on
set(hl,'marker','o','markersize',mksz,'MarkerFaceColor',collistfL, ...
    'MarkerEdgeColor',collistfL,'LineStyle','none')

xlim([0 ntotsol+1])

set(gca,'YScale','log','YMinorTick','on')
set(gca,'XTick',pts_idx,'XTickLabel',xtlab)
xticklabel_rotate([],45,[],'Fontsize',ftsz,'fontname',ftnm);

yll = ylabel('{\itN}_{FO}/{\itN}_{H}','fontsize',ftsz,'fontname',ftnm);
set(yll,'Position',[-0.13 0.5 0],'VerticalAlignment','middle')

grid

SS_mk_label_panel(lab_panel{isub-1},xlim,'lin',ylim,'log',ftsz,ftnm);

set(gca,'position',[0.51034 0.73834 0.29526 0.18666],'fontsize',ftsz,'fontname',ftnm)

%% eta_F and eta_H: methane abundance in FO chamber and Herriott cell
subplot(nrow,ncol,isub); isub = isub + 1;
ileg = 1;

hp = plot(pts_idx,etaF_list); hold on
set(hp,'marker','o','MarkerSize',mksz,'MarkerFaceColor',collistFO, ...
    'MarkerEdgeColor',collistFO,'LineStyle','none')
legname{ileg} = 'FO chamber: \eta_{FO}'; ileg = ileg + 1;

hp = plot(pts_idx,etaH_list,'sb'); hold on
set(hp,'marker','o','MarkerSize',mksz,'MarkerFaceColor',collistHC, ...
    'MarkerEdgeColor',collistHC,'LineStyle','none')
legname{ileg} = 'Herriott cell: \eta_{H}';

hleg = legend(legname,'fontsize',ftsz-2,'fontname',ftnm);
set(hleg,'position',[0.22642 0.63361 0.10615 0.040625])

xlim([0 ntotsol+1])

ylim([0.9999e-1 1e5])

set(gca,'YScale','log','YMinorTick','on','YTick',10.^(-1:5))

set(gca,'XTick',pts_idx,'XTickLabel',xtlab)
xticklabel_rotate([],45,[],'Fontsize',ftsz,'fontname',ftnm);

set(gca,'YMinorTick','on')
yll = ylabel('CH_4 abundance [ppbv]','fontsize',ftsz,'fontname',ftnm);
set(yll,'Position',[-0.13 0.5 0],'VerticalAlignment','middle')

grid

SS_mk_label_panel(lab_panel{isub-1},xlim,'lin',ylim,'log',ftsz,ftnm);

set(gca,'fontsize',ftsz,'fontname',ftnm,'box','on')

%% Ratio eta_FO/eta_HC
subplot(nrow,ncol,isub); isub = isub + 1;

hl = plot(pts_idx,etaF_list./etaH_list,'color',collistfL); hold on
set(hl,'marker','o','markersize',mksz,'MarkerFaceColor',collistfL, ...
    'MarkerEdgeColor',collistfL,'LineStyle','none')

xlim([0 ntotsol+1])

set(gca,'YScale','log','YMinorTick','on')

set(gca,'XTick',pts_idx,'XTickLabel',xtlab)
xticklabel_rotate([],45,[],'Fontsize',ftsz,'fontname',ftnm);

yll = ylabel('{\it\eta}_{FO}/{\it\eta}_{H}','fontsize',ftsz,'fontname',ftnm);
set(yll,'Position',[-0.13 0.5 0],'VerticalAlignment','middle')

grid

SS_mk_label_panel(lab_panel{isub-1},xlim,'lin',ylim,'log',ftsz,ftnm);

set(gca,'position',[0.51034 0.43871 0.29526 0.18666],'fontsize',ftsz,'fontname',ftnm)

%% Mean FO and HCell pressures [Pa]
subplot(nrow,ncol,isub); isub = isub + 1;
set(gca,'fontsize',ftsz,'fontname',ftnm)
ileg = 1;

plot([0 19],[0 0],'--k','handlevisibility','off'); hold on

hp = plot(pts_idx,F_prs,'or'); hold on
set(hp,'marker','o','MarkerSize',mksz,'MarkerFaceColor',collistFO, ...
    'MarkerEdgeColor',collistFO,'LineStyle','none')
legname{ileg} = 'FO chamber: {\itp}_{FO}'; ileg = ileg + 1;

hp = plot(pts_idx,H_prs,'or'); hold on
set(hp,'marker','o','MarkerSize',mksz,'MarkerFaceColor',collistHC, ...
    'MarkerEdgeColor',collistHC,'LineStyle','none')
legname{ileg} = 'Herriott cell: {\it{p}}_{H}';

hleg = legend(legname,'fontsize',ftsz-2,'fontname',ftnm);
set(hleg,'position',[0.22642 0.33361 0.10615 0.040625])

xlim([0 ntotsol+1])
ylim([0 1200])
yl = ylim; yl(1) = -yl(2)/15; ylim(yl)

set(gca,'XTick',pts_idx,'XTickLabel',xtlab)
xticklabel_rotate([],45,[],'Fontsize',ftsz,'fontname',ftnm);

grid

set(gca,'YMinorTick','on')
yll = ylabel('mean pressure {\itp} [Pa]','fontsize',ftsz,'fontname',ftnm);
set(yll,'Position',[-0.15 0.5 0],'VerticalAlignment','middle')

SS_mk_label_panel(lab_panel{isub-1},xlim,'lin',ylim,'lin',ftsz,ftnm);

set(gca,'fontsize',ftsz,'fontname',ftnm)

%% leaking factor f
subplot(nrow,ncol,isub); isub = isub + 1;

plot([0 19],[0 0],'--k','handlevisibility','off'); hold on

hp = plot(pts_idx,f_list); hold on
set(hp,'marker','o','MarkerSize',mksz,'MarkerFaceColor',collistfL, ...
    'MarkerEdgeColor',collistfL,'LineStyle','none')

xlim([0 ntotsol+1])
set(gca,'YScale','log')

set(gca,'XTick',pts_idx,'XTickLabel',xtlab)
xticklabel_rotate([],45,[],'Fontsize',ftsz,'fontname',ftnm);

grid

set(gca,'YMinorTick','on')
yll = ylabel('leaking fraction {\itf}_L [%]','fontsize',ftsz,'fontname',ftnm);
set(yll,'Position',[-0.13 0.5 0],'VerticalAlignment','middle')

SS_mk_label_panel(lab_panel{isub-1},xlim,'lin',ylim,'log',ftsz,ftnm);

set(gca,'position',[0.51034 0.13908 0.29526 0.18666])
set(gca,'fontsize',ftsz,'fontname',ftnm)

%% Save figure
if (savefig == 1)
    sfpath     = '../Figures/';
    sffilename = 'Figure_01.png';
    sfflnm     = fullfile(sfpath,sffilename);
    exportgraphics(gcf,sfflnm,'Resolution',400)
end

toc