% This script makes a figure showing the spectra from the Sol-2627
% experiment.
% 
% The (level 1a) data are taken from NASA's Planetary Data System (PDS)
% (http://pds-183 geosciences.wustl.edu/missions/msl/sam.htm). 
%
%   Elapsed time: 11 seconds
% 
% Author: sebastien.viscardy@aeronomie.be
% 
%% Link to functions
addpath('../functions/');

%%
clearvars

tic

%% Parameters of experiment
SAM_exp = 25556; sol_index = 2626; t_exp = {'full cell runs' 'empty cell runs'};

%% Infos on figure
figid     = 1;
sfigtype  = 'png';
ftsz      = 10;
ftnm      = 'times';
savefig   = 1; % = 1: save Figure
savedata  = 1; % = 1: save Data (Table)
collistF  = [0.8 0 0.8];
collistE  = [0.2 0.8 0.6];
lab_panel = {'(a)' '(b)' '(c)'};
klab      = 1;
yl        = [-0.01 0.015];

%% Load Level 1a: TLS housekeeping data
ldpath     = '../MSL_data';
ldfilename = ['*f',num2str(sol_index,'%04i'),'rdr1a_*_tls_ch4_*.tab'];
ldflnm     = fullfile(ldpath,ldfilename);
list       = dir(ldflnm);
nt         = length(list);
xl         = [3057.65 3057.79];

%% Store spectra
nspec      = 26;
ngr        = 2;

wvnb = zeros(1024,nspec,ngr);
tran = zeros(1024,nspec,ngr);
harm = zeros(1024,nspec,ngr);

for igr = 1:ngr
    for isp = 1:nspec
        it              = (igr-1)*nspec + isp;
        ldfilename      = list(it).name;
        ldflnm          = fullfile(ldpath,ldfilename);
        tmp             = importdata(ldflnm);
        data            = tmp.data;
        
        wvnb(:,isp,igr) = data(:,2);
        tran(:,isp,igr) = data(:,3);
        harm(:,isp,igr) = data(:,4);
    end
end

%% Interpolation on fixed grid
wvnb_grid = xl(1):1e-5:xl(2); ngrid = length(wvnb_grid);
tran_grid  = zeros(ngrid,nspec,ngr);
harm_grid = zeros(ngrid,nspec,ngr);

for igr = 1:ngr
    for isp = 1:nspec
        tran_grid(:,isp,igr) = interp1(squeeze(wvnb(:,isp,igr)),squeeze(tran(:,isp,igr)),wvnb_grid);
        harm_grid(:,isp,igr) = interp1(squeeze(wvnb(:,isp,igr)),squeeze(harm(:,isp,igr)),wvnb_grid);
    end
end

%% Averaged spectra (mean + standard deviation)
meantran = squeeze(mean(tran_grid,2));
stdtran  = squeeze(std(tran_grid,1,2));
meanharm = squeeze(mean(harm_grid,2));
stdharm  = squeeze(std(harm_grid,1,2));

%% Make figure
figure(figid)
delete(get(gcf,'children'))
set(gcf,'paperpositionmode','auto')
set(gcf,'position',[400 -100 600 1000])
nrow = ngr+1; ncol = 1; isub = 1;

%% All spectra vs. averaged spectra
for igr = 1:ngr
    subplot(nrow,ncol,isub); isub = isub + 1;
    
    % Color
    collistavg  = [0.8 0 0.8; 0.2 0.8 0.6];
    
    collist = copper(nspec);
    for isp = 1:nspec
        plot(wvnb_grid,squeeze(harm_grid(:,isp,igr)),'color',collist(isp,:)); hold on
    end
    plot(wvnb_grid,meanharm(:,igr),'color',collistavg(igr,:),'linewidth',3); hold on
    
    xlim(xl)
    ylim(yl)
    
    xlabel('wavenumber [cm^{-1}]','fontsize',ftsz,'fontname',ftnm)
    ylabel('2f signal','fontsize',ftsz,'fontname',ftnm)
    
    SS_mk_label_panel(lab_panel{klab},xlim,'lin',ylim,'lin',ftsz,ftnm)
    klab = klab + 1;
    
    title(t_exp{igr},'fontsize',ftsz,'fontname',ftnm)
    
    set(gca,'XMinorTick','on','YMinorTick','on')
    set(gca,'fontsize',ftsz,'fontname',ftnm,'Box','on')
end

%% Averaged spectra
subplot(nrow,ncol,isub); isub = isub + 1;

collist  = [0.8 0 0.8; 0.2 0.8 0.6];
set(gca,'fontsize',ftsz,'fontname',ftnm)
clear legname

for igr = 1:ngr
    hs = shadedErrorBar(wvnb_grid,meanharm(:,igr),stdharm(:,igr)); hold on
    hs.mainLine.Color     = collist(igr,:);
    hs.mainLine.LineWidth = 3;
    hs.patch.FaceColor    = collist(igr,:);
    hs.patch.FaceAlpha    = 0.2;
    legname{igr}          = t_exp{igr};
end

legend(legname)
xlim(xl)
ylim(yl)

xlabel('wavenumber [cm^{-1}]','fontsize',ftsz,'fontname',ftnm)
ylabel('2f signal','fontsize',ftsz,'fontname',ftnm)

SS_mk_label_panel(lab_panel{klab},xlim,'lin',ylim,'lin',ftsz,ftnm)
klab = klab + 1;

htt1 = text(3057.762,-0.0075,'e','HorizontalAlignment','center','fontsize',ftsz+2,'fontname',ftnm);
htt2 = text(3057.727,-0.0075,'f','HorizontalAlignment','center','fontsize',ftsz+2,'fontname',ftnm);
htt3 = text(3057.687,-0.0075,'g','HorizontalAlignment','center','fontsize',ftsz+2,'fontname',ftnm);

set(gca,'XMinorTick','on','YMinorTick','on')
set(gca,'fontsize',ftsz,'fontname',ftnm,'Box','on')

%% Save figure
if (savefig == 1)
    sfpath     = '../Figures/';
    sffilename = 'Figure_S01.png';
    sfflnm     = fullfile(sfpath,sffilename);
    exportgraphics(gcf,sfflnm,'Resolution',400)
end

toc