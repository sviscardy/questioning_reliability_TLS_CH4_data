% Normal distributions (right column)
% 
% Author: sebastien.viscardy@aeronomie.be
%%
switch t_run
    case 'full'
        y       = F_lines_CH4(:,i_hist);
        ymean   = F_lines_mean(i_hist);
        ystd    = F_lines_std(i_hist);
        npts    = nFpts;
        collist = collistF;
    case 'empty'
        y       = E_lines_CH4(:,i_hist);
        ymean   = E_lines_mean(i_hist);
        ystd    = E_lines_std(i_hist);
        npts    = nEpts;
        collist = collistE;
end

%% Normal distribution
minhist = ymean - 5*ystd;
maxhist = ymean + 5*ystd;

dx      = ystd/3;
dxc     = ystd/30;
xc      = (minhist:dxc:maxhist)+dxc/2;
yc      = npts/(sqrt(2*pi)*ystd)*exp(-0.5*(xc-ymean).^2/ystd^2)*dx;

plot(xc,yc,'color',collist*0.7,'linewidth',2,'HandleVisibility','off'); hold on

%% Axes, labels
ylim(y_hist); yl = ylim;

xlabel([lab_data{i_hist},' CH_4 vmr [ppbv]'],'fontsize',ftsz,'fontname',ftnm)
ylabel('count','fontsize',ftsz,'fontname',ftnm)

set(gca,'YTick',yl(1):1:yl(2))
yt  = get(gca,'YTick');
ytl = get(gca,'YTickLabel');

nyt = length(yt);
for iyt = 2:2:nyt
    ytl{iyt} = ' '; 
end
set(gca,'YTickLabel',ytl)

set(gca,'XMinorTick','on','fontname',ftnm,'fontsize',ftsz)
