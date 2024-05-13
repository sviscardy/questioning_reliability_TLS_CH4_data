% Histograms (right column)
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
        lab     = 'F';
    case 'empty'
        y       = E_lines_CH4(:,i_hist);
        ymean   = E_lines_mean(i_hist);
        ystd    = E_lines_std(i_hist);
        npts    = nEpts;
        collist = collistE;
        lab     = 'E';
end

%% Histogram
minhist = ymean - 5*ystd;
maxhist = ymean + 5*ystd;

dx      = ystd/3;
hh      = histogram(y,minhist:dx:maxhist); hold on
set(hh,'FaceColor',collist,'EdgeColor',collist,'EdgeAlpha',0.2)

maxcount = max([maxcount max(hh.Values)]);

if ( i_hist <= 3 )
    legname{ileg} = ['$$\overline{\eta}_{\rm{',lab_line{i_hist},',',lab, ...
        '}} \pm ','\overline{\sigma}_{\rm{',lab_line{i_hist},',',lab, ...
        '}} = ',num2str(ymean,'%2.2f'),' \pm ', ...
        num2str(ystd,'%2.2f'),'\ \rm{ppbv}$$']; ileg = ileg + 1;
elseif ( i_hist == 4 )
    legname{ileg} = ['$$\overline{\eta}_{\rm{',lab,'}} \pm ', ...
        '\overline{\sigma}_{\rm{',lab,'}} = ', ...
        num2str(ymean,'%2.2f'),' \pm ', ...
        num2str(ystd,'%2.2f'),'\ \rm{ppbv}$$']; ileg = ileg + 1;
end