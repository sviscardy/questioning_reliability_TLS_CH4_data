% Time series of data points (left column)
% 
% Author: sebastien.viscardy@aeronomie.be
%%
switch t_run
    case 'full'
        x       = F_elapsed_t/3600;
        y       = F_lines_CH4(:,i_ts);
        collist = collistF;
    case 'empty'
        x       = E_elapsed_t/3600;
        y       = E_lines_CH4(:,i_ts);
        collist = collistE;
end

%% Data points
hp1 = plot(x,y); hold on
set(hp1,'Marker','o','MarkerSize',mksz,'MarkerFaceColor',collist, ...
    'MarkerEdgeColor',collist,'color',collist,'LineStyle','-', ...
    'HandleVisibility','off'); hold on

%% linear fit
X      = [ones(length(x),1) x];
b      = X\y;
y_lin  = b(2)*x+b(1);
R1     = corrcoef(x,y); % Correlation coefficient
plot(x,y_lin,'color',collist*0.7,'linewidth',2); hold on

%% Legend
switch t_run, case 'full', yvarrun = 'F'; case 'empty', yvarrun = 'E'; end

if (i_ts <= 3)
    yvar = ['\eta_{\rm{',lab_line{i_ts},',',yvarrun,'}}(t)'];
else
    yvar = ['\eta_{\rm{',yvarrun,'}}(t)'];
end
legname{ileg} = SS_mk_eqn_text(b,'t',yvar,'latex'); ileg = ileg + 1;

hl = legend(legname,'fontsize',ftsz,'fontname',ftnm,'interpreter','latex');

xxx = 0.215;
if (i_ts == 1)
    set(hl','position',[xxx 0.9265 0.16413 0.036111])
elseif (i_ts == 2)
    set(hl','position',[xxx 0.7075 0.16413 0.036111])
elseif (i_ts == 3)
    set(hl','position',[xxx 0.4880 0.16413 0.036111])
elseif (i_ts == 4)
    set(hl','position',[xxx 0.2690 0.16413 0.036111])
end

%% axes, labels
xlim(x_ts)
ylim(y_ts)

if (i_ts == 4)
    xlabel('elapsed time {\itt} [h]','fontsize',ftsz,'fontname',ftnm)
end
ylabel([lab_data{i_ts},' vmr \eta [ppbv]'],'fontsize',ftsz,'fontname',ftnm)

set(gca,'XMinorTick','on','YMinorTick','on')
set(gca,'fontsize',ftsz,'fontname',ftnm,'box','on')