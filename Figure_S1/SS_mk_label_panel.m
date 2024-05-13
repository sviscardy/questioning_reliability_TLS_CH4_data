% SS_mk_label_panel(lab,xl,xsc,yl,ysc,ftsz,ftnm)
% 
% This function makes the label of a panel.
% 
% INPUT:
%   - lab  : label (ex. '(a)')
%   - xl   : x-axis range
%   - xsc  : x-axis scale ('lin' or 'log')
%   - yl   : y-axis range
%   - ysc  : y-axis scale ('lin' or 'log')
%   - ftsz : fontsize (ex. 12)
%   - ftnm : fontname (ex. 'times')
% 
% Author: sebastien.viscardy@aeronomie.be
% 
%%
function SS_mk_label_panel(lab,xl,xsc,yl,ysc,ftsz,ftnm)

switch xsc
    case 'lin'
        x1   = xl(1); x2 = xl(2);
        xlab = x1+(x2-x1)/20;
    case 'log'
        x1   = log10(xl(1)); x2 = log10(xl(2));
        xlab = 10^(x1+(x2-x1)/20);
end

switch ysc
    case 'lin'
        y1   = yl(1); y2 = yl(2);
        ylab = y2+(y2-y1)/8;
    case 'log'
        y1   = log10(yl(1)); y2 = log10(yl(2));
        ylab = 10^(y2+(y2-y1)/8);
end

ht = text(xlab,ylab,lab);

set(ht,'fontsize',ftsz+2,'FontWeight','bold','fontname',ftnm,'HorizontalAlignment','center')
end