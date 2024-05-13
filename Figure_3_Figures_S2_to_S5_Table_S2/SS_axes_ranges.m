% Xlim and Ylim (time series and histograms)
% 
% Author: sebastien.viscardy@aeronomie.be
%% 1.a) Time series: Xlim
x_ts = [min([min(E_elapsed_t) min(F_elapsed_t)]) ...
    max([max(E_elapsed_t) max(F_elapsed_t)])]/3600;
x_ts = [sign(x_ts(1))*abs(x_ts(1))*0.95 x_ts(2)*1.05];

%% 1.b) Time series: Ylim
y_ts = [min([min(min(E_lines_CH4)) min(min(F_lines_CH4))]) ...
    max([max(max(E_lines_CH4)) max(max(F_lines_CH4))])];
y_ts = [min([y_ts(1)*1.2 0]) y_ts(2)*1.1];

%% 2.a) Histograms: Xlim
clear xl_distr
if     (sol_index == 2442), xl_distr = [-100 700];
elseif (sol_index == 2446), xl_distr = [-50 150];
elseif (sol_index == 2615), xl_distr = [-100 200];
elseif (sol_index == 2627), xl_distr = [-25 150];
elseif (sol_index == 2644), xl_distr = [-100 200];
end

%% 2.b) histograms: Ylim
y_hist = [0 13];