% Means and standard deviations (each triplet line separately + Wefg)
% 
% 'F_...': full-cell runs
% 'E_...': empty-cell runs
% 
% Author: sebastien.viscardy@aeronomie.be
%%
nEpts = length(E_e_line_CH4); % number of empty data
nFpts = length(F_e_line_CH4); % number of full data

%% Store data (e-line, f-line, g-line, Wefg)
E_lines_CH4      = zeros(nEpts,4);
F_lines_CH4      = zeros(nFpts,4);

E_lines_CH4(:,1) = E_e_line_CH4;
E_lines_CH4(:,2) = E_f_line_CH4;
E_lines_CH4(:,3) = E_g_line_CH4;
E_lines_CH4(:,4) = E_Wefg_CH4;

F_lines_CH4(:,1) = F_e_line_CH4;
F_lines_CH4(:,2) = F_f_line_CH4;
F_lines_CH4(:,3) = F_g_line_CH4;
F_lines_CH4(:,4) = F_Wefg_CH4;

%% means and standard deviations
E_lines_mean = zeros(1,4);
F_lines_mean = zeros(1,4);
E_lines_std  = zeros(1,4);
F_lines_std  = zeros(1,4);
for jline = 1:4
    % Full-cell runs (where x = 'e', 'f', 'g', 'H')
    F_lines_mean(jline) = mean(F_lines_CH4(:,jline)); % 'eta_x_F' 
    F_lines_std(jline)  = std(F_lines_CH4(:,jline));  % 'sig_x_F'
    
    % Empty-cell runs (where x = 'e', 'f', 'g', 'H')
    E_lines_mean(jline) = mean(E_lines_CH4(:,jline)); % 'eta_x_E'
    E_lines_std(jline)  = std(E_lines_CH4(:,jline));  % 'sig_x_E'
end