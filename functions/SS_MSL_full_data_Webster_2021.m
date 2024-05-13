% This script loads (and corrects) TLS data (Sols >= 2442)
% 
% Notes:
%   - type of experiments: 'D' or 'E' stand for 'direct-ingest' and enrichment'
%   - 'sol_index' = Sol on which the TLS experiment was conducted.
%   - 'FO' and HC' stand for 'foreoptics chamber' and Herriott cell,
%   respectively.
% 
% Tables taken from:
% 
% Webster et al. (2021), Day-night differences in Mars methane suggest
% nighttime containment at Gale crater, A&A 650, A166,
% doi:10.1051/0004-6361/202040030
% 
% Author: sebastien.viscardy@aeronomie.be
%
%% enrichment factor
switch t_exp
    case 'D', enr_fct = 1;
    case 'E', enr_fct = 25;
end

%% Load definitions of data
ldpath     = '/home/sebv/mars/matlab/main_data/MSL/';
ldfilename = 'definitions_data_Webster_2021.txt';
ldflnm     = [ldpath,ldfilename];
def        = importdata(ldflnm);

%% 1) EMPTY CELL: Load TLS data
ldfilename = [t_exp,'_MSL_data_sol_',num2str(sol_index),'_empty_cell.dat'];
ldflnm     = [ldpath,ldfilename];
sol_data   = sol_index;

if ~exist(ldflnm,'file')
    ldfilename = [t_exp,'_MSL_data_sol_',num2str(sol_index-1),'_empty_cell.dat'];
    ldflnm     = [ldpath,ldfilename];
    sol_data   = sol_index - 1;
end
if ~exist(ldflnm,'file')
    ldfilename = [t_exp,'_MSL_data_sol_',num2str(sol_index+1),'_empty_cell.dat'];
    ldflnm     = [ldpath,ldfilename];
    sol_data   = sol_index + 1;
end

if ~exist(ldflnm,'file')
    disp([t_exp,' on Sol ',num2str(sol_index),' not found'])
    ind_file = 0;
    return
else
    ind_file = 1;
end
tmp        = load(ldflnm);

%% 1.a) Reorganize data
E_idx        = tmp(:,1);
E_elapsed_t  = tmp(:,2);
E_FO_prs     = tmp(:,3)*100;
E_laser_temp = tmp(:,4) + 273.15;
E_FO_temp    = tmp(:,5) + 273.15;
E_ref_temp   = tmp(:,6) + 273.15;
E_detec_temp = tmp(:,7) + 273.15;
E_HC_prs     = tmp(:,8)*100;
E_e_line_CH4 = tmp(:,9);
E_f_line_CH4 = tmp(:,10);
E_g_line_CH4 = tmp(:,11);
E_Wefg_CH4   = tmp(:,12);

%% 1.b) Correction: 'e line' column (error in Tables from Webster et al., 2021)
ttt          = 4*E_Wefg_CH4 - E_f_line_CH4 - 2*E_g_line_CH4;
E_e_line_CH4(2:end,:) = ttt(2:end);

%% 2) FULL CELL: Load TLS data
ldfilename = [t_exp,'_MSL_data_sol_',num2str(sol_data),'_full_cell.dat'];
ldflnm     = [ldpath,ldfilename];
tmp        = load(ldflnm);

%% 2.a) Reorganize data:
F_idx        = tmp(:,1);
F_elapsed_t  = tmp(:,2);
F_FO_prs     = tmp(:,3)*100;
F_laser_temp = tmp(:,4) + 273.15;
F_FO_temp    = tmp(:,5) + 273.15;
F_ref_temp   = tmp(:,6) + 273.15;
F_detec_temp = tmp(:,7) + 273.15;
F_HC_prs     = tmp(:,8)*100;
F_e_line_CH4 = tmp(:,9);
F_f_line_CH4 = tmp(:,10);
F_g_line_CH4 = tmp(:,11);
F_Wefg_CH4   = tmp(:,12);

%% 3) Definitions of data
def_data = { ...
    'Index' ...
    'Elapsed time [s]' ...
    'Foreoptics pressure [Pa]' ...
    'Laser plate temperature [K]' ...
    'Foreoptics temperature [K]' ...
    'Ref cell temperature [K]' ...
    'Science detector temperature [K]' ...
    'HCell pressure [Pa]' ...
    'e line in situ CH4 [ppbv]' ...
    'f line in situ CH4 [ppbv]' ...
    'g line in situ CH4 [ppbv]' ...
    'Wefg in situ CH4 [ppbv]'};