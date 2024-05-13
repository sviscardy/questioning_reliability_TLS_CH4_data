% This script loads (and corrects) TLS data (Sols < 2442)
% 
% Notes:
%   - type of experiments: 'D' or 'E' stand for 'direct-ingest' and enrichment'
%   - 'sol_index' = Sol on which the TLS experiment was conducted.
%   - 'FO' and HC' stand for 'foreoptics chamber' and Herriott cell,
%   respectively.
% 
% Tables taken from:
% 
% Webster et al. (2015), Mars methane detection and variability at Gale
% crater, Science 347, 415-417, doi:10.1126/science.1261713
% 
% Author: sebastien.viscardy@aeronomie.be
%
%% enrichment factor
switch t_exp
    case 'D', enr_fct = 1;
    case 'E', enr_fct = 25;
end

%% Load definitions of data
ldpath     = '../MSL_data';
ldfilename = 'definitions_data_Webster_2015.txt';
ldflnm     = fullfile(ldpath,ldfilename);
def        = importdata(ldflnm);

%% 1) EMPTY CELL: Load MSL data
ldfilename = [t_exp,'_MSL_data_sol_',num2str(sol_index),'_empty_cell.dat'];
ldflnm     = fullfile(ldpath,ldfilename);
sol_data   = sol_index;

if ~exist(ldflnm,'file')
    ldfilename = [t_exp,'_MSL_data_sol_',num2str(sol_index-1),'_empty_cell.dat'];
    ldflnm     = fullfile(ldpath,ldfilename);
    sol_data   = sol_index - 1;
end
if ~exist(ldflnm,'file')
    ldfilename = [t_exp,'_MSL_data_sol_',num2str(sol_index+1),'_empty_cell.dat'];
    ldflnm     = fullfile(ldpath,ldfilename);
    sol_data   = sol_index + 1;
end

if ~exist(ldflnm,'file')
    disp([t_exp,' on Sol ',num2str(sol_index),' not found'])
    disp(['    ',ldflnm])
    ind_file = 0;
    return
else
    ind_file = 1;
end
tmp        = load(ldflnm);

%% 1.a) Reorganize data
E_idx        = tmp(:,1);
E_elapsed_t  = tmp(:,2);
E_FO_prs     = (tmp(:,3)-3)*100; % correction (see Webster et al., 2015, SI, p. 26)
E_laser_temp = tmp(:,4) + 273.15;
E_FO_temp    = tmp(:,5) + 273.15;
E_ref_temp   = tmp(:,6) + 273.15;
E_detec_temp = tmp(:,7) + 273.15;
E_HC_prs     = tmp(:,8)*100;
E_Wefg_CH4   = tmp(:,9);
E_diff_CH4   = tmp(:,10);

%% 2) FULL CELL: Load MSL data
ldfilename = [t_exp,'_MSL_data_sol_',num2str(sol_data),'_full_cell.dat'];
ldflnm     = fullfile(ldpath,ldfilename);
tmp        = load(ldflnm);

%% 2.a) Reorganize data
F_idx        = tmp(:,1);
F_elapsed_t  = tmp(:,2);
F_FO_prs     = (tmp(:,3)-3)*100; % correction (see Webster et al., 2015, SI, p. 26)
F_laser_temp = tmp(:,4) + 273.15;
F_FO_temp    = tmp(:,5) + 273.15;
F_ref_temp   = tmp(:,6) + 273.15;
F_detec_temp = tmp(:,7) + 273.15;
F_HC_prs     = tmp(:,8)*100;
F_Wefg_CH4   = tmp(:,9);
F_diff_CH4   = tmp(:,10);

%% Correction: Sol-684 experiment (enrichment)
% The temperatures in Table (Column 5), as documented by Webster et al.
% (2015, pp. 42-43), are not correct. Since the temperature remained
% reasonably stable during the experiments, the values have been replaced
% with those measured during the full-cell runs.
if ( strcmpi(t_exp,'E') && sol_index == 684 )
    E_FO_temp = F_FO_temp;
end

%% Definitions of data
def_data = { ...
    'Index' ...
    'Elapsed time [s]' ...
    'Foreoptics pressure [Pa]' ...
    'Laser plate temperature [K]' ...
    'Foreoptics temperature [K]' ...
    'Ref cell temperature [K]' ...
    'Science detector temperature [K]' ...
    'HCell pressure [Pa]' ...
    'Wefg in situ CH4 [ppbv]' ...
    'Diff. with mean(empty) [ppbv]'};