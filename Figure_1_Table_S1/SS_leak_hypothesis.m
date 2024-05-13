% Leak hypothesis: estimate of the fraction 'f' (and other parameters),
% where 'f' is defined as the fraction air inside the foreoptics (FO)
% chamber that needs to leak into the Herriott cell (HC) to obtain the CH4
% vmr measured in this cell and reported in Webster et al. (2015, 2021).
% 
% Author: sebastien.viscardy@aeronomie.be
%
%%
disp(' ')
disp(['Leak into HCell: implications: Sol ',num2str(sol_index)])
disp('=============================')
disp(' ')

%% Main parameters
l_HC   = 16.8;             % HC optical length [m]
l_FO   = 0.09;             % FO optical length [m]
V_HC   = 405;              % HCell volume      [cm3]
V_FO   = 988.2;            % FO volume         [cm3]
T      = mean(F_FO_temp);  % FO temperature    [K]
etaH   = eta_H*1e-9;       % CH4 vmr in HCell  [/]

%% Number of CH4 molecules in FO chamber [moles]
N_FO  = mean(F_HC_prs)/(kboltz*mean(F_FO_temp))*V_FO*1e-6*l_HC/l_FO*mean(E_Wefg_CH4)*1e-9;

%% Before leak
F_FO_prs = F_FO_prs - min(0,F_FO_prs(1));

% 1) FO chamber
pF(1)   = mean(F_FO_prs);         % FO pressure            [Pa]
NF(1)   = N_FO;                     % FO CH4 molecules       [molec.]
NaF(1)  = pF(1)*V_FO/(kboltz*T)*1e-6; % FO air molecules       [molec.]
etaF(1) = NF(1)/NaF(1);             % FO CH4 vmr             [/]

% 2) HCell 
pH(1)   = mean(F_HC_prs);        % HCell pressure         [Pa]
NH(1)   = 0;                        % HCell CH4 molecules    [molec.]
NaH(1)  = pH(1)*V_HC/(kboltz*T)*1e-6; % HCell air molecules    [molec.]

%% Estimate of the fraction 'f' of air molecules leaking from FO into HCell
f       = etaH*NaH(1)/((etaF(1) - etaH)*NaF(1));

%% After leak

% 1) FO chamber
NaF(2)  = NaF(1) - f*NaF(1);         % FO air molecules      [molec.]
pF(2)   = NaF(2)*kboltz*T/(V_FO*1e-6); % FO pressure           [Pa]
NF(2)   = NF(1) - f*NF(1);           % FO CH4 molecules      [molec.]
etaF(2) = NF(2)/NaF(2);              % FO CH4 vmr            [/]

% 2) HCell
NaH(2)  = NaH(1) + f*NaF(1);         % HCell air molecules   [molec.]
pH(2)   = NaH(2)*kboltz*T/(V_HC*1e-6); % HCell pressure        [Pa]
NH(2)   = NH(1) + f*NF(1);           % HCell CH4 molecules   [molec.]

%% Print key results
disp('Before leak')
disp('-----------')
disp(' ')

disp(['FO pressure     = ',num2str(pF(1),'%2.3f'),' Pa'])
disp(['FO air moles    = ',num2str(NaF(1)/NA*1e9,'%2.3e'),' nmol'])
disp(['FO CH4 moles    = ',num2str(NF(1)/NA*1e9,'%2.3e'),' nmol'])
disp(['FO CH4 vmr      = ',num2str(etaF(1)*1e6,'%2.3f'),' ppmv'])

disp(' ')
disp(['HCell pressure  = ',num2str(pH(1),'%2.3f'),' Pa'])
disp(['HCell air moles = ',num2str(NaH(1)/NA*1e9,'%2.3e'),' nmol'])
disp(['HCell CH4 moles = ',num2str(NH(1)/NA*1e9,'%2.3e'),' nmol'])
disp(['HCell CH4 vmr   = ',num2str(NH(1)/NaH(1)*1e9,'%2.3f'),' ppbv'])

disp(' ')
disp('After leak')
disp('-----------')
disp(' ')

disp(['FO pressure     = ',num2str(pF(2),'%2.3f'),' Pa'])
disp(['FO air moles    = ',num2str(NaF(2)/NA*1e9,'%2.3e'),' nmol'])
disp(['FO CH4 moles    = ',num2str(NF(2)/NA*1e9,'%2.3e'),' nmol'])
disp(['FO CH4 vmr      = ',num2str(etaF(2)*1e6,'%2.3f'),' ppmv'])

disp(' ')
disp(['HCell pressure  = ',num2str(pH(2),'%2.3f'),' Pa'])
disp(['HCell air moles = ',num2str(NaH(2)/NA*1e9,'%2.3e'),' nmol'])
disp(['HCell CH4 moles = ',num2str(NH(2)/NA*1e9,'%2.3e'),' nmol'])
disp(['HCell CH4 vmr   = ',num2str(NH(2)/NaH(2)*1e9,'%2.3f'),' ppbv'])

disp(' ')

%% Changes
disp('Differences')
disp('-----------')
disp(' ')

disp(['delta pF   = ',num2str(pF(2)-pF(1),'%2.3f'),' Pa'])
disp(['delta pH   = ',num2str(pH(2)-pH(1),'%2.3f'),' Pa'])
disp(['fraction f = ',num2str(f*100,'%2.3f'),'%'])
