% Run all scripts producing:
%   - Figures 1-4 and S1-S5
%   - Tables S1 and S2
%   - estimates of probabilities presented in Section 5.
% 
% Author: sebastien.viscardy@aeronomie.be
% 
%%
clearvars

%% Figure 1 and Table 1
cd Figure_1_Table_S1/
MS_Figure_1_Table_S1
cd ..

%% Figure 2
cd Figure_2/
MS_Figure_2
cd ..

%% Figures 3, S2-S5, and Table S2
cd Figure_3_Figures_S2_to_S5_Table_S2/
MS_Figure_3_Figures_S2_to_S5_Table_S2
cd ..

%% Figure 4
cd Figure_4/
MS_Figure_4
cd ..

%% Figure S1
cd Figure_S1/
MS_Figure_S1
cd ..

%% Section 5: Uncertainties associated with the retrieval process
cd stat_inconsistency_3_triplet_lines/
MS_nb_exp_consist_3_lines
MS_probability_SIM_CH4_vmr_3_triplet_lines
cd ..

disp(' ')
disp('   --- END ---')