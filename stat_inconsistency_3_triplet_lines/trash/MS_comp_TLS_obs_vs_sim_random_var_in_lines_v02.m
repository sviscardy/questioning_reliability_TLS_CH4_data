clearvars

tic

% marscst

disp(' ')
disp('Simulation of experiments')
disp('=========================')
disp(' ')

%% Simulations: parameter
nsim   = 1e3;
ER_mu  = 1;
ER_sig = 0;

%% Type of experiment
E_sol_list = [2442 2446 2615 2627 2644];
nsol       = length(E_sol_list);
sc_list    = {'exCH4' 'leak' 'noCH4'}; nsc = length(sc_list);
t_exp      = 'E';

%% Infos on figure
figid     = 514;
sfigtype  = 'png';
ftsz      = 11;
ftnm      = 'times';
savefig   = 0;
savedata  = 0;
collistF  = [0.8 0 0.8];
collistE  = [0.2 0.8 0.6];
mksz      = 6;
lab_data  = {'e line' 'f line' 'g line' 'Wefg'};
lab_line  = {'e' 'f' 'g'};
lab_panel = {'(a)' '(b)' '(c)' '(d)' '(e)' '(f)' '(g)' '(h)' '(i)' '(j)' '(k)' '(l)'};

%% Make figure
myfig(figid)
set(gcf,'position',[200 -500 1400 1400])
set(gca,'fontsize',ftsz,'fontname',ftnm,'box','on')
nrow = nsol; ncol = 3; isub = 1;

xminall = NaN;
xmaxall = NaN;
yminall = NaN;
ymaxall = NaN;

%% Loop over Sols
for isol = 1:nsol
    sol_index = E_sol_list(isol);
    
    disp(' ')
    disp('=================')
    disp(['   Sol ',num2str(sol_index)])
    disp('=================')
    disp(' ')
    
    %% Loop over scenarios
    for isc = 1:nsc
        scenario = sc_list{isc};
        
        %% Load TLS obs. data + calculate stat (mean + standard error)
        SS_load_data_and_stat_TLS_obs
        
        [F_sigonmean,F_sigonsig] = ...
            SS_calc_sigonmean_sigonsig(F_lines_mean(1:3),F_lines_std(1:3));
        [E_sigonmean,E_sigonsig] = ...
            SS_calc_sigonmean_sigonsig(E_lines_mean(1:3),E_lines_std(1:3));
        
        return
        
        F_sigonmean = zeros(1,3);
        %         F_sigonsig  = zeros(1,3);
        %         F_sigonmean = F_sigonmean * sqrt(nFpts);
        F_sigonsig  = F_sigonsig * sqrt(nFpts);
        %         F_sigonsig  = zeros(1,3);
        E_sigonmean = zeros(1,3);
        %         E_sigonsig  = zeros(1,3);
        %         E_sigonmean = E_sigonmean * sqrt(nEpts);
        E_sigonsig  = E_sigonsig * sqrt(nEpts);
        
        %% Main parameters for simulations
        switch scenario
            case 'exCH4', isCH4  = 1; isleak = 0;
            case 'leak' , isCH4  = 0; isleak = 1;
            case 'noCH4', isCH4  = 0; isleak = 0;
        end
        
        if (isCH4 == 1) && (isleak == 0)
            disp('initial methane abundance in HCell')
            HC_mu0   = TLS_eta*ER_mu;
            F_Dmu0   = 0;
            leak_vmr = HC_mu0;
        elseif (isCH4 == 0) && (isleak == 1)
            disp('leaks of methane into HCell')
            HC_mu0   = 0;
            F_Dmu0   = 2*TLS_eta*ER_mu;
            leak_vmr = F_Dmu0;
        elseif (isCH4 == 0) && (isleak == 0)
            disp('no methane and no leaks into HCell')
            HC_mu0   = 0;
            F_Dmu0   = 0;
            leak_vmr = 0;
        else, error('wrong parameters')
        end
        
        %% Determine parameters and variables for simulations
        dt        = 2.7*60;               % time between two meas.   [s] (Webster2021, p. 3)
        Wefg_mean = E_lines_mean(4);      % mean CH4 vmr (basis)     [ppbv]
        F_mu      = Wefg_mean+HC_mu0;     % mean of data             [ppbv]
        E_mu      = Wefg_mean;            % mean of data             [ppbv]
        F_sigma   = F_lines_std(1:3);     % std err. in 3 lines      [ppbv]
        E_sigma   = E_lines_std(1:3);     % std err. in 3 lines      [ppbv]
        w_fct     = [1 1 2];              % weights of the 3 lines   [/]
        F_Dmu     = F_Dmu0;               % full : delta(mu)         [ppbv]
        E_Dmu     = 0;                    % empty: delta(mu)         [ppbv]
        F_alpha   = F_Dmu/((nFpts-1)*dt); % full : trend of the mean [ppbv s-1]
        E_alpha   = E_Dmu/((nEpts-1)*dt); % empty: trend of the mean [ppbv s-1]
        nb_sig    = 1;                    % criterion pos. detection [/]
        
        %% Print key parameters
        disp(' ')
        disp('Key parameters')
        disp('--------------')
        disp(' ')
        disp(['nb sim  = ',num2str(nsim)])
        disp(['dt      = ',num2str(dt/60),' min.'])
        disp(['mu(F)   = ',num2str(F_mu),' ppbv'])
        disp(['mu(E)   = ',num2str(E_mu),' ppbv'])
        disp(['F_Dmu   = ',num2str(F_Dmu),' ppbv'])
        disp(['E_Dmu   = ',num2str(E_Dmu),' ppbv'])
        disp(['F_alpha = ',num2str(F_alpha*3600),' ppbv h-1'])
        disp(['E_alpha = ',num2str(E_alpha*3600),' ppbv h-1'])
        disp(['nb_sig  = ',num2str(nb_sig)])
        disp(' ')
        
        %% Initialization of tensors
        F_data = zeros(nsim,nFpts,3);
        E_data = zeros(nsim,nEpts,3);
        Fx_idx = zeros(nsim,nFpts);
        Ex_idx = zeros(nsim,nEpts);
        
        %% Generate data for each line
        for ip = 1:nFpts, Fx_idx(:,ip) = ip; end
        for ip = 1:nEpts, Ex_idx(:,ip) = ip; end
        
        % initial mean + trend + variability
        for iline = 1:3
            F_data(:,:,iline) = ...
                (F_mu+F_sigonmean(iline)*randn(nsim,1)) + ...
                F_alpha*(Fx_idx-1)*dt + ...
                (F_sigma(iline)+F_sigonsig(iline)*randn(nsim,1)).*randn(nsim,nFpts);
            
            E_data(:,:,iline) = ...
                (E_mu+E_sigonmean(iline)*randn(nsim,1)) + ...
                E_alpha*(Ex_idx-1)*dt + ...
                (E_sigma(iline)+E_sigonsig(iline)*randn(nsim,1)).*randn(nsim,nEpts);
        end
        
        %% Statistical analysis of simulated data
        SS_stat_sim_data
        
        %% Calculate pressure and air density changes in both the FO chamber and HCell
        SS_pressure_changes_in_FO_HCell
        
        %% Make panel: histogram2
        SS_GS_comp_eta_vs_sigma_stat_model_vs_TLS_measurements2
        
        %% Save data
        if (savedata == 1)
            sdpath = 'data4/';
            sdfilename = ['Sol_',num2str(sol_index),'_',scenario,'_3_lines_data.mat'];
            sdflnm     = [sdpath,sdfilename];
            save(sdflnm,'HC_mu0','F_Dmu0','leak_vmr', ...
                'F_mu','E_mu','F_Dmu','E_Dmu', ...
                'nsim','nFpts','nEpts','dt','w_fct', ...
                'F_sigma','E_sigma','Fx_idx','Ex_idx','F_data','E_data','-mat')
        end
    end
end

%% make labels + use same limits for x and y axes
for isub = 1:nrow*ncol
    subplot(nrow,ncol,isub)
    xlim([xminall xmaxall])
    ylim([yminall ymaxall])
    SS_mk_label_panel(lab_panel{isub},xlim,'lin',ylim,'lin',ftsz,ftnm)
end

%% Save figure
if (savefig == 1)
    sfpath     = 'Figures/';
    sffilename = 'leak_vs_no_leak11.png';
    sfflnm     = fullfile(sfpath,sffilename);
    %             saveas(gcf,sfflnm,'png');
    exportgraphics(gcf,sfflnm,'Resolution',600)
end

toc