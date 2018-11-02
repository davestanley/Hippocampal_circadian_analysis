
function run_allrats_ergodic
    global freq_band freq_listing

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%% Plotting parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% %%%%%% Set up paths %%%%%% 
    
    % Setup display params
    format compact
    
    % Get home folder location
    if ispc; userdir= getenv('USERPROFILE');
    else userdir= getenv('HOME');
    end
    
    % Setup paths to submodules and supporting functions
    restoredefaultpath
    addpath(genpath(fullfile('..','submodules','lib_MAScPhD_Matlab')));
    addpath(genpath(fullfile('..','submodules','lib_dav')));
    addpath(genpath(fullfile('..','submodules','SigProc-Plott')));
    addpath(genpath(fullfile('..','submodules','chronux')));
    addpath(fullfile('.','cosinor_dav'));
    addpath(fullfile('.','funcs_supporting'));
    addpath(fullfile('.','funcs_supporting_demo'));
    
    

    %% %%%%%% Plot preprocessing parameters %%%%%% 
    
    recalc = 0;     % Re-extract all sinusoids for all frequency bands;
                    % otherwise, load from saved file

    cosinor_mode = 1;       % Setting to 1 will return the theta band power (or "extract_band" if specified)
                            % Setting to 2 will return the power in just theta epochs
                            % Setting to 3 will return the power in non-theta epochs
                            % Setting to 0 will return the theta epoch probability (data_binary)
    rat2plot = 2;

    % Cosinor analysis
    compare_amp = 0;
        take_log = 0;
        normalize_amps = 1;
            percent_difference = 0;
    
    % Raw data analysis
    plot_Amp_vs_rawamp = 0;
    analyze_rawamp = 0;     % Use rawamp values in place of all amp data
    
    analyze_only_seizing = 1; % 1 for seizing group, 0 for non-seizing group, 2 for controls - only works for some codes!  
    prepostchronic = 0;     % 0 for all; 1 for pre; 2 for latency; 3 for chronic
    smoothmode_ts_main = 0;     % This is only used for extracting the initial Amp and phi arrays
                                % Should be left as default (0) unless want to test effects of normalizing
    smoothmode_ts = 0 ;          % Used for plot_timeseries and princomp
                % smoothmode_ts - Determines filter to apply to ts data
                %    = -2 : derivative of unwrapped angle as determined by Hilbert transform
                %    = -1 : 0 + estimate amplitude by taking variance
                %    = 0 : Pre smoothing + baseline subtraction (default value of r.ctrl.d)
                %    = 1 : No pre-smoothing
                %    = 2 : Baseline signal: 1 day 90% aveage
                %    = 3 : Pre smoothing - 6-hr 90% moving average
                %    = 4 : 3 + Pull out amplitude based on cosine fit - cosfit not working, so for now is just taking variance (same as -1)
                %    = 5 : 3 + Pull out phase based on cosine fit
                %    = 11 : From scratch, baseline subtract (should be same as 0)
                %    = 12 : From scratch, baseline subtract and normalize (percent change)

    
    %% %%%%%% Plot Choices %%%%%%     
    
        % Time series plots
    plot_timeseries = 0;
        os.shift = 5;
        
        % Figs specifically for paper
    plot_phasefits_shiftbands = 0;
    plot_phasefits_allbands = 0;
    

        % Cosinor plots
    plot_on_imagesc_individual = 0;
    plot_amp_phase_correlation2 = 0;
    ergodic_mode = 0;
        ratrange0 = -1;             % Rat range - which rats to use. -1 = all rats
        compare_theta_delta = 0;

        % Correlation Plots
    plot_corrcoef_EEG_vs_EEG = 0; % IT's okay, this is fast.
    plot_princomp = 0;              % Standard PCA analysis
    plot_princomp2 = 1;             % Look at PCA across stages (pre, post, chronic)
    plot_corr_phaseshift = 0;       % Measure how well each frequency band correlates with the band that phase shifts; plots as a subset of plot_corrcoef_EEG_vs_EEG
    plot_princomp_amp_vs_ampraw = 0; % Correlation coefficient
    plot_corrcoef_all = 0;              % Look for correlations in everything - EEG amplitude, baseline, and 
    plot_EEG_vs_baseline_corr = 0;      % Compare changes in basline to circadian amp. Need to set smoothmode_ts to -1

    
    
        % Not used stuff
    plot_EEG_vs_spks2 = 0;  % New version. Takes into account smoothemode_ts and also has some better plotting options
        reload_spks = 0;
    plot_movingphase = 0;   % CAn delete
        show_amp = 0;
        
   
    
    clear r1 r4 r9 r10 r8 r5 r11 r12

    FS_axis_sm = 14;
    do_gamma_range = 0;
    do_theta_range = 0;
    do_all = 1;
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%% Plotting presets for demo %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if plot_on_imagesc_individual || plot_amp_phase_correlation2
        normalize_amps = 1;
            percent_difference = 0;
    end
    
    if ergodic_mode
        normalize_amps = 0;
            percent_difference = 0;
    end
    
    if plot_corrcoef_EEG_vs_EEG
        rat2plot = 1;
    end
    
    if plot_princomp
        prepostchronic = 3;
    end
    
    if plot_corrcoef_all
        prepostchronic = 0;
    end
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%% Load data %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        
    freq_band = [1 5; 5 10; 12 25; 25 50; 50 90; 90 140; 140 200; 200 500];      % Gamma as defined in Buzsaki, Wang 2012 Ann Rev - Mechanisms of gamma oscillations; Theta 6-10 as optimized

    if recalc
        
        if cosinor_mode == 1
            load Rallrats_major_freqbands_theta5_10_redo_R8_shortacute.mat
        elseif cosinor_mode == 2
            load Rallrats_major_freqbands_theta5_10_th_redo_R8_shortacute.mat
        elseif cosinor_mode == 3
            load Rallrats_major_freqbands_theta5_10_delta_redo_R8_shortacute.mat
        end
        
        [r1arr] = ratscript_FFT_thetadelta2_arr (1, 2,[5 10],[1 5],[],1,cosinor_mode,freq_band);
        [r4arr] = ratscript_FFT_thetadelta2_arr (4, 2,[5 10],[1 5],[],1,cosinor_mode,freq_band);
        [r9arr] = ratscript_FFT_thetadelta2_arr (9, 2,[5 10],[1 5],[],1,cosinor_mode,freq_band);
        [r10arr] = ratscript_FFT_thetadelta2_arr (10, 2,[5 10],[1 5],[],1,cosinor_mode,freq_band);
        
        [r6arr] = ratscript_FFT_thetadelta2_arr (6, 2,[5 10],[1 5],[],1,cosinor_mode,freq_band);
        [r7arr] = ratscript_FFT_thetadelta2_arr (7, 2,[5 10],[1 5],[],1,cosinor_mode,freq_band);
        
        [r5arr] = ratscript_FFT_thetadelta2_arr (5, 2,[5 10],[1 5],[],1,cosinor_mode,freq_band);
        [r8arr] = ratscript_FFT_thetadelta2_arr (8, 2,[5 10],[1 5],[],1,cosinor_mode,freq_band);
        [r11arr] = ratscript_FFT_thetadelta2_arr (11, 2,[5 10],[1 5],[],1,cosinor_mode,freq_band);
        [r12arr] = ratscript_FFT_thetadelta2_arr (12, 2,[5 10],[1 5],[],1,cosinor_mode,freq_band);
        
        r1 = cellpack (r1arr);
        r4 = cellpack (r4arr);
        r9 = cellpack (r9arr);
        r10 = cellpack (r10arr);
        
        r6 = cellpack (r6arr);
        r7 = cellpack (r7arr);
        
        r5 = cellpack (r5arr);
        r8 = cellpack (r8arr);
        r11 = cellpack (r11arr);
        r12 = cellpack (r12arr);

        if cosinor_mode == 1
            save ('Rallrats_major_freqbands_theta5_10_redo_R8_shortacute.mat','r1','r4','r9','r10','r5','r8','r11','r12','r6','r7')
        elseif cosinor_mode == 2
            save ('Rallrats_major_freqbands_theta5_10_th_redo_R8_shortacute.mat','r1','r4','r9','r10','r5','r8','r11','r12','r6','r7')
        elseif cosinor_mode == 3
            save ('Rallrats_major_freqbands_theta5_10_delta_redo_R8_shortacute.mat','r1','r4','r9','r10','r5','r8','r11','r12','r6','r7')
        end
    else
        if cosinor_mode == 1
            load Rallrats_major_freqbands_theta5_10_redo_R8_shortacute.mat
        elseif cosinor_mode == 2
            load Rallrats_major_freqbands_theta5_10_th_redo_R8_shortacute.mat
        elseif cosinor_mode == 3
            load Rallrats_major_freqbands_theta5_10_delta_redo_R8_shortacute.mat
        end
        
    end
    
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%% Cosinor analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Pull out data from cosinor analysis on each rat, frequency band, and
    % stage
    % Index format: (freq band, rat number, pre/acute/chronic)
    w = 2*pi/1.0;
    alpha = .05;
    
    for i = 1:length(r1)
        % Seizing
        if analyze_only_seizing == 1
            r = r4;
            sout{i,1,1} = cosinor_struct(r{i}.ctrl.t,r{i}.ctrl.d,w,alpha);
            sout{i,1,2} = cosinor_struct(r{i}.acute.t,r{i}.acute.d,w,alpha);
            sout{i,1,3} = cosinor_struct(r{i}.chr.t,r{i}.chr.d,w,alpha);

            r = r9;
            sout{i,2,1} = cosinor_struct(r{i}.ctrl.t,r{i}.ctrl.d,w,alpha);
            sout{i,2,2} = cosinor_struct(r{i}.acute.t,r{i}.acute.d,w,alpha);
            sout{i,2,3} = cosinor_struct(r{i}.chr.t,r{i}.chr.d,w,alpha);

            r = r10;
            sout{i,3,1} = cosinor_struct(r{i}.ctrl.t,r{i}.ctrl.d,w,alpha);
            sout{i,3,2} = cosinor_struct(r{i}.acute.t,r{i}.acute.d,w,alpha);
            sout{i,3,3} = cosinor_struct(r{i}.chr.t,r{i}.chr.d,w,alpha);

            r = r1;
            sout{i,4,1} = cosinor_struct(r{i}.ctrl.t,r{i}.ctrl.d,w,alpha);
            sout{i,4,2} = cosinor_struct(r{i}.acute.t,r{i}.acute.d,w,alpha);
            sout{i,4,3} = cosinor_struct(r{i}.chr.t,r{i}.chr.d,w,alpha);
            
        elseif analyze_only_seizing == 2
            % Control
            r = r6;
            sout{i,1,1} = cosinor_struct(r{i}.ctrl.t,r{i}.ctrl.d,w,alpha);
            sout{i,1,2} = cosinor_struct(r{i}.acute.t,r{i}.acute.d,w,alpha);
            sout{i,1,3} = [];

            r = r7;
            sout{i,2,1} = cosinor_struct(r{i}.ctrl.t,r{i}.ctrl.d,w,alpha);
            sout{i,2,2} = cosinor_struct(r{i}.acute.t,r{i}.acute.d,w,alpha);
            sout{i,2,3} = [];
        else
            % Non-seizing
            r = r8;
            sout{i,1,1} = cosinor_struct(r{i}.ctrl.t,r{i}.ctrl.d,w,alpha);
            sout{i,1,2} = cosinor_struct(r{i}.acute.t,r{i}.acute.d,w,alpha);
            sout{i,1,3} = [];

            r = r11;
            sout{i,2,1} = cosinor_struct(r{i}.ctrl.t,r{i}.ctrl.d,w,alpha);
            sout{i,2,2} = cosinor_struct(r{i}.acute.t,r{i}.acute.d,w,alpha);
            sout{i,2,3} = [];

        end
    end
    

    % Pack all rat data into a single cell array
    r_seiz{1} = r4;
    r_seiz{2} = r9;
    r_seiz{3} = r10;
    r_seiz{4} = r1;
    r_ctrl{1} = r6;
    r_ctrl{2} = r7;
    r_ns{1} = r8;
    r_ns{2} = r11;
    %r_ns{3} = r5;
    %r_ns{4} = r12;
    
    
    % Pull values out of cosinor structure
    phi=[]; CI_phi_min=[]; CI_phi_max=[]; p_3a=[]; M=[]; beta=[]; gamma=[];
    for i = 1:size(sout,1)
        for j = 1:size(sout,2)
            for k = 1:size(sout,3)
                if ~isempty(sout{i,j,k})
                    phi(i,j,k) = sout{i,j,k}.phi;
                    CI_phi_min(i,j,k) = sout{i,j,k}.CI_phi_min;
                    CI_phi_max(i,j,k) = sout{i,j,k}.CI_phi_max;
                    p_3a(i,j,k) = sout{i,j,k}.p_3a;

                    RNE = sout{i,j,k}.RNE;
                    M(i,j,k) = RNE(1,4); beta(i,j,k) = RNE(2,4); gamma(i,j,k) = RNE(3,4);
                    Amp_old(i,j,k) = sqrt(beta(i,j,k)^2 + gamma(i,j,k)^2);
                    Amp(i,j,k) = sout{i,j,k}.Amp;
                else
                    fprintf('Is empty for i=%d, j=%d, k=%d \n',i,j,k);
                    phi(i,j,k) = NaN;
                    CI_phi_min(i,j,k) = NaN;
                    CI_phi_max(i,j,k) = NaN;
                    p_3a(i,j,k) = NaN;

                    RNE = NaN;
                    M(i,j,k) = NaN; beta(i,j,k) = NaN; gamma(i,j,k) = NaN;
                    Amp_old(i,j,k) = NaN;
                    Amp(i,j,k) = NaN;
                end
            end
        end
    end
    
    
    % Circle shifting of phases and confidence intervals.
    % This code is extremely messed up and dated  ... don't look! (particularly the
    % shift_thresholds_confid function)
    thresh = 0.20;
    if analyze_only_seizing
        CI_phi_min = shift_thresholds_confid(CI_phi_min,phi,thresh);
        CI_phi_max = shift_thresholds_confid(CI_phi_max,phi,thresh);
        phi = shift_thresholds(phi,thresh);
    else
        CI_phi_min = shift_thresholds_confid(CI_phi_min,phi,thresh);
        CI_phi_max = shift_thresholds_confid(CI_phi_max,phi,thresh);
        phi = shift_thresholds(phi,thresh);
    end
    
    % Converts confidence values into distance from phi, and then takes
    % average for display purposes
    confid1 = abs(CI_phi_min - phi);
    confid2 = abs(CI_phi_max - phi);
    confid = (confid1 + confid2) / 2;
    
    % Shoot out some phase values for Rat 10. Argue that it experiences a
    % minor phase shift
    if analyze_only_seizing == 1
        phi(5,3,1)*24
        confid(5,3,1)*24

        phi(5,3,2)*24
        confid(5,3,2)*24

        phi(5,3,3)*24
        confid(5,3,3)*24
    end
    
    % Zero amplitude fails
    za_fail = p_3a > 0.05;
    
    
    if take_log
        Amp = log(Amp);
    end
    
    if normalize_amps
        if ~percent_difference
            Amp = Amp ./ (repmat(Amp(:,:,1),[1 1 3]));
        else
            Amp = (Amp - repmat(Amp(:,:,1),[1 1 3])) ./ abs(repmat(Amp(:,:,1),[1 1 3]));
        end

    end
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%% Plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    % Plot images for individual rats
    if plot_on_imagesc_individual
        plot_imagesc_individual(phi,Amp,compare_amp,analyze_only_seizing,za_fail,normalize_amps,rat2plot)
    end
    
    % Correlations
    if plot_amp_phase_correlation2
        func_plot_correlation(phi, Amp);
    end
    
    
    if plot_timeseries || plot_movingphase
        %prepostchronic = 0;
        i=1;r=r4; [T{i} X{i}] = extract_ts(r,0, smoothmode_ts); 
        i=2;r=r9; [T{i} X{i}] = extract_ts(r,0, smoothmode_ts); 
        i=3;r=r10; [T{i} X{i}] = extract_ts(r,0, smoothmode_ts); 
        i=4;r=r1; [T{i} X{i}] = extract_ts(r,0, smoothmode_ts); 
        
        i=5;r=r5; [T{i} X{i}] = extract_ts(r,0, smoothmode_ts); 
        i=6;r=r8; [T{i} X{i}] = extract_ts(r,0, smoothmode_ts);
            T{i} = T{i}(1:end-4,:); X{i} = X{i}(1:end-4,:); % Don't display trailing data points to save space.
        i=7;r=r11; [T{i} X{i}] = extract_ts(r,0, smoothmode_ts); 
        
        i=8;r=r6; [T{i} X{i}] = extract_ts(r,0, smoothmode_ts); 
        i=9;r=r7; [T{i} X{i}] = extract_ts(r,0, smoothmode_ts);
      

        if plot_timeseries
            if analyze_only_seizing == 1
                switch rat2plot
                    case 1; i=1;r=r4;j=4; figure('Color','w','Position',[80 2 788   732]); hold on; plot_matrix2(T{i},X{i},os,'.-','MarkerSize',6);[h1, h2]=add_stimseiz2(j); legend([h1(1) h2(1)],'Stimulation','Seizure')
                        text(2,41,'Healthy','FontSize',20); text(15,41,'Latent','FontSize',20); text(50,41,'Seizing','FontSize',20)
                    case 2; i=2;r=r9;j=9; figure('Color','w','Position',[80 2 788   732]); hold on; plot_matrix2(T{i},X{i},os,'.-','MarkerSize',6);[h1, h2]=add_stimseiz2(j); legend([h1(1) h2(1)],'Stimulation','Seizure')
                        text(2,41,'Healthy','FontSize',20); text(15,41,'Latent','FontSize',20); text(29,41,'Seizing','FontSize',20)
                    case 3; i=3;r=r10;j=10; figure('Color','w','Position',[80 2 788   732]); hold on; plot_matrix2(T{i},X{i},os,'.-','MarkerSize',6);h=add_stimseiz2(j); legend('1','2','3','4','5','6','7','8')
                    case 4; i=4;r=r1;j=1; figure('Color','w','Position',[80 2 788   732]); hold on; plot_matrix2(T{i},X{i},os,'.-','MarkerSize',6);h=add_stimseiz2(j); legend('1','2','3','4','5','6','7','8')
%                 i=2;r=r9;j=9; hold on; plot_matrix2(T{i},X{i},os,'m.','LineWidth',1,'MarkerSize',6);add_stimseiz2(j)
                end
            elseif analyze_only_seizing == 2
                    i=8;r=r6;j=6; figure('Color','w','Position',[80 2 788   732]); hold on; plot_matrix2(T{i},X{i},os,'.','MarkerSize',6);add_stimseiz2(j); legend('1','2','3','4','5','6','7','8')
                    i=9;r=r7;j=7; figure('Color','w','Position',[80 2 788   732]); hold on; plot_matrix2(T{i},X{i},os,'.','MarkerSize',6);add_stimseiz2(j); legend('1','2','3','4','5','6','7','8')
            else
                %i=5;r=r5;j=5; figure('Color','w','Position',[80 2 788   732]); hold on; plot_matrix2(T{i},X{i},os,'k.','LineWidth',1,'MarkerSize',6);add_stimseiz2(j); legend('1','2','3','4','5','6','7','8')
                i=6;r=r8;j=8; figure('Color','w','Position',[80 2 788   732]); hold on; plot_matrix2(T{i},X{i},os,'.','LineWidth',1,'MarkerSize',6);add_stimseiz2(j); legend('1','2','3','4','5','6','7','8')
                i=7;r=r11;j=11; figure('Color','w','Position',[80 2 788   732]); hold on; plot_matrix2(T{i},X{i},os,'.','LineWidth',1,'MarkerSize',6);add_stimseiz2(j); legend('1','2','3','4','5','6','7','8')
            end
            xlabel('Time (days)');
        end
        
        if plot_movingphase
%             fhandle = @(x,y) xcorr(x,y,0,'coeff');      % Function handel for cross correlation with zero lag
%             i=1; cross_ts_func(T{i}, X{i},2,0.9,0.5,fhandle);
            
            fhandle = @(t,y) get_phase_only(t,y,2*pi/1.0,0.2,show_amp);      % Function handel for cross correlation with zero lag
                if analyze_only_seizing
                    i=1;j=4; [Tout Xout] = calc_ts_func(T{i}, X{i},2,0.5,0.5,fhandle); figure(i); set(gcf,'Color','w','Position',[80 2 788   732]); plot_ts_func(T{i},X{i},Tout,Xout,j);
                    i=2;j=9; [Tout Xout] = calc_ts_func(T{i}, X{i},2,0.5,0.5,fhandle); figure(i); set(gcf,'Color','w','Position',[80 2 788   732]); plot_ts_func(T{i},X{i},Tout,Xout,j);
                    i=3;j=10; [Tout Xout] = calc_ts_func(T{i}, X{i},2,0.5,0.5,fhandle); figure(i); set(gcf,'Color','w','Position',[80 2 788   732]); plot_ts_func(T{i},X{i},Tout,Xout,j);
                    i=4;j=1; [Tout Xout] = calc_ts_func(T{i}, X{i},2,0.5,0.5,fhandle); figure(i); set(gcf,'Color','w','Position',[80 2 788   732]); plot_ts_func(T{i},X{i},Tout,Xout,j);
                else
%                     i=5;j=5; [Tout Xout] = calc_ts_func(T{i}, X{i},2,0.5,0.5,fhandle); figure(i); set(gcf,'Color','w','Position',[80 2 788   732]); plot_ts_func(T{i},X{i},Tout,Xout,j);
                    i=6;j=8; [Tout Xout] = calc_ts_func(T{i}, X{i},2,0.5,0.5,fhandle); figure(i); set(gcf,'Color','w','Position',[80 2 788   732]); plot_ts_func(T{i},X{i},Tout,Xout,j);
                    i=7;j=11; [Tout Xout] = calc_ts_func(T{i}, X{i},2,0.5,0.5,fhandle); figure(i); set(gcf,'Color','w','Position',[80 2 788   732]); plot_ts_func(T{i},X{i},Tout,Xout,j);
                end
        end
       
    end
    
    

    if plot_corrcoef_EEG_vs_EEG
 
        %func_plot_corrcoef_EEG_bands_old(T,X,phi,prepostchronic);
            % It is necessary to pass prepostchronic in order to pull out
            % the correct portion of the phi array
        func_plot_corrcoef_EEG_bands(r_seiz,r_ns,r_ctrl,analyze_only_seizing,prepostchronic,smoothmode_ts,phi,rat2plot)
        
    end
    
    if plot_princomp

        
        if analyze_only_seizing == 1
            rcell_curr = r_seiz;
        elseif analyze_only_seizing == 2
            rcell_curr = r_ctrl;
        else
            rcell_curr = r_ns;
        end
        
        func_plot_PCA(rcell_curr, prepostchronic, smoothmode_ts)
        
        
    end
    
    if plot_princomp2
        
        
        % Index format: (rat number, pre/acute/chronic)
        clear T X
        for k=1:3 % 0 for all; 1 for pre; 2 for latency; 3 for chronic
            i=1;r=r4; [T{i,k} X{i,k}] = extract_ts(r, k); 
            i=2;r=r9; [T{i,k} X{i,k}] = extract_ts(r, k); 
            i=3;r=r10; [T{i,k} X{i,k}] = extract_ts(r, k); 
            i=4;r=r1; [T{i,k} X{i,k}] = extract_ts(r, k); 
        end
        
        
        Mcoef = zeros([size(X,1),size(X,2),length(r),length(r)]);
        Mlatent = zeros(size(X,1),size(X,2),length(r));
        princphi = zeros(size(X,1),size(X,2),length(r));
        for i = 1:size(X,1)
            for j = 1:size(X,2)
                [Mcoef(i,j,:,:), Sscore{i,j}, Mlatent(i,j,:)] = pca(X{i,j});
                for k = 1:size(Sscore{i,j},2)
                    dat{i,j,k} = cosinor_struct(T{i,j}(:,1),Sscore{i,j}(:,k),w,alpha,0);
                    princphi(i,j,k) = dat{i,j,k}.phi;
                    princpass(i,j,k) = dat{i,j,k}.p_3a;
                end
            end
        end
        
        
        princphi = shift_thresholds(princphi,thresh);
        for j = 1:3
           figure; imagesc(squeeze(princphi(:,j,1:3))); colormap winter; colorbar; set(gca,'YDir','normal');ylabel('Rat num');xlabel('Eigenmode');
        end
        
        
        Mcoef_temp = (Mcoef).^2;
        for k=1:2   % Mode
            figure
            for j=1:3 % Stage
                subplot(1,3,j); plot(squeeze(Mcoef_temp(:,j,:,k))')
                hold on; plot(squeeze(mean(Mcoef_temp(:,j,:,k),1)),'r','LineWidth',2)
            end
            
        end
        
    end
    
    
    
    if plot_princomp_amp_vs_ampraw
        
        
        if analyze_only_seizing == 1
            func_plot_princomp_amp_vs_ampraw(r_seiz,prepostchronic,smoothmode_ts);
        elseif analyze_only_seizing == 2
            func_plot_princomp_amp_vs_ampraw(r_ctrl,prepostchronic,smoothmode_ts);
        else
            func_plot_princomp_amp_vs_ampraw(r_ns,prepostchronic,smoothmode_ts);
        end
        
        %func_plot_princomp_amp_vs_ampraw_old(r4,r9,r10,r1,r8,r11,r6,r7,prepostchronic,smoothmode_ts,analyze_only_seizing);
    end
    
    
    if plot_corrcoef_all
        
        if analyze_only_seizing == 1
            func_plot_corrcoef_all(r_seiz,prepostchronic,smoothmode_ts,reload_spks);
        elseif analyze_only_seizing == 2
            func_plot_corrcoef_all(r_ctrl,prepostchronic,smoothmode_ts,reload_spks);
        else
            func_plot_corrcoef_all(r_ns,prepostchronic,smoothmode_ts,reload_spks);
        end

    end
    

    
    if plot_EEG_vs_spks2
        func_plot_EEG_vs_spks2(r_seiz, prepostchronic,smoothmode_ts,reload_spks);
        
    end

    
    
    
    if plot_EEG_vs_baseline_corr
        
        if analyze_only_seizing == 1
            rcell_curr = r_seiz;
        elseif analyze_only_seizing == 2
            rcell_curr = r_ctrl;
        else
            rcell_curr = r_ns;
        end
        
        func_plot_EEG_vs_baseline_corr(rcell_curr,prepostchronic, smoothmode_ts)
    end
    
    
    
    if plot_corr_phaseshift

        
        func_plot_corr_phaseshift(r_seiz,r_ns,r_ctrl,analyze_only_seizing,smoothmode_ts);
        
    end
    
    if plot_phasefits_allbands || plot_phasefits_shiftbands
        
        if analyze_only_seizing == 1
            rcell_curr = r_seiz(rat2plot);
        elseif analyze_only_seizing == 2
            rcell_curr = r_ctrl;
        else
            rcell_curr = r_ns;
        end
        
        if plot_phasefits_shiftbands
            figure('Color','w','Position',[ 622   544   768   240]);
            freq_range = [4];
            subplot(131); xlabel('Latent')
        else
            figure('Color','w','Position',[1500 10  757 800]);
            freq_range = [1:3 8];
        end
        func_plot_phasefits(rcell_curr,analyze_only_seizing,freq_range) 
        if plot_phasefits_allbands
            annotation(gcf,'textbox',...
                [0.0131730515191546 0.48625 0.0889392338177015 0.08375],...
                'String',{'Freq bands'},...
                'FontSize',16,...
                'FitBoxToText','off');
            
            % Create arrow
            annotation(gcf,'arrow',[0.0568031704095112 0.0581241743725231],...
                [0.58025 0.82125]);

            % Create arrow
            annotation(gcf,'arrow',[0.0568031704095112 0.0568031704095112],...
                [0.4715 0.24125]);

        end
    end

    if ergodic_mode
        [erg] = ergodic_extract(r1,r4,r9,r10,r5,r8,r11,analyze_only_seizing);
        [phicell Ampcell] = erg_build (erg, ratrange0,analyze_only_seizing);
        
        phimean = 24*meanCell(phicell); phiste = 24*stdCell(phicell) ./ countCell(phicell); phistd = 24*stdCell(phicell); phiconfid = 24*confidCell(phicell,0.05);
        figure('color','w','Position', [90   272   934   382]);
        %figure(1); hold on;
        bar(phimean);
        colormap(gray)
        legend('Pre','Post-L','Post-SS')
        set(gca,'FontSize',30,'Box','off')
        set(gca,'YTick',[0:6:30])
        ylim([0 26]);
        
        offset = 0.22;
        hold on; errorbar([(1:size(phimean,1))-offset; (1:size(phimean,1))+0; (1:size(phimean,1))+offset]',phimean, phiconfid ,'Color',[0.0 0.0 0.0],'Marker','none','LineStyle','none','LineWidth',1)
        for i = 1:size(phicell,1)
            for j = 1:size(phicell,2)
                %hold on; plot( repmat(i+(j-2)*offset,length(phicell{i,j}),1) ,24*phicell{i,j}','.')
            end
        end


        % Phase shift in each frequency band
        N = size(phi,1);
        harr = [zeros(1,N)];
        hrank = [zeros(1,N)];
        parr = [zeros(1,N)];
        prank = [zeros(1,N)];
        for i = 1:N
            [harr(i) parr(i)] = ttest2(phicell{i,1},phicell{i,2});
            [prank(i) hrank(i)] = ranksum(phicell{i,1},phicell{i,2});
%             [prank(i) hrank(i)] = signrank(abs(phi(i,:,1)-phi(i,:,2)));
            if hrank(i)
                text(i-0.07,phimean(i,2)+phiconfid(i,2)+1,'*','FontSize',30)
            end
        end
        harr
        hrank
        parr
        prank
        clear N
        
        if analyze_only_seizing
            N = size(phi,1);
            harr = [zeros(1,N)];
            hrank = [zeros(1,N)];
            parr = [zeros(1,N)];
            prank = [zeros(1,N)];
            for i = 1:N
                [harr(i) parr(i)] = ttest2(phicell{i,1},phicell{i,3});
                [prank(i) hrank(i)] = ranksum(phicell{i,1},phicell{i,3});
    %             [prank(i) hrank(i)] = signrank(abs(phi(i,:,1)-phi(i,:,2)));
                if hrank(i)
                    text(i-0.07+offset,phimean(i,2)+phiconfid(i,2)+1,'*','FontSize',30)
                end
            end
            harr
            hrank
            parr
            prank
            clear N
        end

        % Phase shift relative to theta band (for control only)
        N = size(phi,1);
        harr = [zeros(1,N)];
        hrank = [zeros(1,N)];
        parr = [zeros(1,N)];
        prank = [zeros(1,N)];
        for i = 1:N
            [harr(i) parr(i)] = ttest2(phicell{2,1},phicell{i,1});
            [prank(i) hrank(i)] = ranksum(phicell{2,1},phicell{i,1});
            %[prank(i) hrank(i)] = signrank(abs(phi(2,:,1)-phi(i,:,1)));
                if hrank(i)
                    %text(i-0.07-offset,phimean(i,1)+phiconfid(i,1)+1,'#','FontSize',20)
                end
        end
        harr
        hrank
        parr
        prank
        clear N

                Ampmean = meanCell(Ampcell); Ampste = stdCell(Ampcell) ./ countCell(Ampcell); Ampstd = stdCell(Ampcell); Ampconfid = confidCell(Ampcell,0.05);
                
                
        if normalize_amps
            if ~percent_difference
                Ampmean = Ampmean ./ repmat(Ampmean(:,1),[1 3]);
                Ampconfid = Ampconfid ./ repmat(Ampconfid(:,1),[1 3]);
            else
                Ampmean = (Ampmean - repmat(Ampmean(:,1),[1 3])) ./ abs(repmat(Ampmean(:,1),[1 3]));
                Ampconfid = (Ampconfid) ./ abs(repmat(Ampmean(:,1),[1 3]));
            end
        end
        
        
        
        figure('color','w','Position',[ 445   407   943   292]);
        %figure(3); hold on;
        range=1:8; Ampmean_pl = Ampmean(range,:,:); Ampconfid_pl = Ampconfid(range,:,:); Ampcell_pl = Ampcell(range,:);
        %Ampmean_pl = Ampmean; Ampconfid_pl = Ampconfid; Ampcell_pl = Ampcell;
        for i = 1:size(Ampmean_pl,1)
            subplot(1,size(Ampmean_pl,1),i);hold on; bar(([Ampmean_pl(i,:); [0 0 0]])); xlim([0.5 1.5]); colormap(gray)
        end
        for i = 1:size(Ampmean_pl,1);
            offset=0.22
            subplot(1,size(Ampmean_pl,1),i); errorbar([1-offset 1 1+offset],(Ampmean_pl(i,:)), (Ampconfid_pl(i,:)),'Color',[0.0 0.0 0.0],'Marker','none','LineStyle','none','LineWidth',1);
            [ptemp htemp] = ranksum(Ampcell_pl{i,1},Ampcell_pl{i,2}); if htemp; text(1-0.03,(Ampmean_pl(i,2)+Ampconfid_pl(i,2))*1.1,'*','FontSize',30); end
            [ptemp htemp] = ranksum(Ampcell_pl{i,1},Ampcell_pl{i,3}); if htemp; text(1+offset-0.03,(Ampmean_pl(i,3)+Ampconfid_pl(i,3))*1.1,'*','FontSize',30); end
            set(gca,'XTickLabel',{},'FontSize',20)
        end
        
        for i = 1:size(Ampcell_pl,1)
            for j = 1:size(Ampcell_pl,2)
                %hold on; subplot(1,size(Ampmean_pl,1),i); plot( repmat(j,length(Ampcell_pl{i,j}),1) ,(Ampcell_pl{i,j}'),'.')
            end
        end
        


%         allphis = [];
%         for i = 1:size(phicell,1)
%             for j = 1:size(phicell,2)
%                 allphis = [allphis phicell{i,j}];
%             end
%         end
%         figure; hist(allphis)
%         figure; hist(phicell{1,3},50)   %% Looks Gaussian to me!

    
    end

    
    
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Supporting functions %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
function erg = ergodic_extract(r1,r4,r9,r10,r5,r8,r11,analyze_only_seizing)

    fhandle = @(t,y) cosinor_struct(t,y,2*pi/1.0,0.2);      % Function handel for cross correlation with zero lag

    for i = 1:length(r4)

        binsize = 2;
        binoverlap = 0.0;

        if analyze_only_seizing
            for j = 1:4
                erg{j}{i,1} = [];
                erg{j}{i,2} = [];
                erg{j}{i,3} = [];
            end

            j=1; erg{j} = calc_erg(erg{j},r4,i,binsize,binoverlap,fhandle);
            j=2; erg{j} = calc_erg(erg{j},r9,i,binsize,binoverlap,fhandle);
            j=3; erg{j} = calc_erg(erg{j},r10,i,binsize,binoverlap,fhandle);
            j=4; erg{j} = calc_erg(erg{j},r1,i,binsize,binoverlap,fhandle);
        else
            for j = 1:2
                erg{j}{i,1} = [];
                erg{j}{i,2} = [];
                erg{j}{i,3} = [];
            end
            j=1; erg{j} = calc_erg_nochr(erg{j},r8,i,binsize,binoverlap,fhandle);
            j=2; erg{j} = calc_erg_nochr(erg{j},r11,i,binsize,binoverlap,fhandle);
            %j=3; erg{j} = calc_erg_nochr(erg{j},r5,i,binsize,binoverlap,fhandle);
            %erg = calc_erg_nochr(erg,r12,i,binsize,binoverlap,fhandle);
        end

%             
%             r=r1;
%             [Tout temp] = calc_ts_func(r{i}.ctrl.t(:),r{i}.ctrl.d(:),2,0.5,0.5,fhandle); temp = reshape(temp,1,1,length(temp)); erg{i,1} = cat(3,erg{i,1},temp)
%             [Tout temp] = calc_ts_func(r{i}.acute.t(:),r{i}.acute.d(:),2,0.5,0.5,fhandle); temp = reshape(temp,1,1,length(temp)); erg{i,2} = cat(3,erg{i,2},temp)
%             [Tout temp] = calc_ts_func(r{i}.chr.t(:),r{i}.chr.d(:),2,0.5,0.5,fhandle); temp = reshape(temp,1,1,length(temp)); erg{i,3} = cat(3,erg{i,3},temp)

    end
    
end

function [phicell Ampcell] = erg_build (erg, ratrange0,analyze_only_seizing)

    %Pull values out of cosinor structure using cell/matrix operations
    ergphi = []; ergAmp = [];
    ratrange = ratrange0;
    if ratrange0 < 0; ratrange = 1:length(erg); end
    for j = ratrange
        [ergphi_temp ergAmp_temp] = erg_unpack (erg{j});

        thresh = 0.2;
        if (analyze_only_seizing==0 && j==3);
            thresh = 0.0;
        end 
        ergphi_temp = shift_thresholds(ergphi_temp,thresh);

        ergphi = cat(2,ergphi,ergphi_temp);
        ergAmp = cat(2,ergAmp,ergAmp_temp);
    end



    % Plot bargraph stuff
    phicell = remove_NaN (ergphi);
    Ampcell = remove_NaN (ergAmp);
end


function [Amp phi] = getamptheta(beta, gamma)

    %Calculate amplitude and acrophase from beta and gamma
    Amp = sqrt(beta.^2 + gamma.^2);
    theta = atan(abs(gamma./beta));

    % Calculate acrophase (phi) and convert from radians to degrees
    a = sign(beta);
    b = sign(gamma);
    
    phi = theta;
    index = (a == 1 | a == 0) & b == 1;
    phi(index) = -theta(index);
    
    index = a == -1 & (b == 1 | b == 0);
    phi(index) = -pi + theta(index);
    
    index = (a == -1 | a == 0) & b == -1;
    phi(index) = -pi - theta(index);
    
    index = a == 1 & (b == -1 | b == 0);
    phi(index) = -2*pi + theta(index);
    
end



function CI = shift_thresholds_confid(CI,phi,thresh);

    CI = abs(CI/2/pi);
    phi = abs(phi/2/pi);
    CI(CI < thresh) = CI(CI < thresh) + 1;
    
end

function plot_confid_int(r)

    w = 2*pi/1.0;
    alpha = .05
    
    for i = 1:length(r)
        figure; hold on;
        cosinor_struct_plot(r{i}.ctrl.t,r{i}.ctrl.d,w,alpha,'b');
        cosinor_struct_plot(r{i}.acute.t,r{i}.acute.d,w,alpha,'r');
        cosinor_struct_plot(r{i}.chr.t,r{i}.chr.d,w,alpha,'y');
    end
end





function erg = calc_erg(erg,r,i,binsize,binoverlap,fhandle)
    struct_empty = struct('phi',NaN,'CI_phi_min',NaN,'CI_phi_max',NaN,'RNE',NaN,'p_3a',NaN,'Amp',NaN);
    struct_empty.RNE = repmat(NaN,3,4);

    [Tout temp1] = calc_ts_func(r{i}.ctrl.t(:),r{i}.ctrl.d(:),binsize,binoverlap,0.5,fhandle); temp1 = reshape(temp1,1,1,length(temp1));
    [Tout temp2] = calc_ts_func(r{i}.acute.t(:),r{i}.acute.d(:),binsize,binoverlap,0.5,fhandle); temp2 = reshape(temp2,1,1,length(temp2));
    [Tout temp3] = calc_ts_func(r{i}.chr.t(:),r{i}.chr.d(:),binsize,binoverlap,0.5,fhandle); temp3 = reshape(temp3,1,1,length(temp3));
    N = max([length(temp1) length(temp2) length(temp3)]);
    temp1 = cat(3,temp1,repmat(struct_empty,[1,1,N-length(temp1)])); erg{i,1} = cat(3,erg{i,1},temp1);
    temp2 = cat(3,temp2,repmat(struct_empty,[1,1,N-length(temp2)])); erg{i,2} = cat(3,erg{i,2},temp2);
    temp3 = cat(3,temp3,repmat(struct_empty,[1,1,N-length(temp3)])); erg{i,3} = cat(3,erg{i,3},temp3);            
end


function erg = calc_erg_nochr(erg,r,i,binsize,binoverlap,fhandle)
    struct_empty = struct('phi',NaN,'CI_phi_min',NaN,'CI_phi_max',NaN,'RNE',NaN,'p_3a',NaN,'Amp',NaN);
    struct_empty.RNE = repmat(NaN,3,4);

    [Tout temp1] = calc_ts_func(r{i}.ctrl.t(:),r{i}.ctrl.d(:),binsize,binoverlap,0.5,fhandle); temp1 = reshape(temp1,1,1,length(temp1));
    [Tout temp2] = calc_ts_func(r{i}.acute.t(:),r{i}.acute.d(:),binsize,binoverlap,0.5,fhandle); temp2 = reshape(temp2,1,1,length(temp2));
    
    N = max([length(temp1) length(temp2)]);
    temp1 = cat(3,temp1,repmat(struct_empty,[1,1,N-length(temp1)])); erg{i,1} = cat(3,erg{i,1},temp1);
    temp2 = cat(3,temp2,repmat(struct_empty,[1,1,N-length(temp2)])); erg{i,2} = cat(3,erg{i,2},temp2);
    temp3 = repmat(struct_empty,[1,1,N]); erg{i,3} = cat(3,erg{i,3},temp3);
    
end



function [Tout Xout] = cross_ts_func(T,X,time_bin,fract_overlap,fract_maxgap,fhandle)
    
    plot_debug =1;

    fract_shift = 1.0 - fract_overlap;
    shift = time_bin*fract_shift;
    
    t = T(:,1);
    bin_min = min(t):shift:(max(t)-time_bin);
    bin_max = bin_min + time_bin;

    dt = median(diff(t));
    
    Tout = zeros(length(bin_min),size(T,2));
    Xout = zeros(length(bin_min),size(T,2));
    good_data = logical(zeros(length(bin_min),1));
    for i = 1:length(bin_min)
        index = (t >= bin_min(i) & (t < bin_max(i)));
        
        t_temp = t(index);
        fract_gap = 1.0 - length(t_temp)*dt ./ (time_bin-dt);
        
        if fract_gap <= fract_maxgap
            good_data(i) = 1;
%             for j = 1:size(T,2)-1
%                 x = X(index,j);
%                 y = X(index,j+1);
%                 %Xout(i,j) = fhandle(x,y);
%                 Xout(i,j) = fhandle(x-mean(x),y-mean(y));
%                 Tout(i,j) = mean([bin_min(i) bin_max(i)]);
%             end
            for j = 1:size(T,2)
                x = X(index,2);
                y = X(index,j);
                %Xout(i,j) = fhandle(x,y);
                Xout(i,j) = fhandle(x-mean(x),y-mean(y));
                Tout(i,j) = mean([bin_min(i) bin_max(i)]);
            end
        else
            good_data(i) = 0;
        end        
    end
    
    Xout = Xout(good_data,:);
    Tout = Tout(good_data,:);
    
    if plot_debug 
%         figure;
%         shift = 1.0;
%         for i = 1:size(Xout,2); Xout(:,i) = Xout(:,i) + (shift*(i-1));end
%         plot(Tout,Xout)
        
%         figure;
%         for i = 1:size(X,2)-1
%             subplot(ceil((size(X,2)-1)/2),2,i);
%             hold on; plot(T(:,i),X(:,i),'b');
%             hold on; plot(T(:,i+1),X(:,i+1),'r');
%             hold on; plot(Tout(:,i),Xout(:,i),'k','LineWidth',2);
%             title(['Freq band' num2str(i)]);
%         end
        
        figure;
        for i = 1:size(X,2)
            subplot(ceil((size(X,2)-1)/2),2,i);
            temp = max(abs([X(:,2); X(:,i)]));
            hold on; plot(T(:,i),X(:,2)/temp,'b.','MarkerSize',15);
            hold on; plot(T(:,i),X(:,i)/temp,'r.','MarkerSize',15);
            hold on; plot(Tout(:,i),Xout(:,i),'k.','LineWidth',2,'MarkerSize',25);
            title(['Freq band' num2str(i)]);
            add_stimseiz(9,'k','LineWidth',2);
        end
        
    end
end





function plot_ts_func(T,X,Tout,Xout,ratN)

    global freq_band

    for i = 1:size(X,2)
        subplot(ceil((size(X,2)-1)/2),2,i);
        temp = max(abs([X(:,i)]));
        
        plotemp1_h = @(t,y) plot(t,y,'r.','MarkerSize',7);
%         plotemp2 = @(t,y) plot(t,y,'k.','MarkerSize',25);
        plotemp2_h = @(t,y) plotemp2(t,y);
        
        hold on; [ax ah1 ah2] = plotyy(T(:,i),X(:,i)/temp,Tout(:,i),Xout(:,i),plotemp1_h,plotemp2_h);
        %hold on; [ax ah1 ah2] = plotyy(Tout(:,i),Xout(:,i),T(:,i),X(:,i)/temp,plotemp2_h,plotemp1_h);
        xtemp = get(gca,'XLim'); ytemp = get(gca,'YLim');
        set(gca,'FontSize',20)
        set(gca,'YTick',-1:0.5:1.0)
        set(gca,'Box','off')
        %axis([xtemp(1) xtemp(2) ytemp(1)+0.2 ytemp(2)+0.2]);
        axes(ax(2));
        axis([xtemp(1) xtemp(2) ytemp(1)+0.2 ytemp(2)+0.2]);
        
        
%         hold on; plot(T(:,i),X(:,i)/temp,'r.','MarkerSize',15);
%         hold on; plot(Tout(:,i),Xout(:,i),'k.','LineWidth',2,'MarkerSize',25);
        title([num2str(freq_band(i,1)) '-' num2str(freq_band(i,2)) 'Hz']);
        %add_stimseiz(ratN,'k','LineWidth',2);
        
        add_stimseiz2(ratN)
    end
    
    function h = plotemp2 (t,y) 
        h = plot(t,y,'kx','MarkerSize',10,'LineWidth',2);
        %set(gca,'YLim',[-1.0 1.0])
        ylim([-1 1.0])
        set(gca,'FontSize',20,'Box','off')
        
    end 
  
end

function [h1, h2] = add_stimseiz2(ratN)
    [stim seiz] = get_stimtime2(ratN);
    h1 = add_vert_bars([stim],'k','LineWidth',2)

    h2 = [];
    ratnum = num2str(ratN, '%6.3d');
    if (ratN == 4) || (ratN == 9) || (ratN == 10) || (ratN == 1) 
        sz_days = get_seizure_times (ratN,['../data/Disk_inventory/Seizures/R' ratnum 'FileInfo.mat']);
        if ratN == 10; sz_days = sz_days(sz_days < 46); end
        h2 = add_vert_bars([sz_days],'r','LineWidth',2);
    end
end

function h = add_vert_bars(x,varargin)
        % Make x column
    x = x(:);
    x = x';     
    temp = get(gca,'YLim');
    xmin=temp(1);xmax=temp(2);
    hold on; h = plot([x; x],repmat([xmin; xmax],1,length(x)),varargin{:});
    %for i = min(min(T)):max(max(T)); hold on; plot([i i],[xmin xmax],'k:','LineWidth',2); end
end



function r = cellpack (rarr)
    for i = 1:size(rarr.all.d,2)
        r{i} = rarr;
        r{i}.ctrl.d = rarr.ctrl.d(:,i)';
        r{i}.acute.d = rarr.acute.d(:,i)';
        if ~isempty(rarr.chr.d) r{i}.chr.d = rarr.chr.d(:,i)';
        else r{i}.chr.d = []; end
        r{i}.all.d = rarr.all.d(:,i)';
        
        r{i}.ctrl.t = rarr.ctrl.t';
        r{i}.acute.t = rarr.acute.t';
        if ~isempty(rarr.chr.t) r{i}.chr.t = rarr.chr.t';
        else r{i}.chr.t = []; end
        r{i}.all.t = rarr.all.t';
    end
end



function [ergphi ergAmp] = erg_unpack (erg);
    %Pull values out of cosinor structure using cell/matrix operations
    ergmat = cell2mat(erg);
    ergmat = permute(ergmat,[1 3 2]);
    ergcell = struct2cell(ergmat);
    ergphi = squeeze(cell2mat(ergcell(1,:,:,:)));
    ergCI_phi_min = squeeze(cell2mat(ergcell(2,:,:,:)));
    ergCI_phi_max = squeeze(cell2mat(ergcell(3,:,:,:)));

        ergcell_temp = reshape(ergcell, [1 size(ergcell,1) size(ergcell,2) size(ergcell,3) size(ergcell,4)]);   % Add an extra dimension to ergcell
    ergRNE = (cell2mat(ergcell_temp(1,4,:,:,:)));                                                          % Drop RNE matrix into these dimensions
    ergM = squeeze(ergRNE(1,4,:,:,:)); ergbeta = squeeze(ergRNE(2,4,:,:,:)); erggamma = squeeze(ergRNE(3,4,:,:,:));       % Extract

    ergp_3a = squeeze(cell2mat(ergcell(5,:,:,:)));
    ergAmp = squeeze(cell2mat(ergcell(6,:,:,:)));

end



function y = remove_NaN (x)

    N = size(x,1);
    M = size(x,3);
    
    y = {[]};
    y=repmat(y,N,M);
    
    for i = 1:N
        for j = 1:M
            xtemp = x(i,:,j);
            y{i,j} = xtemp(~isnan(xtemp));
        end
    end
end

function y = meanCell(x)
    
    N = size(x,1);
    M = size(x,2);
    
    for i = 1:N
        for j = 1:M
            y(i,j) = mean(x{i,j});
        end
    end

end



function y = medianCell(x)
    
    N = size(x,1);
    M = size(x,2);
    
    for i = 1:N
        for j = 1:M
            y(i,j) = median(x{i,j});
        end
    end

end


function y = stdCell(x)
    
    N = size(x,1);
    M = size(x,2);
    
    for i = 1:N
        for j = 1:M
            y(i,j) = std(x{i,j});
        end
    end

end


function y = confidCell(x,alpha)
    
    N = size(x,1);
    M = size(x,2);
    
    for i = 1:N
        for j = 1:M
            y(i,j) = confidence(x{i,j},alpha);
        end
    end

end

function y = countCell(x)
    
    N = size(x,1);
    M = size(x,2);
    
    for i = 1:N
        for j = 1:M
            y(i,j) = length(x{i,j});
        end
    end

end



function Amp_raw = calculate_Ampraw(r)

    Amp_raw = [];
    for prepostchronic = 1:3
        [T X] = extract_ts(r, prepostchronic,1,0); 
        temp = mean(X);
        Amp_raw = cat(3,Amp_raw, temp(:));
    end

end


function X = normalize_dimension (X,dim)

    index = zeros(size(X));
    Xtemp = X()

end



function r = deposit_TX (r,T,X)
    
    plot_debug = 0;
    
    if plot_debug
        os.shift = 3;
        figure('Color','w','Position',[80 2 788   732]); hold on;
        Xdef = []; Tdef = [];
        for prepostchronic = 1:3
            [Ttemp Xtemp] = extract_ts (r, prepostchronic, 0, 1);
            Tdef = [Tdef; Ttemp];
            Xdef = [Xdef; Xtemp];
        end
        plot_matrix2(Tdef,Xdef,os,'k','LineWidth',1,'MarkerSize',6);
    end

    for i = 1:length(r)
        r{i}.ctrl.t = T{1}(:,i)'; r{i}.ctrl.d = X{1}(:,i)';
        r{i}.acute.t = T{2}(:,i)'; r{i}.acute.d = X{2}(:,i)';
        r{i}.chr.t = T{3}(:,i)'; r{i}.chr.d = X{3}(:,i)';
    end
    
    if plot_debug
        Xdef = []; Tdef = [];
        for prepostchronic = 1:3
            [Ttemp Xtemp] = extract_ts (r, prepostchronic, 0, 1);
            Tdef = [Tdef; Ttemp];
            Xdef = [Xdef; Xtemp];
        end
        plot_matrix2(Tdef,Xdef,os,'r--','MarkerSize',6);
        legend('Original','Repacked');
    end

end



function y = dummy_STD(t,varargin)
    y = std(varargin{:});
end
