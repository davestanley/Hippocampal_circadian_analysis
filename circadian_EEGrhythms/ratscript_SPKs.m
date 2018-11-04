

function [phi_means data_struct] = ratscript_SPKs (ratN, chanN, pattern,chosen_simversion,theta_or_nontheta)

    if (~exist('ratN','var')); ratN=1; end
    if (~exist('chanN','var')); chanN=2; end
    if (~exist('pattern','var')); pattern=1; end                            % Pattern=1 for downwards, 2 for upwards; although this is swapped for Rat 004
    if (~exist('chosen_simversion','var')); chosen_simversion=3; end        % Version4 = all spikes merged; Version3 = up and downwards separated
    if (~exist('theta_or_nontheta','var')); theta_or_nontheta=0; end        % =1 for only SPKs in theta epochs; =2 for only SPKs in non-theta epochs; <=0 for no separation


    addpath('./cosinor_dav')
    global ratN_gl chanN_gl pattern_gl
    ratN_gl = ratN;
    chanN_gl = chanN;
    pattern_gl = pattern;

    %close all

    FS_axis = 12;
    FS_axis_sm = 12;
    FS_axis_timeseries = 30;
    start_from_scratch=1;
    spect_from_scratch = 0;
    filter_out_bad_data = 1;
    pre_smooth_data = 1;
    use_tcell_search = 0;
    fract_maxgap0 = 0.5;
    merge_type1_type2 = 0;
    plot_cosinor = 1;
    use_normalized_SPK_in_theta_rates = 0;
    short_acute_period = 0;

    plot_on=0;
        shift_1st_seizure_time_to_zero = 1;         % Makes the time of first seizure correspond to t=0 on SPK plot

    basepath{4} = '~/Anco/Evol/EEG/RawOneKSpike/';
    basepath{3} = '/Volumes/MYPASSPORT/data/SortedDave_all/';
    basepath{2} = '/Volumes/MYPASSPORT/data/SortedDave2/';
    basepath{1} = '/Volumes/MYPASSPORT/data/SpikeDetectionAndSorting/Output/';
    basepath = basepath{chosen_simversion};
    
    basepath_sachintimelogs = ('~/Nex/PhD/Evol_to_epil/Data/Disk_inventory/Sachin_Timelogs/');
    %path_ratlog = ('~/Anco/Evol/EEG/RatData_out');
    % path_ratlog='../../RatData_out';
    path_ratlog='../Disk_inventory/RatData_out';
    
    %path_ratlog='../../RatData_out';
    path_savephase = './Phase_Roll';
    
    if ispc; userdir= getenv('USERPROFILE'); 
    else; userdir= getenv('HOME'); 
    end

    clc
    format compact


    ratnum = num2str(ratN, '%6.3d');
    channum = num2str(chanN, '%6.2d');

    if ratnum == '009'; 
        ii = 1;
        for i=[1:32]; channum_arr{ii}=num2str(i,'%05.2d'); ii=ii+1; end
        stim_time = 28;
        first_seizure = 73;
    end
    if ratnum == '004';
        ii = 1;
        for i=[1:32]; channum_arr{ii}=num2str(i,'%05.2d'); ii=ii+1; end
        stim_time = 28;
        first_seizure = 99;
    end
    if ratnum == '006';
        ii = 1;
        for i=[1:32]; channum_arr{ii}=num2str(i,'%05.2d'); ii=ii+1; end
        stim_time = 15;
        first_seizure = 46;
    end
    if ratnum == '007';
        ii = 1;
        for i=[1:32]; channum_arr{ii}=num2str(i,'%05.2d'); ii=ii+1; end
        stim_time = 22;
        first_seizure = 106;
    end
    if ratnum == '001';
        ii = 1;
        for i=[1:8 9 10:16]; channum_arr{ii}=num2str(i,'%05.2d'); ii=ii+1; end
        % if ~strcmp(plot_type, 'skew'); clear channum_arr; ii = 1; for i=[1:8 10:16]; channum_arr{ii}=num2str(i,'%05.2d'); ii=ii+1; end; end % We're missing channel 9 for stdev and kurt for rat 001
        stim_time = 78;
        first_seizure = 109;
    end
    if ratnum == '013'; 
        ii = 1;
            % We're missing channels 25 and 28
        for i=[1:24 26:27 29:32]; channum_arr{ii}=num2str(i,'%05.2d'); ii=ii+1; end
        stim_time = 38;
        first_seizure = 75; % Note Rat 13 actually doesn't exhibit seizures; just stuck this in to make code work.
    end
    if ratnum == '010'; 
        ii = 1;
        for i=[1:32]; channum_arr{ii}=num2str(i,'%05.2d'); ii=ii+1; end
        stim_time = 43;
        first_seizure = 75;
    end
    if ratnum == '005';
        ii = 1;
        for i=[1:32]; channum_arr{ii}=num2str(i,'%05.2d'); ii=ii+1; end
        stim_time = 11;
        %first_seizure = 113;
        first_seizure = 40; % Data becomes crappy at this point
    end
    if ratnum == '008';
        ii = 1;
        for i=[1:32]; channum_arr{ii}=num2str(i,'%05.2d'); ii=ii+1; end
        stim_time = 15;
        first_seizure = 128;
        first_seizure = 130;
    end
    if ratnum == '011';
        ii = 1;
        for i=[1:32]; channum_arr{ii}=num2str(i,'%05.2d'); ii=ii+1; end
        stim_time = 24;
        first_seizure = 87;
    end
    if ratnum == '012';
        ii = 1;
        for i=[1:32]; channum_arr{ii}=num2str(i,'%05.2d'); ii=ii+1; end
        stim_time = 30;
        first_seizure = 76;
        first_seizure = 78;
    end

    
    outlog_path = ['./Ratoutmat_FFT'];
    if start_from_scratch
        [outlog_path outlog_name] = save_log(path_ratlog, ratnum, outlog_path);   % Save log file
    else
        outlog_name = strcat('Rat',ratnum,'log');
    end
    load ([outlog_path '/' outlog_name])    % Load log file

    for ii = 1:length(fileNames);
        fileNums(ii) = str2num(fileNames{ii}(1:4));
    end
    %tabs_old = tabs;
    
%     temp = load ([basepath_sachintimelogs 'Rat' ratnum 'File_TimeLog.dat']);
%     fileNums = temp(:,1);
%     tabs = temp(:,2) / 3600 / 24 + floor(tabs_old(1));


    if spect_from_scratch
        if theta_or_nontheta > 0
            FFT_path = ['../RatFFT/03_theta_delta_2sec_bins'];
            ratio_thresh0 = 2.0;
            %load ([FFT_path '/' 'Rat' ratnum 'Ch' channum 'Rat001Ch02full_file2D.mat']);
            load ([FFT_path '/' 'Rat' ratnum 'Ch' channum 'full_file2D_th6_10_del0_6.mat']);
            ttheta = torig;
            
            ratio_serial = theta_serial ./ delta_serial;
            if theta_or_nontheta == 1; data_binary = ratio_serial > ratio_thresh0;
            elseif theta_or_nontheta == 2; data_binary = ratio_serial <= ratio_thresh0;
            else fprintf('Error, incorrect theta_or_nontheta value \n'); data_binary = []; ttheta = [];
            end
        else
            data_binary =[];
            ttheta = [];
            theta_serial = [];
            delta_serial = [];
        end
        
        if chosen_simversion ~= 4
            fprintf (['Preparing to save ' 'Rat' ratnum 'Ch' channum 'Spk' num2str(pattern) 'Ver' num2str(chosen_simversion) 'ThetaVer' num2str(theta_or_nontheta) '.mat' '\n']);
            datapath = [basepath 'Rat' ratnum '/PreStim'];
            filenames = ['/ch' channum '_Final*.mat'];
            [spike_times_arr filename_arr is_upward_spike] = get_spike_times (datapath, filenames,pattern);
            [tpre,rpre, rpre_norm] = get_num_spikes_per_hour (spike_times_arr,0.9,0.4,[ttheta(:) data_binary(:)], theta_serial, delta_serial);

            datapath = [basepath 'Rat' ratnum];
            filenames = ['/ch' channum '_Final*.mat'];
            [spike_times_arr filename_arr is_upward_spike] = get_spike_times (datapath, filenames,pattern);
            [tpost,rpost, rpost_norm] = get_num_spikes_per_hour (spike_times_arr,0.9,0.4,[ttheta(:) data_binary(:)], theta_serial, delta_serial);

            torig = [tpre(:); tpost(:)] / 24;
            dfull = [rpre(:); rpost(:)];
            dfull_norm = [rpre_norm(:); rpost_norm(:)];
            
%             save(['Rat' ratnum 'Ch' channum 'Spk' num2str(pattern) 'Ver' num2str(chosen_simversion) 'ThetaVer' num2str(theta_or_nontheta) '.mat'],'torig','dfull','dfull_norm');
        else
            fprintf ([basepath 'Rat' ratnum '-OneKSpike' '/*ch' channum '*' '\n']);
            datapath = [basepath 'Rat' ratnum '-OneKSpike'];
            filenames = ['/*ch' channum '*'];
            timeslog = [fileNums(:) (tabs(:) - min(floor(tabs))) ];
            [spike_times_arr filename_arr waveforms] = get_spike_times_unsorted (datapath, filenames,timeslog);
            [torig, dfull] = get_num_spikes_per_hour_unsorted (spike_times_arr,0.9,0.4,[ttheta(:) data_binary(:)]);
            torig = torig / 24;
            
%             save(['Rat' ratnum 'Ch' channum 'SpkAll' 'Ver' num2str(chosen_simversion) 'ThetaVer' num2str(theta_or_nontheta) '.mat'],'torig','dfull');
        end
        
        if use_normalized_SPK_in_theta_rates; dfull = dfull_norm; end

    else
        
        if chosen_simversion ~= 4
            
            path_spk_circ = fullfile(userdir,'/Crystalized/Anco/Evol/EEG/SPK_circadian_analysis_dave');
            
            load ([path_spk_circ '/Rat' ratnum 'Ch' channum 'Spk' num2str(pattern) 'Ver' num2str(chosen_simversion) 'ThetaVer' num2str(theta_or_nontheta) '.mat']); if use_normalized_SPK_in_theta_rates; dfull = dfull_norm; end

            if merge_type1_type2
                load ([path_spk_circ '/Rat' ratnum 'Ch' channum 'Spk' num2str(1) 'Ver' num2str(chosen_simversion) 'ThetaVer' num2str(theta_or_nontheta) '.mat']); if use_normalized_SPK_in_theta_rates; dfull = dfull_norm; end
                dfull1 = dfull; torig1 = torig;
                load ([path_spk_circ '/Rat' ratnum 'Ch' channum 'Spk' num2str(2) 'Ver' num2str(chosen_simversion) 'ThetaVer' num2str(theta_or_nontheta) '.mat']); if use_normalized_SPK_in_theta_rates; dfull = dfull_norm; end
                dfull2 = dfull; torig2 = torig;
        %         if max(abs(torig1-torig2)) > 0.1; fprintf('Torig difference too large! \n'); return; end
                if length(dfull1) >= length(dfull2)
                    dfull = dfull1(1:length(dfull2))+dfull2;
                    torig = torig2;
                else
                    dfull = dfull1+dfull2(1:length(dfull1));
                    torig = torig1;
                end
                if ratN == 4
                    figure; plot(torig1,dfull1,'r.'); hold on; plot(torig2,dfull2,'g.'); hold on; plot(torig,dfull); legend('up','down','total');
                else
                    figure; plot(torig2,dfull2,'r.'); hold on; plot(torig1,dfull1,'g.'); hold on; plot(torig,dfull); legend('up','down','total');
                end
            end
        else
            %load(['Rat' ratnum 'Ch' channum 'SpkAll' 'Ver' num2str(chosen_simversion) '.mat']);
            load([path_spk_circ '/Rat' ratnum 'Ch' channum 'SpkAll' 'Ver' num2str(chosen_simversion) 'ThetaVer' num2str(theta_or_nontheta) '.mat']); if use_normalized_SPK_in_theta_rates; dfull = dfull_norm; end
        end
    end
    
%     data_struct = [];

    
    if filter_out_bad_data

        index = find((dfull >= 1) .* (torig >= 0) );
%         if ratN == 10 & theta_or_nontheta == 0
%             index = find((dfull >= 50) .* (torig >= 0) ); end
        if ratnum == '010' & pattern == 1 & theta_or_nontheta ~= 1; index = find((dfull >= 40) .* (torig >= 0) ); end
        if ratnum == '011' & pattern == 1; index = find((dfull >= 50) .* (torig >= 0) ); end
            
        tfilt2 = torig(index); dfilt2 = dfull(index);
        
        if plot_on
            figure; plot(torig, dfull,'b.');
            hold on; plot(tfilt2, dfilt2,'r.');
        end
        
        torig = torig(index); dfull = dfull(index); % Remove outliers!
        
%         diff_filt_thresh = 1;
%         dfull_p = abs(diff(dfull));
%         index = find((dfull_p <= diff_filt_thresh*std(dfull_p)) .* (torig(1:end-1) >= 0) );
%         tfilt = torig(index); dfilt = dfull(index);
%         if plot_on
%             figure; plot(torig, dfull,'b.');
%             hold on; plot(tfilt, dfilt,'r.');
%             hold on; plot(torig(1:end-1), dfull_p,'g.');
%             hold on; plot([torig(1) torig(end)], [diff_filt_thresh*std(dfull_p) diff_filt_thresh*std(dfull_p)], 'k:');
%         end
    end
    
    
    
    if filter_out_bad_data
        if ratnum == '004'; index = find ( ~((torig>=23.5) .* (torig<=27)) );
            %torig = torig(index); theta_serial = theta_serial(index); delta_serial = delta_serial(index); % Remove outliers!
        end
        if ratnum == '010'; index = find ( ~((torig<3.5)) .* ~((torig>11) .* (torig<12.5) ) .* ~((torig>18) .* (torig<20) ) .* ~((torig>23) .* (torig<24)) );
%             torig = torig(index); dfull = dfull(index); % Remove outliers!
        end
    end
    
    
%     ratio_serial = theta_serial ./ delta_serial;
%     data_binary = ratio_serial > 1.0;
    
    
    
    if ~isempty(find(diff(fileNums)==0)); fprintf('Warning: Overlapping filenums. Check log! \n'); end
    stim_index = find(fileNums >= stim_time, 1, 'first'); stim_tabs=(tabs(stim_index)-floor(tabs(1)));
    seiz_index = find(fileNums >= first_seizure, 1, 'first'); seiz_tabs=(tabs(seiz_index)-floor(tabs(1)));
    tabs_ctrl_1 = 0;
    tabs_ctrl_2 = stim_tabs;
    tabs_acute_1 = stim_tabs;
    tabs_acute_2 = seiz_tabs;
    tabs_chr_1 = seiz_tabs;
    tabs_chr_2 = 100;
    
    
    if ratnum == '009'; tabs_acute_1 = stim_tabs + 2; end
    if ratnum == '004'; tabs_acute_1 = 11; tabs_acute_2 = 30.5; end
    %if ratnum == '010'; tabs_ctrl_1 = 4; end
    if ratnum == '010'; tabs_acute_1 = 19; end
    %if ratnum == '010' & pattern == 1; tabs_chr_2 = 50; end
    if ratnum == '005'; tabs_acute_1 = 5; end
    %if ratnum == '004'; tabs_acute_1 = 16.5; end
%     if ratnum == '001'; tabs_acute_1 = 21.5; end
%     if ratnum == '001'; tabs_acute_2 = 26.5; end
    if ratnum == '008' tabs_acute_1 = 5; tabs_ctrl_2 = 5; end
    
%     if ratnum == '010'; tabs_chr_1 = 40; end        %Final day = 50
%     if ratnum == '009'; tabs_chr_1 = 36; end        %Final day = 46
%     if ratnum == '004'; tabs_chr_1 = 63; end        %Final day = 73
    
    if short_acute_period
        tabs_acute_2 = tabs_acute_1 + 5;
        if ratN == 4;
            %tabs_acute_2 = tabs_acute_1 + 15;
        end
    end

    % Calculate power in theta and delta bands as a function of pre- and
    % latent periods
    index_pre = find ( (tabs_ctrl_1 <= torig) .* (torig <= tabs_ctrl_2) );
    index_acute = find ( (tabs_acute_1 <= torig) .* (torig <= tabs_acute_2) );
    index_chronic = find ( (tabs_chr_1 <= torig) .* (torig <= tabs_chr_2) );
    
    data_struct.rates.SPKs(1) = mean(dfull(index_pre));
    data_struct.rates.SPKs(2) = mean(dfull(index_acute));
    data_struct.rates.SPKs(3) = mean(dfull(index_chronic));
    
    
    
    bin_size = 1/24;
    fract_overlap = 0.0;
    fract_maxgap = fract_maxgap0;
    [times_serial data_serial] = daveMVAVG_bin (torig, dfull, bin_size, fract_overlap, fract_maxgap,use_tcell_search); % Take average time spent in theta state.
%     [times_serial data_serial] = daveMVAVG_bin (torig_theta, theta_in_theta, bin_size, fract_overlap, fract_maxgap,use_tcell_search); % Take average time spent in theta state.
%     [times_serial data_serial] = daveMVAVG_bin (torig_delta, delta_in_delta, bin_size, fract_overlap, fract_maxgap,use_tcell_search); % Take average time spent in theta state.
%     [times_serial theta_serial_filt] = daveMVAVG_bin (torig, theta_serial, bin_size, fract_overlap, fract_maxgap,use_tcell_search); % Take average time spent in theta state.
%     [times_serial delta_serial_filt] = daveMVAVG_bin (torig, delta_serial, bin_size, fract_overlap, fract_maxgap,use_tcell_search); % Take average time spent in theta state.
    times_serial = torig;
    data_serial = dfull;


    ds = abs(downsample(data_serial, 1));
    ts = downsample(times_serial, 1);
    if plot_on
        figure; set (gcf,'Color','w'); 
        hold on; plot(ts, ds);
        xlabel ('Time (days)','FontSize',FS_axis);
        ylabel('SPK rate','FontSize',FS_axis);
        set(gca,'FontSize',FS_axis);
    end

    bin_size = 6/24;
    fract_overlap = 0.9;
%     bin_size = 2/24;
%     fract_overlap = 1-15/120;
    fract_maxgap = fract_maxgap0;
    if pre_smooth_data
        [t_sm dat_sm] = daveMVAVG_bin (ts, ds, bin_size, fract_overlap, fract_maxgap,use_tcell_search);
    else
        t_sm = ts;
        dat_sm = ds;
    end

    bin_size = 1;
    fract_overlap = 0.5;
    fract_maxgap = fract_maxgap0;
    [t_base dat_base] = daveMVAVG_bin (ts, ds, bin_size, fract_overlap, fract_maxgap0,use_tcell_search);


    
    
    if plot_on; hold on; plot(t_sm, dat_sm, 'k.','MarkerSize',18); end

    dat_basei = interp1(t_base,dat_base, t_sm);
    dat_sub = dat_sm - dat_basei;
    if plot_on
        hold on; plot(t_sm, dat_basei, 'g')
        hold on; plot(t_base, dat_base, 'k:');
        hold on; plot([tabs_ctrl_1 tabs_ctrl_1],[0 max(dat_sm)],'k','LineWidth',2);
        hold on; plot([tabs_ctrl_2 tabs_ctrl_2],[0 max(dat_sm)],'k','LineWidth',2);
        hold on; plot([tabs_acute_1 tabs_acute_1],[0 max(dat_sm)],'m','LineWidth',2);
        hold on; plot([tabs_acute_2 tabs_acute_2],[0 max(dat_sm)],'m','LineWidth',2);
        hold on; plot([tabs_chr_1 tabs_chr_1],[0 max(dat_sm)],'r','LineWidth',2);
        hold on; plot([tabs_chr_2 tabs_chr_2],[0 max(dat_sm)],'r','LineWidth',2);
        for kk = 1:max(t_sm);
            hold on; plot([kk kk], [min(dat_sm) max(dat_sm)*2],'k-.');
        end
        legend('original','smoothed','baseline');
        set(gca,'FontSize',FS_axis_timeseries);
    end
    
    if shift_1st_seizure_time_to_zero && plot_on
        first_seiz = 17;
        figure; set (gcf,'Color','w'); 
        hold on; plot(ts-first_seiz, ds);
        xlabel ('Time (days)','FontSize',FS_axis);
        ylabel('SPK rate','FontSize',FS_axis);
        set(gca,'FontSize',FS_axis);
        hold on; plot(t_sm-first_seiz, dat_sm, 'k.','MarkerSize',18); 
        hold on; plot(t_sm-first_seiz, dat_basei, 'g')
        hold on; plot(t_base-first_seiz, dat_base, 'k:');
        hold on; plot([tabs_ctrl_1 tabs_ctrl_1]-first_seiz,[0 max(dat_sm)],'k','LineWidth',2);
        hold on; plot([tabs_ctrl_2 tabs_ctrl_2]-first_seiz,[0 max(dat_sm)],'k','LineWidth',2);
        hold on; plot([tabs_acute_1 tabs_acute_1]-first_seiz,[0 max(dat_sm)],'m','LineWidth',2);
        hold on; plot([tabs_acute_2 tabs_acute_2]-first_seiz,[0 max(dat_sm)],'m','LineWidth',2);
        hold on; plot([tabs_chr_1 tabs_chr_1]-first_seiz,[0 max(dat_sm)],'r','LineWidth',2);
        hold on; plot([tabs_chr_2 tabs_chr_2]-first_seiz,[0 max(dat_sm)],'r','LineWidth',2);
        for kk = 1:max(t_sm);
            hold on; plot([kk kk]-first_seiz, [min(dat_sm) max(dat_sm)*2],'k-.');
        end
        legend('original','smoothed','baseline');
        set(gca,'FontSize',FS_axis_timeseries);
    end

    if plot_on
        figure; plot(t_sm, dat_sub); set (gcf,'Color','w'); 
        hold on; plot([tabs_ctrl_1 tabs_ctrl_1],[0 max(dat_sub)],'k','LineWidth',2);
        hold on; plot([tabs_ctrl_2 tabs_ctrl_2],[0 max(dat_sub)],'k','LineWidth',2);
        hold on; plot([tabs_acute_1 tabs_acute_1],[0 max(dat_sub)],'m','LineWidth',2);
        hold on; plot([tabs_acute_2 tabs_acute_2],[0 max(dat_sub)],'m','LineWidth',2);
        hold on; plot([tabs_chr_1 tabs_chr_1],[0 max(dat_sub)],'r','LineWidth',2);
        hold on; plot([tabs_chr_2 tabs_chr_2],[0 max(dat_sub)],'r','LineWidth',2);
        for kk = 1:max(t_sm);
            hold on; plot([kk kk], [min(dat_sub) max(dat_sub)*2],'k:');
        end
        xlabel('time, (days)','FontSize',FS_axis); ylabel('� EMD Amp','FontSize',FS_axis);
    end

    %     [A phi t_sin] = fit_sinusoids (t_sm, dat_sub, 1.0);
    %     figure;
    %     subplot(211); plot(t_sin, A);
    %     subplot(212); plot(t_sin, phi);



    % % % % %     Average everything
    %     t_smm= mod(t_sm, 1);
    %     [A phi t_sin] = fit_sinusoids (t_smm, dat_sub, 1.0);

    A_means=[];
    phi_means=[];
    A_errs=[];
    phi_errs=[];
    use_rollback = 1;
    % % % % %   Control
    %index = find (t_sm <= stim_tabs);
    index = find ( (tabs_ctrl_1 <= t_sm) .* (t_sm <= tabs_ctrl_2) );
    i=1;
    if ~isempty(index)
        t_section = t_sm(index); dat_section = dat_sub(index);
    %         figure; plot(t_section, dat_section,'o');
        [A_means(i) phi_means(i) A_errs(i) phi_errs(i)] = get_phase (t_section, dat_section, 1.0, use_rollback);
        data_struct.ctrl.t = t_section;
        data_struct.ctrl.d = dat_section;
        if plot_cosinor
            w = 2*pi/1.0;
            alpha = .05;
            index = ~isnan(dat_section);
            [cosphi(i) cosminphi(i) cosmaxphi(i)] = cosinor(t_section(index), dat_section(index),w,alpha);
            %[cMESOR,cAMPL,cPH,cP] = cosinor2(t_section(index), dat_section(index),1.0)
        end
    else
        A_means(i)=0; phi_means(i)=0; A_errs(i)=0; phi_errs(i)=0; cosphi(i)=0; cosminphi(i)=0; cosmaxphi(i)=0;
    end

    % % % % %   Acute (Latent)
    %index = find ( (stim_resume_tabs <= t_sm) .* (t_sm <= seiz_tabs) );
    index = find ( (tabs_acute_1 <= t_sm) .* (t_sm <= tabs_acute_2) );
%     if ratnum == '004'; index = find ( ((tabs_acute_1 <= t_sm) .* (t_sm < 23.5)) | ((t_sm <= tabs_acute_2) .* (t_sm > 27)) ); end
%     if ratnum == '010'; index = find ( (tabs_acute_1 <= t_sm) .* (t_sm <= tabs_acute_2) .* ~((t_sm>18) .* (t_sm<20) ) .* ~((t_sm>23) .* (t_sm<24)) ); end  % Discount isolated data points
    %if ratnum == '001'; index = find ( (21 <= t_sm) .* (t_sm <= tabs_acute_2) ); end
    i=2;
    if ~isempty(index)
        t_section = t_sm(index); dat_section = dat_sub(index);
    %         figure; plot(t_section, dat_section,'o');
        [A_means(i) phi_means(i) A_errs(i) phi_errs(i)] = get_phase (t_section, dat_section, 1.0, use_rollback);
        data_struct.acute.t = t_section;
        data_struct.acute.d = dat_section;
        if plot_cosinor
            w = 2*pi/1.0;
            alpha = .05;
            index = ~isnan(dat_section);
            [cosphi(i) cosminphi(i) cosmaxphi(i)] = cosinor(t_section(index), dat_section(index),w,alpha);
            %[cMESOR,cAMPL,cPH,cP] = cosinor2(t_section(index), dat_section(index),1.0)
        end
    else
        A_means(i)=0; phi_means(i)=0; A_errs(i)=0; phi_errs(i)=0; cosphi(i)=0; cosminphi(i)=0; cosmaxphi(i)=0;
    end



    % % % % %   Chronic
    %index = find (seiz_tabs <= t_sm);        
    index = find ( (tabs_chr_1 <= t_sm) .* (t_sm <= tabs_chr_2) );
    i=3;
    if ~isempty(index)
        t_section = t_sm(index); dat_section = dat_sub(index);
    %         figure; plot(t_section, dat_section,'o');
        [A_means(i) phi_means(i) A_errs(i) phi_errs(i)] = get_phase (t_section, dat_section, 1.0, use_rollback);
        data_struct.chr.t = t_section;
        data_struct.chr.d = dat_section;
        if plot_cosinor
            w = 2*pi/1.0;
            alpha = .05;
            index = ~isnan(dat_section);
            [cosphi(i) cosminphi(i) cosmaxphi(i)] = cosinor(t_section(index), dat_section(index),w,alpha);
        end
    else
        A_means(i)=0; phi_means(i)=0; A_errs(i)=0; phi_errs(i)=0; cosphi(i)=0; cosminphi(i)=0; cosmaxphi(i)=0;
    end
    
    data_struct.all.t = times_serial;
    data_struct.all.d = data_serial;
    data_struct.all.phi = phi_means;    
    data_struct.all.cosphi = cosphi;    
    data_struct.all.cosminphi = cosminphi;    
    data_struct.all.cosmaxphi = cosmaxphi;    
    data_struct.all.A_means = A_means;
    data_struct.all.A_errs = A_errs;
    data_struct.all.phi_errs = phi_errs;
    data_struct.all.stim_tabs = stim_tabs;
    data_struct.all.seiz_tabs = seiz_tabs;
    
    
    if plot_on
        phi_means = phi_means * 24;
        phi_errs = phi_errs * 24;
        for i = 1:length(phi_means)
            if phi_means(i) < 12
                phi_means(i) = phi_means(i) + 24;
            end
        end
        phi_means= phi_means/24;
        phi_errs = phi_errs/24;
        figure;
        colormap (flipud(wkeep(copper(20),[10 3])));
        bar(phi_means);
        %hold on; errorbar(1:length(phi_means),phi_means, phi_errs,'ko')
        ylabel('Circadian phase (days)','FontSize',FS_axis_sm);
        set(gca,'FontSize',FS_axis_sm);
        xlabel_arr = {'PreInjury','Latent','Chronic'};
        set(gca,'XTick',1:length(xlabel_arr))
        set(gca,'XTickLabel',xlabel_arr);


    %         figure;bar(A_means,'w');
    %         hold on; errorbar(1:length(A_means),A_means, A_errs,'ko')
    end

    path_savephase_roll = [path_savephase '_' num2str(use_rollback)];
        [s,mess,messid] = mkdir(path_savephase_roll);
    out_path=strcat(path_savephase_roll,'/','Rat',ratnum);
        [s,mess,messid] = mkdir(out_path);
    out_name = strcat('DataRat',ratnum,'ch',channum,'_','.mat');
    savephase = [out_path '/' out_name];


    save (savephase, 't_sm','dat_sub','stim_tabs','seiz_tabs','A_means','phi_means');
    clear t_sm dat_sub stim_tabs seiz_tabs;

end



function [tabs Fout] = extract_FFT_band (file,freq_band)
    
%     f = struct2matrix(file,'f',[],0,0);
%     F = struct2matrix(file,'F',[],0,0);
%     tabs = struct2matrix(file,'tabs',[],0,0);
    
    f = file.f;
    F = file.F;
    tabs = file.tabs;

    freq_band = freq_band(:);
    fband = zeros(1,size(F,2));
    Fband = zeros(1,size(F,2));
    
    for k=1:size(F,2)
        [fband(k) Fout(k)] = psd_freqbins (f(:,k), F(:,k), freq_band);
    end

end






%%% New functions for SPK analysis

function [spike_times_arr filename_arr is_upward_spike] = get_spike_times (datapath, filenames,spike_num)
    global ratN_gl
    plot_on = 0;
    plot_waveforms = 1;
    
    file_list = dir([datapath filenames]);
    
    is_upward_arr = [];
    spike_times_arr = [];
    
    if isempty(file_list)
        spike_times_arr = [];
        filename_arr = [];
        is_upward_spike = [];
        return;
    end

    for i = 1:length(file_list)
        fprintf (['Loading file ' [datapath '/' file_list(i).name] '\n']);
        load ([datapath '/' file_list(i).name]);


        if length(spikes.Shape) == 2
            temp=spikes.Shape(spike_num).waveforms';
            if ~isempty(temp)
                maxes = max(temp(130:170,:)); mins = min(temp(130:170,:));

                is_upwards = (abs(maxes) > abs(mins));
                if plot_waveforms; figure; plot(temp(:,1:50)); title(['Is upwards = ' num2str(mean(is_upwards))]);end
                is_upward_arr = [is_upward_arr mean(is_upwards)];

                curr_spiketimes = spikes.Shape(spike_num).spiketimes;
                spike_times_arr = [spike_times_arr (curr_spiketimes(:))'];
            end
            filename_arr{i} = file_list(i).name;
        elseif length(spikes.Shape) == 1 && ratN_gl == 8         % Guess SPK type. Only do this for Rat 008 since it's missing so much data
            
            % Calculate direction of the mean trace 2 different ways in
            % order to be sure we have the correct direction estimate.
            temp2 = spikes.Shape(1).Mean; maxes = max(temp2(130:170)); mins = min(temp2(130:170));
            is_upwards_mean = (abs(maxes) > abs(mins));
            
            temp=spikes.Shape(1).waveforms'; maxes = max(temp(130:170,:)); mins = min(temp(130:170,:));
            is_upwards = (abs(maxes) > abs(mins));
            mean_is_upwards = (mean(is_upwards) > 0.5);
            
            
            
            if is_upwards_mean == mean_is_upwards       % Make sure both our estimates match; if not skip.
                if (is_upwards_mean == 1 && spike_num == 2) || ((is_upwards_mean == 0 && spike_num == 1))       % If direction of SPK matches desired spike pattern
                    is_upwards = (abs(maxes) > abs(mins));

                    if plot_waveforms; figure; plot(temp(:,1:50)); title(['Is upwards = ' num2str(mean(is_upwards))]);end
                    is_upward_arr = [is_upward_arr mean(is_upwards)];

                    curr_spiketimes = spikes.Shape(spike_num).spiketimes;
                    spike_times_arr = [spike_times_arr (curr_spiketimes(:))'];
                    filename_arr{i} = file_list(i).name;
                end
            else
                fprintf('Spike shape estimate mismatch. Skipping...\n');
            end
            
        else
            fprintf('Incorrect number of shapes \n');
        end
        

    end
    
    if plot_on
        figure; plot(is_upward_arr); end
    
    is_upward_spike = (mean(is_upward_arr) > 0.5);
    spike_times = 1;
    
end


function [t_center, r, r_normalized] = get_num_spikes_per_hour (spike_times_arr,binsize,overlap,theta_data,theta_serial,delta_serial)

    global ratN_gl chanN_gl pattern_gl
    ratN = ratN_gl;
    chanN = chanN_gl;
    pattern = pattern_gl;
    

    %binsize in hours
    %overlap in fractional
    
    plot_on = 0;
    recalculate_FFT_ratio_from_scratch = 0; % Don't need to do this! The algorithm using pre-extracted data works well
    plot_individual_spks = 0;
    
    if ~isempty(theta_data)
        ttheta = theta_data(:,1);
        data_binary = theta_data(:,2);
        ttheta = ttheta * 24;   % Convert to hours
        ttheta_width = 2.1 / 3600; % 2 seconds width of theta array
    end
    
    spk_hrs = spike_times_arr / 3600;
    
    starttime = min(spk_hrs);
    endtime = max(spk_hrs);
    
    t_starts = [starttime:(binsize*(1-overlap)):(endtime-binsize)];
    t_stops = t_starts + binsize;
    
    if isempty(t_starts)
        t_center = [];
        r = [];
        r_normalized = [];
        return
    end
    
    if isempty(theta_data)
        for i = 1:length(t_starts)
           goods = ( spk_hrs > t_starts(i) & spk_hrs < t_stops(i) );
           t_center(i) = mean([t_starts(i) t_stops(i)]);
           r(i) = sum(goods);
        end
        r_normalized = r;
    else
        for i = 1:length(t_starts)
            fprintf(['Tallying SPKs in hour ' num2str(t_starts(i)) ' to ' num2str(t_stops(i)) '\n']);
            
            ttheta_temp = ttheta(ttheta >= t_starts(i) & ttheta < t_stops(i) );
            data_binary_temp = data_binary(ttheta >= t_starts(i) & ttheta < t_stops(i) );
            theta_serial_temp = theta_serial(ttheta >= t_starts(i) & ttheta < t_stops(i) );
            delta_serial_temp = delta_serial(ttheta >= t_starts(i) & ttheta < t_stops(i) );
            
            
            goods = ( spk_hrs > t_starts(i) & spk_hrs < t_stops(i) ); % Identify spks in time range
            t_goodspks = spk_hrs(goods);
            is_in_theta = [];
            
            for ii = 1:length(t_goodspks)
                theta_index = find( (t_goodspks(ii) >= ttheta_temp - ttheta_width/2) & (t_goodspks(ii) < ttheta_temp + ttheta_width/2),1,'first');
                theta_index = wkeep(theta_index,1,'c');     % Keep only center entry
                if ~isempty(theta_index); is_in_theta(ii) = data_binary_temp(theta_index);
                else is_in_theta(ii) = -1;
                end
                
                if plot_individual_spks || recalculate_FFT_ratio_from_scratch
                    theta_band = [6 10];
                    delta_band = [0 6];
                    %x = Find_original_times (ratN,chanN,t_goodspks(ii)*3600);         % Takes raw data from 2-seconds interval with SPK centered at 1 second
                    x = Find_original_times (ratN,chanN,ttheta_temp(theta_index)*3600);     % Takes raw data from the original 2-second theta interval into which the SPK falls
                    Fs=12207.03125/12;
                    t = (((0:length(x)-1)) ) /Fs;
                    x=x(:);
                    xfilt = qif(t,x,[0 theta_band(1); theta_band(2) Inf]);
                    [f F] = generate_FFT_file (t,x,Fs);
                    index = find(f <= 100); f = f(index); F = F(index);
                    itheta = find((f >= theta_band(1)) .* (f < theta_band(2)));
                    idelta = find((f >= delta_band(1)) .* (f < delta_band(2)));
                    theta_F2 = mean(F(itheta));
                    delta_F2 = mean(F(idelta));
                    df = f(2)-f(1);
%                     theta_F2 = sum(F(itheta)) * df;
%                     delta_F2 = sum(F(idelta)) * df;
                end
                
                if recalculate_FFT_ratio_from_scratch; is_in_theta(ii) = (theta_F2 / delta_F2) > 2.0; end
                
                if plot_individual_spks && is_in_theta(ii)

                    figure;
                    subplot(211); plot(f,F);
                        hold on; plot(f(itheta), F(itheta),'r');
                        hold on; plot(f(idelta), F(idelta),'m');
%                         title(['Theta est vs saved serial ' num2str(theta_F2) ' ' num2str(theta_serial_temp(theta_index))]);
%                         xlabel(['Delta est vs saved serial ' num2str(delta_F2) ' ' num2str(delta_serial_temp(theta_index))]);
                        title(['Ratio est vs saved serial ' num2str(theta_F2/delta_F2) ' ' num2str(theta_serial_temp(theta_index)/delta_serial_temp(theta_index))]);
                    
                    subplot(212); hold on; plot(t, x,'g');
                    subplot(212); hold on; plot(t, xfilt,'r'); 
                    title(['Is in theta ' num2str(is_in_theta(ii))])
                end
            end
            
            t_goodspks = t_goodspks(is_in_theta>=0);
            is_in_theta = is_in_theta(is_in_theta>=0);
            
            if plot_on
                figure; plot(ttheta_temp,data_binary_temp,'b');
                hold on; plot(t_goodspks,is_in_theta,'ko');
                hold on; plot(t_goodspks(is_in_theta == 1),ones(1,sum(is_in_theta)),'ro');
            end
            
            fract_intheta = sum(data_binary) / length(data_binary);
            
            t_center(i) = mean([t_starts(i) t_stops(i)]);
            r(i) = sum(is_in_theta);
            r_normalized(i) = sum(is_in_theta) / fract_intheta; % normalize r based on the amount of theta actually present in the 1-hour bin
            
        end
    end
end


function [spike_times_arr filename_arr waveforms] = get_spike_times_unsorted (datapath, filenames,times_log)
    plot_on = 0;
    plot_waveforms = 0;
    
    file_list = dir([datapath filenames]);
    
    spike_times_arr = [];
    waveforms = [];
    
    if isempty(file_list)
        spike_times_arr = [];
        filename_arr = [];
        waveforms = [];
        return;
    end

     for i = 1:length(file_list)
%    for i = 33:35
        %fprintf (['Loading file ' [datapath '/' file_list(i).name] '\n']);
        load ([datapath '/' file_list(i).name]);
        
        curr_filenum = str2num(file_list(i).name(12:15));
        index = find(times_log(:,1) == curr_filenum,1,'first');
        i
        curr_spiketimes = spikes.spiketimes + times_log(index,2)*24*3600;


        spike_times_arr = [spike_times_arr (curr_spiketimes(:))'];
        filename_arr{i} = file_list(i).name;
       % waveforms = [waveforms; spikes.PeakAdjustwaveforms];

    end
    
    %clusters_SPKs(waveforms');
    
end

function [t_center,r] = get_num_spikes_per_hour_unsorted (spike_times_arr,binsize,overlap)

    %binsize in hours
    %overlap in fractional
    
    spk_hrs = spike_times_arr / 3600;
    
    starttime = min(spk_hrs);
    endtime = max(spk_hrs);
    
    t_starts = [starttime:(binsize*(1-overlap)):(endtime-binsize)];
    t_stops = t_starts + binsize;
    
    if isempty(t_starts)
        t_center = [];
        r = [];
        return
    end
    
    for i = 1:length(t_starts)
       goods = ( spk_hrs > t_starts(i) & spk_hrs < t_stops(i) );
       t_center(i) = mean([t_starts(i) t_stops(i)]);
       r(i) = sum(goods);
    end
    
    
end

function T = clusters_SPKs(units_arr)
    plot_on = 1;

    Fs = round(12207 / 10);
    [coef score,latent] = princomp(units_arr);
    unit_t = (1:size(units_arr,1)) / Fs;
    
    lx = coef(:,1); ly = coef(:,2);lz = coef(:,3);
    X = [lx ly lz];
    num_clusters = 4;
    T = clusterdata(X,'linkage','ward','savememory','on','maxclust',num_clusters);

%     figure; scatter(X(:,1),X(:,2),100,T,'filled')
    figure; set(gcf,'Color',[1 1 1])
    colourarr = get(gca,'ColorOrder');
    set(gcf,'Color',[1 1 1]); scatter3(X(:,1),X(:,2),X(:,3),100,T,'filled')
    colormap(colourarr(1:num_clusters,:))
    xlabel ('PCA1','FontSize',16); ylabel('PCA2','FontSize',16); zlabel('PCA3','FontSize',16); set(gca,'FontSize',16);
    
    if plot_on
        figure; set(gcf,'Color',[1 1 1])
        for i = 1:num_clusters
            subplot(num_clusters,1,i);
            index = find(T == i);
            plot(unit_t,units_arr(:,index),'Color',colourarr(i,:));
            hold on; plot(unit_t,mean(units_arr(:,index)'),'k','LineWidth',2);
            legend(['Type ' num2str(i)]);
            xlabel ('time(ms)','FontSize',16); ylabel('LFP (mV)','FontSize',16);set(gca,'FontSize',16);
        end
    end
end



function x = Find_original_times (ratN,chanN,time)

    Fs=12207.03125/12;
    if (~exist('ratN','var')); ratN=1; end
    if (~exist('chanN','var')); chanN=2; end
    if (~exist('time','var')); time=1; end
    
    ratnum = num2str(ratN, '%6.3d');
    channum = num2str(chanN, '%6.2d');

    
    
    plot_on = 0;
    
    
    path_ratlog = ('~/Anco/Evol/EEG/RatData_out');
    raw_path = ('/Volumes/MYPASSPORT/data/RawOneK');
    
    
    outlog_path = './Ratlog_out';
    [outlog_path outlog_name] = save_log(path_ratlog, ratnum, outlog_path);   % Save log file
    load ([outlog_path '/' outlog_name])    % Load log file
    times_dat = [fileNums(:) (tabs(:) - floor(min(tabs)))*24*3600];
    
    index = find(times_dat(:,2) <= time,1,'last');
    correct_file_num = times_dat(index ,1);
    file_start_time = times_dat(index , 2);
    dt = time - file_start_time;
    correct_file_num_str = num2str(correct_file_num, '%6.4d');
    
    
    fprintf(['Loading ...' raw_path '/Rat' ratnum 'OneK/Rat' ratnum 'ch' channum 'F' correct_file_num_str '_DownSampled_Dec.bin \n']);
    len = Fs * 2;
    pos = round(dt * Fs) - round(Fs*1);
    x = readbin([raw_path '/Rat' ratnum 'OneK/Rat' ratnum 'ch' channum 'F' correct_file_num_str '_DownSampled_Dec.bin'],'int16',pos,len);
    
    
    if plot_on
        
        %figure; plot((0:(length(x)-1))/Fs + time - 0,x,'g');
%         x2 = readbin([raw_path '/Rat' ratnum 'OneK/Rat' ratnum 'ch' channum 'F' correct_file_num_str '_DownSampled_Dec.bin'],'int16',0,pos + 20*Fs);
%         figure;
%         plot(1:length(x2),x2);
%         hold on; plot(pos + (1:length(x)) - 1, x,'g');
        
        figure;
%         plot((0:length(x2)-1) / Fs + file_start_time,x2);
        hold on; plot((((0:length(x)-1)) ) /Fs , x','g');
        
    end
    

end



function [f F] = generate_FFT_file (t,x, fs)

    if (~exist('fs','var')); fs=round(12207/12); end
    
    plot_on = 0;
    dt = 1/fs;

    %[f F] = daveFFT(t,x,1);
    [f F]  = dave_binoverlap_FFT(t, x, 10); % 10 second bin size.
    %[f F] = dave_welch2_FFT(t,x,10);    

    N = length(f);
    T = dt*length(t);
    df = 1/T;
    dw = 2*pi*df;

    f = f(1:round(N/2)); F = F(1:round(N/2));  % Take only positive half
    F = abs(F).^2 * T / (2*pi)^2 *2;   % Get into mV^2 / Hz (multiply x2 to compensate for -ve half)
    if plot_on;
        subplot(211); hold on; plot(t/3600, x,'r'); xlabel('time (h)');
        subplot(212); loglog(f, F,'g'); hold on; 
        pow_time = sum(abs(x).^2)/length(x)
        pow_freq = sum(F) * df
        differen = pow_time - pow_freq

    end
    
end



