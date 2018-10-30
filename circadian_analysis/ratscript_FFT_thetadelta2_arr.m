%

function data_struct =  ratscript_FFT_thetadelta2_arr (ratN, chanN,theta_band,delta_band,timerange0,extract_circadian_info,cosinor_mode,extract_band)

    global tabs_ctrl_1 tabs_ctrl_2 tabs_acute_1 tabs_acute_2 tabs_chr_1 tabs_chr_2 tabs_dark_1 tabs_dark_2

    if (~exist('ratN','var')); ratN=6; end
    if (~exist('chanN','var')); chanN=2; end
    if (~exist('theta_band','var')); theta_band=[5 10]; end
    if (~exist('delta_band','var')); delta_band=[1 5]; end
    if (~exist('timerange0','var')); timerange0=[]; end
    if (~exist('extract_circadian_info','var')); extract_circadian_info=1; end  % Setting to 0 will skip calculating phases and just return powers. Used for generating data for "GetFigs"
    
    if (~exist('cosinor_mode','var')); cosinor_mode=1; end                      % Setting to 1 or 2 will prevent from running "get_phase" command
                                                                                % Setting to 1 will return the theta band power (or "extract_band" if specified)
                                                                                % Setting to 2 will return the power in just theta epochs
                                                                                % Setting to 3 will return the power in non-theta epochs
                                                                                % Setting to 0 will return the theta epoch probability (data_binary)
%     if (~exist('extract_band','var')); extract_band=[1 5; 5 10; 12 30; 30 70; 70 100; 100 140; 140 300; 300 500]; end
    if (~exist('extract_band','var')); extract_band=[5 10]; end
    if isempty(timerange0); timerange0 = 0; end

    
    
    ratio_thresh0 = 1.75;

    FS_axis = 20;
    FS_axis_sm = 14;
    FS_axis_timeseries = 30;
    start_from_scratch=1;
    smart_filter = 1;                   % If use smart filter, don't use the manual filters below.
        spect_from_scratch = 0;         % If set to 1, filter redo filtering of spectra
    pre_smooth_data = 1;
    use_tcell_search = 0;
    fract_maxgap0 = 0.9;

    plot_on=0;
        plot_filterdata = 0;
        shift_1st_seizure_time_to_zero = 1;
        
    
   basepath = '/Users/davestanley/Nex/PhD/Evol_to_epil/Data/RawOneK';
%      basepath = '/data/SpikeStudyEpileptogenesis';
  path_ratlog='../Disk_inventory/RatData_out';
%      path_ratlog=[basepath '/Logs'];
    path_savephase = './Phase_Roll';

    clc
    format compact
    
%     if low_mid_high == 1; freq_band=[1 5];
%     elseif low_mid_high == 2; freq_band=[5 10];     % As defined in Buzsaki, Wang 2012 Ann Rev - Mechanisms of gamma oscillations
%     elseif low_mid_high == 3; freq_band=[12 25];
%     elseif low_mid_high == 4; freq_band=[25 50];
%     elseif low_mid_high == 5; freq_band=[50 90];
%     elseif low_mid_high == 6; freq_band=[90 140];
%     elseif low_mid_high == 7; freq_band=[140 300];
%     elseif low_mid_high == 8; freq_band=[300 500];
%     end


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
    if ratnum == '005';
        ii = 1;
        for i=[1:32]; channum_arr{ii}=num2str(i,'%05.2d'); ii=ii+1; end
        stim_time = 11;
        first_seizure = 113;
    end
    if ratnum == '008';
        ii = 1;
        for i=[1:32]; channum_arr{ii}=num2str(i,'%05.2d'); ii=ii+1; end
        stim_time = 18;
        first_seizure = 128;
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
    end



    outlog_path = ['./Ratoutmat_FFT'];
    if start_from_scratch
        [outlog_path outlog_name] = save_log(path_ratlog, ratnum, outlog_path);   % Save log file
    else
        outlog_name = strcat('Rat',ratnum,'log');
    end
    load ([outlog_path '/' outlog_name])    % Load log file

    if smart_filter
        load (['../14a_def_chronux/wrkspc_Rat_temp2_' ratnum 'Ch' channum 'file2D.mat'])
        file2D = clean_file2D (file2D, file2D.tabs < 0 ); % Remove timepoints that are < 0
        
        if spect_from_scratch
            [t_temp d_temp] = extract_FFT_band(file2D,[0 100]);
            bad_indices = smartfilter_files (t_temp,d_temp,file2D.fnum,0,ratN);
            bad_indices2 = smartfilter_files (t_temp,d_temp,file2D.fnum,1,ratN);
            bad_indices3 = bad_indices | bad_indices2;
%             if ratN == 1;           
%                 [t_temp2 d_temp2] = extract_FFT_band(file2D,[100 500]);     % Remove high frequency noisy points as well
%                 d_temp2(bad_indices3) = mean(d_temp2(~bad_indices3));
%                 bad_indices3b = smartfilter_files (t_temp2,d_temp2,file2D.fnum,0,ratN);
%                 bad_indices3 = [bad_indices3 | bad_indices3b]; clear bad_indices3b t_temp2 d_temp2; end
            
            
            bad_indices4 = smartfilter_datapoints(t_temp(~bad_indices3),d_temp(~bad_indices3),0);
            if ratN == 1; bad_indices4b = smartfilter_datapoints(t_temp(~bad_indices3),d_temp(~bad_indices3),1); bad_indices4 = [bad_indices4 | bad_indices4b]; clear bad_indices4b; end
            if ratN == 11; bad_indices4b = (t_temp(~bad_indices3)>=1.92 & t_temp(~bad_indices3)<=2.5); bad_indices4 = [bad_indices4 | bad_indices4b]; clear bad_indices4b; end
%             if ratN == 1;           
%                 [t_temp2 d_temp2] = extract_FFT_band(file2D,[100 500]);     % Remove high frequency noisy points as well
%                 bad_indices4c = smartfilter_datapoints_r1(t_temp2(~bad_indices3),d_temp2(~bad_indices3));
%                 bad_indices4 = [bad_indices4 | bad_indices4c]; clear bad_indices4c t_temp2 d_temp2; end

%             if ratN == 10
%                 %fnum_temp = file2D.fnum(~bad_indices3);
%                 tabs_temp = file2D.tabs(~bad_indices3);
%                 bad_indices4 = bad_indices4 | (tabs_temp >= 18 & tabs_temp <= 20) | (tabs_temp >= 23 & tabs_temp <= 24);
%                 clear tabs_temp
%             end
            clear t_temp d_temp;
            
            save(['wrkspc_filt_Rat' ratnum 'Ch' channum '.mat'],'bad_indices','bad_indices2','bad_indices4');  % Save just the indices to save space
        else
            load (['wrkspc_filt_Rat' ratnum 'Ch' channum '.mat']);
        end
         
        %%% Apply filters to invoke indices.
        file2D_cleaned = clean_file2D (file2D,(bad_indices | bad_indices2));
        file2D_cleaned = clean_file2D (file2D_cleaned,bad_indices4);
        
        if plot_filterdata
            d_temp = sum(file2D.F);
            t_temp = file2D.tabs;
            bad_indices3 = bad_indices | bad_indices2;
            t_temp = t_temp(~bad_indices3);
            d_temp = d_temp(~bad_indices3);
            figure;
            %subplot(211);
            plot(file2D.tabs,sum(file2D.F),'.')
            hold on; plot(file2D_cleaned.tabs,sum(file2D_cleaned.F),'r.')
            hold on; plot(t_temp(bad_indices4),d_temp(bad_indices4),'g.')
            legend('original','good data','removed')
            xlabel('time days');
            
%             t_temp = file2D.fnum;
%             t_temp = t_temp(~bad_indices3);
%             subplot(212); plot(file2D.fnum,sum(file2D.F),'.')
%             hold on; plot(file2D_cleaned.fnum,sum(file2D_cleaned.F),'r.')
%             bad_indices3 = bad_indices | bad_indices2;
%             hold on; plot(t_temp(bad_indices4),d_temp(bad_indices4),'g.')
%             legend('original','good data','removed')
%             xlabel('file number');
        end
        file2D = file2D_cleaned;
        clear file2D_cleaned
        
    else
        load (['./wrkspc_Rat_temp2_' ratnum 'Ch' channum 'file2D.mat'])
        
%         if ratN == 10; 
%             bad_files = [10 11 18 21 33:35 43 45:51 61 62 5:9 52 53 57 64 44];
%             bad_files = [ 44 45:51 52 61:63 64];
%             fnum =file2D.fnum;
%             bad_indices = logical(zeros(1,length(fnum)));
%             for i = 1:length(bad_files)
%                 bad_indices = bad_indices | (fnum == bad_files(i));
%             end
%             file_orig = file2D
%             file2D = clean_file2D (file2D,bad_indices);
% 
%             
%         end
        
    end
    
    

    [torig dfull] = extract_FFT_band (file2D,[0 100]);
    [torig theta_serial] = extract_FFT_band (file2D,theta_band);
    [torig delta_serial] = extract_FFT_band (file2D,delta_band);
    if cosinor_mode ~= 0;
        extract_serial = zeros(size(extract_band,1),length(torig));
        for i = 1:size(extract_band,1)
            [torig extract_serial(i,:)] = extract_FFT_band (file2D,extract_band(i,:));
        end
    else
        extract_serial = theta_serial;
    end
    
    ratio_serial = theta_serial ./ delta_serial;
    data_binary = ratio_serial > ratio_thresh0;
    
    


    for ii = 1:length(fileNames);
        fileNums(ii) = str2num(fileNames{ii}(1:4));
    end
    if ~isempty(find(diff(fileNums)==0)); fprintf('Warning: Overlapping filenums. Check log! \n'); end
    stim_index = find(fileNums >= stim_time, 1, 'first'); stim_tabs=(tabs(stim_index)-floor(tabs(1)));
    seiz_index = find(fileNums >= first_seizure, 1, 'first'); seiz_tabs=(tabs(seiz_index)-floor(tabs(1)));
    if ratnum == '006'; seiz_tabs = 16; end
    tabs_ctrl_1 = 0;
    tabs_ctrl_2 = stim_tabs;
    tabs_acute_1 = stim_tabs;
    tabs_acute_2 = seiz_tabs;
    tabs_chr_1 = seiz_tabs;
    tabs_chr_2 = 100;
    
    if ratnum == '001'
        dark_index = find(fileNums >= 188, 1, 'first'); dark_tabs =(tabs(dark_index)-floor(tabs(1))); 
        tabs_chr_2 = dark_tabs;
        tabs_dark_1 = dark_tabs;
        tabs_dark_2 = 100;
    end
    
    
    if ratnum == '009'; tabs_acute_1 = stim_tabs + 2; end
    if ratnum == '009'; tabs_ctrl_2 = 9.0; end
    if ratnum == '004'; tabs_acute_1 = 9; tabs_acute_2 = 30.5; tabs_chr_1 = 47.5; end
    if ratnum == '010'; tabs_chr_2 = 46.0; end
    if ratnum == '010'; tabs_acute_1 = 20; end
    %if ratnum == '010'; tabs_ctrl_1 = 3.5; end      % Remove discontinuity at starting
    %if ratnum == '001'; tabs_chr_2 = 37.0; end
    %if ratnum == '004'; tabs_acute_1 = 16.5; end
    if ratnum == '001'; tabs_acute_1 = 20.0; end
%     if ratnum == '001'; tabs_acute_2 = 26.5; end
    %if ratnum == '001'; tabs_acute_2 = 31.46; tabs_chr_1 = tabs_acute_2; end        % Enforce min acute time of 2 weeks.
    if ratnum == '008' tabs_acute_1 = 5; tabs_ctrl_2 = 5; tabs_acute_2 = 17; end        % I chopped the ending of the acute period off for R8, since it seems to trail off.
    if ratnum == '005'; tabs_ctrl_2 = 3.193; tabs_acute_1 = 6.0; end

    data_struct.power_timerange.theta = calculate_powers (torig,theta_serial,timerange0);
    data_struct.power_timerange.delta = calculate_powers (torig,delta_serial,timerange0);
    data_struct.theta_fract = calculate_powers (torig,data_binary,timerange0);

    torig_theta = torig(data_binary);
    torig_delta = torig(~data_binary);
    theta_in_theta = theta_serial(data_binary); % Theta power in only theta epochs
    delta_in_delta = delta_serial(~data_binary); % Delta power in only delta epochs
    
    data_struct.power_epochs_timerange.theta = calculate_powers (torig_theta,theta_in_theta,timerange0);
    data_struct.power_epochs_timerange.delta = calculate_powers (torig_delta,delta_in_delta,timerange0);
    
    
    
    theta_temp1 = extract_serial(repmat(data_binary,size(extract_serial,1),1));
    extract_band_theta = reshape(theta_temp1,size(extract_serial,1),sum(data_binary));
    
    theta_temp1 = extract_serial(repmat(~data_binary,size(extract_serial,1),1));
    extract_band_delta = reshape(theta_temp1,size(extract_serial,1),sum(~data_binary));
    
    data_struct.power_ergodic.all = calculate_powers_ergodic(torig,extract_serial,timerange0);
    data_struct.power_ergodic.intheta = calculate_powers_ergodic(torig_theta,extract_band_theta,timerange0);
    data_struct.power_ergodic.indelta = calculate_powers_ergodic(torig_delta,extract_band_delta,timerange0);
    
    data_struct.power_allbands.all = calculate_powers2(torig,extract_serial,timerange0);
    data_struct.power_allbands.intheta = calculate_powers2(torig_theta,extract_band_theta,timerange0);
    data_struct.power_allbands.indelta = calculate_powers2(torig_delta,extract_band_delta,timerange0);

    
    data_struct.param.ratio_thresh0 = ratio_thresh0;
    data_struct.param.timerange = timerange0;
    data_struct.param.theta_band = theta_band;
    data_struct.param.delta_band = delta_band;
    
    if extract_circadian_info
    
        bin_size = 1/24;
        fract_overlap = 0.0;
        fract_maxgap = fract_maxgap0;
%         [times_serial data_serial] = daveMVAVG_bin (torig, data_binary, bin_size, fract_overlap, fract_maxgap,use_tcell_search); % Take average time spent in theta state.
%         [times_serial data_serial] = daveMVAVG_bin (torig_theta, theta_in_theta, bin_size, fract_overlap, fract_maxgap,use_tcell_search); % Take average time spent in theta state.
%         [times_serial data_serial] = daveMVAVG_bin (torig_delta, delta_in_delta, bin_size, fract_overlap, fract_maxgap,use_tcell_search); % Take average time spent in theta state.
%         [times_serial theta_serial_filt] = daveMVAVG_bin (torig, theta_serial, bin_size, fract_overlap, fract_maxgap,use_tcell_search); % Take average time spent in theta state.
%         [times_serial delta_serial_filt] = daveMVAVG_bin (torig, delta_serial, bin_size, fract_overlap, fract_maxgap,use_tcell_search); % Take average time spent in theta state.
        
        
        if cosinor_mode == 1
            [times_serial data_serial] = daveMVAVG_MAT (torig, extract_serial', bin_size, fract_overlap, fract_maxgap); % Take average time spent in theta state.
            %[times_serial data_serial] = daveMVAVG_bin (torig, extract_serial(3,:), bin_size, fract_overlap, fract_maxgap,use_tcell_search); % Take average time spent in theta state.
        elseif cosinor_mode == 2
            theta_temp1 = extract_serial(repmat(data_binary,size(extract_serial,1),1));
            theta_temp2 = reshape(theta_temp1,size(extract_serial,1),sum(data_binary));
            [times_serial data_serial] = daveMVAVG_MAT (torig(data_binary), theta_temp2', bin_size, fract_overlap, fract_maxgap); % Take average time spent in theta state.
            clear theta_temp1 theta_temp2
        elseif cosinor_mode == 3
            theta_temp1 = extract_serial(repmat(~data_binary,size(extract_serial,1),1));
            theta_temp2 = reshape(theta_temp1,size(extract_serial,1),sum(~data_binary));
            [times_serial data_serial] = daveMVAVG_MAT (torig(~data_binary), theta_temp2', bin_size, fract_overlap, fract_maxgap); % Take average time spent in theta state.
            clear theta_temp1 theta_temp2
        else
            [times_serial data_serial] = daveMVAVG_MAT (torig, data_binary', bin_size, fract_overlap, fract_maxgap); % Take average time spent in theta state.
%             [times_serial data_serial] = daveMVAVG_bin (torig_theta, theta_in_theta, bin_size, fract_overlap, fract_maxgap,use_tcell_search); % Take average time spent in theta state.
%             [times_serial data_serial] = daveMVAVG_bin (torig_delta, delta_in_delta, bin_size, fract_overlap, fract_maxgap,use_tcell_search); % Take average time spent in theta state.
%             [times_serial data_serial] = daveMVAVG_bin (torig, delta_serial, bin_size, fract_overlap, fract_maxgap,use_tcell_search); % Take average time spent in theta state.
%             [times_serial data_serial] = daveMVAVG_bin (torig, theta_serial, bin_size, fract_overlap, fract_maxgap,use_tcell_search); % Take average time spent in theta state.

        end
        


%         data_serial = data_serial / mean(data_serial);
        ds = abs(downsample(data_serial, 1));
        ts = downsample(times_serial, 1);
        if plot_on
            figure; set(gcf,'Color','w')
            %plot(torig_theta,theta_in_theta/mean(theta_in_theta))
            hold on; plot(ts, ds,'.b');
%             xlabel ('Time (days)','FontSize',FS_axis);
%             ylabel('� EMD Amp','FontSize',FS_axis);
        end

        bin_size = 6/24;
        fract_overlap = 0.9;
        fract_maxgap = fract_maxgap0;
        if pre_smooth_data
            [t_sm dat_sm] = daveMVAVG_MAT (ts, ds, bin_size, fract_overlap, fract_maxgap);
        else
            t_sm = ts;
            dat_sm = ds;
        end

        bin_size = 1;
        fract_overlap = 0.9;
        fract_maxgap = fract_maxgap0;
        [t_base dat_base] = daveMVAVG_MAT (ts, ds, bin_size, fract_overlap, fract_maxgap);




        if plot_on; hold on; plot(t_sm, dat_sm, 'k.','MarkerSize',18); end

        dat_basei = interp1(t_base,dat_base, t_sm);
        dat_sub = dat_sm - dat_basei;
        if plot_on
            %hold on; plot(t_sm, dat_basei, 'k')
            hold on; plot(t_base, dat_base, 'Color',[0 0.5 0],'LineWidth',2);
            hold on; plot([tabs_ctrl_1 tabs_ctrl_1],[0 max(max(dat_sm))],'k','LineWidth',2);
            hold on; plot([tabs_ctrl_2 tabs_ctrl_2],[0 max(max(dat_sm))],'k','LineWidth',2);
            hold on; plot([tabs_acute_1 tabs_acute_1],[0 max(max(dat_sm))],'m','LineWidth',2);
            hold on; plot([tabs_acute_2 tabs_acute_2],[0 max(max(dat_sm))],'m','LineWidth',2);
            hold on; plot([tabs_chr_1 tabs_chr_1],[0 max(max(dat_sm))],'r','LineWidth',2);
            hold on; plot([tabs_chr_2 tabs_chr_2],[0 max(max(dat_sm))],'r','LineWidth',2);
            if ratnum == '001';
                hold on; plot([tabs_dark_1 tabs_dark_1],[0 max(max(dat_sm))],'b','LineWidth',2);
                hold on; plot([tabs_dark_2 tabs_dark_2],[0 max(max(dat_sm))],'b','LineWidth',2);
            end
            for kk = 1:max(t_sm);
                hold on; plot([kk kk], [min(min(dat_sm)) max(max(dat_sm))*2],'k-.');
            end
            legend('original','smoothed','baseline');
            set(gca,'FontSize',FS_axis_timeseries);
        end

        if shift_1st_seizure_time_to_zero && plot_on
            day_of_SE = floor(stim_tabs);

            figure; set(gcf,'Color','w')
            %hold on; plot(ts-day_of_SE, ds);
%             xlabel ('Time (days)','FontSize',FS_axis);
%             ylabel('� EMD Amp','FontSize',FS_axis);
            set(gca,'FontSize',FS_axis);
            hold on; plot(t_sm-day_of_SE, dat_sm, 'k.','MarkerSize',18); 
            %hold on; plot(t_sm-day_of_SE, dat_basei, 'k')
            hold on; plot(t_base-day_of_SE, dat_base, 'Color',[0 0.5 0],'LineWidth',2);
            hold on; plot([tabs_ctrl_1 tabs_ctrl_1]-day_of_SE,[0 max(max(dat_sm))],'k','LineWidth',2);
            hold on; plot([tabs_ctrl_2 tabs_ctrl_2]-day_of_SE,[0 max(max(dat_sm))],'k','LineWidth',2);
            hold on; plot([tabs_acute_1 tabs_acute_1]-day_of_SE,[0 max(max(dat_sm))],'m','LineWidth',2);
            hold on; plot([tabs_acute_2 tabs_acute_2]-day_of_SE,[0 max(max(dat_sm))],'m','LineWidth',2);
            hold on; plot([tabs_chr_1 tabs_chr_1]-day_of_SE,[0 max(max(dat_sm))],'r','LineWidth',2);
            hold on; plot([tabs_chr_2 tabs_chr_2]-day_of_SE,[0 max(max(dat_sm))],'r','LineWidth',2);
            if ratnum == '001';
                hold on; plot([tabs_dark_1 tabs_dark_1]-day_of_SE,[0 max(max(dat_sm))],'b','LineWidth',2);
                hold on; plot([tabs_dark_2 tabs_dark_2]-day_of_SE,[0 max(max(dat_sm))],'b','LineWidth',2);
            end
            for kk = 1:max(t_sm);
                hold on; plot([kk kk]-day_of_SE, [min(min(dat_sm)) max(max(dat_sm))*2],'k-.');
            end
            %legend('original','smoothed','baseline');
            set(gca,'FontSize',FS_axis_timeseries);
        end

        if plot_on
            figure; plot(t_sm, dat_sub);
            hold on; plot([tabs_ctrl_1 tabs_ctrl_1],[min(min(dat_sub)) max(max(dat_sub))],'k');
            hold on; plot([tabs_ctrl_2 tabs_ctrl_2],[min(min(dat_sub)) max(max(dat_sub))],'k');
            hold on; plot([tabs_acute_1 tabs_acute_1],[min(min(dat_sub)) max(max(dat_sub))],'m');
            hold on; plot([tabs_acute_2 tabs_acute_2],[min(min(dat_sub)) max(max(dat_sub))],'m');
            hold on; plot([tabs_chr_1 tabs_chr_1],[min(min(dat_sub)) max(max(dat_sub))],'r');
            hold on; plot([tabs_chr_2 tabs_chr_2],[min(min(dat_sub)) max(max(dat_sub))],'r');
            if ratnum == '001';
                hold on; plot([tabs_dark_1 tabs_dark_1],[min(min(dat_sub)) max(max(dat_sub))],'b','LineWidth',2);
                hold on; plot([tabs_dark_2 tabs_dark_2],[min(min(dat_sub)) max(max(dat_sub))],'b','LineWidth',2);
            end
            for kk = 1:max(t_sm);
                hold on; plot([kk kk], [min(min(dat_sub)) max(max(dat_sub))*2],'k-.');
            end
            xlabel('time, (days)','FontSize',FS_axis); ylabel('� EMD Amp','FontSize',FS_axis);
            set(gca,'FontSize',FS_axis_timeseries);
        end

        %     [A phi t_sin] = fit_sinusoids (t_sm, dat_sub, 1.0);
        %     figure;
        %     subplot(211); plot(t_sin, A);
        %     subplot(212); plot(t_sin, phi);



        % % % % %     Average everything
        %     t_smm= mod(t_sm, 1);
        %     [A phi t_sin] = fit_sinusoids (t_smm, dat_sub, 1.0);
        %cosinor_mode = 0;
        use_rollback = 1;
        % % % % %   Control
        %index = find (t_sm <= stim_tabs);
        index = find ( (tabs_ctrl_1 <= t_sm) .* (t_sm <= tabs_ctrl_2) );
        i=1;
        if ~isempty(index)
            t_section = t_sm(index); dat_section = dat_sub(index,:);
        %         figure; plot(t_section, dat_section,'o');
            
            data_struct.ctrl.t = t_section;
            data_struct.ctrl.d = dat_section;
        else
            data_struct.chr.t = [];
            data_struct.chr.d = [];
        end

        % % % % %   Acute (Latent)
        %index = find ( (stim_resume_tabs <= t_sm) .* (t_sm <= seiz_tabs) );
        index = find ( (tabs_acute_1 <= t_sm) .* (t_sm <= tabs_acute_2) );
    %     if ratnum == '004'; index = find ( ((tabs_acute_1 <= t_sm) .* (t_sm < 23.5)) | ((t_sm <= tabs_acute_2) .* (t_sm > 27)) ); end
    %     if ratnum == '010'; index = find ( (tabs_acute_1 <= t_sm) .* (t_sm <= tabs_acute_2) .* ~((t_sm>18) .* (t_sm<20) ) .* ~((t_sm>23) .* (t_sm<24)) ); end  % Discount isolated data points
        %if ratnum == '001'; index = find ( (21 <= t_sm) .* (t_sm <= tabs_acute_2) ); end
        i=2;
        if ~isempty(index)
            t_section = t_sm(index); dat_section = dat_sub(index,:);
        %         figure; plot(t_section, dat_section,'o');
          
            data_struct.acute.t = t_section;
            data_struct.acute.d = dat_section;
        else
            data_struct.chr.t = [];
            data_struct.chr.d = [];
        end


        % % % % %   Chronic
        %index = find (seiz_tabs <= t_sm);        
        index = find ( (tabs_chr_1 <= t_sm) .* (t_sm <= tabs_chr_2) );
        i=3;
        if ~isempty(index)
            t_section = t_sm(index); dat_section = dat_sub(index,:);
        %         figure; plot(t_section, dat_section,'o');
         
            data_struct.chr.t = t_section;
            data_struct.chr.d = dat_section;
        else
            data_struct.chr.t = [];
            data_struct.chr.d = [];
        end
        
        
        % % % % %   Darkness
        if ratnum == '001'
            %index = find (seiz_tabs <= t_sm);        
            index = find ( (tabs_dark_1 <= t_sm) .* (t_sm <= tabs_dark_2) );
            i=4;
            if ~isempty(index)
                t_section = t_sm(index); dat_section = dat_sub(index,:);
            %         figure; plot(t_section, dat_section,'o');

                data_struct.dark.t = t_section;
                data_struct.dark.d = dat_section;
            else
                data_struct.dark.t = [];
                data_struct.dark.d = [];
            end
        end

        data_struct.all.t = times_serial;
        data_struct.all.d = data_serial;
        data_struct.all.stim_tabs = stim_tabs;
        data_struct.all.seiz_tabs = seiz_tabs;
        data_struct.all.tabs_ctrl_1 = tabs_ctrl_1;
        data_struct.all.tabs_ctrl_2 = tabs_ctrl_2;
        data_struct.all.tabs_acute_1 = tabs_acute_1;
        data_struct.all.tabs_acute_2 = tabs_acute_2;
        data_struct.all.tabs_chr_1 = tabs_chr_1;
        data_struct.all.tabs_chr_2 = tabs_chr_2;
        data_struct.all.tabs_dark_1 = tabs_dark_1;
        data_struct.all.tabs_dark_2 = tabs_dark_2;

        clear t_sm dat_sub stim_tabs seiz_tabs;
    
    end

end


function out = calculate_powers (t,x,timerange)

    global tabs_ctrl_1 tabs_ctrl_2 tabs_acute_1 tabs_acute_2 tabs_chr_1 tabs_chr_2

    timerange;
    if timerange == 1; min_timerange = 0.9375; max_timerange = 0.0625; cross_midnight=1;              %11:30pm to 1:30am
    elseif timerange == 2; min_timerange = 0.0625; max_timerange = 0.1875; cross_midnight=0;          %1:30am to 4:30am
    elseif timerange == 3; min_timerange = 0.1875; max_timerange = 0.3125; cross_midnight=0;          %...
    elseif timerange == 4; min_timerange = 0.3125; max_timerange = 0.4375; cross_midnight=0;          %...
    elseif timerange == 5; min_timerange = 0.4375; max_timerange = 0.5625; cross_midnight=0;          %...
    elseif timerange == 6; min_timerange = 0.5625; max_timerange = 0.6875; cross_midnight=0;          %...
    elseif timerange == 7; min_timerange = 0.6875; max_timerange = 0.8125; cross_midnight=0;          %...
    elseif timerange == 8; min_timerange = 0.8125; max_timerange = 0.9375; cross_midnight=0;          %7:30pm to 10:30pm
    elseif timerange == 0; min_timerange = -0.1; max_timerange = 1.1; cross_midnight=0;               %All times
    end
    
    if cross_midnight
        ind = (mod(t,1) < max_timerange) | (mod(t,1) >= min_timerange);
    else
        ind = (mod(t,1) < max_timerange) & (mod(t,1) >= min_timerange);
    end


    index_pre =  (tabs_ctrl_1 <= t) & (t <= tabs_ctrl_2) .* (ind) ;
    index_acute =  (tabs_acute_1 <= t) & (t <= tabs_acute_2) .* (ind)  ;
    index_chronic =  (tabs_chr_1 <= t) & (t <= tabs_chr_2) .* (ind)  ;
    out(1) = mean(x(index_pre));
    out(2) = mean(x(index_acute));
    out(3) = mean(x(index_chronic));
end

function out = calculate_powers2 (t,x,timerange)

    global tabs_ctrl_1 tabs_ctrl_2 tabs_acute_1 tabs_acute_2 tabs_chr_1 tabs_chr_2

    timerange;
    if timerange == 1; min_timerange = 0.9375; max_timerange = 0.0625; cross_midnight=1;              %11:30pm to 1:30am
    elseif timerange == 2; min_timerange = 0.0625; max_timerange = 0.1875; cross_midnight=0;          %1:30am to 4:30am
    elseif timerange == 3; min_timerange = 0.1875; max_timerange = 0.3125; cross_midnight=0;          %...
    elseif timerange == 4; min_timerange = 0.3125; max_timerange = 0.4375; cross_midnight=0;          %...
    elseif timerange == 5; min_timerange = 0.4375; max_timerange = 0.5625; cross_midnight=0;          %...
    elseif timerange == 6; min_timerange = 0.5625; max_timerange = 0.6875; cross_midnight=0;          %...
    elseif timerange == 7; min_timerange = 0.6875; max_timerange = 0.8125; cross_midnight=0;          %...
    elseif timerange == 8; min_timerange = 0.8125; max_timerange = 0.9375; cross_midnight=0;          %7:30pm to 10:30pm
    elseif timerange == 0; min_timerange = -0.1; max_timerange = 1.1; cross_midnight=0;               %All times
    end
    
    if cross_midnight
        ind = (mod(t,1) < max_timerange) | (mod(t,1) >= min_timerange);
    else
        ind = (mod(t,1) < max_timerange) & (mod(t,1) >= min_timerange);
    end

    index_pre =  (tabs_ctrl_1 <= t) & (t <= tabs_ctrl_2) .* (ind) ;
    index_acute =  (tabs_acute_1 <= t) & (t <= tabs_acute_2) .* (ind)  ;
    index_chronic =  (tabs_chr_1 <= t) & (t <= tabs_chr_2) .* (ind)  ;
    
    theta_temp1 = x(repmat(index_pre,size(x,1),1));
    x_sect = reshape(theta_temp1,size(x,1),sum(index_pre)); % Pull out section of x matrix    
    out{1} = mean(x_sect');
    
    theta_temp1 = x(repmat(index_acute,size(x,1),1));
    x_sect = reshape(theta_temp1,size(x,1),sum(index_acute)); % Pull out section of x matrix    
    out{2} = mean(x_sect');
    
    theta_temp1 = x(repmat(index_chronic,size(x,1),1));
    x_sect = reshape(theta_temp1,size(x,1),sum(index_chronic)); % Pull out section of x matrix    
    out{3} = mean(x_sect');
end

function out = calculate_powers_ergodic (t,x,timerange)

    global tabs_ctrl_1 tabs_ctrl_2 tabs_acute_1 tabs_acute_2 tabs_chr_1 tabs_chr_2

    timerange;
    if timerange == 1; min_timerange = 0.9375; max_timerange = 0.0625; cross_midnight=1;              %11:30pm to 1:30am
    elseif timerange == 2; min_timerange = 0.0625; max_timerange = 0.1875; cross_midnight=0;          %1:30am to 4:30am
    elseif timerange == 3; min_timerange = 0.1875; max_timerange = 0.3125; cross_midnight=0;          %...
    elseif timerange == 4; min_timerange = 0.3125; max_timerange = 0.4375; cross_midnight=0;          %...
    elseif timerange == 5; min_timerange = 0.4375; max_timerange = 0.5625; cross_midnight=0;          %...
    elseif timerange == 6; min_timerange = 0.5625; max_timerange = 0.6875; cross_midnight=0;          %...
    elseif timerange == 7; min_timerange = 0.6875; max_timerange = 0.8125; cross_midnight=0;          %...
    elseif timerange == 8; min_timerange = 0.8125; max_timerange = 0.9375; cross_midnight=0;          %7:30pm to 10:30pm
    elseif timerange == 0; min_timerange = -0.1; max_timerange = 1.1; cross_midnight=0;               %All times
    end
    
    if cross_midnight
        ind = (mod(t,1) < max_timerange) | (mod(t,1) >= min_timerange);
    else
        ind = (mod(t,1) < max_timerange) & (mod(t,1) >= min_timerange);
    end

    out{1} = [];
    out{2} = [];
    out{3} = [];
    
    tmin = floor(min(t));
    tmax = floor(max(t));
    
    if cross_midnight
        tmin = tmin - 0.5; tmax = tmax + 0.5;
    end
    
    for tbin = tmin:(tmax-1)

        index = (ind) & (tbin < t) & (t <= (tbin+1));   % Extract section of time series for this 'bin'
        tcurr = t(index);

        if sum(index) > 120         % Must have > 2 minutes data in each time bin 

            theta_temp1 = x(repmat(index,size(x,1),1));
            x_sect = reshape(theta_temp1,size(x,1),sum(index)); % Pull out section of x matrix

            if sum( (tcurr >= tabs_ctrl_1) & ( tcurr <= tabs_ctrl_2 )) > 0
                out{1} = [out{1}; mean(x_sect')];
            elseif sum( (tcurr >= tabs_acute_1) & ( tcurr <= tabs_acute_2 )) > 0
                out{2} = [out{2}; mean(x_sect')];
            elseif sum( (tcurr >= tabs_chr_1) & ( tcurr <= tabs_chr_2 )) > 0
                out{3} = [out{3}; mean(x_sect')];
            else
                %fprintf('Somethinng wrong! \n');
            end
        end
    end
    
    
end


function [tabs Fout] = extract_FFT_band (file,freq_band)
    
%     f = struct2matrix(file,'f',[],0,0);
%     F = struct2matrix(file,'F',[],0,0);
%     tabs = struct2matrix(file,'tabs',[],0,0);
    
    f = file.f;
    F = file.F;
    tabs = file.tabs;

    freq_band = freq_band(:);
    
    index = find((f(:,1) >= freq_band(1)) & (f(:,1) < freq_band(2)));
    f_width = abs(diff(freq_band));
    
    Fout = sum(F( index,: ),1);
    Fout = Fout / f_width;
    
    Fout = Fout / 3276700^2 * (1e3)^2; % Add TDT scaling factor

end





function [fout Fout] = psd_freqbins (f, F, f_bins, f_centers, f_width)
    

    filter_out_60Hz = 1;
    plot_on = 0;
    colourarr = 'bgrmcykbgrmcykbgrmcykbgrmcykbgrmcykbgrmcykbgrmcyk';
    
    df = f(2) - f(1);

    if isempty(f_bins)
        f_bins = [];
        f_centers = (f_centers(:))';
        f_bins = [f_centers-f_width/2; f_centers+f_width/2];
        f_width_arr = diff(f_bins);
    else
        f_centers = mean(f_bins);
        f_width_arr = diff(f_bins);
    end
    
    if plot_on;
        figure;
        hold on; plot(f,F);
    end
    Fout = zeros(1,length(f_centers));
    for i=1:length(f_centers)
        index = find((f >= f_bins(1,i)) .* (f < f_bins(2,i)));
        if filter_out_60Hz
            index = find((f >= f_bins(1,i)) .* (f < f_bins(2,i)) .* ~((f > 58) .* (f < 62)) );
        end
        
        Fout(i) = sum(F(index));
        Fout(i) = Fout(i) / f_width_arr(i);
        
        if plot_on
            hold on; plot(f(index),F(index),[colourarr(mod(i,length(colourarr))) '.']);
            hold on; plot(f_centers(i),Fout(i),[colourarr(mod(i+1,length(colourarr))) 'o']);
        end
    end
    
    fout = f_centers;

end



function file2D_cleaned = clean_file2D (file2D,bad_indices)
    file2D_cleaned = file2D;
    file2D_cleaned.f = file2D_cleaned.f(:,~bad_indices);
    file2D_cleaned.F = file2D_cleaned.F(:,~bad_indices);
    file2D_cleaned.tcent = file2D_cleaned.tcent(:,~bad_indices);
    file2D_cleaned.fnum = file2D_cleaned.fnum(:,~bad_indices);
    file2D_cleaned.fstart_tabs = file2D_cleaned.fstart_tabs(:,~bad_indices);
    file2D_cleaned.tabs = file2D_cleaned.tabs(:,~bad_indices);

end


function bad_indices = smartfilter_files (t,d0,fnum,invert,ratN)
    
    if ~invert
        plot_on = 1;
        bin_size = 5;
        prefilter_threshold = 10;
        envelope_threshold = 20;
        max_allowed_badfiles = 10;
        if ratN == 7
            prefilter_threshold = 2.5;
        end
    else
        plot_on = 1;
        bin_size = 5;
        prefilter_threshold = 0.2;
        envelope_threshold = 15;
        max_allowed_badfiles = 10;
        invertmax = 1/5;
        invertmax = 1/5 * 3276700^2 / (1e3)^2;
        if ratN == 1
            % Use this for more strict filtering. - instead of doing this to remove entire files,
            % below I just filter out the low-amp data points point-by-point
%             invertmax = 1/100;
%             envelope_threshold = 10;
            max_allowed_badfiles = 500;
        end
        if ratN == 7
            envelope_threshold = 5;
        end
    end

    if invert
        d = 1./d0;
        d(d>invertmax) = invertmax;
    else
        d = d0;
    end
    mrksize = 20;
    ind = find(t > bin_size,1,'first');
    ind2 = find(t > t(end) - bin_size,1,'first');
    d_temp = [fliplr(d(1:ind)) d fliplr(d(ind2:end))];
    t_temp = [-1*fliplr(t(1:ind)) t t(end) + -1*(fliplr(t(ind2:end))-t(end)) ];
    ind = d_temp < (mean(d_temp) + std(d_temp)*prefilter_threshold);
    [t_sm d_sm d_std] = daveMVAVG_bin (t_temp(ind), d_temp(ind), bin_size, 0.5, 0.9,0); % Take average time spent in theta state.

    if plot_on
        figure; plot(t,d,'b.','MarkerSize',mrksize)
        hold on; plot(t_temp(ind), d_temp(ind),'k.','MarkerSize',mrksize)
        if ~invert hold on; errorbar(t_sm,d_sm,d_std*envelope_threshold,'r.','MarkerSize',mrksize)
        else hold on; errorbar(t_sm,d_sm,d_sm*envelope_threshold,'r.','MarkerSize',mrksize); end
        legend('original','prefiltered','envelope');
        xlabel('time days');
        
        %hold on; plot(t_sm,d_sm+ mean(d_std)*5,'y.')
    end
    clear ind ind2
    
    temp = interp1(t_sm,d_sm,t);
    %if ~invert envelope = temp + temp*envelope_threshold;
    if ~invert envelope = temp + interp1(t_sm,d_std,t)*envelope_threshold;
    else envelope = temp + temp*envelope_threshold; end
    ind = d > envelope;
    bad_files_candidates = fnum(ind);
    
    filebins = unique(fnum);
    [nfiles] = hist(fnum,filebins);
    [nbad] = hist(bad_files_candidates,filebins);
    
    if plot_on
%         figure;
%         subplot(311); bar(filebins,nfiles)
%         subplot(312); bar(filebins,nbad)
%         subplot(313); bar(filebins,nbad./nfiles)
    end
    
    %ind = (nbad./nfiles) > 0.001;
    ind = (nbad) >= max_allowed_badfiles;
    bad_files = (filebins(ind));
    if ratN == 1; bad_files = [bad_files 188]; end
    if ratN == 1; bad_files = [bad_files 187 151 152 153 201 202 218 250 251 252 253 264 267 115 129 141 149 151:154]; end
    if ratN == 10; bad_files = [bad_files 5:11]; end        % Added to remove discontinuity during start.
    if ratN == 7; bad_files = [bad_files 108:111]; end        % High frequency data explodes near the end. Remove this


    %%% Remove bad files
    bad_indices = logical(zeros(1,length(fnum)));
    for i = 1:length(bad_files)
        bad_indices = bad_indices | (fnum == bad_files(i));
    end


    if plot_on
       figure; plot(t,d,'b.');
       hold on; plot(t(~bad_indices),d(~bad_indices),'k.')
       figure; plot(fnum,d,'b.');
       hold on; plot(fnum(~bad_indices),d(~bad_indices),'k.')
       
%        figure; plot(t,log(d),'b.');
%        hold on; plot(t(~bad_indices),log(d(~bad_indices)),'k.')
    end


end

function bad_indices = smartfilter_datapoints(t,d0,invert)
    plot_on = 1;
    if ~invert
        bin_size = 5;
        envelope_threshold = 20;
    else
        plot_on = 0;
        bin_size = 5;
        envelope_threshold = 4;
        invertmax = 1/5;
    end
    
    if invert
        d = 1./d0;
        d(d>invertmax) = invertmax;
    else
        d = d0;
    end
    
    
    ind = find(t > bin_size,1,'first');
    ind2 = find(t > t(end) - bin_size,1,'first');
    d_temp = [fliplr(d(1:ind)) d fliplr(d(ind2:end))];
    t_temp = [-1*fliplr(t(1:ind)) t t(end) + -1*(fliplr(t(ind2:end))-t(end)) ];
    [t_sm d_sm d_std] = daveMVAVG_bin (t_temp, d_temp, bin_size, 0.5, 0.9,0); % Take average time spent in theta state.

    clear ind ind2
    
    temp = interp1(t_sm,d_sm,t);
    if ~invert envelope = temp + interp1(t_sm,d_std,t)*envelope_threshold;
    else envelope = temp + temp*envelope_threshold; end
    %envelope = temp + interp1(t_sm,d_std,t)*envelope_threshold;
    
    bad_indices = d > envelope;
    
    
    
    
    if plot_on
        figure; plot(t,d,'b.')
        
        if ~invert hold on; errorbar(t_sm,d_sm,d_std*envelope_threshold,'r.')
        else hold on; hold on; errorbar(t_sm,d_sm,d_sm*envelope_threshold,'r.'); end
        %hold on; errorbar(t_sm,d_sm,d_std*envelope_threshold,'r.')
        
        hold on; plot(t(bad_indices),d(bad_indices),'r.')
        legend('original data','envelope','good data')
        xlabel('time days');
    end

end
% 
% 
% function bad_indices = smartfilter_datapoints_r1(t,d0)
%     plot_on = 0;
%     bin_size = 5;
%     envelope_threshold = 10;
% 
%     d = d0;
% 
%     ind = find(t > bin_size,1,'first');
%     ind2 = find(t > t(end) - bin_size,1,'first');
%     d_temp = [fliplr(d(1:ind)) d fliplr(d(ind2:end))];
%     t_temp = [-1*fliplr(t(1:ind)) t t(end) + -1*(fliplr(t(ind2:end))-t(end)) ];
%     [t_sm d_sm d_std] = daveMVAVG_bin (t_temp, d_temp, bin_size, 0.5, 0.9,0); % Take average time spent in theta state.
% 
%     clear ind ind2
% 
%     temp = interp1(t_sm,d_sm,t);
%     envelope = temp + temp*envelope_threshold;
%     bad_indices = d > envelope;
% 
% 
% 
% 
%     if plot_on
%         figure; plot(t,d,'b.')
% 
% 
%         hold on; hold on; errorbar(t_sm,d_sm,d_sm*envelope_threshold,'r.');
%         %hold on; errorbar(t_sm,d_sm,d_std*envelope_threshold,'r.')
% 
%         hold on; plot(t(bad_indices),d(bad_indices),'r.')
%         legend('original data','envelope','good data')
%         xlabel('time days');
%     end
% 
% end
