

function generate_FFT_spect (ratN, chanN, fname_output)

    if (~exist('ratN','var')); ratN=9; end
    if (~exist('chanN','var')); chanN=2; end
    if (~exist('fname_output','var')); fname_output='blah.mat'; end

    addpath(genpath('/data/David/chronux'));

    ratnum = num2str(ratN, '%6.3d');
    channum = num2str(chanN, '%6.2d');

    
%     basepath = '/Users/davestanley/Nex/PhD/Evol_to_epil/Data/RawOneK';
    basepath = '/data/SpikeStudyEpileptogenesis';
%     path_ratlog='../../RatData_out';
    path_ratlog=[basepath '/Logs'];
    ratpath = ['/Rat' ratnum 'OneK'];
    
	filenames = dir ([basepath ratpath '/*ch' channum '*'])
    
     for k = 1:length(filenames)
%    for k = 1:2
        
        curr_file = filenames(k).name;
        if isempty(strfind(curr_file,'.bin')); continue; end
        fprintf (['Loading file ' [ratpath filenames(k).name] '\n']);


        [fout Fout tcent] = generate_smoothed_FFT ([basepath ratpath '/' curr_file]);
        file{k}.path = [basepath ratpath];
        file{k}.name = [curr_file];
        file{k}.fout = fout;
        file{k}.Fout = Fout;
        file{k}.tcent = tcent;
    end

    save(['wrkspc_Rat_temp2b_' ratnum 'Ch' channum 'filestruct.mat'],'file');
    
    file2D = rearrange_file (file, ratnum,path_ratlog);

    save(['wrkspc_Rat_temp2_' ratnum 'Ch' channum 'file2D.mat'],'file2D');
end


function file_out = rearrange_file (file,ratnum,path_ratlog)

    start_from_scratch = 1;

    outlog_path = ['./Ratoutmat_FFT'];
    if start_from_scratch
        [outlog_path outlog_name] = save_log(path_ratlog, ratnum, outlog_path);   % Save log file
    else
        outlog_name = strcat('Rat',ratnum,'log');
    end
    
    load ([outlog_path '/' outlog_name])    % Load log file
    fileNums = filenames2numbers(fileNames);
    
    tabs_zeroed = tabs - min(floor(tabs));
    k=0;
    for i = 1:length(file)
        if ~isempty (file{i})
            for j = 1:length(file{i}.fout)
                k=k+1;
                file_out.f(:,k) = file{i}.fout{j};
                file_out.F(:,k) = file{i}.Fout{j};
                file_out.tcent(k) = file{i}.tcent{j};
                file_out.fnum(k) = str2num(file{i}.name(12:15));
                index = find( fileNums == file_out.fnum(k),1,'first');
                if ~isempty(index)
	                file_out.fstart_tabs(k) = tabs_zeroed(index);
    	            file_out.tabs(k) = tabs_zeroed(index) + (file_out.tcent(k))/3600/24;
    	        else
    	        	file_out.fstart_tabs(k) = -1;
    	            file_out.tabs(k) = -1;
   	            end
            end
        end
    end
end

function [fout Fout tcent] = generate_smoothed_FFT (fname)

    plot_on=0;
    len=round(12207/12)*65;
    len=Inf;
    
    [f F t x tcent] = psd_timebins (fname, 2,len, round(12207)/12); % Get spectrogram for each hour of the dataset

    for i = 1:length(f)
        [fout{i} Fout{i}] = psd_freqbins (f{i},F{i}, [1 5; 5 10; 12 25; 25 50; 50 90; 90 140; 140 200; 200 500]'); % Bin spectrogram
    end

    if plot_on
        for i = 1:length(f)
            figure; 
            subplot(211); hold on; plot(t{i},x{i});
            subplot(212); hold on; 
            % plot(f{i}, F{i})
            plot(fout{i}, Fout{i},'g')
        end
    end
end


function [f_bins F_bins t_bins x_bins tcent_bins] = psd_timebins (fname, timebins, len, fs)

    if (~exist('fname','var')); fname='blah.bin'; end
    if (~exist('timebins','var')); timebins=3600; end % Duration of FFT bins to extract (def 1 hour)
    if (~exist('len','var')); len=Inf; end % Total length of data to read from specified file
    if (~exist('fs','var')); fs=round(12207/12); end
    
    plot_on = 0;
%     fs = 12207;
%     if is_downsampled; 
%         fs = fs/12;   % We are using downsampled data
%         len = round(len/12);
%     end

    dt = 1/fs;

    [x]= readbin(fname, 'int16', 0, len);
    x = x - mean(x);
    t = [0:length(x)-1];
    t = t * dt;
%     x=1*cos(2*pi*40*t);

    length_bin = round(timebins / dt);
    if (length_bin > length(x))
        fprintf('Timebin size is longer than dataset. Decreasing bin size\n');
        length_bin = length(x);
    end

    
    curr_index = 1;
    nbins = 0;
    
    while (curr_index + length_bin - 1 <= length(x))
        x_bin = x(curr_index:curr_index+length_bin-1);
        x_bin = x_bin - mean(x_bin);
%         t_bin = t(curr_index:curr_index+length_bin-1);
        t_bin = (0:length(x_bin)-1)*dt;

%         [f_bin F_bin] = daveFFT(t_bin,x_bin,1);
%         [f_bin F_bin]  = dave_binoverlap_FFT(t_bin, x_bin, 10); % 10 second bin size.
%         [f_bin F_bin] = dave_welch2_FFT(t_bin,x_bin,10);    
        
%         N = length(f_bin);
%          %T = dt*length(t_bin); %This is wrong because the actual "bin size" is 10 or length of t_bin
%          %df = 1/T;
%          df = f_bin(2) - f_bin(1);
%          T = 1/df;
%          dw = 2*pi*df;
%
%          f_bin = f_bin(1:round(N/2)); F_bin = F_bin(1:round(N/2));  % Take only positive half
%          F_bin = abs(F_bin).^2 * T / (2*pi)^2 *2;   % Get into mV^2 / Hz (multiply x2 to compensate for -ve half)
        
%         [Pxx,f] = pmtm(x_bin,[2],length(x_bin),fs);
%         f_bin = f;
%         F_bin = Pxx;
        
%         mtparams.tapers = [2 3];
%         mtparams.Fs = fs;
        mtparams.tapers = [1 2 1];
        mtparams.Fs = fs;
        [F_bin,f_bin]=mtspectrumc(x_bin,mtparams);

        if plot_on;
            
            subplot(211); hold on; plot(t_bin/3600, x_bin,'r'); xlabel('time (h)');
            subplot(212); loglog(f_bin, F_bin*2,'g'); hold on; 
%             set(gca,'XScale','log')
%             set(gca,'YScale','log')

%             pow_time = sum(abs(x_bin).^2)/length(x_bin)
%             pow_freq = sum(F_bin) * df
%             differen = pow_time - pow_freq

            

        end;
        nbins = nbins + 1;
        tcent_bins{nbins} = (curr_index + length_bin/2) * dt;
        f_bins{nbins} = f_bin;
        F_bins{nbins} = F_bin;
        t_bins{nbins} = t_bin;
        x_bins{nbins} = x_bin;
        curr_index = curr_index + length_bin;
        
    end
    
    fprintf (['Number of timebins ' num2str(nbins) '.\n']);
    

end


function fileNums = filenames2numbers(fileNames)

    fileNums=[];
    for i = 1:length(fileNames)
        temp = fileNames{i};
        fileNums(i) = str2num(fileNames{i}(1:4));
    end

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
        
        Fout(i) = sum(F(index)) * df;
        
        if plot_on
            hold on; plot(f(index),F(index),[colourarr(mod(i,length(colourarr))) '.']);
            hold on; plot(f_centers(i),Fout(i),[colourarr(mod(i+1,length(colourarr))) 'o']);
        end
    end
    
    fout = f_centers;

end

