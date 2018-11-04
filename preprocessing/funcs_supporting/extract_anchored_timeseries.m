
function [t,x] = extract_anchored_timeseries(ratN,chanN,fileN,offset,len,skip)
    % % Return values % %
    % x - data values
    % t - time in days
    plot_on = 0;
    if (~exist('ratN','var')); ratN=4; end
    if (~exist('chanN','var')); chanN=2; end
    if (~exist('fileN','var')); fileN=26; end
    if (~exist('offset','var')); offset=0; end
    if (~exist('len','var')); len=26; end       % Length of data to read in samples
    if (~exist('skip','var')); skip=0; end      % Number of datapoints to skip in samples
    
    
    %path_ratlog='../RatData_out';
    path_ratlog = ('../data/Disk_inventory/RatData_out');

    ratnum = num2str(ratN, '%6.3d');
    channum = num2str(chanN, '%6.2d');
    filenum = num2str(fileN,'%6.4d');

    fname = ['../data/Raw_ln/Evol/RawOneK/Rat' ratnum 'OneK/Rat' ratnum 'ch' channum 'F' filenum '_DownSampled_Dec.bin'];              % If mounted data
    if ~exist(fname,'file')
        fname = ['../data/Raw/Evol/RawOneK_partial/Rat004OneK/Rat' ratnum 'ch' channum 'F' filenum '_DownSampled_Dec.bin'];         % For local copy of data
    end
    
    
    fs=round(12207/12);
    
    
    outlog_path = ['./Ratoutmat_FFT'];
    [outlog_path outlog_name] = save_log(path_ratlog, ratnum, outlog_path);   % Save log file
    
    load ([outlog_path '/' outlog_name])    % Load log file
    
    for ii = 1:length(fileNames);
        fileNums(ii) = str2num(fileNames{ii}(1:4));
    end
    
    index = find(fileN == fileNums,1,'first');
    tabs_days = tabs;
    tstart = tabs_days(index) - floor(tabs_days(1));

    
    
%     dt = 1/fs;
%     [x]= readbin(fname, 'int16', 0, 10000,0) / 3276700 * 1e3;
%     x = x - mean(x);
%     t = [0:length(x)-1];
%     t = t * dt;
%     
%     ds = 10;
%     x = downsample(x,ds);
%     t = downsample(t,ds);
%     
% 
%     figure; plot(t/3600,x,'b');

    
    
    %skip = 9;
    dt = 1/fs*(skip+1);
    [x]= readbin(fname, 'int16', offset, len,skip) / 3276700 * 1e3;
    x = x - mean(x);
    t = [0:length(x)-1];
    t = t * dt;
    t = t /3600/24;
    t = t + tstart;
    
%     ds = 100;
%     x = downsample(x,ds);
%     t = downsample(t,ds);
    
    if plot_on
        figure; hold on; plot(t,x,'r');
    end

%     tstart

end

