


function [seizure_times] = Seizure_summary

    addpath('./cosinor');

    sz_r4 = get_sz_tod (4,'./R004FileInfo.mat');
    sz_r9 = get_sz_tod (9,'./R009FileInfo.mat');
    sz_r10 = get_sz_tod (10,'./R010FileInfo.mat');
    sz_r1 = get_sz_tod (4,'./R001FileInfo.mat');
    
    all_sz = [sz_r4;sz_r9;sz_r10;sz_r1];
    
    [n, xout] = hist(all_sz,6);
    figure; bar(xout, n);
    xlabel('time (h)','FontSize',20);
    ylabel('number of seizures','FontSize',20);
    set(gca,'FontSize',20);
    set(gcf,'Color','w');
    
    seizure_times.r4 = sz_r4;
    seizure_times.r9 = sz_r9;
    seizure_times.r10 = sz_r10;
    seizure_times.r1 = sz_r1;

end


function [sz_tod] = get_sz_tod (ratN, fileinfo)

    if (~exist('ratN','var')); ratN=10; end
    if (~exist('chanN','var')); chanN=2; end
    if (~exist('pattern','var')); pattern=1; end                            % Pattern=1 for downwards, 2 for upwards; although this is swapped for Rat 004
    if (~exist('chosen_simversion','var')); chosen_simversion=3; end        % Version4 = all spikes merged; Version3 = up and downwards separated
    if (~exist('theta_or_nontheta','var')); theta_or_nontheta=2; end        % =1 for only SPKs in theta epochs; =2 for only SPKs in non-theta epochs; <=0 for no separation

    
    start_from_scratch = 1;
    path_ratlog = ('./logs_corrected');
    
    path_savephase = './Phase_Roll';


    %ratN=10;
    ratnum = num2str(ratN, '%6.3d');
    channum = num2str(chanN, '%6.2d');
    
    outlog_path = ['./Ratoutmat'];
    if start_from_scratch
        [outlog_path outlog_name] = save_log_pres(path_ratlog, ratnum, outlog_path);   % Save log file
    else
        outlog_name = strcat('Rat',ratnum,'log');
    end
    load ([outlog_path '/' outlog_name])    % Load log file

    for ii = 1:length(fileNames);
        fileNums(ii) = str2num(fileNames{ii}(1:4));
    end
    
    tabs = tabs - floor(tabs(1));
    tabs = tabs * 24;
    
    load(fileinfo)
    
    sz_filestarts = [];
    sz_offsets = [];
    for i = 1:size(Sz,1)
        index = find(fileNums == Sz(i,1),1,'first');
        sz_filestarts = [sz_filestarts; tabs(index)];
        sz_offsets = [sz_offsets; Sz(i,2)/3600];
    end

    sz_tod = sz_filestarts + sz_offsets;
    sz_filestarts = mod(sz_filestarts,24);
    sz_tod = mod(sz_tod, 24);
    
end

