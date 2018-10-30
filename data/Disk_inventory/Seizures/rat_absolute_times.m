
% Returns the "aboslute" times for each rat file.
% This is based on the content of log files 'RatXYZEEG_best.log' in the
% logs_corrected folder.
% INPUT
% ratN = input rat number, either 1, 4, 5, 6, 7, 8, 9, 10, 11, or 12
% 
% OUTPUT
% fileNames = cell array of file names
% fileNums = corresponding matrix of file numbers
% tabs = corresponding matrix of start times (in days), with day "0" corresponding to the starting time of the 1st file.
% 
% EXAMPLE CODE
% [fileNames, fileNums, tabs] = rat_absolute_times(4)
% 
% 

function [fileNames, fileNums, tabs] = rat_absolute_times (ratN)

    if (~exist('ratN','var')); ratN=4; end
    
    start_from_scratch = 1;
    path_ratlog = ('./logs_corrected');
    
    ratnum = num2str(ratN, '%6.3d');
    
    
    outlog_path = ['./temp'];
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
    
end