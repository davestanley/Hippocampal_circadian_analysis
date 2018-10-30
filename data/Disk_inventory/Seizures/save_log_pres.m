
% This file reads in the log file and saves the time and other information
% to be read in later


function [outlog_path outlog_name]= save_log_pres(path_ratlog, ratNum, outlog_path)
%Ratnum must be the string with the rat number like '004' and '008'

logName = strcat(path_ratlog,'/','Rat',ratNum,'EEG_best.log');
fid = fopen(logName,'r');
if fid == -1
    fprintf (['Error, log file ' logName ' not found. Exiting! \n'])
    return
end
format = '%s %s %f %f %f %f %f %f %s %s %s %s %s %s %d32 %s %s';
format = strcat('Rat',ratNum,'ch01F',format);
A=textscan(fid,format,'delimiter',':- \b\t','multipledelimsasone',1,'Headerlines',1);
fclose(fid);

fileNames=A{1};
dates=A{2};
hrStart=A{3};
minStart=A{4};
secStart=A{5};
hrEnd = A{6};
minEnd = A{7};
secEnd = A{8};
sampFreq = A{12};
numSamples = A{15};
[numFiles,n]=size(numSamples);

for i = 1:length(dates);
    temp=dates{i};
    year(i)=str2num(temp(7:8));
    month(i)=str2num(temp(1:2));
    day(i)=str2num(temp(4:5));
    tabs(i) = datenum([(year(i)) (month(i)) (day(i)) (hrStart(i)) (minStart(i)) (secStart(i))]);
    trel(i) = datenum([0 0 0 (hrStart(i)) (minStart(i)) (secStart(i))]);
    freq(i) = str2num(sampFreq{i});
end
clear temp


for ii = 1:length(fileNames);
    fileNums(ii) = str2num(fileNames{ii}(1:4));
end


% outlog_path = './RatoutMat';
[s,mess,messid] = mkdir(outlog_path);
outlog_name = strcat('Rat',ratNum,'log');
save([outlog_path '/' outlog_name],'fileNames','dates','hrStart','minStart','secStart','hrEnd','minEnd','secEnd','freq', 'numSamples','numFiles','year','month','day','tabs','trel','fileNums')



%cd ~/rat/RatData_out/Rat004/bin1000/
%read in the same file

%fid = fopen('Rat004ch01F0012_stdev.bin','rb');   %open file
%data1 = fread(fid, inf, 'double');  %read in the data
%fclose(fid);  %close file
%fid = fopen('Rat004ch01F0013_stdev.bin','rb');   %open file
%data2 = fread(fid, inf, 'double');  %read in the data
%fclose(fid);   %close file
%cd ~/Documents/Ma'tlab Files'/