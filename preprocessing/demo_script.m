


%% Setup plotting defaults and path
set(0, 'DefaultFigureColor', 'White', ...
    'DefaultTextFontSize', 10, ...
    'DefaultAxesFontSize', 10);
set(0,'DefaultFigurePosition',[190   286   715   520]) % Above command window
clear all
emd_from_scratch = 0;

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

addpath(genpath('funcs_supporting'));         % Supporting processing functions
addpath(genpath('funcs_plotting'));           % Plotting functions
addpath('rParabEmd__L');                      % EMD code


%% Load and plot a singe file (7 hours raw data)

% Chosen rat, channel, and file
ratN = 4;
chanN = 2;
fnum = 22;
offset = 0;
len = Inf;

% Load highly decimated data
skip = 9;                                                                           % Extract every tenth data point to save space (fs = 101.7 Hz)                                       
traw=[];xraw=[];
[traw, xraw] = extract_anchored_timeseries(ratN,chanN,fnum,offset,len,skip);

% Load data with denser sampling
skip = 1;                                                                           % Extract every second data point (fs = 508.5 Hz). Note, this will introduce some aliasing
[traw2, xraw2] = extract_anchored_timeseries(ratN,chanN,fnum,offset,len,skip);
% dt = mode(diff(traw2))*24*3600; df = 1/dt

% Plot full dataset
figure; plot_rawdata2 (traw,xraw);
keyboardnavigate on

%% Plot subset of data showing delta-theta transitions
% (See start of green traces)

% Specify region of interest
tstart = 6.876; tjump = 0.005; Nblocks = 2; tskip = 0.005;
TM = buildTimeMatrix(Nblocks,tstart,tjump,tskip);

% Plot showing region of interest
figure; plot_rawdata2 (traw,xraw);
hold on; plotTM(traw,xraw,TM);
keyboardnavigate on

% Spectrograms
subplotsTM_specTiled(traw2*24*3600,xraw2,TM*24*3600);


%% Plot subset of data showing delta-theta transitions
% (See green traces)

% Specify region of interest
tstart = 6.767; tjump = 0.001; Nblocks = 4; tskip=0;
TM = buildTimeMatrix(Nblocks,tstart,tjump,tskip);

% Plot showing region of interest
figure; plot_rawdata2 (traw,xraw);
hold on; plotTM(traw,xraw,TM);
keyboardnavigate on

% Spectrograms
subplotsTM_specTiled(traw2*24*3600,xraw2,TM*24*3600);



%% Load a bad file (showing 60 Hz corruption)
% Bad files for rat9 are: [31    73    96   142];

% Chosen rat, channel, and file
ratN = 9;
chanN = 2;
fnum = 31;
offset = 0;
len = Inf;

% Load highly decimated data
skip = 9;
[traw, xraw] = extract_anchored_timeseries(ratN,chanN,fnum,offset,len,skip);

% Load data with denser sampling
skip = 4;                                                                       % Extract every fifth data point to save space (fs = 203.4 Hz). Note, this will introduce some aliasing
[traw2, xraw2] = extract_anchored_timeseries(ratN,chanN,fnum,offset,len,skip);
% dt = mode(diff(traw2))*24*3600; df = 1/dt

% Plot full dataset
figure; plot_rawdata2 (traw,xraw);
keyboardnavigate on

%% Plot a bad region of data

% Specify region of interest
tstart = 11.255; tjump = 0.01; Nblocks = 3;
TM = buildTimeMatrix(Nblocks,tstart,tjump);

% Plot showing region of interest
figure; plot_rawdata2 (traw,xraw);
hold on; plotTM(traw,xraw,TM);
keyboardnavigate on

% Spectrograms
subplotsTM_specTiled(traw2*24*3600,xraw2,TM*24*3600);


%% Plot a bad region of data

% Specify region of interest
tstart = 11.34; tjump = 0.003; Nblocks = 3; tskip = 0.007;
TM = buildTimeMatrix(Nblocks,tstart,tjump,tskip);

% Plot showing region of interest
figure; plot_rawdata2 (traw,xraw);
hold on; plotTM(traw,xraw,TM);
keyboardnavigate on

% Spectrograms
subplotsTM_specTiled(traw2*24*3600,xraw2,TM*24*3600);


%% Load and plot a bad file - Seizure
% Bad files for rat9 are: [31    73    96   142];

% Chosen rat, channel, and file
ratN = 9;
chanN = 2;
fnum = 96;
offset = 0;
len = Inf;

% Load highly decimated data
skip = 9;                                                                       % Extract every tenth data point to save space (fs = 101.7 Hz)
[traw, xraw] = extract_anchored_timeseries(ratN,chanN,fnum,offset,len,skip);

% Load data with denser sampling
skip = 1;                                                                       % Extract every second data point (fs = 508.5 Hz). Note, this will introduce some aliasing
[traw2, xraw2] = extract_anchored_timeseries(ratN,chanN,fnum,offset,len,skip);
% dt = mode(diff(traw2))*24*3600; df = 1/dt

% Plot full dataset
figure; plot_rawdata2 (traw,xraw);
keyboardnavigate on


%% Plot 1st seizure

% Specify region of interest
tstart = 30.5115; tjump = 0.001; Nblocks = 4;
TM = buildTimeMatrix(Nblocks,tstart,tjump);

% Plot showing region of interest
figure; plot_rawdata2 (traw,xraw);
hold on; plotTM(traw,xraw,TM);
keyboardnavigate on

% Spectrograms
subplotsTM_specTiled(traw2*24*3600,xraw2,TM*24*3600);


%% Plot 2nd seizure

% Specify region of interest
tstart = 30.607; tjump = 0.001; Nblocks = 5;
TM = buildTimeMatrix(Nblocks,tstart,tjump);

% Plot showing region of interest
figure; plot_rawdata2 (traw,xraw);
hold on; plotTM(traw,xraw,TM);
keyboardnavigate on

% Spectrograms
subplotsTM_specTiled(traw2*24*3600,xraw2,TM*24*3600);



%% Plot good file with lots of SPWs
% Bad files for rat9 are: [31    73    96   142];

% Chosen rat, channel, and file
ratN = 9;
chanN = 2;
fnum = 80;
offset = 0;
len = Inf;

% Load highly decimated data
skip = 9;                                                                       % Extract every tenth data point to save space (fs = 101.7 Hz)
[traw, xraw] = extract_anchored_timeseries(ratN,chanN,fnum,offset,len,skip);

% Load data with denser sampling
skip = 1;                                                                       % Extract every second data point (fs = 508.5 Hz). Note, this will introduce some aliasing
[traw2, xraw2] = extract_anchored_timeseries(ratN,chanN,fnum,offset,len,skip);
% dt = mode(diff(traw2))*24*3600; df = 1/dt

% Plot full dataset
figure; plot_rawdata2 (traw,xraw);
keyboardnavigate on


%% Plot SPWs

tstart = 26.374; tjump = 0.002; Nblocks = 4;
TM = buildTimeMatrix(Nblocks,tstart,tjump);

% Plot showing region of interest
figure; plot_rawdata2 (traw,xraw);
hold on; plotTM(traw,xraw,TM);
keyboardnavigate on

% Spectrograms
subplotsTM_specTiled(traw2*24*3600,xraw2,TM*24*3600);


%% Load data and plot frequency decomposition

% Load data
ratN = 4;
chanN = 2;
fs=round(12207/12);
offset = round(30*fs); len = round(12*fs);
skip = 0;                                % Sampling rate of extracted data
[t, x] = extract_anchored_timeseries(ratN,chanN,22,offset,len,skip);
t = t - t(1); t = t*24*3600;            % Convert from days to seconds, starts at t=0
dt = mode(diff(t));
fs = 1/dt;

% Plot decomposition
figure; xarr = plot_frequencydecomp(t,x);


%% Plot EMD

if emd_from_scratch
    rParabEmd = rParabEmd__L (x, 40, 40, 1);
    save('emd.mat','rParabEmd')
else
    load emd.mat
end

% Plot EMD decomposition
figure; plot_matrix3D(rParabEmd,'active_dim',3,'do_shift',.3);
emdh = hilbert(rParabEmd);
emda = unwrap(angle(emdh))/2/pi;
instf = (diff(emda,[],1) / dt);
instt = repmat(t(1:end-1)',1,size(instf,2));

% pause
%%
% Plot EMD spectra
figure; plott_ani(rParabEmd,'fs',fs,'fname',{@plott_fs, @plott_spect},'plotargs',{'axis_lims',[0 200]})

figure; plott_ani_pairs(fliplr(xarr(:,1:8)),@plott_spect,rParabEmd(:,1:8),@plott_spect,'fsubplot',@subplotcols,'plotargs',{'fs',1017,'axis_lims',[0 260]});


%% Theta state autodetection

% Plot decomposition
figure; plot_frequencydecomp(t,x); pause

% Highlight theta
theta_delta_ratio = 1.5;        % Theta-over-delta power ratio
hold on; plot_thetahighlight(t,x,theta_delta_ratio);

% Plot spectrogram
%figure('Position',[260 334 1104 450]); plot_all(t,x,'psd_on',0,'axis_lims',[0 20])

%% Theta state manual detection

path_theta = (fullfile('..','data','Anco','Evol','EEG','RatFFT','03_theta_delta_2sec_bins','Manual_compare_theta'));

% Load a set of 2-second epochs
st = load(fullfile(path_theta,'RAll4_CT10_14.mat'));
tep = st.r09_pre.twin;
xep = zscore(st.r09_pre.xwin);
dt = mode(diff(tep(:,1)));
fs = 1/dt;

% Plot the first 10 entries
%figure; plot_matrix3D(tep(:,1),xep(:,1:10),'active_dim',3,'do_shift',3,'fs',fs); pause

% My manual classification (yes, no, yes, no)
figure('Position',[92 384 1364 356]); plot_ani(tep(:,1),xep(:,1:10),'fs',fs,'fname',@plot_all, 'plotargs',{'axis_lims',[0 20],'Nwind',0.5*fs,'fsubplot',@subplotcols});
% figure('Position',[92 384 1364 356]); plot_ani(tep(:,1),xep(:,1:10),'fs',fs,'fname',{@plot_fs, @plot_spect}, 'plotargs',{'axis_lims',[0 20]});


%% RoC
roc = load ('allrats_indiv_d1_5_t5_10_c');
plot_roc_graph(roc.predict_all_r1,roc.actual_all_r1,roc.theta_ratio_ranges0)





