


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

% Change to correct working directory
cd(fullfile(userdir,'GDrive','from_Dropbox','Career','2018','Insight Data Science Boston','demo'));

% Setup paths
restoredefaultpath
addpath(genpath(fullfile(userdir,'GDrive','from_Dropbox','MATLAB')))
addpath(genpath('~/src/ds_kb3/funcs_general'))
addpath(genpath('~/src/chronux'))
addpath(genpath('./funcs_supporting'));
addpath(genpath('./funcs_plotting'));


%% Load and plot 7 hours raw data

ratN = 4;
chanN = 2;
offset = 0;
len = 10000000;

sampling = 100;
traw=[];xraw=[];
[traw, xraw] = extract_anchored_timeseries(ratN,chanN,21,offset,len,sampling-1);

figure; plot_rawdata (traw,xraw);



%% Load data and plot frequency decomposition

% Load data
ratN = 4;
chanN = 2;
fs=round(12207/12);
offset = round(30*fs); len = round(12*fs);
sampling = 1;
[t x] = extract_anchored_timeseries(ratN,chanN,22,offset,len,sampling-1);
t = t - t(1); t = t*24*3600;
dt = mode(diff(t));
fs = 1/dt;

% Plot decomposition
figure; xarr = plot_frequencydecomp(t,x);





%% Plot EMD
addpath('./rParabEmd__L');

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
hold on; plot_thetahighlight(t,x,theta_delta_ratio); pause

% Plot spectrogram
%figure('Position',[260 334 1104 450]); plot_all(t,x,'psd_on',0,'axis_lims',[0 20])

%% Theta state manual detection

path_theta = (fullfile(userdir,'Crystalized','Anco','Evol','EEG','RatFFT','03_theta_delta_2sec_bins','Manual_compare_theta'));

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





