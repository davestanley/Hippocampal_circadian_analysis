
function h = figsIM_spectTiled(traw,xraw,IM)
% Tiles time series and spectrogram plots across the screen, each
% plot corresponding to a row of indices in IM
% 
% Part of a series of functions designed to operate on index matrix, IM.
% IM is a matrix whereby each row is an index vector, used to index input
% data (traw, xraw)
%

% Make sure IM are aligned as rows
if isvector(IM)
    IM = IM(:)';
end

colourarr = get_colourarr;

h = cell(1,size(IM,1));

for i = 1:size(IM,1)
    ind = IM(i,:);
    
    tseg = traw(ind); tseg = tseg - tseg(1); xseg = xraw(ind);      % Pull out data and align to t=0;
    fs = 1/mode(diff(tseg));
    h{i} = figure('Position', [ 50+300*(i-1), 60, 582, 746]); 
    subplotrows(2,1);plot(tseg,xseg,colourarr(i));
    subplotrows(2,2);plott_spect(tseg,xseg,'logplot',1,'mode',2,'Nwind',floor(2*fs),'fract_overlap',0.5,'axis_lims',[0 100]);
    title('Log Power');
    xlabel('Time (s)'); ylabel('Freq (Hz)');


end

    
% 
% ind = traw >= 11.254 & traw <= 11.264;
% tseg = traw(ind)*24*3600; tseg = tseg - tseg(1); xseg = xraw(ind);
% fs = 1/mode(diff(tseg));
% figure('Position', [ 50, 60, 582, 746]); 
% subplotrows(2,1); plot(tseg,xseg);
% subplotrows(2,2);plott_spect(tseg,xseg,'logplot',1,'mode',2,'Nwind',floor(10*fs),'fract_overlap',0.5);
% xlabel('Time (s)'); ylabel('Freq (Hz)');
% 
% 
% ind = traw >= 11.264 & traw <= 11.274;
% tseg = traw(ind)*24*3600; tseg = tseg - tseg(1); xseg = xraw(ind);
% fs = 1/mode(diff(tseg));
% figure('Position', [ 450, 60, 582, 746]); 
% subplotrows(2,1); plot(tseg,xseg);
% subplotrows(2,2);plott_spect(tseg,xseg,'logplot',1,'mode',2,'Nwind',floor(10*fs),'fract_overlap',0.5);
% xlabel('Time (s)'); ylabel('Freq (Hz)');
% 
% ind = traw >= 11.274 & traw <= 11.284;
% tseg = traw(ind)*24*3600; tseg = tseg - tseg(1); xseg = xraw(ind);
% fs = 1/mode(diff(tseg));
% figure('Position', [ 850, 60, 582, 746]); 
% subplotrows(2,1); plot(tseg,xseg);
% subplotrows(2,2);plott_spect(tseg,xseg,'logplot',1,'mode',2,'Nwind',floor(10*fs),'fract_overlap',0.5);
% xlabel('Time (s)'); ylabel('Freq (Hz)');




end

