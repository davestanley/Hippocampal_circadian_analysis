
function [h, hxp] = subplotsIM_spectTiled(traw,xraw,IM)
% Like figsIM_spectTiled, but tiles subplots across the screen, instead of
% separate figuresa
%

% Make sure IM are aligned as rows
if isvector(IM)
    IM = IM(:)';
end

colourarr = get_colourarr;


h = figure('units','normalized','outerposition',[0.00,0.0,0.8,0.6],'Color','w');
set(h,'Visible','off');   % Hide figure until it's fully plotted
Ncols = size(IM,1);
hxp = subplot_grid(2,Ncols);
for i = 1:Ncols
    ind = IM(i,:);
    tseg = traw(ind); tseg = tseg - tseg(1); xseg = xraw(ind);      % Pull out data and align to t=0;
    fs = 1/mode(diff(tseg));
    hxp.set_gca(1,i); plot(tseg,xseg,colourarr(i)); set(gca,'YTickLabels',{});
    hxp.set_gca(2,i); plott_spect(tseg,xseg,'logplot',1,'mode',2,'Nwind',floor(2*fs),'fract_overlap',0.5,'axis_lims',[0 100]);
    title('Log Power');
    xlabel('Time (s)'); ylabel('');
    
    % Sync axes
    hxp.sync_axes([hxp.hax(1,i),hxp.hax(2,i)],'x');
    
    % Set y label
    if i == 1
        %ylabel('Freq (Hz)');
    end
end

%hcurr.figtitle('Time series vs Log Power Spectrogram');
hxp.rowtitles({'','Freq (Hz)'});
set(h,'Visible','on');
    
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

