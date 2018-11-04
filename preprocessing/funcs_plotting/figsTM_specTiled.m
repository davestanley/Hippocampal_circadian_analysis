
function h = figsTM_specTiled(traw,xraw,TM)
% Tiles time series and spectrogram plots across the screen, each
% plot corresponding to a row of time ranges in TM
% 
% Part of a series of functions designed to operate on time range matrix, TM.
% TM is a matrix whereby each row contains a pair of values representing
% time ranges to plot between. Applied to the input time series
% data (traw, xraw)
%

% Convert TM time range to IM (index matrix)
IM = false(size(TM,1),length(traw));
for i = 1:size(TM,1)
    % Calculate IM matrix based on time ranges
    IM(i,:) = traw >= TM(i,1) & traw <= TM(i,2);
end


h = figsIM_spectTiled(traw,xraw,IM);

end
