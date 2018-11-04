function h = plotTM (traw,xraw,TM)
% Plots time series based on time range matrix TM.
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

h = plotIM (traw,xraw,IM);

end