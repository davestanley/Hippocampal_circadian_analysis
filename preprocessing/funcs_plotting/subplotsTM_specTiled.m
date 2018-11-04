
function [h, hxp] = subplotsTM_specTiled(traw,xraw,TM)
% Like figsTM_specTiled, but tiles subplots across the screen, instead of
% separate figuresa
%

% Convert TM time range to IM (index matrix)
IM = false(size(TM,1),length(traw));
for i = 1:size(TM,1)
    % Calculate IM matrix based on time ranges
    IM(i,:) = traw >= TM(i,1) & traw <= TM(i,2);
end


[h, hxp] = subplotsIM_spectTiled(traw,xraw,IM);

end
