
function bounds = get_bounds(r)
    %bounds = zeros(3,2);
    
    % Use new mode to determine bounds - baesd on values stored in the structure
    if isfield(r.all,'tabs_acute_1')
        bounds(1,:) = [r.all.tabs_ctrl_1 r.all.tabs_ctrl_2];
        bounds(2,:) = [r.all.tabs_acute_1 r.all.tabs_acute_2];
        if ~isempty(r.chr.t) bounds(3,:) = [r.all.tabs_chr_1 r.all.tabs_chr_2]; end
        if isfield(r,'dark') bounds(4,:) = [r.all.tabs_dark_1 r.all.tabs_dark_2]; end
    else
%         If the structure is an old version, extract the bounds based on the time values.
        bounds(1,:) = [r.ctrl.t(1) r.ctrl.t(end)];
        bounds(2,:) = [r.acute.t(1) r.acute.t(end)];
        if ~isempty(r.chr.t)
            bounds(3,:) = [r.chr.t(1) r.chr.t(end)];
        else
            %bounds(3,:) = [r.acute.t(end) r.acute.t(end)];
        end
        if isfield(r,'dark')
            bounds(4,:) = [r.dark.t(end) r.dark.t(end)];
        else
            %bounds(4,:) = [r.acute.t(end) r.acute.t(end)];
        end
    end
end

