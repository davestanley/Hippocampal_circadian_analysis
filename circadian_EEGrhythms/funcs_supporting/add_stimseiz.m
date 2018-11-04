

function add_stimseiz(ratN,varargin)
    temp = get(gca,'YLim');
    ymin=temp(1);ymax=temp(2);
    temp = get(gca,'XLim');
    xmin=temp(1);xmax=temp(2);
    [stim seiz] = get_stimtime2(ratN);    
    hold on; plot([stim stim],[ymin ymax],varargin{:})
    hold on; plot([seiz seiz],[ymin ymax],varargin{:})
    for i = xmin:xmax; hold on; plot([i i],[ymin ymax],'k:','LineWidth',1); end

end
