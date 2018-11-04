

function plot_rawdata (traw,xraw,axis_on,color)

    if nargin < 3
        axis_on = 1;
    end
    
    if nargin < 4;
        color = 'k'; 
    end
    
    set(gcf,'Color','w');
    plot(traw,xraw,color);
    
    % Get axis limits
    xl = xlim;
    %yl=ylim;
    yl = [min(xraw) max(xraw)];
    
    % We'll work in units of 1/7th the axis length
    scalex = (xl(2) - xl(1))/7;     % Take 1/7th of the x-axis limits. This should be about 1 hour
    scaley = (yl(2) - yl(1))/7;
    
    % Expand axis limits by 10% in bottom-left direction, to create space
    % for scalebars
    %xlim([xl(1)-scalex, xl(2)]);
    ylim([yl(1)-scaley, yl(2)]);
    xl = xlim; yl=ylim;
    
    
    ox=xl(1)+scalex/2; oy=yl(1)+scaley/2; lenx=1/24; leny=(scaley);  % x-axis is in units of days
    plot_scale([(ox+lenx/2) oy], [1], [lenx],'k','sec','h');
    plot_scale([ox (oy+leny/2)], [1], [leny],'k','mV','v');
    xlims = get(gca,'XLim');
    ylims = get(gca,'YLim');
    text(ox,oy-scaley*0.1,'1 hour')
    text(ox-scalex*0.1,oy,'0.25 mV','Rotation',90)
    
    if axis_on
        set(gca,'YTickLabels',{});
        xlabel('Time (days)');
    else
        set(gca,'Box','off');
        set(gca,'Visible','off');
    end
    

end