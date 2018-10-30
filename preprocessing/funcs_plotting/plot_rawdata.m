

function plot_rawdata (traw,xraw)

    set(gcf,'Color','w');
    plot(traw,xraw,'k');
    xl = xlim; yl=ylim;
    %ox=6.45; oy=-2.2; 
    ox=xl(1); oy=yl(1); lenx=1/24; leny=0.5; 
    plot_scale([(ox+lenx/2) oy], [1], [lenx],'k','sec','h');
    plot_scale([ox (oy+leny/2)], [1], [leny],'k','mV','v');
    xlims = get(gca,'XLim');
    ylims = get(gca,'YLim');
    set(gca,'YLim',[-2.5 2]);
    set(gca,'Box','off');
    set(gca,'Visible','off');
    text(ox,oy-0.1,'1 hour')
    text(ox-0.01,oy,'0.25 mV','Rotation',90)
    

end