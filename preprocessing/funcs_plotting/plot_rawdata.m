

function plot_rawdata (traw,xraw)

    set(gcf,'Color','w');
    plot(traw,xraw,'k');
    ox=6.45; oy=-2.2; lenx=1/24; leny=0.5; 
    plot_scale([(ox+lenx/2) oy], [1], [lenx],'k','sec','h');
    plot_scale([ox (oy+leny/2)], [1], [leny],'k','mV','v');
    xlims = get(gca,'XLim');
    ylims = get(gca,'YLim');
    set(gca,'YLim',[-2.5 2]);
    set(gca,'Box','off');
    set(gca,'Visible','off');
    text(6.45,-2.3,'1 hour')
    text(6.44,-2.2,'0.25 mV','Rotation',90)

end