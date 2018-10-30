function xarr = plot_frequencydecomp(t,x)
freq_bins = [1 5; 5 10; 12 25; 25 50; 50 90; 90 140; 140 200; 200 500];
[xarr] = filter_freqbands (t,x,freq_bins);

xarr = [xarr x(:)];
xarr(:,1) = xarr(:,1) - 0.5;
xarr(:,2) = xarr(:,2) - 0.25;
xarr(:,end) = xarr(:,end) + 0.5;
xarr(:,7:8) = xarr(:,7:8)*5;

set(gcf,'Color','w','Position',[ 131         291        1192         457]); hold on;
%plot_matrix2(t,xarr,os,'k');
plot_matrix3D(t,xarr,'showErrorbars',0,'do_shift',0.5,'active_dim',3,'LineSpec',{'k'})

ox=0.5; oy=-1.75; lenx=0.5; leny=0.5; 
plot_scale([(ox+lenx/2) oy], [1], [lenx],'k','sec','h');
plot_scale([ox (oy+leny/2)], [1], [leny],'k','mV','v');
xlims = get(gca,'XLim');
ylims = get(gca,'YLim');
set(gca,'YLim',[-2 6]);
set(gca,'XLim',[0 14]);
set(gca,'Box','off');
set(gca,'Visible','off');
text(.5,-2.0,'0.5 sec');
%text(.2,-1.6,'0.5 sec)','Rotation',90);

textstrings = arrayfunu(@(x) num2str(x), freq_bins);
addtext(12.3,mean(xarr(:,1:end-1)),0.5,textstrings);




end


function addtext(xpos,yoffsets,spacing,textstrings)
    ypos = 0:size(textstrings,1)-1;
    ypos = (ypos*spacing+yoffsets);
    
    
    for i = 1:length(ypos)
        text(xpos,ypos(i),[textstrings{i,1} '-' textstrings{i,2} 'Hz']);
    end
end