    

function plot_cosinor (t,y,f,w)

    apply_modulus = 1;


%         xlabel('x-axis');
%         ylabel('y-axis');

    if apply_modulus
        period = 2*pi/w;
        t = mod(t,period)';
        Y = [t(:) f(:)];
        Y = sortrows(Y,1);
        plot(t*24,y./1,'k.','MarkerSize',12); hold on;
        plot(Y(:,1)*24,Y(:,2)./1,'r','LineWidth',2); hold on;
    else
        plot(t,y./1,'k.','MarkerSize',12); hold on;
        plot(t,f/1,'r','LineWidth',2); hold on;
    end

        %legend('Original', 'Cosinor');
        xlim([min(t) max(t)*24]);
        %set(gca,'FontSize',30)
        set(gcf,'Color','w');
        set(gca,'Box','off');
        set(gca,'XTick',[0 6 12 18 24]);
        set(gca,'XLim',[0 24])
end