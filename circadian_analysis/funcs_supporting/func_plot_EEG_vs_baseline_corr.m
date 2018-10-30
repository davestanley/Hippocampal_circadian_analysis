function func_plot_EEG_vs_baseline_corr(rcell,prepostchronic, smoothmode_ts)

        plottest_timeseries = 0;        
        
        [T X] = rat_2_ts (rcell,prepostchronic, smoothmode_ts);
        [Tb Xb] = rat_2_ts (rcell,prepostchronic, 2);               % Get baseline
        
        Nrats = length(X);
        Nbands=size(X{1},2);
   
        correlations = zeros(Nrats,Nbands,Nbands);
        for i = 1:Nrats
            if i == 1; rrat = 4; elseif i == 2; rrat = 9; elseif i == 3; rrat = 10; elseif i == 4; rrat = 1; end
            for j = 1:Nbands
                for k = 1:Nbands

                    T1 = T{i}(:,j);
                    X1 = X{i}(:,j);
                    T2 = Tb{i}(:,k);
                    X2 = Xb{i}(:,k);
                    
                    
                    if plottest_timeseries % && j==k
                        figure;
                        hold on; plot(T1,X1);hold on; plot(T2,X2,'r')
                        add_stimseiz(rrat,'k','LineWidth',2')
                    end

                    [T1,T2,X1,X2] = interpmin(T1,T2,X1,X2);
                    
                    if plottest_timeseries && j==k
                        hold on;
                        plot(T1,X1,'b.');hold on; plot(T2,X2,'r.')
                        add_stimseiz(rrat,'k','LineWidth',2')
                    end
                    

%                     [T1 X1] = calc_ts_func(T1,X1,2,0.5,0.5,fhandle_AMP);
%                     if plottest_timeseries; hold on; plot(T1, X1,'b.'); end
%                     [T2 X2] = calc_ts_func(T2,X2,2,0.5,0.5,fhandle_PHASE);
%                     if plottest_timeseries; hold on; plot(T2, X2,'r.'); end
% 
                    % Remove NaN's
                    index = ~isnan(X1) & ~isnan(X2); 
                    X1 = X1(index); X2 = X2(index); T1 = T1(index); T2 = T2(index);

                    correlations(j,k,i) = xcorr(X1-mean(X1),X2-mean(X2),0,'coeff');
                    correlations(j,k,i)
                end
            end
        end
        
        for ii = 1:Nrats
            figure;
            imagesc((squeeze(correlations(:,:,ii))));
            colormap(winter);colorbar; set(gca,'YDir','normal','Visible','on','FontSize',20,'XTick',[]); 
            xlabel('Baseline Freq');ylabel('Circ Freq');
        end
end        





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



function [T1,T2,X1,X2] = interpmin(T1,T2,X1,X2)

    dt1 = mode(diff(T1));
    dt2 = mode(diff(T2));
    
    if dt1 < dt2
        X2 = interp1(T2,X2,T1);
        T2=T1;
    else
        X1 = interp1(T1,X1,T2);
        T1 = T2;
    end

end