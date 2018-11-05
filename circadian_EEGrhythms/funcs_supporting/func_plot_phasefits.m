
function func_plot_phasefits(rcell_curr,analyze_only_seizing,freq_range)

    if ~exist('freq_range','var'); freq_range = []; end

    w = 2*pi/1.0;
    alpha = .05;
    
    %N=5
    fs_subplot = 12;
    
    Nfreqs = length(rcell_curr{1});
    if isempty(freq_range); freq_range=1:Nfreqs; end
    N=length(freq_range);
    

    for j = 1:length(rcell_curr)
        r = rcell_curr{j};
        
        % Open new figure
        if N > 1
            figure('Color','w','Position',[1500 10  757 800]);
        else
            figure('Color','w','Position',[ 622   544   768   240]);
        end
        
        % Subplot for each frequency range
        for i = 1:length(freq_range)
            [temp pd] = cosinor_struct(r{freq_range(i)}.ctrl.t,r{freq_range(i)}.ctrl.d,w,alpha,0); ax1 = subplot(N,3,(N*3+1)-((3)*(i-1)+3)); title(num2str(3*(i-1)+1)); plot_cosinor(pd.t,pd.y,pd.f,w); set(gca,'FontSize',fs_subplot)
            if (N*3+1)-((3)*(i-1)+3) == N*3-2
                if N>1
                    xlabel('Healthy','Color','b');
                else
                    legend('Healthy');
                end
            end        
            [temp pd] = cosinor_struct(r{freq_range(i)}.acute.t,r{freq_range(i)}.acute.d,w,alpha,0); ax2 = subplot(N,3,(N*3+1)-((3)*(i-1)+2)); title(num2str(3*(i-1)+2)); plot_cosinor(pd.t,pd.y,pd.f,w); set(gca,'FontSize',fs_subplot)
            if (N*3+1)-((3)*(i-1)+3) == N*3-2
                if N > 1
                    xlabel('Latent','Color','m');
                else
                    legend('Latent');
                end
            end
            ax3 = [];
            if analyze_only_seizing == 1; [temp pd] = cosinor_struct(r{freq_range(i)}.chr.t,r{freq_range(i)}.chr.d,w,alpha,0); ax3 = subplot(N,3,(N*3+1)-((3)*(i-1)+1)); title(num2str(3*(i-1)+3)); plot_cosinor(pd.t,pd.y,pd.f,w); set(gca,'FontSize',fs_subplot);
                if (N*3+1)-((3)*(i-1)+3) == N*3-2
                    if N>1
                        xlabel('Seizing','Color','r');
                    else
                        legend('Seizing');
                    end
                end
            end
            linkaxes([ax1,ax2,ax3]);
        end
    end

end