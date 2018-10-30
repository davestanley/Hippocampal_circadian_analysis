
function func_plot_correlation(phi, Amp, colour_by_freqband)

    markerarr = 'x+od';
    colour_by_freqband = 1;
    
%     if uselog Amp = log(Amp);
%     else Amp = Amp; end
%     
%     Amp = Amp ./ repmat(Amp(:,:,1),[1 1 3]);
%     
    figure;
    freq_range = 1:size(Amp,1);
    animal_range = 1:size(Amp,2);
    %animal_range = 4;
    stages_2_plot = 2:3;

    if colour_by_freqband
        colourarr = summer;
        index = linspace(1,size(colourarr,1),size(Amp,1));
        colourarr = colourarr(index(:),:);
        % Plot the symbols
        for i = freq_range
            for j = animal_range
                amp_post = (Amp(i,j,stages_2_plot)); phi_post = phi(i,j,stages_2_plot);
                index = ~isnan(amp_post); amp_post = amp_post(index); phi_post = phi_post(index);
                amp_post = amp_post(:); phi_post = phi_post(:);
                hold on; plot(amp_post,phi_post,[markerarr(j)],'Color',colourarr(i,:),'LineWidth',2,'MarkerSize',9);
%                 if j == 4 % For Rat 4, highlight latency in red
%                     amp_post = (Amp(i,j,2)); phi_post = phi(i,j,2);
%                     index = ~isnan(amp_post); amp_post = amp_post(index); phi_post = phi_post(index);
%                     amp_post = amp_post(:); phi_post = phi_post(:);
%                     hold on; plot(amp_post,phi_post,[markerarr(j)],'Color','r','LineWidth',2,'MarkerSize',9);
%                 end
            end
        end
        legend('R1','R2','R3','R4','Location','SouthEast');
        % Plot the numbers
        for i = freq_range
            for j = animal_range
                amp_post = (Amp(i,j,stages_2_plot)); phi_post = phi(i,j,stages_2_plot);
                index = ~isnan(amp_post); amp_post = amp_post(index); phi_post = phi_post(index);
                amp_post = amp_post(:); phi_post = phi_post(:);
                for k = 1:length(amp_post)
                    percent_shift = 0.01; xshift = get(gca,'XLim'); yshift = get(gca,'YLim'); xshift = diff(xshift)*percent_shift; yshift = diff(yshift)*percent_shift;
                    hold on; text(amp_post(k)+xshift,phi_post(k)+yshift,num2str(i),'Color','k','FontSize',12);
                end
            end
        end
    else
        for j = animal_range
            amp_latent = (Amp(freq_range,j,2)); phi_latent = phi(freq_range,j,2);
            index = ~isnan(amp_latent); amp_latent = amp_latent(index); phi_latent = phi_latent(index);

            amp_chr = (Amp(freq_range,j,3)); phi_chr = phi(freq_range,j,3);
            index = ~isnan(amp_chr); amp_chr = amp_chr(index); phi_chr = phi_chr(index);

    %             figure; hold on; plot(amp_temp(index));
    %             hold on; plot(phi_temp(index),'r');
    %             title(['Corr = ' num2str(temp(1,2))])

            %hold on; plot(amp_temp,phi_temp,[colourarr(j) '.']);
            hold on; plot(amp_chr,phi_chr,['k' markerarr(j)],'LineWidth',2,'MarkerSize',9);
            hold on; plot(amp_latent,phi_latent,[markerarr(j)],'Color',[0.3 0.3 0.3],'LineWidth',2,'MarkerSize',9);

        end
        
    end

    amp_temp = (Amp(freq_range,animal_range,stages_2_plot));
    phi_temp = phi(freq_range,animal_range,stages_2_plot);
    index = ~isnan(amp_temp(:)); amp_temp = amp_temp(index); phi_temp = phi_temp(index);

    p = polyfit(amp_temp,phi_temp,1);
    xplotvals = [min(amp_temp) max(amp_temp)];
    hold on; plot(xplotvals, polyval(p,xplotvals),'k','LineWidth',2);
    [cR,cP,cRLO,cRUP] = corrcoef(amp_temp,phi_temp);
    cR = cR(1,2); cP = cP(1,2); cRLO = cRLO(1,2); cRUP = cRUP(1,2); cCI = (cRUP - cRLO)/2;
    title(['r = ' num2str(cR) '+/-' num2str(cCI) ', p = ' num2str(cP)]);
    set(gca,'FontSize',20)
    
    xlims = get(gca,'XLim'); ylims = get(gca,'YLim');
    hold on; plot([0 0],[ylims],'k');
    hold on; plot([xlims],[0 0],'k');
    
    xlabel('Amp (normalized) ');ylabel('Phi');

%     figure; plot(zscore(amp_temp(index)),'b');
%     hold on; plot(zscore(phi_temp(index)),'g');                
%     hold on; plot((amp_temp),(phi_temp),'r.');
    
end

