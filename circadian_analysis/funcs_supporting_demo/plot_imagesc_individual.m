function plot_imagesc_individual(phi,Amp,compare_amp,analyze_only_seizing,za_fail,normalize_amps,rat2plot)
    global freq_listing
    
    if ~compare_amp
        for j = rat2plot
                figure; set(gcf,'Color','w');
                if analyze_only_seizing == 1
                    phi_plot = phi;
                    za_fail_pl = za_fail;
                else
                    set(gcf,'Position',[ 440   356   560*2/3   420])
                    phi_plot = phi(:,:,1:2);
                    za_fail_pl = za_fail(:,:,1:2);
                end
                imagesc(squeeze(phi_plot(:,j,:))*24);colormap(winter);colorbar; set(gca,'YDir','normal','Visible','on','FontSize',20,'XTick',[]); 
                h = findall(gca,'Type','image');
                set (h(1), 'AlphaData', squeeze(~za_fail_pl(:,j,:)));                 % Set transparency to account for error.

                if size(phi_plot,1) > 10
                    set(gca,'YTick',2:3:size(phi_plot,1))
                    iindex=0;
                    for itemp = 2:3:size(phi_plot,1)
                        iindex=iindex+1;
                        ylabelcells{iindex} = num2str(freq_listing(itemp,3),'%6.3g');
                    end

                    set(gca,'YTickLabel',ylabelcells);
                    clear itemp ylabelcells iindex
                end
            xlabel('Healthy / Latent / Seizing');
            ylabel('Freq Band');
        end
    else

        for j = rat2plot
            figure; set(gcf,'Color','w'); 
            
            Amp_pl = Amp;
            za_fail_pl = za_fail;
            
            if normalize_amps
                clims = [0 2];
            else
                clims = [min(min((Amp_pl(:,j,:)))) max(max((Amp_pl(:,j,:))))];
            end
            
            hottemp = hot;
            imagesc(squeeze(Amp_pl(:,j,:)), clims);
            %colormap(hottemp(round(end*1/8):round(end*7/8),:));
            colormap(hottemp(1:round(end*4/5),:));
            %colormap(gray)
            colorbar; set(gca,'YDir','normal','Visible','on','FontSize',20,'XTick',[]); 
%                     if normalize_amps; figure; set(gcf,'Color','w'); imagesc(squeeze(Amp_pl(:,j,:)), clims);colormap(hot);colorbar; set(gca,'YDir','normal'); else
%                     figure; set(gcf,'Color','w'); imagesc(squeeze(Amp_pl(:,j,:)));colormap(hot);colorbar; set(gca,'YDir','normal');  end
            h = findall(gca,'Type','image');
            
            if normalize_amps
                set (h(1), 'AlphaData', squeeze(~isnan(Amp_pl(:,j,:))));                 % Set transparency to account for error.
            else
                set (h(1), 'AlphaData', squeeze(~za_fail_pl(:,j,:)));                 % Set transparency to account for error.
            end
            %if normalize_amps; hold on; image(zeros(size(Amp,1),1,3)); end


            if size(Amp_pl,1) > 10
                set(gca,'YTick',2:3:size(Amp_pl,1))
                iindex=0;
                for itemp = 2:3:size(Amp_pl,1)
                    iindex=iindex+1;
                    ylabelcells{iindex} = num2str(freq_listing(itemp,3),'%6.3g');
                end

                set(gca,'YTickLabel',ylabelcells);
                clear itemp clear ylabelcells
            end
            xlabel('Pre / Latent / Chronic');
            ylabel('Freq Band');
        end
        
    end

end
