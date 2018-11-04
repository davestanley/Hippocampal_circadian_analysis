function func_plot_PCA(rcell,prepostchronic, smoothmode_ts)


        freq_bands_sort= 0;  % Separate frequency bands based on dominant mode
        plot_ICA=0;
        plot_biplot=0;
        plot_first_mode=0;
        
        [T X] = rat_2_ts (rcell,prepostchronic, smoothmode_ts);
        
        
        
        
        % End of get corrcof from princcomp
        
        % Calculate percent power contribution of each eigenvalue
        figure('Color','w');
        hold on; latent_arr = [];
        Nrats = length(X);
        for i=1:Nrats; [coef, score, latent] = ordered_princomp(T{i},X{i});
            plot(cumsum(latent)/sum(latent),'ko','MarkerSize',2,'LineWidth',3); latent_arr = [latent_arr latent];
        end
        xlabel('Principle Components','FontSize',16);
        ylabel('Power Fraction','FontSize',16);
        
        ylim([0 1]); if length(latent) > 10; xlim([0 8]); end
        hold on; plot([0 8], [0.95 0.95],'k--');
        hold on; plot(cumsum(mean(latent_arr'))/sum(mean(latent_arr')),'rx','MarkerSize',20,'LineWidth',2);
        set(gca,'FontSize',15)

 
        % Plot eigenvectors for the first 3 principle components
        figure('Color','w','Position',[20 20  560   420*1.5]); 
        colourarr = 'brgm';
        coefarr = [];
        for i = 1:Nrats
            [coef, score, latent] = ordered_princomp(T{i},X{i}); for j = 1:3; subplot(3,1,j); hold on; plot(coef(:,j),'k');end % Plot first 3 eigenvectors
            coefarr = cat(3,coefarr,coef);
        end
        subplot(3,1,1); title('Eigenvectors of first three principle components');
        
        % Overlay mean of first 3 eigenvectors
        for j = 1:3; subplot(3,1,j); hold on; plot(squeeze(mean(coefarr(:,j,:),3)),'r','LineWidth',4); end 
        
        
        if size(coef,1) > 10
            for j = 1:3; subplot(3,1,j); set(gca,'XTick',2:3:size(coef,1),'XTickLabel',{},'Box','off','FontSize',16); end
            iindex=0;
            for itemp = 2:3:size(coef,1)
                iindex=iindex+1;
                ylabelcells{iindex} = num2str(freq_listing(itemp,3),'%0.0f');
            end

            set(gca,'XTickLabel',ylabelcells,'Box','off');
            clear itemp ylabelcells iindex
        else
            for j = 1:2; subplot(3,1,j); set(gca,'XTickLabel',{},'Box','off','FontSize',16); end
            j=3; subplot(3,1,j); set(gca,'Box','off','FontSize',16); 
        end
        xlabel('Freq Band','FontSize',16);
        
        % Plot phase-wrapped score for rat 4
        %figure('Color','w','Position',[20 20 300 800]); % Plot first 3 transformed modes
        figure('Color','w','Position',[ 445   378   460   406]); % Plot first 3 transformed modes
        colourarr = 'brgy';
        ratrange=  1:1;
        for i = ratrange
            subplot(length(ratrange),1,i);
            [coef, score, latent] = ordered_princomp(T{i},X{i});
            score_norm = score ./ repmat(var(score(:,1)),size(score,1),size(score,2));
            for j = 1:3
                if smoothmode_ts ~= 0  % Take mod only for baseline case
                    hold on; plot((T{i}(:,1)),score_norm(:,j),[colourarr(j) '.']);
                else
                    hold on; plot(mod(T{i}(:,1),1)*24,score_norm(:,j),[colourarr(j) 'o'],'MarkerSize',2,'LineWidth',3);
                end
                
            end
            set(gca,'FontSize',14);
            ylim([-1 1]);
            zero_amp_arr =[];
            w = 2*pi/1.0;
            alpha = .05;
            for j = 1:3; dat = cosinor_struct(T{i}(:,1),score_norm(:,j),w,alpha,0); zero_amp_arr = [zero_amp_arr dat.p_3a]; end
            zero_amp_arr 
        end
        xlabel('Time (hours)','FontSize',14);
        ylabel('Score','FontSize',14);
        legend('PC1','PC2','PC3')
%         
%         figure('Color','w'); % Plot first 3 transformed modes
%         colourarr = 'brgy';
%         ratrange=  1;
%         for i = ratrange
%             subplot(length(ratrange),1,i);
%             [coef, score, latent] = ordered_princomp(T{i},X{i});
%             score_norm = score ./ repmat(var(score(:,1)),size(score,1),size(score,2));
%             hold on; plot3(score_norm(:,1),score_norm(:,2),score_norm(:,3),'.')
%             set(gca,'FontSize',30);
%         end
%         
        
        
        
%         i=1; [coef, score, latent] = ordered_princomp(T{i},X{i});
%         coefinv= inv(coef);
%         %figure; plot(var(score(:,1)*coefinv(1,:)) ./ var(X{i}))
%         figure; plot(coef(:,1:3).^2 .* repmat(latent(1:3)',8,1))

        %Plot fraction of variance in each of the frequency bands by keeping the first X modes
        figure('Color','w'); % Plot first 3 transformed modes
        colourarr = 'bgr';
        ratrange= 1:1;
        for i = ratrange
            [coef, score, latent] = ordered_princomp(T{i},X{i});
            subplot(length(ratrange),1,i)
            for j = 1:3
                range = 1:j;
                coefset = coef(:,range).^2
                latentset = repmat(latent(range)',size(coef,1),1)
                hold on; plot( sum(coefset .* latentset,2),colourarr(j)); ylim([0 1]);
            end
            if i==1; legend('PC1','PC2','PC3','Location','southeast'); end
            if i==1; ylabel('Cumulative Power Contribution','FontSize',20); end
            
            %princomp_addticks(coef);
            
        end
        set(gca,'FontSize',16);
        xlabel('Freq Band','FontSize',20)
        
        % % Separate frequency bands based on dominant modess
        if freq_bands_sort
            colourarr = 'bgr';
            ratrange= 1:4;
            for i = ratrange
                modes_contrib_arr = [];
                [coef, score, latent] = ordered_princomp(T{i},X{i});
                figure;set(gcf,'Position',[  278         136        1120         415],'Color','w');
                subplot(121);
                for j = 1:3
                    modescontrib = coef(:,j).^2*latent(j);
                    modescontrib = coef(:,j);

                    hold on; plot( modescontrib,colourarr(j))
                    modes_contrib_arr = [modes_contrib_arr modescontrib];
                end
                princomp_addticks(coef);

                modes_ratio = (modes_contrib_arr(:,1) > modes_contrib_arr(:,2));
    %             figure
    %             subplot(2,1,1)
    %             hold on; plot(score(:,1),'k','LineWidth',2)
    %             hold on; plot(X{i}(:,modes_ratio),'b')
    %             
    %             
    %             subplot(2,1,2)
    %             hold on; plot(score(:,2),'k','LineWidth',2)
    %             hold on; plot(X{i}(:,~modes_ratio),'r')
    % % 
                  subplot(122);
                  shiftval = 2;
                  Xshift = X{i} + repmat((0:size(X{i},2)-1)*shiftval,size(X{i},1),1);
                  hold on; plot(T{i}(:,modes_ratio),Xshift(:,modes_ratio),'b')
                  hold on; plot(T{i}(:,~modes_ratio),Xshift(:,~modes_ratio),'r')


    %               figure; shiftval = 0;
    %               Xshift = X{i} + repmat((0:size(X{i},2)-1)*shiftval,size(X{i},1),1);
    %               hold on; plot(mod(T{i}(:,modes_ratio),1),Xshift(:,modes_ratio),'b.')
    %               hold on; plot(mod(T{i}(:,~modes_ratio),1),Xshift(:,~modes_ratio),'r.')

            end
        end
        
        
        % % ICA
        if plot_ICA
            figure;
            dim2plot=1:2;
            coefsetall = [];
            addpath ('./ICA/PCA_and_ICA')
            Nmodes = 2;
            ratrange=1:4;
            for i = ratrange
                [z_ic A T mean_z] = myICA(X{i}',Nmodes);
                z = X{i}';
                z_LD = T \ pinv(A) * z_ic + repmat(mean_z,1,size(z,2));
                z_ic=z_ic'; z_LD = z_LD';


                coefset = zeros(size(X{i},2),Nmodes);
                for k = 1:size(X{i},2)
                    for j = 1:Nmodes
                        coefset(k,j) = corr(z_ic(:,j),X{i}(:,k));
                    end
                end

                %figure; plot(z_ic)

                %hold on; plot(coefset.^2,'LineWidth',2);legend('m1','m2','m3');


                coefset2 = (T \ pinv(A) + repmat(mean_z,1,size(A,1))).^2;
                figure; hold on; plot(coefset2(:,dim2plot))

                coefsetall = cat(3,coefsetall,coefset2);

            end
        
        
            meancoefs = mean(coefsetall,3);
            hold on; plot(meancoefs(:,dim2plot),'Linewidth',2);
        
        end
        
        if plot_biplot
            for i = ratrange
                [coef, score, latent] = ordered_princomp(T{i},X{i});
                vbls = {'1','2','3','4','5','6','7','8'};
                figure; biplot(coef(:,1:2),'scores',score(:,1:2),...
                    'varlabels',vbls);
            end
        end
        
        if plot_first_mode
            ratrange=1;
            for i = ratrange
                [coef, score, latent] = ordered_princomp(T{i},X{i});
                coefinv=  inv(coef);
                mode1 = score(:,1)*coefinv(1,:);
                figure; plot(mode1);
                mode1b = score(:,1)*coef(:,1)';
                hold on; plot(mode1b)
                mode1_pow1 = latent(1).*coef(:,1).^2;
                mode1_pow1b = var(mode1b)';
                mode1_pow1 - mode1_pow1b;
            end
        end
        
        
end



function [coef, score, latent] = ordered_princomp(T,X)
    use_ordered_princomp = 1;
    w=2*pi/1;
    alpha=0.05;

    [coef, score, latent] = pca(X);
    % Note: score = X * coef (columns of coef are the loadings for each
    % principal component)
    
    if use_ordered_princomp
        for j = 1:size(coef,2)

    %         p=polyfit(1:size(coef,2),coef(:,j)',1);
    %         if p(1) < 0
    %             coef(:,j) = -coef(:,j);
    %             score(:,j) = -score(:,j);
    %         end

    %         if mean(coef(1:2,j)) < 0
    %             coef(:,j) = -coef(:,j);
    %             score(:,j) = -score(:,j);
    %         end

    %         s = cosinor_struct(T,score(:,j),w,alpha); phi = abs(s.phi)/2/pi;
    %         if (j==1 && (abs(phi-0.5)<0.25)) || (j==2 && (abs(phi-0.5)>=0.25))
    %             coef(:,j) = -coef(:,j);
    %             score(:,j) = -score(:,j);
    %         end
            % Not sure what this is doing, but I think it's shifting the
            % signs of the first and 2nd principle components based on the
            % estimated phase of the 1st component.
            s = cosinor_struct(T,score(:,j),w,alpha); phi = abs(s.phi)/2/pi;
            if (j==1 && (abs(phi-0.5)<0.25))
                coef(:,j) = -coef(:,j); coef(:,j+1) = -coef(:,j+1);
                score(:,j) = -score(:,j); coef(:,j+1) = -coef(:,j+1);       % Shouldn't this be the score?
            end
        end
    end
end


function princomp_addticks(coef)
    global freq_listing

    if size(coef,1) > 10
        set(gca,'XTick',2:3:size(coef,1),'XTickLabel',{},'Box','off','FontSize',20);
        iindex=0;
        for itemp = 2:3:size(coef,1)
            iindex=iindex+1;
            ylabelcells{iindex} = num2str(freq_listing(itemp,3),'%0.0f');
        end

        set(gca,'XTickLabel',ylabelcells,'Box','off');
        clear itemp ylabelcells iindex
    else
        set(gca,'XTickLabel',{},'Box','off','FontSize',20);
        set(gca,'Box','off','FontSize',20); 
    end
end
