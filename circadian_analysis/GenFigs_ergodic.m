

function GenFigs

    % Based on script previously called: JPHY_theta_epochs
    
    
    fig4_theta_power_vs_circadian_ergodic = 0;
    fig4_theta_power_vs_circadian_ergodic = 1;
    
    analyze_only_seizing = 1;
    band_number = 2;
    ergodic_mode = 0;
    theta_mode = 1;
    

    basepath = './Rats_theta_delta_epochs_power2/';
        circadian_filename = 'All4T1.mat';
    if analyze_only_seizing
        basepath_power = './t3_seizing_ratio2/';
    else
        %basepath_power = './t2_ratio_1_5_erg_nonseizing/';
        basepath_power = './t4_nonseizing_ratio2/';
    end

    
        %power_filename = 'All4T4_3pm_9pm.mat';
        %power_filename = 'All4_theta_delta_2_4_R01_longacute.mat';
        
%     basepath = '../../Data/RatFFT/03_theta_delta_2sec_bins/Rats_vary_frequrange_theta6_10_delta0_6/';
%         circadian_filename = 'All4T1_9pm_3am_ratio1p5_theta6_10_delta0_6.mat';    
        
    FS_axislabels_def=20;
    FS_axis_def=30;

    font_scaling=1;
    FS_axislabels=FS_axislabels_def*font_scaling;
    FS_axis=FS_axis_def*font_scaling;

    if fig4_theta_power_vs_circadian_ergodic
        
        
        if ~ergodic_mode
        theta_epochs{1} = get_theta_power ([basepath_power 'All4T1.mat'],theta_mode)
        theta_epochs{2} = get_theta_power ([basepath_power 'All4T2.mat'],theta_mode)
        theta_epochs{3} = get_theta_power ([basepath_power 'All4T3.mat'],theta_mode)
        theta_epochs{4} = get_theta_power ([basepath_power 'All4T4.mat'],theta_mode)
        theta_epochs{5} = get_theta_power ([basepath_power 'All4T5.mat'],theta_mode)
        theta_epochs{6} = get_theta_power ([basepath_power 'All4T6.mat'],theta_mode)
        theta_epochs{7} = get_theta_power ([basepath_power 'All4T7.mat'],theta_mode)
        theta_epochs{8} = get_theta_power ([basepath_power 'All4T8.mat'],theta_mode)
        else
        theta_epochs{1} = get_theta_power_ergodic ([basepath_power 'All4T1.mat'],theta_mode)
        theta_epochs{2} = get_theta_power_ergodic ([basepath_power 'All4T2.mat'],theta_mode)
        theta_epochs{3} = get_theta_power_ergodic ([basepath_power 'All4T3.mat'],theta_mode)
        theta_epochs{4} = get_theta_power_ergodic ([basepath_power 'All4T4.mat'],theta_mode)
        theta_epochs{5} = get_theta_power_ergodic ([basepath_power 'All4T5.mat'],theta_mode)
        theta_epochs{6} = get_theta_power_ergodic ([basepath_power 'All4T6.mat'],theta_mode)
        theta_epochs{7} = get_theta_power_ergodic ([basepath_power 'All4T7.mat'],theta_mode)
        theta_epochs{8} = get_theta_power_ergodic ([basepath_power 'All4T8.mat'],theta_mode)
        end
        
        
        if ~analyze_only_seizing
            for itemp = 1:length(theta_epochs)
                theta_epochs_temp{itemp} = theta_epochs{itemp}(1:2);
            end
            theta_epochs = theta_epochs_temp;
            clear itemp theta_epochs_temp
        end
        
        
        theta_all = repackage(theta_epochs);
        theta_merged = merge_rats(theta_all);
        
        analyze_mat = squeeze(theta_merged(band_number,:,:));
        analyze_mat = permute(analyze_mat,[2 1]);
        
        
        % Plot individual rat bargraph
        if ~ergodic_mode
            figure
            ratmat = cell2mat(permute(analyze_mat,[3 1 2]))
            ratmat = permute(ratmat,[2,3,1]);
            for i = 1:length(analyze_mat{1,1})
                subplot(length(analyze_mat{1,1}),1,i);
                %plot(ratmat(:,:,i),'LineWidth',2);
                bar(ratmat(:,:,i));
            end
        end
        
        % Plot mean bargraph stuff
        
        powermean = meanCell(analyze_mat); powerste = stdCell(analyze_mat) ./ countCell(analyze_mat); powerstd = stdCell(analyze_mat); powerconfid = confidCell(analyze_mat,0.05);
        figure('color','w','Position', [90   272   934   382]); bar(powermean); colormap(gray);
        legend('Pre','Post-L','Post-SS')
        set(gca,'FontSize',30,'Box','off')
        set(gca,'XTickLabel',[0 3 6 9 12 15 18 21])
        
        offset = 0.22;
        hold on; errorbar([(1:size(powermean,1))-offset; (1:size(powermean,1))+0; (1:size(powermean,1))+offset]',powermean, powerconfid ,'Color',[0.0 0.0 0.0],'Marker','.','LineStyle','none','LineWidth',1)
%         for i = 1:size(analyze_mat,1)
%             for j = 1:size(analyze_mat,2)
%                 hold on; plot( repmat(i+(j-2)*offset,length(analyze_mat{i,j}),1) ,analyze_mat{i,j}','.')
%             end
%         end


        % Phase shift in each frequency band
        N = size(analyze_mat,1);
        harr = [zeros(1,N)];
        hrank = [zeros(1,N)];
        parr = [zeros(1,N)];
        prank = [zeros(1,N)];
        for i = 1:N
            if ergodic_mode; [harr(i) parr(i)] = ttest2(analyze_mat{i,1},analyze_mat{i,2});
            else [harr(i) parr(i)] = ttest(analyze_mat{i,1},analyze_mat{i,2}); end
            [prank(i) hrank(i)] = ranksum(analyze_mat{i,1},analyze_mat{i,2});
            if hrank(i)
                text(i-0.07+0,(powermean(i,3)+powerconfid(i,3))*1.1,'*','FontSize',30)
            end
        end
        harr
        hrank
        parr
        prank
        clear N
        
        % Phase shift in each frequency band
        N = size(analyze_mat,1);
        harr = [zeros(1,N)];
        hrank = [zeros(1,N)];
        parr = [zeros(1,N)];
        prank = [zeros(1,N)];
        for i = 1:N
            if ergodic_mode; [harr(i) parr(i)] = ttest2(analyze_mat{i,1},analyze_mat{i,3});
            else [harr(i) parr(i)] = ttest(analyze_mat{i,1},analyze_mat{i,3}); end
            [prank(i) hrank(i)] = signrank(analyze_mat{i,1},analyze_mat{i,3});
            if hrank(i)
                text(i-0.07+offset,(powermean(i,3)+powerconfid(i,3))*1.1,'*','FontSize',30)
            end
        end
        harr
        hrank
        parr
        prank
        clear N
        
        

    end
end


function out = get_theta_power_ergodic (power_filename, theta_mode)
    if ~exist('theta_mode'); theta_mode = 0; end
    load ([power_filename]);
    
    if theta_mode == 1
        out{1} = r04_theta.power_ergodic.all;
        out{2} = r09_theta.power_ergodic.all;
        out{3} = r10_theta.power_ergodic.all;
        out{4} = r01_theta.power_ergodic.all;
    elseif theta_mode == 2
        out{1} = r04_theta.power_ergodic.intheta;
        out{2} = r09_theta.power_ergodic.intheta;
        out{3} = r10_theta.power_ergodic.intheta;
        out{4} = r01_theta.power_ergodic.intheta;
    elseif theta_mode == 3
        out{1} = r04_theta.power_ergodic.indelta;
        out{2} = r09_theta.power_ergodic.indelta;
        out{3} = r10_theta.power_ergodic.indelta;
        out{4} = r01_theta.power_ergodic.indelta;
    else
        out = [];
    end
    
end

function out = get_theta_power (power_filename, theta_mode)
    if ~exist('theta_mode'); theta_mode = 0; end
    load ([power_filename]);
    
    if theta_mode == 1
        out{1} = r04_theta.power_allbands.all;
        out{2} = r09_theta.power_allbands.all;
        out{3} = r10_theta.power_allbands.all;
        out{4} = r01_theta.power_allbands.all;
    elseif theta_mode == 2
        out{1} = r04_theta.power_allbands.intheta;
        out{2} = r09_theta.power_allbands.intheta;
        out{3} = r10_theta.power_allbands.intheta;
        out{4} = r01_theta.power_allbands.intheta;
    elseif theta_mode == 3
        out{1} = r04_theta.power_allbands.indelta;
        out{2} = r09_theta.power_allbands.indelta;
        out{3} = r10_theta.power_allbands.indelta;
        out{4} = r01_theta.power_allbands.indelta;
    else
        out = [];
    end
    
end

function theta_all = repackage(theta_epochs)
% 
    for i = 1:size(theta_epochs{1}{1}{1},2)  % Band
        for j = 1:length(theta_epochs{1})          % Rat number
            for k = 1:3      % Pre post etc
                for l = 1:length(theta_epochs)  % Timebin

                    if ~isempty(theta_epochs{l}{j}{k})
                        theta_all{i,j,k,l} = theta_epochs{l}{j}{k}(:,i);
                    else
                        theta_all{i,j,k,l} = [];
                    end
                end
            end
        end
    end

end

function theta_merged = merge_rats(theta_all)
    
    for i = 1:size(theta_all,1)
        for k = 1:size(theta_all,3)
            for l = 1:size(theta_all,4)
                theta_merged{i,k,l} = [];
                for j = 1:size(theta_all,2)
                    theta_merged{i,k,l} = [theta_merged{i,k,l}; theta_all{i,j,k,l}(:)];
                end
            end
        end
    end
        
end

function plot_theta_delta_circ (t,d,A,phi)
    
    t = t*24; t = mod(t,24);
    t = t(find(~isnan(d))); d = d(find(~isnan(d)));
    [t2 d2 d2_std] = daveMVAVG_bin (t, d, 1, 0.1, 0.9,0); % Bin data
    figure;
    hold on; plot(t,d,'k.','MarkerSize',15)
%     hold on; errorbar(t2,d2,d2_std,'k.','LineWidth',2)
    t3 = 0:0.1:24;
    hold on; plot(t3,A*cos(2*pi/24 * (t3-phi*24)),'k','LineWidth',2)
end

function xout = norm_mat_row (x)  % Normalized matrix by values in top-most row
    x_toprow = repmat(x(1,:),size(x,1),1);
    xout = x ./ x_toprow;
end


function y = meanCell(x)
    
    N = size(x,1);
    M = size(x,2);
    
    for i = 1:N
        for j = 1:M
            y(i,j) = mean(x{i,j});
        end
    end

end



function y = medianCell(x)
    
    N = size(x,1);
    M = size(x,2);
    
    for i = 1:N
        for j = 1:M
            y(i,j) = median(x{i,j});
        end
    end

end


function y = stdCell(x)
    
    N = size(x,1);
    M = size(x,2);
    
    for i = 1:N
        for j = 1:M
            y(i,j) = std(x{i,j});
        end
    end

end


function y = confidCell(x,alpha)
    
    N = size(x,1);
    M = size(x,2);
    
    for i = 1:N
        for j = 1:M
            y(i,j) = confidence(x{i,j},alpha);
        end
    end

end

function y = countCell(x)
    
    N = size(x,1);
    M = size(x,2);
    
    for i = 1:N
        for j = 1:M
            y(i,j) = length(x{i,j});
        end
    end

end