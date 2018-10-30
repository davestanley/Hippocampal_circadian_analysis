
function func_plot_princomp_amp_vs_ampraw(rcell,prepostchronic,smoothmode_ts)
    
    plot_individual_rats = 0;

    %prepostchronic = 1;
    ratrange = 1:length(rcell);

    [T X] = rat_2_ts (rcell,prepostchronic, 2);
    corrcoef_baseline = get_corrcoef(X,ratrange,plot_individual_rats);

    [T X] = rat_2_ts (rcell,prepostchronic, smoothmode_ts);
    corrcoef_default = get_corrcoef(X,ratrange,plot_individual_rats);


    figure;

    markerarr = 'x+od';
    colourarr = autumn;
    index = linspace(1,round(size(colourarr,1)),size(corrcoef_default,1));
    colourarr = colourarr(index(:),:);
    [freqs1 freqs2] = get_freq_meshgrid(corrcoef_default);
    % Plot basic to get legend
    for j=1:size(corrcoef_default,3)
        tria_def = pull_out_triangle(corrcoef_default(:,:,j));
        tria_baseline = pull_out_triangle(corrcoef_baseline(:,:,j));
        freqs1_tri = pull_out_triangle(freqs1(:,:,j));
        freqs2_tri = pull_out_triangle(freqs2(:,:,j));
        hold on; plot(tria_def,tria_baseline,[markerarr(j)],'Color',colourarr(1,:),'LineWidth',2,'MarkerSize',9);
        %hold on; plot(freqs1_tri-freqs2_tri,tria_def,[markerarr(j)],'Color',colourarr(1,:),'LineWidth',2,'MarkerSize',9);
    end
    
%     for j=1:size(corrcoef_default,3)
%         figure
%         tria_def = (corrcoef_default(:,:,j)); tria_def = tria_def(:);
%         tria_baseline = (corrcoef_baseline(:,:,j)); tria_baseline = tria_baseline(:);
%         freqs1_tri = (freqs1(:,:,j)); freqs1_tri = freqs1_tri(:);
%         freqs2_tri = (freqs2(:,:,j)); freqs2_tri = freqs2_tri(:);
%         hold on; hist(tria_def,10)
%         
%     end
    
    % Plot to get coloring
    for j=1:size(corrcoef_default,3)
        tria_def = pull_out_triangle(corrcoef_default(:,:,j));
        tria_baseline = pull_out_triangle(corrcoef_baseline(:,:,j));
        freqs1_tri = pull_out_triangle(freqs1(:,:,j));
        freqs2_tri = pull_out_triangle(freqs2(:,:,j));
        for i = 1:length(tria_def)
            hold on; plot(tria_def(i),tria_baseline(i),[markerarr(j)],'Color',colourarr(freqs1_tri(i),:),'LineWidth',2,'MarkerSize',9);
        end
    end

    [freqs1 freqs2] = get_freq_meshgrid(corrcoef_default);
    i=1:size(corrcoef_default,1);
    freqs1 = pull_out_triangle(freqs1,i);
    freqs2 = pull_out_triangle(freqs2,i);
    tria_def = pull_out_triangle(corrcoef_default,i);
    tria_baseline = pull_out_triangle(corrcoef_baseline,i);
    for i = 1:length(tria_def)
       percent_shift = 0.01; xshift = get(gca,'XLim'); yshift = get(gca,'YLim'); xshift = diff(xshift)*percent_shift; yshift = diff(yshift)*percent_shift;
       hold on; text(tria_def(i)+xshift,tria_baseline(i)+yshift,[num2str(freqs1(i)) ',' num2str(freqs2(i))],'Color','k','FontSize',12)
    end

    xlabel('Corrcoef Default - Defined by smoothmode ts'); ylabel('Corrcoef Baseline'); ylim([-1 1]); xlim([-1 1])
    legend('R1','R2','R3','R4')
    
end


function corrcoef_arr = get_corrcoef(X,ratrange,plot_individual_rats)
    %plot_individual_rats=0;
    
    %ratrange = 1:4;
    corrcoef_arr = [];
    for i = ratrange
        [coef, score, latent] = princomp(X{i});
        corrcoef_PC = (coef)*diag(latent)*inv(coef);
        if plot_individual_rats
            figure;
            colourmapval = [summer; flipud(autumn)]; colormap(colourmapval);
            imagesc(corrcoef_PC,[-1 1]); colorbar;set(gca,'YDir','normal')
        end
        corrcoef_arr = cat(3,corrcoef_arr, corrcoef_PC);
    end
    

end



function [freqs1 freqs2] = get_freq_meshgrid(corrcoef_default)
    [freqs1 freqs2] = meshgrid(1:size(corrcoef_default,1),1:size(corrcoef_default,2));
    freqs1 = repmat(freqs1,[1,1,size(corrcoef_default,3)]);
    freqs2 = repmat(freqs2,[1,1,size(corrcoef_default,3)]);
end



function output = pull_out_triangle(input,i)

    if ~exist('i','var')
        i = 1:size(input,2);
    end

    ut = logical(repmat(triu(ones(size(input,1),size(input,2)),1),[1 1 size(input,3)]));  % Generate indices to pull out upper triangle
    ut = ut(:,i,:);
    input = input(:,i,:);
    output = input(ut);
end