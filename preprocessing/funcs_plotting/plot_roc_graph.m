
function plot_roc_graph(predict,actual,theta_ratio_ranges0)


    false_positives = sum((predict == 1) & (actual ~= 1));
    false_negatives = sum((predict == 2) & (actual == 1));
    true_positives = sum((predict == 1) & (actual == 1));
    true_negatives = sum((predict == 2) & (actual ~= 1));
    positives = sum((actual == 1));
    error_rate = (false_positives + false_negatives) ./ (500*ones(1,length(false_positives))) * 100;
    TPR = (true_positives) ./ (true_positives + false_negatives) * 100;
    FPR = (false_positives) ./ (false_positives + true_negatives) * 100;

    figure; plot(FPR,TPR,'k.','MarkerSize',20); set(gcf,'Color','w')
    hold on; plot([0 100], [0 100], 'k--','LineWidth',2)
    xlabel ('1-Specificity (false positive rate)','FontSize',26); ylabel('Sensitivity (true positive rate)','FontSize',26);
    for i = 1:length(FPR)
        hold on; text(FPR(i),TPR(i),[num2str(theta_ratio_ranges0(i))],'FontSize',20);
    end
    set(gca,'FontSize',20);
    

    
%     actual_old = actual;
%     actual(actual == 2) = 1;
%     predict = predict_all;
%     
%     
%     
%     false_positives = sum((predict == 1) & (actual ~= 1));
%     false_negatives = sum((predict == 2) & (actual == 1));
%     true_positives = sum((predict == 1) & (actual == 1));
%     true_negatives = sum((predict == 2) & (actual ~= 1));
%     positives = sum((actual == 1));
%     error_rate = (false_positives + false_negatives) ./ (500*ones(1,length(false_positives))) * 100;
%     TPR = (true_positives) ./ (true_positives + false_negatives) * 100;
%     FPR = (false_positives) ./ (false_positives + true_negatives) * 100;
% 
%     hold on; plot(FPR,TPR,'.r')
%     xlabel ('FPR or 1-specificity'); ylabel('TPR or sensitivity');
%     for i = 1:length(FPR)
%         hold on; text(FPR(i),TPR(i),[num2str(theta_ratio_ranges0(i))]);
%         title('Let all actual=2 be theta');
%     end

    
    

end

