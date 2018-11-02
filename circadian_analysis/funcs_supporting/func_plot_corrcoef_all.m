
function func_plot_corrcoef_all(rcell,prepostchronic, smoothmode_ts, reload_spks)
    % This analyzes correlation coefficient between the following:
        % Low frequency band 6-hour moving average
        % High freuency band 6-hour moving average
        % LF 2-day
        % HF 2-day
        % SPW rhythm
    
    plot_all_timeseries = 0;
    plot_debug1 = 0;
    plot_debug2 = 1;
    
    chosen_freqbands = [2 5];
    do_average = 1;
        average_freqbands_hf = [5 6 7 8];
    [Tb Xb] = rat_2_ts (rcell,prepostchronic, 2);
    [T X] = rat_2_ts (rcell,prepostchronic, smoothmode_ts);
    
    % Load spike data
    if reload_spks
        [phi_means_rat4 data_struct_rat4] = ratscript_SPKs (4, 2, 2, 3);
        [phi_means_rat9 data_struct_rat9] = ratscript_SPKs (9, 2, 1, 3);
        [phi_means_rat10 data_struct_rat10] = ratscript_SPKs (10, 2, 1, 3);
        [phi_means_rat1 data_struct_rat1] = ratscript_SPKs (1, 2, 2, 3);
        save('ratscript_SPKs.mat','data_struct_rat4','data_struct_rat9','data_struct_rat10','data_struct_rat1');
    else
        load('ratscript_SPKs.mat');
    end

    i=1;r={data_struct_rat4}; [TSPK{i} XSPK{i}] = extract_ts(r, prepostchronic,smoothmode_ts); 
    i=2;r={data_struct_rat9}; [TSPK{i} XSPK{i}] = extract_ts(r, prepostchronic,smoothmode_ts); 
    i=3;r={data_struct_rat10}; [TSPK{i} XSPK{i}] = extract_ts(r, prepostchronic,smoothmode_ts); 
%     i=4;r={data_struct_rat1}; [TSPK{i} XSPK{i}] = extract_ts(r, prepostchronic,smoothmode_ts); 


%     Convert rhythm data to same time window as spikes. Crop to remove NaNs
    Nrats = length(TSPK);
    Nbands = size(T{1},2);
    for i = 1:Nrats
        i
        XR{i} = interp1(T{i}(:,1),X{i}(:,:),TSPK{i});
        Xbi{i} = interp1(Tb{i}(:,1),Xb{i}(:,:),TSPK{i});
        
        index1 = ~isnan(XR{i}(:,1));
        index2 = ~isnan(Xbi{i}(:,1));
        index = index1 & index2;
        
        XR{i} = XR{i}(index,:);
        Xbi{i} = Xbi{i}(index,:);
        TR{i} = repmat(TSPK{i}(index),1,Nbands);
        Tbi{i} = TR{i};

        XS{i} = XSPK{i}(index);
        TS{i} = TSPK{i}(index);
    end
    
    if plot_debug1
        for i = 1:3
            TSPKmat{i} = repmat(TSPK{i},1,8);
            XSPKmat{i} = repmat(XSPK{i},1,8);
            
            TSmat{i} = repmat(TS{i},1,8);
            XSmat{i} = repmat(XS{i},1,8);

            os.shift = 3;
            figure('Color','w','Position',[80 2 788   732]);
            hold on; plot_matrix2(T{i},X{i},os,'-','MarkerSize',6); legend('1','2','3','4','5','6','7','8')
            hold on; plot_matrix2(TSPKmat{i},XSPKmat{i},os,'k','MarkerSize',6); legend('1','2','3','4','5','6','7','8')
            hold on; plot_matrix2(TR{i},XR{i},os,'kx','MarkerSize',6); legend('1','2','3','4','5','6','7','8')
            title('Rhythm interpolations');
            
            os.shift = 3;
            figure('Color','w','Position',[80 2 788   732]);
            hold on; plot_matrix2(T{i},X{i},os,'-','MarkerSize',6); legend('1','2','3','4','5','6','7','8')
            hold on; plot_matrix2(TSPKmat{i},XSPKmat{i},os,'k','MarkerSize',6); legend('1','2','3','4','5','6','7','8')
            hold on; plot_matrix2(TSmat{i},XSmat{i},os,'k.','MarkerSize',6); legend('1','2','3','4','5','6','7','8')
            title('Spike interpolations');
            
        end
    end

    if plot_all_timeseries
        % Test plotting
        for i = 1:3
            TSPKmat{i} = repmat(TSPK{i},1,8);
            XSPKmat{i} = repmat(XSPK{i},1,8);

            os.shift = 3;
            figure('Color','w','Position',[80 2 788   732]);
            hold on; plot_matrix2(T{i},X{i},os,'.-','MarkerSize',6); legend('1','2','3','4','5','6','7','8')
            hold on; plot_matrix2(TSPKmat{i},XSPKmat{i},os,'k.','MarkerSize',6); legend('1','2','3','4','5','6','7','8')


        end
    end

    
    % Merge Data into 1 matrix (not including Rat 4)
    Nbands = length(chosen_freqbands);
    for i = 1:Nrats
        if ~do_average
            Xall{i} = [XR{i}(:,chosen_freqbands) Xbi{i}(:,chosen_freqbands) XS{i}];
            Tall{i} = repmat(TS{i},1,size(Xall{i},2));
        else
            Xall{i} = [XR{i}(:,chosen_freqbands(1)) mean(XR{i}(:,average_freqbands_hf),2) Xbi{i}(:,chosen_freqbands(1)) mean(Xbi{i}(:,average_freqbands_hf),2) XS{i}];
            Tall{i} = repmat(TS{i},1,size(Xall{i},2));
        end
        
    end
    
    
    if plot_debug2
        os.shift = 3;
        for i = 1:Nrats
            figure('Color','w','Position',[80 2 788   732]);
            hold on; plot_matrix2(Tall{i},Xall{i},os,'.-','MarkerSize',6);
        end
    end
    
    phi_stage = zeros(size(Xall{1},2),length(Xall));
    [corrcoef_3D hi] = corrcoef_EEG_bands(Tall,Xall,phi_stage,0,1,1);
    add_axis_labels(hi);
    
end


function add_axis_labels(h)
    for i = 1:length(h)
        set(h(i),'XTick',1:5);
        set(h(i),'XTickLabel',{'A LF','A HF','B LF','B HF','SPW'})
        set(h(i),'YTick',1:5);
        set(h(i),'YTickLabel',{'A LF','A HF','B LF','B HF','SPW'})
    end
end
