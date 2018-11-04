

function func_plot_EEG_vs_spks2(r_seiz, prepostchronic,smoothmode_ts,reload_spks);
        %prepostchronic = 1;     % 0 for all; 1 for pre; 2 for latency; 3 for chronic
        plottest_timeseries = 0;
        plot_all_timeseries = 0;

        i=1;[T{i} X{i}] = extract_ts(r_seiz{i}, prepostchronic,smoothmode_ts); 
        i=2;[T{i} X{i}] = extract_ts(r_seiz{i}, prepostchronic,smoothmode_ts); 
        i=3;[T{i} X{i}] = extract_ts(r_seiz{i}, prepostchronic,smoothmode_ts); 
        i=4;[T{i} X{i}] = extract_ts(r_seiz{i}, prepostchronic,smoothmode_ts); 

%         i=5;r=r5; [T{i} X{i}] = extract_ts(r, prepostchronic); 
%         i=6;r=r8; [T{i} X{i}] = extract_ts(r, prepostchronic); 
%         i=7;r=r11; [T{i} X{i}] = extract_ts(r, prepostchronic); 

        if reload_spks
            [phi_means_rat4 data_struct_rat4] = ratscript_SPKs (4, 2, 2, 3);
            [phi_means_rat9 data_struct_rat9] = ratscript_SPKs (9, 2, 1, 3);
            [phi_means_rat10 data_struct_rat10] = ratscript_SPKs (10, 2, 1, 3);
            [phi_means_rat1 data_struct_rat1] = ratscript_SPKs (1, 2, 2, 3);
            save('ratscript_SPKs.mat','data_struct_rat4','data_struct_rat9','data_struct_rat10','data_struct_rat1');
        else
            load('ratscript_SPKs.mat');
        end

%         i=1;r=data_struct_rat4; [TSPK{i} XSPK{i}] = extract_ts_SPK(r, prepostchronic); 
%         i=2;r=data_struct_rat9; [TSPK{i} XSPK{i}] = extract_ts_SPK(r, prepostchronic); 
%         i=3;r=data_struct_rat10; [TSPK{i} XSPK{i}] = extract_ts_SPK(r, prepostchronic); 
% %         i=4;r=data_struct_rat1; [TSPK{i} XSPK{i}] = extract_ts_SPK(r, prepostchronic); 

        i=1;r={data_struct_rat4}; [TSPK{i} XSPK{i}] = extract_ts(r, prepostchronic,smoothmode_ts); 
        i=2;r={data_struct_rat9}; [TSPK{i} XSPK{i}] = extract_ts(r, prepostchronic,smoothmode_ts); 
        i=3;r={data_struct_rat10}; [TSPK{i} XSPK{i}] = extract_ts(r, prepostchronic,smoothmode_ts); 
%         i=4;r={data_struct_rat1}; [TSPK{i} XSPK{i}] = extract_ts(r, prepostchronic,smoothmode_ts); 

    if plot_all_timeseries
        % Test plotting
        for i = 1:3
            TSPKmat{i} = repmat(TSPK{i},1,8);
            XSPKmat{i} = repmat(XSPK{i},1,8);

            os.shift = 3;
            figure('Color','w','Position',[80 2 788   732]); hold on; plot_matrix2(T{i},X{i},os,'.-','MarkerSize',6); legend('1','2','3','4','5','6','7','8')
            hold on; plot_matrix2(TSPKmat{i},XSPKmat{i},os,'k.-','MarkerSize',6); legend('1','2','3','4','5','6','7','8')

        end

    end

 
 
%         % Princomp with both EEG and SPK rates
%         i=1;
%         XR{i} = interp1(T{i}(:,1),X{i},TSPK{i});
%         index = ~isnan(XR{i}(:,1));
%         for j = 1:size(XR{i},2); XR2(:,j) = XR{i}(index,j); end
%         TR{i} = TSPK{i}(index);
%         XS{i} = XSPK{i}(index);
%         i=1; X_col = [XR2 XS{i}]; t_princ = TR{i};
%         clear XR TR TS index
%         [coef, score, latent] = princomp(X_col);
%         figure; plot(cumsum(latent)/sum(latent))
%         figure;for i = 1:8; subplot(8,1,i);plot(coef(:,i)); end
%         figure; plot(mod(t_princ,1),score(:,1),'.');hold on; plot(mod(t_princ,1),score(:,2),'r.');plot(mod(t_princ,1),score(:,3),'g.')
        

        N_freqbands=size(X{1},2);
        correlations = zeros(N_freqbands,length(TSPK));
        for chosen_band = 1:N_freqbands
            for i = 1:length(TSPK)
                if i == 1; rrat = 4; elseif i == 2; rrat = 9; elseif i == 3; rrat = 10; elseif i == 4; rrat = 1; end

                XR{i} = interp1(T{i}(:,1),X{i}(:,chosen_band),TSPK{i});
                index = ~isnan(XR{i});
                XR{i} = XR{i}(index);
                TR{i} = TSPK{i}(index);

                XS{i} = XSPK{i}(index);


%                 XS{i} = interp1(TSPK{i},XSPK{i},T{i}(:,1));
%                 index = ~isnan(XS{i});
%                 XS{i} = XS{i}(index); 
%                 XR{i} = X{i}(index,chosen_band);
%                 TR{i} = T{i}(index,chosen_band);

                if plottest_timeseries
                    figure;
                    plot(TR{i},XS{i});hold on; 
                    plot(T{i}(:,chosen_band),X{i}(:,chosen_band),'r')
                    plot(TR{i}(:,chosen_band),XR{i}(:,chosen_band),'r.')
                    add_stimseiz(rrat,'k','LineWidth',2')
                end

                %xcorr(XS{i},XR{i},0,'coeff')
                use_EEG_amp = 1;
                fhandle_EEG = @(t,y) get_phase_only(t,y,2*pi/1.0,0.2,use_EEG_amp);      % Function handel for cross correlation with zero lag
                fhandle_SPK = @(t,y) get_phase_only(t,y,2*pi/1.0,0.2,1);      % Always use PHASE for SPKs
%                 fhandle_SPK = @(t,y) dummy_STD(t,y);        % To just look at standard deviation
%                 fhandle_EEG = @(t,y) dummy_STD(t,y);        % To just look at standard deviation
                
                [TSphase{i} XSphase{i}] = calc_ts_func(TR{i}, XS{i},2,0.5,0.5,fhandle_SPK);
                if plottest_timeseries; hold on; plot(TSphase{i}, XSphase{i},'b.'); end
                [TRphase{i} XRphase{i}] = calc_ts_func(TR{i}, XR{i},2,0.5,0.5,fhandle_EEG);
                if plottest_timeseries; hold on; plot(TRphase{i}, XRphase{i},'r.'); end

                index = ~isnan(XSphase{i}) & ~isnan(XRphase{i}); 
                XSphase{i} = XSphase{i}(index); XRphase{i} = XRphase{i}(index); TSphase{i} = TSphase{i}(index); TRphase{i} = TRphase{i}(index);

                %correlations(chosen_band,i) = xcorr(XSphase{i}-mean(XSphase{i}),XRphase{i}-mean(XRphase{i}),0,'coeff');
                correlations(chosen_band,i) = abs(xcorr(zscore(XS{i}),zscore(XR{i}),0,'coeff'));

            end
        end
        figure; subplot(211); bar(correlations');
        subplot(212); bar(mean(correlations'))
end


function [T X] = extract_ts_SPK (r, prepostchronic)

    if (~exist('prepostchronic','var')); prepostchronic=0; end

    plot_debug = 0;
    detrend_ts = 1;

    bounds = get_bounds(r);
    i=1;
    if prepostchronic == 0
        index = (r.all.t >= bounds(1,1) & r.all.t <= bounds(1,2)) | ...
            (r.all.t >= bounds(2,1) & r.all.t <= bounds(2,2)) | ...
            (r.all.t >= bounds(3,1) & r.all.t <= bounds(3,2));
    elseif prepostchronic == 1
        index = (r.all.t >= bounds(1,1) & r.all.t <= bounds(1,2));
    elseif prepostchronic == 2
        index = (r.all.t >= bounds(2,1) & r.all.t <= bounds(2,2));
    elseif prepostchronic == 3
        index = (r.all.t >= bounds(3,1) & r.all.t <= bounds(3,2));
    end
    
    t = r.all.t; d = r.all.d;
    
    if plot_debug
        figure; plot(t,d,'k.');
        hold on; plot(t(index),d(index),'r.');
        s = r;
        hold on; plot(s.ctrl.t,s.ctrl.d,'c.');
        hold on; plot(s.acute.t,s.acute.d,'g.');
        hold on; plot(s.chr.t,s.chr.d,'y.');
    end
    

    if ~detrend_ts
        t = r.all.t(index); d = r.all.d(index);
    else
        
        if prepostchronic == 0
            d = [r.ctrl.d, r.acute.d,r.chr.d];
            t = [r.ctrl.t, r.acute.t,r.chr.t];
        elseif prepostchronic == 1
            d = [r.ctrl.d];
            t = [r.ctrl.t];
        elseif prepostchronic == 2
            d = [r.acute.d];
            t = [r.acute.t];
        elseif prepostchronic == 3
            d = [r.chr.d];
            t = [r.chr.t];
        end
    end
    X = [d(:)];
    T = [t(:)];

    
    index2 = ~isnan(X(:,1));
    X = X(index2,:);
    T = T(index2,:);
    mu = mean(X); sig = std(X);
    X = X - repmat(mu,size(X,1),1);
    X = X ./ repmat(sig,size(X,1),1);

    
    if plot_debug
        figure; plot(T,X,'.');
    end
end

