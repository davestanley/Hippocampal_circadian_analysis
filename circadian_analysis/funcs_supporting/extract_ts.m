


function [T X] = extract_ts (r, prepostchronic, smoothmode_ts, zscore_on)

    if (~exist('prepostchronic','var')); prepostchronic=0; end
%                             0 = all stages
%                             1/2/3 = pre/post/chronic
    if (~exist('smoothmode_ts','var')); smoothmode_ts=0; end
%                 smoothmode_ts - Determines filter to apply to ts data
%                     = -1 : 0 + estimate amplitude by taking variance
%                     = 0 : Pre smoothing + baseline subtraction
%                     = 1 : No pre-smoothing
%                     = 2 : Baseline signal: 1 day 90% aveage
%                     = 3 : Pre smoothing - 6-hr 90% moving average
%                     = 4 : 3 + Pull out amplitude based on cosine fit
%                     = 5 : 3 + Pull out phase based on cosine fit
%                     = 11 : From scratch, baseline subtract (should be same as 0)
%                     = 12 : From scratch, baseline subtract and normalize (percent change)
    
    if (~exist('zscore_on','var')); zscore_on=1; end
    plot_on = 0;
    plot_debug = 0;
    remove_bad_data_r1 = 1; % Remove bad data between 30.6 and 30.9. This corrupts only Rat #1 (4th in our system here)
                            % But there is no easy way to glean this
                            % information
    separate_stages_after = 1;      % Reduce entire dataset into pre/acute/chronic AFTER conducting
                                    % all smoothing, etc. This results in
                                    % slower processing times, but also
                                    % prevents some edges of the dataset
                                    % from getting cut off

    bounds = get_bounds(r{1});
    i=1;
        

    if separate_stages_after
        % This is just a placeholder. The separation occurs at the end of
        % the code
        index = logical(ones(1,length(r{i}.all.t)));
    else
        if prepostchronic == 0
            index = logical(zeros(1,length(r{i}.all.t)));
            for ii = 1:size(bounds,1)
                index = index | (r{i}.all.t >= bounds(ii,1) & r{i}.all.t <= bounds(ii,2));
            end
    %         index = (r{i}.all.t >= bounds(1,1) & r{i}.all.t <= bounds(1,2)) | ...
    %             (r{i}.all.t >= bounds(2,1) & r{i}.all.t <= bounds(2,2)) | ...
    %             (r{i}.all.t >= bounds(3,1) & r{i}.all.t <= bounds(3,2));
        elseif prepostchronic == 1
            index = (r{i}.all.t >= bounds(1,1) & r{i}.all.t <= bounds(1,2));
        elseif prepostchronic == 2
            index = (r{i}.all.t >= bounds(2,1) & r{i}.all.t <= bounds(2,2));
        elseif prepostchronic == 3
            index = (r{i}.all.t >= bounds(3,1) & r{i}.all.t <= bounds(3,2));
        elseif prepostchronic == 4
            index = (r{i}.all.t >= bounds(4,1) & r{i}.all.t <= bounds(4,2));        
        end
    end
    
    
    t = r{i}.all.t; d = r{i}.all.d;
    
    if plot_debug
        figure; plot(t,d,'k.');
        hold on; plot(t(index),d(index),'r.');
        s = r{i};
        hold on; plot(s.ctrl.t,s.ctrl.d,'c.');
        hold on; plot(s.acute.t,s.acute.d,'g.');
        hold on; plot(s.chr.t,s.chr.d,'y.');
    end
    
    
    X=[];
    T=[];
    
    
    if smoothmode_ts >= 1
        % Pull out time series
        for i = 1:length(r)
            t = r{i}.all.t(index); d = r{i}.all.d(index);
            %if i == 4
                if remove_bad_data_r1; index2 = (t<=30.6 | t>=30.9);
                t=t(index2); d=d(index2); end
            %end
            X = [X d(:)];
            T = [T t(:)];
        end
        
%         fhandle = @(t,y) dummy_STD(t,y);      % Get variance
%         [T X] = calc_ts_func(T, X,2,0.9,0.5,fhandle);
        
        % Basline smooth
        if smoothmode_ts == 2
            bin_size = 3;
            fract_overlap = 0.97;
            fract_maxgap = 0.5;
            [T X] = daveMVAVG_MAT (T(:,1), X, bin_size, fract_overlap, fract_maxgap);
            T = repmat(T,1,size(X,2));
        end
            
        % The following options (3 .. Inf) use the basic 6-hour prefilter
        if smoothmode_ts >= 3 && smoothmode_ts <= 10

            bin_size = 6/24;
            fract_overlap = 0.9;
            fract_maxgap = 0.9;
            [T X] = daveMVAVG_MAT (T(:,1), X, bin_size, fract_overlap, fract_maxgap);
            T = repmat(T,1,size(X,2));
                
            if smoothmode_ts == 4
                % The commands with get_phase_only don't work properly now
                %fhandle = @(t,y) get_phase_only(t,y,2*pi/1.0,0.2,1);      % Get sinusoid amplitude
                fhandle = @(t,y) dummy_STD(t,y);      % Get variance
                [T X] = calc_ts_func(T, X,2,0.9,0.5,fhandle);
            elseif smoothmode_ts == 5
                % The commands with get_phase_only don't work properly now
                fhandle = @(t,y) get_phase_only(t,y,2*pi/1.0,0.2,0);      % Get phase
                [T X] = calc_ts_func(T{i}, X{i},2,0.5,0.5,fhandle);
            end
        end
        
        if smoothmode_ts == 11 || smoothmode_ts == 12

            bin_size = 6/24;
            fract_overlap = 0.9;
            fract_maxgap = 0.9;
        
            [t_sm dat_sm] = daveMVAVG_MAT (T(:,1), X, bin_size, fract_overlap, fract_maxgap);
        
            bin_size = 2;
            fract_overlap = 0.9;
            fract_maxgap = 0.2;
            [t_base dat_base] = daveMVAVG_MAT (T(:,1), X, bin_size, fract_overlap, fract_maxgap);

            dat_basei = interp1(t_base,dat_base, t_sm);
            if smoothmode_ts == 11
                dat_sub = dat_sm - dat_basei;
            elseif smoothmode_ts == 12
                dat_sub = (dat_sm - dat_basei) ./ dat_basei;
            end
            
            if plot_on;
                figure;
                subplot(121); hold on;
                plot(T, X, 'b','LineWidth',2);
                plot(t_sm, dat_sm, 'g.','MarkerSize',18);
                plot(t_base, dat_base, 'k.','MarkerSize',18);
                subplot(122)
                
                dat_sub_temp=dat_sub(~isnan(dat_sub));
                dat_sub_temp = reshape(dat_sub_temp,sum(~isnan(dat_sub(:,1))),size(dat_sub,2));
                os.shift=2; plot_matrix2(t_sm(~isnan(dat_sub(:,1))),zscore(dat_sub_temp),os);
            end
            
            T = repmat(t_sm,1,size(X,2));
            X = dat_sub;
            
        end

    elseif smoothmode_ts <= 0 
        for i = 1:length(r)
            if prepostchronic == 0
                d = [r{i}.ctrl.d, r{i}.acute.d,r{i}.chr.d];
                t = [r{i}.ctrl.t, r{i}.acute.t,r{i}.chr.t];
                if remove_bad_data_r1; index2 = (t<=30.2 | t>=31);
                t=t(index2); d=d(index2); end
            elseif prepostchronic == 1
                d = [r{i}.ctrl.d];
                t = [r{i}.ctrl.t];
            elseif prepostchronic == 2
                d = [r{i}.acute.d];
                t = [r{i}.acute.t];
            elseif prepostchronic == 3
                d = [r{i}.chr.d];
                t = [r{i}.chr.t];
                if remove_bad_data_r1; index2 = (t<=30.2 | t>=31);
                t=t(index2); d=d(index2); end
            end
            X = [X d'];
            T = [T t'];
        end
        if smoothmode_ts == -1
            fhandle = @(t,y) dummy_STD(t,y);      % Get variance
            fract_overlap = 0.97;
            fract_maxgap = 0.2;
            [T X] = calc_ts_func(T, X,1,fract_overlap,fract_maxgap,fhandle);
        end
        if smoothmode_ts == -2
            % remove nans first
            index2 = ~isnan(X(:,1));
            X = X(index2,:);
            T = T(index2,:);
            Xi = hilbert(X);
            myangle = unwrap(angle(Xi));
            X = diff(myangle);
            
        end
        end
    
    index2 = ~isnan(X(:,1));
    X = X(index2,:);
    T = T(index2,:);
    
    if separate_stages_after
        if prepostchronic == 0
            index = logical(zeros(1,length(T(:,1))));
            index = index';
            for ii = 1:size(bounds,1)
                index = index | (T(:,1) >= bounds(ii,1) & T(:,1) <= bounds(ii,2));
            end
        elseif prepostchronic == 1
            index = (T(:,1) >= bounds(1,1) & T(:,1) <= bounds(1,2));
        elseif prepostchronic == 2
            index = (T(:,1) >= bounds(2,1) & T(:,1) <= bounds(2,2));
        elseif prepostchronic == 3
            index = (T(:,1) >= bounds(3,1) & T(:,1) <= bounds(3,2));
        elseif prepostchronic == 4
            index = (T(:,1) >= bounds(4,1) & T(:,1) <= bounds(4,2));        
        end
        X = X(index,:);
        T = T(index,:);
    end
    
%     bin_size=3;
%     fract_overlap=0.9;
%     fract_maxgap=0.9;
%     [T X] = daveMVAVG_MAT (T(:,1), X, bin_size, fract_overlap, fract_maxgap); % Take average time spent in theta state.
%     T = repmat(T,1,size(X,2));
    

    if zscore_on
        
        X = zscore(X);
        
%         mu = mean(X); sig = std(X);
%         X = X - repmat(mu,size(X,1),1);
%         X = X ./ repmat(sig,size(X,1),1);
%         
        
    end
%     shift = 5;
%     for i = 1:size(X,2); X(:,i) = X(:,i) + (shift*(i-1));end
    

    if plot_debug
        figure; plot(T,X,'.');
        %[stim seiz] = get_stimtime2(ratN);    
%         hold on; plot([stim stim],[min(min(X)) max(max(X))],'LineWidth',2)
%         hold on; plot([seiz seiz],[min(min(X)) max(max(X))],'LineWidth',2)
%         for i = min(min(T)):max(max(T)); hold on; plot([i i],[min(min(X)) max(max(X))],'k:','LineWidth',2); end
    end
end




function y = dummy_STD(t,varargin)
    y = std(varargin{:});
end
