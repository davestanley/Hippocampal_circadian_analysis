
function func_plot_corr_phaseshift(r_seiz,r_ns,r_ctrl,analyze_only_seizing,smoothmode_ts)
    % These numbers correspond to the freq bands that undergo phase
    % shifting
    
    take_corr_abs = 1;
    
    % So far only works for seizing animals.
    r4ind = 3;
    r9ind = 3;
    r10ind = 4;
    r1ind = 3;
    
    r4ind = 3;
    r9ind = 4;
    r10ind = 4;
    r1ind = 4;
    
    take_abs = 0;
    
    yzpairs = [r4ind r9ind, r10ind, r1ind; 1 2 3 4];
%     yzpairs = [2 2 2 2; 1 2 3 4];
%     yzpairs = [7 7 7 7; 1 2 3 4];
%     yzpairs = [20 20 20 20; 1 2 3 4];

    corrcoef_ppc = [];
    
%     figure;

    prepostchronic_temp = 1;
    [T X] = rat_2_ts (r_seiz,prepostchronic_temp, smoothmode_ts);
    [corrcoef3D] = corrcoef_EEG_bands(T,X,prepostchronic_temp,0,0,0);
    if take_abs corrcoef3D=abs(corrcoef3D); end
    corrcoef2D = mat3Dto2Dyz (corrcoef3D,yzpairs); % Pull out section of correlation matrix corresponding to corr of phase shifting band vs other bands for each rat
    corrcoef_ppc = cat(3,corrcoef_ppc,corrcoef2D);
    
%     for j = 1:size(corrcoef2D,2)
%         subplot(4,1,j); hold on; plot(squeeze(corrcoef2D(:,j)),'k'); ylim([0 1])
%     end
%     
%     figure;
%     subplot(411);plot(corrcoef3D(3,:,1),'k');
%     subplot(412);plot(corrcoef3D(3,:,2),'k');
%     subplot(413);plot(corrcoef3D(4,:,3),'k');
%     subplot(414);plot(corrcoef3D(3,:,4),'k');
    
%     figure;
%     subplot(221);imagesc(corrcoef3D(:,:,1));
%     subplot(222);imagesc(corrcoef3D(:,:,2));
%     subplot(223);imagesc(corrcoef3D(:,:,3));
%     subplot(224);imagesc(corrcoef3D(:,:,4));
    
    
    
    
    
    prepostchronic_temp = 2;
    [T X] = rat_2_ts (r_seiz,prepostchronic_temp, smoothmode_ts);
    [corrcoef3D] = corrcoef_EEG_bands(T,X,prepostchronic_temp,0,0,0);
    if take_abs corrcoef3D=abs(corrcoef3D); end
    corrcoef2D = mat3Dto2Dyz (corrcoef3D,yzpairs); % Pull out section of correlation matrix corresponding to corr of phase shifting band vs other bands for each rat
    corrcoef_ppc = cat(3,corrcoef_ppc,corrcoef2D);
    
    prepostchronic_temp = 3;
    [T X] = rat_2_ts (r_seiz,prepostchronic_temp, smoothmode_ts);
    [corrcoef3D] = corrcoef_EEG_bands(T,X,prepostchronic_temp,0,0,0);
    if take_abs corrcoef3D=abs(corrcoef3D); end
    corrcoef2D = mat3Dto2Dyz (corrcoef3D,yzpairs); % Pull out section of correlation matrix corresponding to corr of phase shifting band vs other bands for each rat
    corrcoef_ppc = cat(3,corrcoef_ppc,corrcoef2D);
    
    if take_corr_abs
        corrcoef_ppc = abs(corrcoef_ppc);
    end
    
    figure;
    for j = 1:size(corrcoef_ppc,2)
        subplot(4,1,j); bar(squeeze(corrcoef_ppc(:,j,:)));
        if take_abs ylim([0 1]);
        else ylim([-1 1]);
        end
    end
    
    addpath(genpath('~/src/ds_kb/funcs_general/lib_dav'))
    
    mu_ccp = mean(corrcoef_ppc,2);
    std_ccp = std(corrcoef_ppc,[],2);
    
    
    %figure; bar(squeeze(mean(corrcoef_ppc,2))); legend('Pre','Post','Chronic'); title('Average Correlations')
    figure; bar_matrix3D(corrcoef_ppc); legend('Pre','Post','Chronic'); title('Average Correlations')
    figure; plot_matrix3D(corrcoef_ppc,'fs',1,'showErrorbars',1); title('Average Correlations')
    xlabel('freq band');ylabel('corrl')
    legend('pre','post','chr')
    %sub2ind(size(corrcoef3D),repmat(repmat(1:8,1,8),1,4),[],[])

    prepostchronic_temp = 1;     % 0 for all; 1 for pre; 2 for latency; 3 for chronic
    plottest_timeseries = 0;
    
end




function mat2D = mat3Dto2Dyz (mat3D,yzpairs)


    [Nx Ny Nz] = size(mat3D);
    %yzpairs = [3 3 4 3; 1 2 3 4];
    Nyz = size(yzpairs,2);
    xind = repmat(1:Nx,1,Nz);
    yind = reshape(repmat(yzpairs(1,:),Nx,1),1,Nyz*Nx);
    zind = reshape(repmat([yzpairs(2,:)],Nx,1),1,Nyz*Nx);
    linind = sub2ind(size(mat3D),xind,yind,zind);

    %Check
    %Atemp = reshape(1:(Nx*Ny*Nz),[Nx Ny Nz]);
    %linind2 = [];
    %for i = 1:Nyz
    %    linind2 = [linind2 Atemp(:,yzpairs(1,i),yzpairs(2,i))'];
    %end
    %sum(abs(linind-linind2))

    mat2D = reshape(mat3D(linind),[Nx Nz]);
end

