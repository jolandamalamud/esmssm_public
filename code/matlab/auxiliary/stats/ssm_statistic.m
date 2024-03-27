function statistics = ssm_statistic(option, stats)

    for s = 1:size(stats, 2)

        netest = stats(s).netest;
        nettrue = stats(s).nettrue;
        
%         aw = (nettrue.A + nettrue.W)^option.gap;
%         nettrue.A = diag(aw); nettrue.W = aw - diag(aw);
%         nettrue.S = option.gap * nettrue.S;
%         nettrue.h = option.gap * nettrue.h;
%         aw = (netest.A + netest.W)^option.gap;
%         netest.A = diag(aw); netest.W = aw - diag(aw);
%         netest.S = option.gap * netest.S;
%         netest.h = option.gap * netest.h;
        
        variables = {'A', 'W'};
        if option.bias; variables{length(variables) + 1} = 'h'; end
        if option.input; variables{length(variables) + 1} = 'C'; end
        if ~option.S_fixed; variables{length(variables) + 1} = 'S'; end
        if ~option.B_fixed && option.B_full; variables{length(variables) + 1} = 'B'; end
        variables{length(variables) + 1} = 'G';
        
        % calculate root mean squared error and correlation between 
        % true and estimated parameter for all parameters
        for p = 1:length(variables)
            try 
                eval(['pred=' sprintf('netest.%s', variables{p}) ';'])
                eval(['tru=' sprintf('nettrue.%s', variables{p}) ';'])
                tmp1 = tru(:); tmp1 = tmp1(tmp1 ~= 0);
                tmp2 = pred(:); tmp2 = tmp2(tmp2 ~= 0);
                mseP(s, p) = nanmean((tmp1' - tmp2').^2);
                corrP(s, p) = corr(tmp1, tmp2);
            catch
                mseP(s, p) = nan;
                corrP(s, p) = nan;
            end
                
        end
%         if any(mseP(s,:) > 2)
%           	keyboard
%         end
        
        % calculate root mean squared error and correlation between 
        % true and estimated states
        if option.state_stats
            idx_nan = isnan(nettrue.x(1,:));
            z_true = nettrue.z(:, ~idx_nan);
            z_est = netest.z(:, ~idx_nan);
            mseZ(s) = mean((z_true(:) - z_est(:)).^2); 
            corrZ(s) = corr(z_true(:), z_est(:));
        end
        
        fp_est(:,s) = ssm_fixed_point(netest);
        fp_true(:,s) = stats(s).nettrue.FP;
        mseP(s, length(variables) + 1) = nanmean((stats(s).nettrue.FP - fp_est(:,s)).^2);
        corrP(s, length(variables) + 1) = corr(stats(s).nettrue.FP, fp_est(:,s));
    end

%     idx_outlier = find(isoutlier(rmseP(:,6)));
%     disp(['deleting ' num2str(length(idx_outlier)) ' outliers: ' num2str(idx_outlier)]);
%     rmseZ(idx_outlier) = [];
%     rmseP(idx_outlier, :) = [];
%     corrZ(idx_outlier) = [];
%     corrP(idx_outlier, :) = [];
    
    %% plot RMSE and correlation

    label = {'dynamics diagonals', 'dynamics off-diagonals'};
    if option.bias; label{length(label) + 1} = 'bias'; end
    if option.input; label{length(label) + 1} = 'input weights'; end
    if ~option.S_fixed; label{length(label) + 1} = 'process noise'; end
    if ~option.B_fixed && option.B_full; label{length(label) + 1} = 'omission weights'; end
    label{length(label) + 1} = 'observation noise';
    label{length(label) + 1} = 'fixed point';

    figure();
    
    if option.state_stats
        n_subplots = size(mseP, 2) + 1;
    else
        n_subplots = size(mseP, 2);
    end
    
    for k=1:size(mseP, 2)
        subplot(3, n_subplots, k);
        errorbar_plot(mseP(:, k), label{k}, 'MSE', false);
        subplot(3, n_subplots, k + n_subplots);
        errorbar_plot(corrP(:, k), label{k}, 'correlation', false);
        ylim([-0.4 1.2])
%         disp([nanmean(mseP(mseP(:, k) ~= 0, k)), nanstd(rmseP(rmseP(:,k) ~= 0, k))]);
    end
    
    statistics.parameter = [mseP, corrP];
    
    if option.state_stats
        subplot(3, n_subplots, n_subplots);
        errorbar_plot(mseZ, 'latent state', 'MSE', false);
        subplot(3, n_subplots, 2 * n_subplots);
        errorbar_plot(corrZ, 'latent state', 'correlation', false);
%         disp([nanmean(rmseZ(rmseZ ~= 0)), nanstd(rmseZ(rmseZ ~= 0))]);
        statistics.state = [mseZ, corrZ];
        
    end
    
    if option.sparse
        sgtitle('MSE between true and estimated sparse system', 'FontSize', 18);
    else
        sgtitle('MSE between true and estimated dense system', 'FontSize', 18);
    end

%% plot all true against all estimated parameter and state values
    datatable1 = []; datatable2 = [];
    nees = horzcat(stats.netest);
    netu = horzcat(stats.nettrue);
    for p = 1:length(variables)
        subplot(3,n_subplots,2*n_subplots+p);
        eval(['x = [netu.' variables{p} '];']); 
        eval(['y = [nees.' variables{p} '];']); 
        x = x(x ~= 0); y = y(y ~= 0);
%         datatable1 = [datatable1; reshape(x,length(x)/1000,1000)];
%         datatable2 = [datatable2; reshape(y,length(y)/1000,1000)];
        plot(x, y, 'o'); 
        xlabel(['true ' label{p}]); ylabel(['estimated ' label{p}]); 
        hold on; ax = axis;
        plot([min(ax), max(ax)],[min(ax), max(ax)], 'r');    
    end
    
    subplot(3,n_subplots,2*n_subplots+length(variables)+1);
    plot(fp_true(:), fp_est(:), 'o'); 
    xlabel('true fixed point'); ylabel('estimated fixed point'); 
    hold on; ax = axis;
    plot([min(ax), max(ax)],[min(ax), max(ax)], 'r'); 

    if option.state_stats
        subplot(3,n_subplots,2*n_subplots+length(variables)+2);
        x = [netu.z]; x = x(x~=0);
        y = [nees.z]; y = y(y~=0);
        plot(x, y, 'o'); 
        xlabel('true latent states'); ylabel('estimated latent states'); 
        hold on; ax = axis;
        plot([min(ax), max(ax)],[min(ax), max(ax)], 'r');   
    end
    set(findall(gcf,'-property','FontSize'),'FontSize',18);
    
%% scatter plot
end