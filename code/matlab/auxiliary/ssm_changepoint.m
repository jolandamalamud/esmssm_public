function icp = ssm_changepoint(X)

    XnonNAN = X(:,~isnan(X(1,:)));
    ipt = findchangepts(XnonNAN);%,'Statistic',stats,'MinThreshold',5);
    loc = find(~isnan(X(1,:)) == 1);
    icp = loc(ipt);
    
%     figure(); 
%     plot(X', 'o'); 
%     xline(icp, 'r', 'linewidth', 3);
%     ylabel('ratings'); xlabel('minutes');

end