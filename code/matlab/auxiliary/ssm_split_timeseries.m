function [splitdata] = ssm_split_timeseries(data, inputs, set_length, option)

    icp = ssm_changepoint(data);
    icp_split = [1, icp, size(data,2)];
    
    if isempty(set_length)
        for k = 1:length(icp_split)-1
            splitdata{k}.data = data(:, icp_split(k):icp_split(k+1));
            if option.input
                splitdata{k}.input = inputs(:, icp_split(k):icp_split(k+1));
            end
        end
    else
        ii1 = find(sum(isnan(data(:, 1:set_length)),1) == 0);
        ii2 = find(sum(isnan(data(:, icp:icp+set_length)),1) == 0);
        splitdata{1}.data = data(:, 1:ii1(end));
        splitdata{2}.data = data(:, icp:icp+ii2(end)-1);
        if option.input
            splitdata{1}.input = inputs(:, 1:ii1(end));
            splitdata{2}.input = inputs(:, icp:icp+ii2(end)-1);
        end
    end

end