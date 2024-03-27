function [cvdata, cvinp] = ssm_cvsplit(data, inp)
    
    nonNANidx = find(~isnan(data(1,:)));
    
    ii = rand(size(data(:,nonNANidx),2), 1) < 0.9;
    cvdata = nan(size(data));
    cvdata(:, nonNANidx(ii)) = data(:, nonNANidx(ii));
    
    if ~isempty(inp)
        cvinp = nan(size(inp));
        cvinp(:, nonNANidx(ii)) = inp(:, nonNANidx(ii));
    else
        cvinp = zeros(2, size(cvdata,2));
    end
end