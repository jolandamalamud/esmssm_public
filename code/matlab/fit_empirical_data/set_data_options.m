function option = set_data_options(option, data)
    
    option.nonNaNidx = ~isnan(data.data(1, :)); % non NaN entries
    option.dimZ = size(data.data, 1); % number hidden states
    option.dimX = size(data.data, 1); % number observations
    option.dimInp = size(data.input, 1); % number observations
    option.timing = find(~isnan(data.data(1,:))); % observation timings
    option.gap = mean(option.timing(2:end) - option.timing(1:end-1));
    option.T = size(data.data,2);
    option.tau = size(data.data(:,option.nonNaNidx), 2);
  
end