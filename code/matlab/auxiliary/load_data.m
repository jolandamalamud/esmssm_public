function fit = load_data(option)

    data = readtable([option.filepath, option.dataset, '/mooddata.csv']);
    columns = data.Properties.VariableNames;
    columns = convertCharsToStrings(columns);
    mood_idx = find(contains(columns, 'mood'));
    ns = length(mood_idx);

    sub_list = unique(data.id(~isnan(data.id)));
    Nsj = length(sub_list);

    if option.input
        inputdata = readtable([option.filepath, option.dataset,'/inputdata.csv']);
    end

    for k = 1:Nsj
       idx = find(data.id == k);
       T = data.timing(idx);
       try fit{k}.id = data.originalid(idx(1)); end
       try fit{k}.group = data.group(idx(1)); end
       fit{k}.timing = T;
       fit{k}.data = nan(ns, T(end)); 
       fit{k}.data(:, T) = table2array(data(idx,mood_idx))';
       if option.input
           input = table2array(inputdata(idx,3:end));
           nonnanidx = ~isnan(input(1,:));
           fit{k}.input = zeros(sum(nonnanidx), T(end)); 
           fit{k}.input(:, T) = input(:,nonnanidx)';
           fit{k}.input(sum(fit{k}.input,2)==0,:) = [];
       else
           fit{k}.input = zeros(2, T(end));
       end
    end
end
