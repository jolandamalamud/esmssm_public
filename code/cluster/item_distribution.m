function hh = item_distribution(data, all_items, variable)
    item_matrix = nan(length(data), length(all_items));
    for k = 1:length(data)
        eval(['dd = data{k}.' variable 'items_included']);
        item_matrix(k,:) = ismember(all_items, dd); 
    end
    hh = nansum(item_matrix,1)';
%     figure(); bar(hh(hh~=0)); xticks([1:length(all_items)]);
%     xticklabels(all_items(hh~=0)); xtickangle(45);
end