option.filepath = '';
option.dataset = 'twindataset';
ask_status = false;
if ~exist([option.filepath '/data_prep_cluster/' option.dataset], 'dir')
   mkdir([option.filepath '/data_prep_cluster/' option.dataset])
end

moodfile = importdata(sprintf('%s/data_prep/%s/mooddata.csv', ...
        option.filepath, option.dataset));

ii = (strfind(moodfile.textdata, 'mood_'));
mood_items = find(not(cellfun('isempty',ii)));%[3:length(moodfile.textdata)-1]; %find(not(cellfun('isempty',ii)));%[3,9,7,11];%[3:length(moodfile.textdata)]; %find(not(cellfun('isempty',ii)));%[3:28];%[3,9,7,12];%
mood_labels = moodfile.textdata(mood_items);
for k = 1:length(mood_labels); mood_labels{k} = mood_labels{k}(6:end); end


inputfile = importdata(sprintf('%s/data_prep/%s/inputdata.csv', ...
    option.filepath, option.dataset));


subjid_column = ismember(moodfile.textdata,'subjid');
subjid_timing = ismember(moodfile.textdata,'timing');
if ask_status
    subjid_status = ismember(moodfile.textdata,'status');
end
subj_list = unique(moodfile.data(:, subjid_column));

for s = 1:length(subj_list) 

    subj_idx = moodfile.data(:, subjid_column) == subj_list(s);

    timing = moodfile.data(subj_idx, subjid_timing);
    [cleaned_subject_mood, ii] = ...
        ssm_clean_variables(moodfile.data(subj_idx, mood_items)');
    timing = timing(~isnan(cleaned_subject_mood(1,:)));
    timing = timing - (timing(1) - 1);

    subject_mood = nan(size(cleaned_subject_mood,1),  max(timing));
    subject_mood(:, timing) = ...
        cleaned_subject_mood(:,~isnan(cleaned_subject_mood(1,:)));

    datafit.id = subj_list(s);
    if ask_status
        datafit.status = nanmean(moodfile.data(subj_idx,subjid_status));
        ss(s) = datafit.status;
    end
    datafit.data  = subject_mood;
    datafit.mooditems_included = mood_labels(ii);


    subject_input = nan(size(inputfile,2), max(timing));
    sub_input_full = inputfile(subj_idx, :)';
    subject_input(:, timing) = ...
        sub_input_full(:, ~isnan(cleaned_subject_mood(1,:)));
    [cleaned_subject_input, ii] = ssm_clean_input(subject_input);
    datafit.input = ssm_input_extension(cleaned_subject_input);
    datafit.inputitems_included = find(ii);

%     save(sprintf('%s/data_prep_cluster/%s/datafit_%s.csv', ...
%         option.filepath, option.dataset, num2str(s)), 'datafit');
   
    all_data{s} = datafit;
end
%%
item_list = ...
item_distribution(all_data, mood_labels, 'mood');
figure(); bar([item_list(item_list~=0)], 'stacked');
xticks([1:length(item_list~=0)]); xticklabels(mood_labels(item_list~=0)); xtickangle(45);

sss = unique(ss);
for k = 1:length(unique(ss))
    item_list{k} = ...
        item_distribution(all_data(ss==sss(k)), mood_labels, 'mood');
    item_list_size(k) = length(item_list{k});
    input_list{k} = item_distribution(all_data(ss==sss(k)), ...
        [1:size(inputfile,2)], 'input');
    

end
figure(); bar([item_list{1}(item_list{1}~=0)], 'stacked');
xticks([1:length(item_list{1}~=0)]); xticklabels(mood_labels(item_list{1}~=0)); xtickangle(45);
ylabel('subjects'); title('study 7.0');
set(findall(gcf,'-property','FontSize'),'FontSize',30);

figure(); bar([input_list], 'stacked');
xlabel('input items');
ylabel('subjects'); title('study 7.0');
set(findall(gcf,'-property','FontSize'),'FontSize',30);


[~,m] = max(item_list_size);
figure(); bar([item_list{1}(item_list{m}~=0),item_list{2}(item_list{m}~=0)], 'stacked');
xticks([1:length(item_list{m}~=0)]); xticklabels(mood_labels(item_list{m}~=0)); xtickangle(45);
ylabel('subjects'); title('study 7.0'); legend('controls', 'patients');
set(findall(gcf,'-property','FontSize'),'FontSize',30);

figure(); bar([input_list{1},input_list{2}], 'stacked');
xlabel('input items');
ylabel('subjects'); title('study 7.0'); legend('controls', 'patients');
set(findall(gcf,'-property','FontSize'),'FontSize',30);

%%


function [reduced_data, ii] = ssm_clean_variables(data)
    
    ii = ~isnan(nanmean(data, 2));%nanvar(data, 0, 2) ~= 0 & ~isnan(nanvar(data, 0, 2));
    reduced_data = data(ii,:);
    
    for t = 1:size(data,2)
        if any(isnan(reduced_data(:,t)))
            reduced_data(:,t) = nan;
        end
    end

end

function [inp, ii] = ssm_clean_input(fullinp)

    % clear inputs with no variance and less than 25%
    percent_inp = nansum(fullinp,2) ./ sum(~isnan(fullinp),2);
    ii = 0 < percent_inp & percent_inp < 1;
    inp = fullinp(ii,:);
end