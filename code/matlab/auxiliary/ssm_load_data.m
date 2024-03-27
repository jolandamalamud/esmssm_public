function datafit = ssm_load_data(option)

    moodfiles = dir(sprintf('%s/data_prep/%s/mooddata*.csv', ...
            option.filepath, option.dataset));

    if option.input
        inputfiles = dir(sprintf('%s/data_prep/%s/inputdata*.csv', ...
            option.filepath, option.dataset));
    end
    
    
    for s = 1:length(moodfiles) 
        % load mooddata per subject
        subject_mood = importdata(sprintf('%s/data_prep/%s/%s', ...
                    option.filepath, option.dataset, moodfiles(s).name));

        timing = find(~isnan(subject_mood(1,:)));
        datafit{s}.data = subject_mood(:, 1:timing(end));

        if option.input
            % load inputdata per subject
            subject_input = importdata(sprintf('%s/data_prep/%s/%s', ...
                    option.filepath, option.dataset, inputfiles(s).name));
            datafit{s}.input = subject_input(:, 1:timing(end));
        else
            datafit{s}.input = zeros(2, timing(end));
        end
    end
end