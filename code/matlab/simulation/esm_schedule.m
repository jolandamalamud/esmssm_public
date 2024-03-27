function time_schedule = esm_schedule(days, start_time, end_time, beeps)

    % semi random esm scheduling ensuring at least 30min inbetween:
    % days = 2;
    % start_time = 9*60;
    % end_time = 22*60;
    % beeps = 10;

    awake_interval = [start_time:end_time];
    sub_intervals = (end_time - start_time) / beeps;
    time_schedule = [];
    
    for d = 1:days
        for b = 1:beeps
            interval = awake_interval((b-1)*sub_intervals+1:b*sub_intervals);
            randomIndex = randi(length(interval), 1);
            timing(b) = interval(randomIndex);
            
            while (b>1) && (timing(b) - timing(b-1)) < 30
                randomIndex = randi(length(interval), 1);
                timing(b) = interval(randomIndex);
            end
        end
        
        time_schedule = [time_schedule;[repmat(d,beeps,1), ...
            timing' + (d-1)*24*60]];
    end
    
end