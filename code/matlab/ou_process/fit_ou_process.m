clear all; close all;

option.filepath = '~/phd/projects/esmssm/data/data_prep/';
option.goalpath = '~/phd/projects/esmssm/results/modelfits';
datasets = {'leemput', 'twindata'};

for d = 1:2
    option.dataset = datasets{d};
    option.input = false;
    data_dataset{d} = load_data(option);
end

% concat data
data = horzcat(data_dataset{1},data_dataset{2});

option.Nsj = length(data);
option.ns = size(data{1}.data,1);
fminoptions = optimoptions(@fminunc,'Display', 'off');

%% 

parfor sj = 1 : option.Nsj

    fprintf('subject %i   \n',sj);
    fprintf('**********************************\n');
    
    subject_mood = data{sj}.data;

    subject_timing = find(~isnan(subject_mood(1,:)));
    subject_time_interval = subject_timing(2:end) - subject_timing(1:end-1);

    mood = subject_mood(:, subject_timing);
    
    ll = []; exit_flag = [];
    for m = 1:option.ns
        [ou_est, ll(m), exit_flag(m)] = fminunc(@(p)loglike_ou(mood(m,:),p(1),p(2),p(3),subject_time_interval),[1,1,1], fminoptions);
        oufit{sj}.parameter(:,m) = ou_est;
    end
    
    if any(exit_flag > 3); ll = nan; end
    
    oufit{sj}.LL = sum(ll);
    oufit{sj}.data = mood;
    oufit{sj}.timing = subject_timing;
    
end

save([option.goalpath, '/oufit.mat'], 'oufit');

%% ou process recovery
maxrun = 10;
parfor sj = 1 : option.Nsj
    
    subject_timing = oufit{sj}.timing;
    subject_time_interval = subject_timing(2:end) - subject_timing(1:end-1);
    
    for r = 1:maxrun
        for m = 1:option.ns

             x = simulate_ou(oufit{sj}.data(m,1), oufit{sj}.parameter(1,m), ...
                  oufit{sj}.parameter(2,m),  oufit{sj}.parameter(3,m), ...
                  oufit{sj}.timing);

             recovery{sj}.est(:,m,r) = fminunc(@(p)loglike_ou(x,p(1),p(2),p(3),...
                 subject_time_interval),[1,1,1], fminoptions);

        end
    end
    
    recovery{sj}.true(:,:) = oufit{sj}.parameter;
    
end

save([option.goalpath, '/oufit_recovery.mat'], 'recovery');

%% plot recovery correlation
for k = 1:option.Nsj
    true_par(:,:,k) = recovery{k}.true(:,:,end); 
    est_par(:,:,k,:) = recovery{k}.est; 
end

for r = 1:maxrun
    for p = 1:3
        t = true_par(p,:,:); 
        e = est_par(p,:,:,r); 
        par_corr(r, p) = corr(t(:), e(:)); 
    end
end

errorbar_plot(par_corr, {'mu', 'sigma', 'lambda'}, ...
    'correlation between true and estimated parameter', true)

