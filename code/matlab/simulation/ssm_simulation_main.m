%% ssm_master
% EM on simulated time series from a linear state space model:
% latent model: z_t = (A+W) * z_(t-1) + h + C * s_t + err1, err1 = N(0, S)
% observation model: x_t = z_t + err2, err2 = N(0, G)

clear all; close all; addpath(genpath(pwd));
maxruns = 10;
% set options
option = ssm_sim_set_options();
option.runs = maxruns;
days_fields = [];

% change length of time-series, i.e. days
for d = [2, 6, 10, 14]
    option.days = d;
    parfor r = 1:option.runs
        rng(r);
        % for different settings same randomizations
        sub=struct;
        suboption = option;
        suboption.consider_gap=false;
        % create ground truth for the ssm
        suboption.empirical = false;
        [nettrue, suboption] = ssm_ground_truth(suboption);

        % simulate time series
        if suboption.sparse
            nettrue = ssm_simulate_sparse(nettrue, suboption);
        else
            nettrue= ssm_simluate(nettrue, suboption);
        end

        suboption.nonNaNidx = ~isnan(nettrue.x(1,:));
        suboption.empirical = true;
        netest = ssm_run_algorithm(suboption, nettrue, []);

        if isempty(netest)
            simulations_fields{r}.netest = nan;
            continue
        else
            simulations_fields{r}.netest = netest;
            simulations_fields{r}.nettrue = nettrue;
            simulations_fields{r}.option = suboption;
        end
    end
    days_fields = [days_fields; simulations_fields];
end

save('simulations_time', 'days_fields','-v7.3');
kf_simulations.days = test_recovery(days_fields, option);

%%

% set options
option = ssm_sim_set_options();
option.runs = maxruns;
noise_fields = [];

% change measurment noise
for d = [0.1, 0.5, 1, 2, 5]
    parfor r = 1:option.runs
        rng(r);
        % for different settings same randomizations
        suboption = option
        suboption.consider_gap=false;
        % create ground truth for the ssm
        suboption.empirical = false;
        [nettrue, suboption] = ssm_ground_truth(suboption);
        
        nettrue.G = d * eye(suboption.dimX);
        
        % simulate time series
        if suboption.sparse
            nettrue = ssm_simulate_sparse(nettrue, suboption);
        else
            nettrue= ssm_simluate(nettrue, suboption);
        end

        suboption.nonNaNidx = ~isnan(nettrue.x(1,:));
        suboption.empirical = true;
        netest = ssm_run_algorithm(suboption, nettrue, []);

        if isempty(netest)
            simulations_fields{r}.netest = nan;
            continue
        else
            simulations_fields{r}.netest = netest;
            simulations_fields{r}.nettrue = nettrue;
            simulations_fields{r}.option = suboption;
        end
    end
    noise_fields = [noise_fields; simulations_fields];
end
save('simulations_noise', 'noise_fields','-v7.3');
kf_simulations.noise = test_recovery(noise_fields, option);

%%

% set options
option = ssm_sim_set_options();
option.runs = maxruns;
dim_fields = [];

% change dimensionality of emotions
for d = [2,4,8,16]
    parfor r = 1:option.runs
        % for different settings same randomizations
        suboption = option;
        suboption.dimX = d;
        suboption.dimZ = d;
        suboption.consider_gap=false;
        % create ground truth for the ssm
        suboption.empirical = false;
        [nettrue, suboption] = ssm_ground_truth(suboption);

        % simulate time series
        if suboption.sparse
            nettrue = ssm_simulate_sparse(nettrue, suboption);
        else
            nettrue= ssm_simluate(nettrue, suboption);
        end

        suboption.nonNaNidx = ~isnan(nettrue.x(1,:));
        suboption.empirical = true;

        netest = ssm_run_algorithm(suboption, nettrue, []);

        if isempty(netest)
            simulations_fields{r}.netest = nan;
            continue
        else
            simulations_fields{r}.netest = netest;
            simulations_fields{r}.nettrue = nettrue;
            simulations_fields{r}.option = suboption;
        end
    end
    dim_fields = [dim_fields; simulations_fields];
end

save('simulations_dim', 'dim_fields','-v7.3');
kf_simulations.dim = test_recovery(dim_fields, option);

%%

% set options
option = ssm_sim_set_options();
option.runs = maxruns;
gap_fields = [];

% compare sparse and dense
for gap = [false, true]
    parfor r = 1:option.runs
        % for different settings same randomizations
        rng(r);
        suboption = option;
        suboption.consider_gap = false;
        suboption.sparse = gap;
        % create ground truth for the ssm
        suboption.empirical = false;
        [nettrue, suboption] = ssm_ground_truth(suboption);

        % simulate time series
        if suboption.sparse
            nettrue = ssm_simulate_sparse(nettrue, suboption);
        else
            nettrue= ssm_simulate(nettrue, suboption);
        end

        suboption.nonNaNidx = ~isnan(nettrue.x(1,:));
        suboption.empirical = true;
        netest = ssm_run_algorithm(suboption, nettrue, []);

        if isempty(netest)
            simulations_fields{r}.netest = nan;
            continue
        else
            simulations_fields{r}.netest = netest;
            simulations_fields{r}.nettrue = nettrue;
            simulations_fields{r}.option = suboption;
        end
    end
    gap_fields = [gap_fields; simulations_fields];
end
%%

% set options
option = ssm_sim_set_options();
option.runs = 1000;
gap_fields = [];

% compare slow and fast processes
for gap = [true, false]
%     simulations_fields = struct();
    parfor r = 1:option.runs
        % for different settings same randomizations
        rng(r);
        suboption = option;
        suboption.consider_gap = gap;
        % create ground truth for the ssm
        suboption.empirical = false;
        [nettrue, suboption] = ssm_ground_truth(suboption);

        % simulate time series
        if suboption.sparse
            nettrue = ssm_simulate_sparse(nettrue, suboption);
        else
            nettrue= ssm_simulate(nettrue, suboption);
        end
%         figure();plot(nettrue.x(1,1:1000),'o'); hold on; plot(nettrue.z(1,1:1000));
        suboption.nonNaNidx = ~isnan(nettrue.x(1,:));
        suboption.empirical = true;
        netest = ssm_run_algorithm(suboption, nettrue, []);

        if isempty(netest)
            simulations_fields{r}.netest = nan;
            continue
        else
            simulations_fields{r}.netest = netest;
            simulations_fields{r}.nettrue = nettrue;
            simulations_fields{r}.option = suboption;
        end
    end
    gap_fields = [gap_fields; simulations_fields];
end
keyboard();

%%
% clear all; close all; addpath(genpath(pwd));
maxruns = 10;
% set options
option = ssm_sim_set_options();
option.input = true;
option.runs = maxruns;
input_fields = [];

for inp = [5, 10, 15, 20]
    option.dimInp = inp;
    option.days  = 100;
    parfor r = 1:option.runs
        rng(r);
        % for different settings same randomizations
        sub=struct;
        suboption = option;
        suboption.consider_gap=false;
        % create ground truth for the ssm
        suboption.empirical = false;
        [nettrue, suboption] = ssm_ground_truth(suboption);

        % simulate time series
        if suboption.sparse
            nettrue = ssm_simulate_sparse(nettrue, suboption);
        else
            nettrue= ssm_simluate(nettrue, suboption);
        end

        suboption.nonNaNidx = ~isnan(nettrue.x(1,:));
        suboption.empirical = true;
        nettrue.input(isnan(nettrue.input)) = 0;
        netest = ssm_run_algorithm(suboption, nettrue, []);

        if isempty(netest)
            simulations_fields{r}.netest = nan;
            continue
        else
            simulations_fields{r}.netest = netest;
            simulations_fields{r}.nettrue = nettrue;
            simulations_fields{r}.option = suboption;
        end
    end
    input_fields = [input_fields; simulations_fields];
end

kf_simulations.inputs = test_recovery(input_fields);
% save('simulations_inputs.mat', input_fields);
%%
clear all; close all; addpath(genpath(pwd));
maxruns = 1000;
% set options
option = ssm_sim_set_options();
option.runs = maxruns;

% change length of time-series, i.e. days
option.days = 14;
parfor r = 1:option.runs
    rng(r);
    % for different settings same randomizations
    sub=struct;
    suboption = option;
    suboption.consider_gap=false;
    % create ground truth for the ssm
    suboption.empirical = false;
    [nettrue, suboption] = ssm_ground_truth(suboption);
    nettrue.G = 0.1 * eye(suboption.dimX);

    % simulate time series
    if suboption.sparse
        nettrue = ssm_simulate_sparse(nettrue, suboption);
    else
        nettrue= ssm_simluate(nettrue, suboption);
    end

    suboption.nonNaNidx = ~isnan(nettrue.x(1,:));
    suboption.empirical = true;
    netest = ssm_run_algorithm(suboption, nettrue, []);

    if isempty(netest)
        simulations_fields{r}.netest = nan;
        continue
    else
        simulations_fields{r}.netest = netest;
        simulations_fields{r}.nettrue = nettrue;
        simulations_fields{r}.option = suboption;
    end
end


%%
% remove outliers
% ex = isoutlier(kf_simulations.noise.mseA(:,2), "quartiles");
% VAR = varm(4,1);
% r = 1;
% for i = 1:option.runs
%     if ex(i) == 0
%         data{r} = noise_fields{1,i};
%         try
%             varest{r} = estimate(VAR, data{r}.nettrue.x');
%         end
%         r=r+1;
%     end
% end
VAR = varm(4,1);
data = days_fields;
for i = 1:length(data)
    data{i} = days_fields{2,i};
    try
        varest{i} = estimate(VAR, data{i}.nettrue.x(:,~isnan(data{i}.nettrue.x(1,:)))');
    end
end


for i =1:length(data)
    Mt = data{i}.nettrue.A+data{i}.nettrue.W; 
    At(:, i) = Mt(:);
    Me = varest{i}.AR{1,1}; 
    var.Ae(:, i) = Me(:);
    
    var.mseA(i) = sqrt(mean((Mt - Me).^2, 'all'));
    var.mseAdiag(i) = sqrt(mean((diag(Mt) - diag(Me)).^2, 'all'));
    var.mseAoff(i) = sqrt(mean(((Mt - diag(Mt)) - (Me - diag(Me))).^2, 'all'));

    vv1 = get_eigvec(data{i}.nettrue.A+data{i}.nettrue.W);
    vv2 = get_eigvec(varest{i}.AR{1,1});

    var.dotprod(:, i) = diag(abs(vv1' * vv2));

end


%%
for i =1:length(data)
    T(:,:,i) = (data{i}.nettrue.A+data{i}.nettrue.W);
    vvt(:,:,i) = get_eigvec(data{i}.nettrue.A+data{i}.nettrue.W);
    TH(:,:,i) = (data{i}.nettrue.A+data{i}.nettrue.W)^60;
    M(:,:,i) = (data{i}.netest.A+data{i}.netest.W);
    vve(:,:,i) = get_eigvec(data{i}.netest.A+data{i}.netest.W);
    MH(:,:,i) = (data{i}.netest.A+data{i}.netest.W)^60;
    AR(:,:,i) = varest{i}.AR{1,1};
    vvvar(:,:,i) = get_eigvec(varest{i}.AR{1,1});
end
%%
k=1;
for i = [0.1,0.3,0.7]
    for j = 1:1000
        y = squeeze(T(4,4,:))*i+0.001*randn(length(data),1); 
        [rt(k,j), pt(k,j)] = corr(squeeze(TH(4,4,:)), y, 'type', 'spearman');
        [rkf(k,j), pkf(k,j)] = corr(squeeze(MH(4,4,:)), y, 'type', 'spearman'); 
        [rvar(k,j), pvar(k,j)] = corr(squeeze(AR(4,4,:)), y, 'type', 'spearman');
        if abs(rkf(k,j) - rvar(k,j)) > 0.1; keyboard(); end
    end
    k = k+1;
end

[nanmean(rt,2),nanmean(rkf,2),nanmean(rvar,2)]

