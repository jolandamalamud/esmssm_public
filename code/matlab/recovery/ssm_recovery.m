% ssm recovery
option.filepath = '~/phd/projects/esmssm/results/modelfits/';
option.dataset = 'leemput';
option = set_universal_options(option);
load([option.filepath, option.dataset, '/ssmfit.mat'], 'fit');
ssmfit = fit;

%% run recovery for all subjects
maxrun = 100;
Nsj = length(ssmfit);

parfor k = 1:Nsj
    
    fprintf('run recovery for subject %i \n', k);
    fprintf('**********************************\n');
    
    nettrue = ssmfit{k}.est;
    nettrue.input = ssmfit{k}.input;
    
    suboption = set_data_options(option, ssmfit{k});
    suboption.empirical = true;
    for r = 1:maxrun
        sim = ssm_simulate_sparse(nettrue, suboption);
        
        nettrue.x = sim.x;
        netest = ssm_run_algorithm(suboption, nettrue);
        recovery{k}.est{r} = netest;
        recovery{k}.true{r} = sim;
    end
    
end

%%
for r = 1:maxrun
    for k = 1:Nsj
        [eig_true(:,:,k,r),~,~] = svd(ctrb(recovery{k}.true{r}.A+recovery{k}.true{r}.W, ...
            recovery{k}.true{r}.C));
        [eig_est(:,:,k,r) ,~,~] = svd(ctrb(recovery{k}.est{r}.A+recovery{k}.est{r}.W, ...
            recovery{k}.est{r}.C));
        dot_product(:,k,r) = diag(eig_true(:,:,k,r)' * eig_est(:,:,k,r));
    end
end

for r = 1:maxrun
    for k = 1:Nsj
        [v, e] = eig(recovery{k}.true{r}.A+recovery{k}.true{r}.W);
        e = diag(e);
        [~,I] = sort(e, 'descend');
        eigval_true(:,k,r) = e(I);
        eigvec_true(:,:,k,r) = v(:,I);
        [v, e] = eig(recovery{k}.est{r}.A+recovery{k}.est{r}.W);
        e = diag(e);
        [~,I] = sort(e, 'descend');
        eigval_est(:,k,r) = e(I);
        eigvec_est(:,:,k,r) = v(:,I);
        
        dot_product(:,k,r) = diag(eigvec_true(:,:,k,r)' * eigvec_est(:,:,k,r));
        A_true(:,:,k,r) = (recovery{k}.true{r}.A+recovery{k}.true{r}.W);
        A_est(:,:,k,r) = (recovery{k}.est{r}.A+recovery{k}.est{r}.W);
        
    end
end

% errorbar_plot(dot_product, {'1st', '2nd', '3rd', '4th'}, ...
%     'unit vectors of controllability Gramian', true)
% xlabel('unit vectors')