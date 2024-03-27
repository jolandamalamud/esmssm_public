option.filepath = '~/phd/projects/esmssm/results/modelfits/';

datasets = {'leemput', 'twindata'};
for d = 1:2
    option.dataset = datasets{d};
    fit = load([option.filepath, option.dataset, '/ssmfit.mat']);
    subfit{d} = fit.fit;
end
ssmfit = horzcat(subfit{1}, subfit{2});
load([option.filepath, '/oufit.mat']);
option = set_universal_options(option);

%% compare SSM, OU process and Gaussian to fit esm data

ns = size(ssmfit{1}.data,1);
num_params = [(ns^2+ns)/2 + ns, 2*ns, 3*ns, ns+ns^2+ns+ns];
Nsj = length(ssmfit);
for s = 1:Nsj
    data = ssmfit{s}.data;
    num_obs(s) = sum(~isnan(ssmfit{s}.data(1,:)));
    
    % Gaussian
    try
        pd = fitgmdist(data', 1, 'CovarianceType', 'full');
        ll_gauss(s) = -pd.NegativeLogLikelihood;
        d = fitgmdist(data', 1, 'CovarianceType', 'diagonal');
        ll_gauss_diag(s) = -d.NegativeLogLikelihood;
    catch
        ll_gauss(s) = nan;
        ll_gauss_diag(s) = nan;
    end
    
    ll_ou(s) = -oufit{s}.LL;
    
    ll_ssm(s) = ssmfit{s}.est.LL(end);
    
    [aic(:,s),bic(:,s)] = aicbic([ll_gauss(s), ll_gauss_diag(s),  ll_ou(s), ll_ssm(s)], ...
    num_params, repmat(num_obs(s),1,length(num_params)));

    bestmodel(s) = find(bic(:,s) == min(bic(:,s)));
    
end



% [aic,bic] = aicbic([sum(ll_gauss), sum(ll_ou), sum(ll_ssm)], ...
%     num_params, [mean(num_obs), mean(num_obs), mean(num_stats)]);
