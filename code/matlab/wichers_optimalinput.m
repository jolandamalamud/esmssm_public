s = 3;
timing = [1:120:16*60, 24*60:120:40*60, 48*60:120:64*60];
T = timing(end);%size(datafit{s}.data,2);
timing = find(~isnan(fit{s}.data(1,1:T)));
input_timings = find(sum(fit{s}.input(:,1:T),1) ~= 0);
% input_category = inputitems(datafit{s}.inputitems_included);
% input_category = table2array(inpLabWichers(:,2));
input_category = {[1:8],[9:16]};
nl = size(fit{s}.data,1);
ni = size(fit{s}.input,1);

A = fit{s}.est.A + fit{s}.est.W;
h = fit{s}.est.h;
C = fit{s}.est.C;
S = fit{s}.est.S;
G = fit{s}.est.G;

rN = [7 7 1 1]';
X = nan(nl, T);
Z = nan(nl, T);
Z(:,1) = fit{s}.data(:,1);
X(:,1) = mvnrnd(Z(:,1), G);
inp = zeros(ni, T);

% for k = 1:2
%     combi{k} = nchoosek([0;input_category{k}'], ...
%         min(3, length(input_category{k})+1));
% end
% 
% [ca, cb] = ndgrid(combi{1},combi{2});
% combs = [ca(:), cb(:)];
% 
% dd = [input_category{1}', zeros(length(input_category{1}), 1); ...
%     zeros(length(input_category{2}), 1), input_category{2}'; combs];
[ca, cb] = ndgrid([0;input_category{1}'],[0;input_category{2}']);
dd = [ca(:), cb(:)];

% for t = 2:T
%      if ismember(t, timing)
%         y = nan(nl, size(dd,1));
%         for k = 1:size(dd,1)
%             u = zeros(ni,1);
%             u(dd(k,dd(k,:) ~= 0)) = 1;
%             y(:,k) = mvnrnd(A * Z(:,t-1) + h + C * u, S);
%         end
%         [m,l] = min(sum((y - rN).^2,1));
% 
%         inp(dd(l,dd(l,:) ~= 0),t) = 1;
% 
%         Z(:,t) = y(:,l);
%         X(:,t) = mvnrnd(Z(:,t), G);
%      else
%         Z(:,t) = mvnrnd(A * Z(:,t-1) + h, S);
%      end
%     clear y zz;    
% end

for t = 2:T-1
     if ismember(t, timing(2:end-1) - 15)
        y = nan(nl, size(dd,1));
        for k = 1:size(dd,1)
            u = zeros(ni,1);
            u(dd(k,dd(k,:) ~= 0)) = 1;
            zz(:,1) = Z(:,t-1);
            for i = 1:30
                zz(:,i+1) = mvnrnd(A * zz(:,i) + h + C * u, S);
            end
            y(:,k) = zz(:,end);
        end
        [m,l] = min(sum((y - rN).^2,1));

        inp(dd(l,dd(l,:) ~= 0),t:t+30) = 1;

        Z(:,t) = y(:,l);
        X(:,t) = mvnrnd(Z(:,t), G); 
    else
        Z(:,t) = mvnrnd(A * Z(:,t-1) + h + C * inp(:,t), S);
    end
    clear y zz;    
end

subplot(211); plot(Z');legend('cheerful','content','anxious','sad');subplot(212); imagesc(inp)

% datafit{s}.steered_state = Z;
% datafit{s}.steered_obs = X;
% datafit{s}.optimalinput = inp;
% 
% input_diff_better = sum(datafit{1}.input(:,timing) ==  datafit{3}.optimalinput(:,timing), 'all') ...
% / numel(datafit{3}.optimalinput(:,timing));