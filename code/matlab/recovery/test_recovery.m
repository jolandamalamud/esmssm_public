function x = test_recovery(fields, option)

    [m,n] = size(fields);
    maxdim = fields{end,end}.option.dimZ;
    var = {'nettrue', 'netest'};
    x.At = nan(maxdim^2,n,m);
    x.Ae = nan(maxdim^2,n,m);
    x.dotprod = nan(maxdim,n,m);
    
    for i = 1:n
        for j = 1:m
            
            Mt = fields{j,i}.nettrue.A+fields{j,i}.nettrue.W; 
            x.At(1:length(Mt(:)), i, j) = Mt(:);
            Me = fields{j,i}.netest.A+fields{j,i}.netest.W; 
            x.Ae(1:length(Me(:)), i, j) = Me(:);
            if option.consider_gap
                    Mt = (fields{j,i}.nettrue.A+fields{j,i}.nettrue.W)^fields{j,i}.option.gap; 
                    x.At(1:length(Mt(:)), i, j) = Mt(:);
                    Me = (fields{j,i}.netest.A+fields{j,i}.netest.W)^fields{j,i}.option.gap; 
                    x.Ae(1:length(Me(:)), i, j) = Me(:);
            end

            x.nrmseA(i,j) = sqrt(mean((Mt - Me).^2, 'all')) / abs(mean(Mt(:)));
            
            for s = 1:2
                M = fields{j,i}.(var{s}).A+fields{j,i}.(var{s}).W; 
                vv.(var{s}) = get_eigvec(M);
            end
            dd = diag(abs(vv.('nettrue')' * vv.('netest')));
            x.dotprod(1:length(dd), i, j) = dd;
            
            x.corrZ(i,j) = mean(diag(corr(fields{j,i}.nettrue.z', fields{j,i}.netest.z')));
            x.nrmseZ(i,j) = sqrt(mean((fields{j,i}.nettrue.z - fields{j,i}.netest.z).^2, 'all'))...
                / abs(mean(fields{j,i}.nettrue.z(:)));
        end
    end

%     for j = 1:m
%         for i = 1:16 
%             corrA(i, j) = corr(squeeze(x.At(i,:,j))',squeeze(x.Ae(i,:,j))');
%         end
%     end
%     
%     diagE = [1,6,11,16];
%     offdiagE = [1:16]; offdiagE(diagE) = [];
%     
%     x.Adiag = corrA(diagE, :);
%     x.Aoff = corrA(offdiagE, :);

    if option.input == true
        for i = 1:n
            for j = 1:m
                 x.C(i, j) =  corr(fields{j,i}.nettrue.C(:), ...
                     fields{j,i}.netest.C(:), 'rows', 'complete');
            
                for s = 1:2
                    [cc.(var{s}),~,~] = svd(ctrb(fields{j,i}.(var{s}).A+ ...
                    fields{j,i}.(var{s}).W, fields{j,i}.(var{s}).C));
                end
                
                dotprodC(:, i, j) = diag(abs(cc.('nettrue')' * cc.('netest')));
                
            end
        end
        x.dotprodC = dotprodC;
    end
    
end