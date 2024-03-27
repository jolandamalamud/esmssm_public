function [inp, ii] = ssm_clean_inputs(fullinp)

    % clear inputs with no variance and less than 25%
    percent_inp = nansum(fullinp,2) / sum(~isnan(fullinp(1,:)));
    ii = 0 > percent_inp < 1;
    inp = fullinp(ii,:);


end