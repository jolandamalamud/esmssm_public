function [z_swapped,state_err] = ssm_swap_latents(Ez, Z)

    comb = perms(1:size(Z,1));
    for k = 1:size(comb,1)      
          Ezs = Ez(comb(k,:),:);
          state_err(k) = nanmean((Ezs(:) - Z(:)).^2);
    end
    
    z_swapped = Ez(comb(state_err == min(state_err),:), :);
    if ~isequal(z_swapped, Ez)
        disp('swapped')
    end

end