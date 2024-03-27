function new_inp = ssm_input_extension(inp)

    tt = find(~isnan(inp(1,:)));
    new_inp = zeros(size(inp));
    
    new_inp(:,1:tt(1)+15) = repmat(inp(:,tt(1)),1,16);
    for k = 2:length(tt)-1
        new_inp(:,tt(k)-15:tt(k)+15) = repmat(inp(:,tt(k)),1,31);
    end
    new_inp(:,tt(end)-15:end) = repmat(inp(:,tt(end)),1,16);
end