function FP = ssm_fixed_point(parameter)
    
    nl = size(parameter.A,1);
    FP = (eye(nl) - (parameter.A + parameter.W)) \ parameter.h;

end