function V = ssm_vector_transformation(v, dim)
    
    V = zeros(dim, dim*length(v));

    for k = 1:dim
        V(k,(k-1)*length(v)+1:(k-1)*length(v)+length(v)) = v'; 
    end
    
end