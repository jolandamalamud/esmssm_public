function B = ssm_vector_transformation(A, dim)
    
    B = zeros(size(A,1), size(A,1)*size(A,2));

    for k = 1:size(A,1)
        B(k,(k-1)*size(A,1)+1:(k-1)*size(A,1)+size(A,1)) = A(:); 
    end
    
end