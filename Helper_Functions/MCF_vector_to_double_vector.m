
function result = MCF_vector_to_double_vector(vector)
    N = length(vector)
    result = zeros(1, N)
    for i = 1:N
        result(i) = double(vector(i))
    end
end