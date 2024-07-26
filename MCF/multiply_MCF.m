
function product = multiply_MCF(MCF1, MCF2, max_terms)
    N = length(MCF2);
    product = 0;
    for i = 1:N
        product = sum_MCF(MCF1*MCF2(i), product, max_terms);
    end
end