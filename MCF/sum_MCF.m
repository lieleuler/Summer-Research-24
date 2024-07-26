
function sum = sum_MCF(MCF1, MCF2, max_terms)
    if ~isscalar(MCF1)
        sum_terms = MCF1;
        for i = 1:length(MCF2)
            sum_terms = sum_MCF_and_float(sum_terms, MCF2(i), max_terms);
        end
    else
        sum_terms = MCF2;
        for i = 1:length(MCF1)
            sum_terms = sum_MCF_and_float(sum_terms, MCF1(i), max_terms);
        end
    end
    sum = MCF(sum_terms);
end