
function sum = sum_MCF_and_float(MCF1, a, max_terms)
    % Sum
    N = length(MCF1);
    sum = zeros(1, N + 1);
    q = a;
    for i = 1:N
        curr_sum = sum_floats_to_MCF(q, MCF1(i));
        q = curr_sum(1);
        sum(i) = curr_sum(2);
    end
    sum(N + 1) = q;
    sum = unique(sum);

    % Normalize
    N = length(sum) - 1;
    if N >= max_terms
        s = sum(N + 1);
        t = zeros(1, N);
        for i = N:-1:1
            curr_sum = sum_floats_to_MCF(sum(i), s);
            s = curr_sum(1);
            t(i) = curr_sum(2);
        end
    
        k = 1;
        sum = zeros(1, N);
        for i = 1:N
            curr_sum = sum_floats_to_MCF(s, t(i));
            if curr_sum(2) ~= 0
                sum(k) = curr_sum(1);
                s = curr_sum(2);
                k = k + 1;
            end
        end
    end
end
