
function sum = sum_floats_to_MCF(a, b)
    float_sum = a + b;
    b_proj = float_sum - a;
    a_proj = float_sum - b_proj;
    a_roundoff = a - a_proj;
    b_roundoff = b - b_proj;
    sum = [float_sum, a_roundoff + b_roundoff];
end