
function dist = dist_H(p1, p2)
    dist = 2*atanh(abs((p1 - p2)/(p1 - conj(p2))));
end