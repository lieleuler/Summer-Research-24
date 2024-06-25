
function dist = dist_H(p1, p2)
    x1 = real(p1);
    x2 = real(p2);
    if x1 == x2
        dist = abs(log(abs(imag(p1)/imag(p2))));
    else
        y1 = imag(p1);
        y2 = imag(p2);
        
        c = (x1^2 + y1^2 - x2^2 - y2^2) / (2 * (x1 - x2));
        r = sqrt((x1 - c)^2 + y1^2);


        u = c - r;
        v = c + r;

        dist = abs(log(abs((p1-u) * (p2-v))) - log(abs((p2-u) * (p1-v))));
    end
end