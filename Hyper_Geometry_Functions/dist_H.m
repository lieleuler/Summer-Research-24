
function dist = dist_H(p1, p2)
    dist = acosh(1 + ((real(p1) - real(p2))^2 + (imag(p1) - imag(p2))^2)/(2*imag(p1)*imag(p2)));
end

% p1 = 1i;
% p2 = 10^3 + 10^-10i;
% dist = 2*atanh(abs((p1 - p2)/(p1 - conj(p2))))
