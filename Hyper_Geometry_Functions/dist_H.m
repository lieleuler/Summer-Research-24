
function dist = dist_H(p1, p2)
    dist = acosh(1 + (double(real(p1) - real(p2))^2 + double(imag(p1) - imag(p2))^2)/(2*imag(double(p1))*imag(double(p2))));
end

% dist = acosh(1 + (MCF_to_float(sum_MCF(real(p1), - real(p2)))^2 + MCF_to_float(sum_MCF(imag(p1), -imag(p2)))^2)...
%         /...
%         (2*MCF_to_float(imag(p1))*MCF_to_float(imag(p2))));