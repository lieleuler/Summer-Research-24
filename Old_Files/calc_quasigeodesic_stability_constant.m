
function R = calc_quasigeodesic_stability_constant(lambda, epsilon, delta)
    %k1 = lambda * (lambda + epsilon);
    %k2 = (2*lambda*(lambda + epsilon) + 3)*(lambda + epsilon);

    k1 = lambda * epsilon;
    k2 = (lambda * epsilon^2 + 3*epsilon) / 2;

    a = log(2)/delta;
    b = -exp(a) * (6*k1 + 2);
    c = exp(a)*k2;

    w = a/b * exp(a*c/b);
    D_0 = c/b - lambertw(-1, w)/a;

    R_prime = D_0*(1 + k1) + k2/2;
    R = R_prime + lambda + epsilon;
    