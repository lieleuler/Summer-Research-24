% calculation for R in THM 1.7

del = delta;
lam = lambda;
eps = epsilon;

k_1 = lam * (lam + eps);
k_2 = (2*lam * (lam + eps) + 3) * (lam + eps);

a = (1/del) * log(2);
b = (-6*k_1 - 2) * exp(log(2)/del);
c = k_2 * exp(log(2)/del);

% use the W_-1 branch of the Lamber-W function
w = lambertw(-1, (a/b) * exp(a*c/b));

D_0 = c/b - w/a;
R = D_0 * (1 + k_1) + k_2 / 2;
