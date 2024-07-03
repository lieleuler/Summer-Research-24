

eps_vals = linspace(1, 1000, 500);
lambda = 2;

y_vals = [];
for i=1:500
    R = get_R_value(log(1 + sqrt(2)), lambda, eps_vals(i));
    y_vals = [y_vals, R];
end

plot(eps_vals, y_vals)
