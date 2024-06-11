
% equation for the extended s-neighborhood
% used global constant "lambda" and "step_size"
function s_Ngbh = get_S_ngbh(geodesic, ti, tn, lambda, eps, step_size)
    fat = (1 - 1/lambda) * step_size;
    di = getLowerBd(ti, tn, lambda, eps);
    s_Ngbh = getLowerBdEq(geodesic,di+fat); % fatten lower bound
end


% express a circle as a function f(x)
function circleFunc = getCircleFunc(center, radius)
    sym s
    circleFunc = dist_H(s,center)-radius;
end


% get the lower bound d_i based on the quasi-geodesic condition
% used global constant "lambda" and "epsilon" 
function lowerBd = getLowerBd(ti, tn, lambda, eps)
    lowerBd = (1/lambda) * abs(tn-ti) - eps; % quasi-geodesic condition
end

% equation for the lower bound
function lowerBdEq = getLowerBdEq(geodesic,di)
    lowerBdEq = @(tr) geodesic.dist_from_point(tr) <= di; % Eq that parameterizes the lower bound on a line
end
