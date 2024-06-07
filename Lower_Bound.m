% Define the center and hyperbolic radius of the circle
t = 1i;   % Center at t = i
D = 1.0;      % Hyperbolic radius

% Extract real and imaginary parts of the center
x_t = real(t);
y_t = imag(t);

% equation of a circle
function circleEq = getCircleEq(center, radius)
    % Define symbolic variable
    sym s
   
    % Construct the symbolic equation of the circle
    circleEq = dist_H(s,center) == radius;
end

% express a circle as a function f(x)
function circleFunc = getCircleFunc(center, radius)
    sym s
    circleFunc = dist_H(s,center)-radius;
end


% get the lower bound d_i based on the quasi-geodesic condition
% used global constant "lambda" and "epsilon" 
function lowerBd = getLowerBd(ti,tn)
    lowerBd = (1/lambda) * abs(tn-ti) - eps; % quasi-geodesic condition
end

% equation for the lower bound
function lowerBdEq = getLowerBdEq(geodesic,di)
    sym tr

    % calling functions from file "GeodesicSegment"
    lowerBdEq = dist_from_point(geodesic,tr) == di; % Eq that parameterizes the lower bound on a line
end

% equation for the extended s-neighborhood
% used global constant "lambda" and "step_size"
function s_Ngbh = getS_ngbh(geodesic)
    fat = (1 - 1/lambda) * step_size;
    s_Ngbh = getLowerBdEq(geodesic,di+fat); % fatten lower bound
end
