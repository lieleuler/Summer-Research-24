% express a circle as a function f(x)
function circleFunc = getCircleFunc(center, radius)
    syms s
    circleFunc = dist_H(s,center)-radius;
end