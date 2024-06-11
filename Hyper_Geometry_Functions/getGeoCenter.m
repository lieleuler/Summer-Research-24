% calculate the center of the great geodesic circle given two points on it
function geoCenter = getGeoCenter(s,t)
    x1 = real(s);
    y1 = imag(s);
    x2 = real(t);
    y2 = imag(t);
        
    geoCenter = (x1^2 + y1^2 - x2^2 - y2^2) / (2 * (x1 - x2));
end
