
for k = 1:1000000
    geodesic = GeodesicSegment(rand(), rand());
    new_z = geodesic.get(2);
end
x = "Done!";