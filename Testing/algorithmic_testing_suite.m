

% =============================
% Intersections + Getting Range
% =============================

% ========
% TEST 1-1
% ========
z = -1.6 + 3.36i;
a1 = 1.6i;
a2 = 3.73i;
point_radius = 0.1;
segment_radius = 0.4;
ints = intersections_of_point_and_segment_ngbhs(z, a1, a2, point_radius, segment_radius, 1, 0, 0, 1);z = -1.6 + 3.36i;
angle_range = getRangeTn(z, ints, GeodesicSegment(a1, a2), point_radius, segment_radius, true);
expected_ints = [-1.495 + 3.697i; -1.3141 + 3.199i];
expected_angle_range = [5.9952, 10.5269];
"TEST 1-1-1: " + (round(abs(ints - expected_ints), 3))
"TEST 1-1-2: " + (round(abs(angle_range - expected_angle_range), 3))
 
% ========
% TEST 1-2
% ========
z = 0.76 + 5.55i;
a1 = 1.6i;
a2 = 3.73i;
point_radius = 0.1;
segment_radius = 0.4;
ints = intersections_of_point_and_segment_ngbhs(z, a1, a2, point_radius, segment_radius, 1, 0, 0, 1);
angle_range = getRangeTn(z, ints, GeodesicSegment(a1, a2), point_radius, segment_radius, true);
expected_ints = [0.205 + 5.551i; 1.078 + 5.121i];
expected_angle_range = [3.8092, 7.8026];
"TEST 1-2-1: " + (round(abs(ints - expected_ints), 3))
"TEST 1-2-2: " + (round(abs(angle_range - expected_angle_range), 3))

% ========
% TEST 1-3
% ========
z = -0.2196 + 0.048i;
a1 = 0.19i;
a2 = 0.48i;
point_radius = 1.3;
segment_radius = 2.2;
ints = intersections_of_point_and_segment_ngbhs(z, a1, a2, point_radius, segment_radius, 1, 0, 0, 1);
angle_range = getRangeTn(z, ints, GeodesicSegment(a1, a2), point_radius, segment_radius, true);
expected_ints = [-0.2961 + 0.06643i; -0.16219 + 0.03673i];
expected_angle_range = [0.74506, 5.1141];
"TEST 1-3-1: " + (round(abs(ints - expected_ints), 3))
"TEST 1-3-2: " + (round(abs(angle_range - expected_angle_range), 3))

% ========
% TEST 1-4
% ========
z = 0.69 + 0.09i;
a1 = 0.02i;
a2 = 1.3i;
point_radius = 2.5;
segment_radius = 0.35;
ints = intersections_of_point_and_segment_ngbhs(z, a1, a2, point_radius, segment_radius, 1, 0, 0, 1);
angle_range = getRangeTn(z, ints, GeodesicSegment(a1, a2), point_radius, segment_radius, true);
expected_ints = [0.3487 + 0.9762i; 0.1571 + 0.4399i];
expected_angle_range = [0.2016, 6.3410];
"TEST 1-4-1: " + (round(abs(ints - expected_ints), 3))
"TEST 1-4-2: " + (round(abs(angle_range - expected_angle_range), 3))

% ========
% TEST 1-5
% ========
z = -0.69 + 0.09i;
a1 = 0.02i;
a2 = 1.3i;
point_radius = 2.5;
segment_radius = 0.35;
ints = intersections_of_point_and_segment_ngbhs(z, a1, a2, point_radius, segment_radius, 1, 0, 0, 1);
angle_range = getRangeTn(z, ints, GeodesicSegment(a1, a2), point_radius, segment_radius, true);
expected_ints = [-0.3487 + 0.9762i; -0.1571 + 0.4399i];
expected_angle_range = [4*pi - 6.3410, 4*pi - 0.2016];
"TEST 1-5-1: " + (round(abs(ints - expected_ints), 3))
"TEST 1-5-2: " + (round(abs(angle_range - expected_angle_range), 3))

% ========
% TEST 1-6
% ========
z = 0.56 + 2.18i;
a1 = 3.9i;
a2 = 6.1i;
point_radius = 0.6;
segment_radius = 0.6;
ints = intersections_of_point_and_segment_ngbhs(z, a1, a2, point_radius, segment_radius, 1, 0, 0, 1);
angle_range = getRangeTn(z, ints, GeodesicSegment(a1, a2), point_radius, segment_radius, true);
expected_ints = [1.8827 + 3.0046i; -0.7918 + 2.2700i];
expected_angle_range = [1.2094 , 5.5196];
"TEST 1-6-1: " + (round(abs(ints - expected_ints), 3))
"TEST 1-6-2: " + (round(abs(angle_range - expected_angle_range), 3))

% ========
% TEST 1-7
% ========
z = -2.71 + 4.85i;
a1 = 2.2i;
a2 = 8.8i;
point_radius = 0.23;
segment_radius = 0.54;
ints = intersections_of_point_and_segment_ngbhs(z, a1, a2, point_radius, segment_radius, 1, 0, 0, 1);
angle_range = getRangeTn(z, ints, GeodesicSegment(a1, a2), point_radius, segment_radius, true);
expected_ints = [-3.3469 + 5.9067i; -2.2414 + 3.9557i];
expected_angle_range = [0.4833, 3.6774];
"TEST 1-7-1: " + (round(abs(ints - expected_ints), 3))
"TEST 1-7-2: " + (round(abs(angle_range - expected_angle_range), 3))

% ========
% TEST 1-8
% ========
z = -2.55 + 4.85i;
a1 = 2.2i;
a2 = 8.8i;
point_radius = 0.23;
segment_radius = 0.54;
ints = intersections_of_point_and_segment_ngbhs(z, a1, a2, point_radius, segment_radius, 1, 0, 0, 1);
angle_range = getRangeTn(z, ints, GeodesicSegment(a1, a2), point_radius, segment_radius, true);
expected_ints = [-3.2977 + 5.8199i; -2.2128 + 3.9052i];
expected_angle_range = [0.5868, 3.5229];
"TEST 1-8-1: " + (round(abs(ints - expected_ints), 3))
"TEST 1-8-2: " + (round(abs(angle_range - expected_angle_range), 3))

% ========
% TEST 1-9
% ========
z = -0.83 + 8.84i;
a1 = 2.2i;
a2 = 8.6i;
point_radius = 0.2;
segment_radius = 0.2;
ints = intersections_of_point_and_segment_ngbhs(z, a1, a2, point_radius, segment_radius, 1, 0, 0, 1);
angle_range = getRangeTn(z, ints, GeodesicSegment(a1, a2), point_radius, segment_radius, true);
expected_ints = [0.1596 + 10.4967i; -1.4821 + 7.3613i];
expected_angle_range = [5.7958, 8.9692];
"TEST 1-9-1: " + (round(abs(ints - expected_ints), 3))
"TEST 1-9-2: " + (round(abs(angle_range - expected_angle_range), 3))

% =========
% TEST 1-10
% =========
z = -0.62 + 8.74i;
a1 = 2.2i;
a2 = 8.6i;
point_radius = 0.2;
segment_radius = 0.2;
ints = intersections_of_point_and_segment_ngbhs(z, a1, a2, point_radius, segment_radius, 1, 0, 0, 1);
angle_range = getRangeTn(z, ints, GeodesicSegment(a1, a2), point_radius, segment_radius, true);
expected_ints = [0.1502 + 10.4975i; -1.4866 + 7.3839i];
expected_angle_range = [5.9102, 2*pi + 2.5192];
"TEST 1-10-1: " + (round(abs(ints - expected_ints), 3))
"TEST 1-10-2: " + (round(abs(angle_range - expected_angle_range), 3))

% =========
% TEST 1-11
% =========
z = -0.024815 + 0.32431i;
a1 = -3.7539e-17 + 0.36555i;
a2 = -3.7539e-17 + 0.33076i;
point_radius = 0.1;
segment_radius = 5.5511e-17;
ints = intersections_of_point_and_segment_ngbhs(z, a1, a2, point_radius, segment_radius, 1, 0, 0, 1)
angle_range = getRangeTn(z, ints, GeodesicSegment(a1, a2), point_radius, segment_radius, true);
expected_ints = [0.1502 + 10.4975i; -1.4866 + 7.3839i]; % 0.8951 + 0.0000i doubled
expected_angle_range = [5.9102, 2*pi + 2.5192];
"TEST 1-11-1: " + (round(abs(ints - expected_ints), 3))
"TEST 1-11-2: " + (round(abs(angle_range - expected_angle_range), 3))




