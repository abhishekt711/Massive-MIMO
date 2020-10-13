function okay = checkHexagonal(points,radius)


%Extract distances and angle
angles = angle(points);
distances = abs(points);

%Symmetry allows us to rotate all angles to lie in the area 0, pi/3
angles_modulus = mod(angles,pi/3);

%Extract the Cartesian coordinates for the rotated points
x = distances .* cos(angles_modulus);
y = distances .* sin(angles_modulus);

%Check if the points are in the hexagon, in an area limited by three lines
okay = (x<radius) & ( y<radius*(sqrt(3)/2)) & (x < radius - y/sqrt(3));
