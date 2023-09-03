// Create geometry
Point(1) = {0, 0, 0, 0.04};
Point(2) = {2.2, 0, 0, 0.04};
Point(3) = {2.2, 0.41, 0, 0.04};
Point(4) = {0, 0.41, 0, 0.04};
Point(5) = {0.2, 0.2, 0};
Point(6) = {0.2 + 0.05, 0.2, 0, 0.005};
Point(7) = {0.2, 0.2 + 0.05, 0, 0.005};
Point(8) = {0.2 - 0.05, 0.2, 0, 0.005};
Point(9) = {0.2, 0.2 - 0.05, 0, 0.005};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Circle(5) = {6, 5, 7};
Circle(6) = {7, 5, 8};
Circle(7) = {8, 5, 9};
Circle(8) = {9, 5, 6};


Line Loop(6) = {1, 2, 3, 4};
Curve Loop(7) = {5, 6, 7, 8};
Plane Surface(8) = {6, 7};

Physical Curve(9) = {4};
Physical Curve(10) = {1};
Physical Curve(11) = {2};
Physical Curve(12) = {3};
Physical Curve(13) = {5, 6, 7, 8};
Physical Surface(14) = {8};

// Generate mesh
Mesh 2;

