//SetFactory("OpenCASCADE");
// 0 for triangles, 1 for quads
tri_or_quad = 0;

Point(1) = {0, 0, 0, 1.0};
Point(2) = {4, 0, 0, 1.0};
Point(3) = {4, 1, 0, 1.0};
Point(4) = {0, 1, 0, 1.0};


Point(5) = {0.5+0.15, 0.5-0.15, 0, 0.5};
Point(6) = {0.5+0.15, 0.5+0.15, 0, 0.5};
Point(7) = {0.5-0.15, 0.5+0.15, 0, 0.5};
Point(8) = {0.5-0.15, 0.5-0.15, 0, 0.5};

Point(9) = {0.5+0.125, 0.5-0.125, 0, 0.5};
Point(10) = {0.5+0.125, 0.5+0.125, 0, 0.5};
Point(11) = {0.5-0.125, 0.5+0.125, 0, 0.5};
Point(12) = {0.5-0.125, 0.5-0.125, 0, 0.5};

Point(13) = {0.5+0.2, 0.5-0.2, 0, 0.5};
Point(14) = {0.5+0.2, 0.5+0.2, 0, 0.5};
Point(15) = {0.5-0.2, 0.5+0.2, 0, 0.5};
Point(16) = {0.5-0.2, 0.5-0.2, 0, 0.5};

Characteristic Length {1} = 0.010;
Characteristic Length {2} = 0.010;
Characteristic Length {3} = 0.010;
Characteristic Length {4} = 0.010;
Characteristic Length {5} = 0.001;
Characteristic Length {6} = 0.001;
Characteristic Length {7} = 0.001;
Characteristic Length {8} = 0.001;
Characteristic Length {9} = 0.014;
Characteristic Length {10} = 0.014;
Characteristic Length {11} = 0.014;
Characteristic Length {12} = 0.014;
Characteristic Length {13} = 0.010;
Characteristic Length {14} = 0.010;
Characteristic Length {15} = 0.010;
Characteristic Length {16} = 0.010;

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};

Line(9) = {9, 10};
Line(10) = {10, 11};
Line(11) = {11, 12};
Line(12) = {12, 9};

Line(13) = {13, 14};
Line(14) = {14, 15};
Line(15) = {15, 16};
Line(16) = {16, 13};


Curve Loop(1) = {3, 4, 1, 2};
Plane Surface(1) = {1};
Point{5} In Surface{1};
Point{6} In Surface{1};
Point{7} In Surface{1};
Point{8} In Surface{1};
Point{9} In Surface{1};
Point{10} In Surface{1};
Point{11} In Surface{1};
Point{12} In Surface{1};
Point{13} In Surface{1};
Point{14} In Surface{1};
Point{15} In Surface{1};
Point{16} In Surface{1};
Line{5} In Surface{1};
Line{6} In Surface{1};
Line{7} In Surface{1};
Line{8} In Surface{1};
Line{9} In Surface{1};
Line{10} In Surface{1};
Line{11} In Surface{1};
Line{12} In Surface{1};
Line{13} In Surface{1};
Line{14} In Surface{1};
Line{15} In Surface{1};
Line{16} In Surface{1};


If (tri_or_quad == 1)
   Recombine Surface {1};
EndIf
Physical Curve(1) = {1};
Physical Curve(2) = {2};
Physical Curve(3) = {3};
Physical Curve(4) = {4};
Physical Surface(1) = {1};


//Mesh.MshFileVersion = 2.2;
Mesh 2;


//+



