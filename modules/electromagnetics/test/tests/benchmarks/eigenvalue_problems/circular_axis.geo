// Gmsh project created on Tue Dec 17 08:52:16 2024


//+
Point(1) = {0, 0, 0, 0.1};
//+
Point(2) = {0, -1, 0, 0.1};
//+
Point(3) = {1, 0, 0, 0.1};
//+
Point(4) = {0, 1, 0, 0.1};
//+
Point(5) = {-1, 0, 0, 0.1};
//+
Circle(1) = {2, 1, 3};
//+
Circle(2) = {3, 1, 4};
//+
Circle(3) = {4, 1, 5};
//+
Circle(4) = {5, 1, 2};
//+
Line(5) = {2, 1};
//+
Line(6) = {1, 4};
//+
Physical Curve("wall", 7) = {2, 3, 4, 1};
//+
Curve Loop(1) = {5, 6, -2, -1};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {5, 6, 3, 4};
//+
Plane Surface(2) = {2};
//+
Physical Surface("port", 8) = {1, 2};
