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
Curve Loop(1) = {1, 2, 3, 4};
//+
Plane Surface(1) = {1};
//+
Physical Curve("wall", 5) = {2, 3, 4, 1};
//+
Physical Surface("port", 6) = {1};
