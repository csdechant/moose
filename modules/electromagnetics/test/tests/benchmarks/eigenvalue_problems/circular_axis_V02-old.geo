// Gmsh project created on Tue Dec 17 08:52:16 2024


//+
Point(1) = {0, 1, 0, 0.1};
//+
Point(2) = {0, 0, 0, 0.1};
//+
Point(3) = {1, 1, 0, 0.1};
//+
Point(4) = {0, 2, 0, 0.1};
//+
Point(5) = {-1, 1, 0, 0.1};
//+
Circle(1) = {2, 1, 3};
//+
Circle(2) = {3, 1, 4};
//+
Circle(3) = {4, 1, 5};
//+
Circle(4) = {5, 1, 2};

//+
Line(5) = {5, 1};
//+
Line(6) = {1, 4};
//+
Line(7) = {1, 3};
//+
Curve Loop(1) = {7, 2, -6};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {6, 3, 5};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {5, 7, -1, -4};
//+
Plane Surface(3) = {3};
//+
Physical Curve("wall", 8) = {2, 3, 4, 1};
//+
Physical Surface("port", 9) = {1, 2, 3};
