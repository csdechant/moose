// Gmsh project created on Mon Mar 17 14:05:44 2025
//+
Point(1) = {0, 0, 0, 0.1};
//+
Point(2) = {1, 0, 0, 0.1};
//+
Point(3) = {1, 1, 0, 0.1};
//+
Point(4) = {0, 1, 0, 0.1};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Curve Loop(1) = {1, 2, 3, 4};
//+
Surface(1) = {1};
//+
Transfinite Line {1, -3} = 10 Using Progression 1;
Transfinite Line {2, -4} = 10 Using Progression 1;
Transfinite Surface {1};
//+
Physical Curve("boundary", 5) = {1, 2, 3, 4};
//+
Physical Surface("surface", 9) = {1};
