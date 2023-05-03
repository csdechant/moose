Mesh.ElementOrder = 2;
Mesh.SecondOrderLinear = 0;

uniform_refine = .1e-2;
// Dimensions
//-----------------------------------------------------

offset = 4.4e-3;

Ceramic_Diameter = 175e-3;

Ceramic_Thickness = 9e-3;

Ceramic_Column_Diameter = 47.3e-3;

Ceramic_Column_Height = 35e-3;

Ceramic_Metal_Ring_Thickness = 16e-3;

Ceramic_Metal_Ring_Height = 10e-3;

Chamber_Diameter = 195e-3;

Chamber_Height = 240e-3;

Metal_Plate_Thickness = 8e-3;

Ressonator_Pin_Diameter = 3.6e-3;

Ressonator_Pin_Height = 22e-3;

Ressontor_Body_Height = 47e-3;

Ressontor_Body_Width = 52e-3;
//Points
//-----------------------------------------------------

// Axis of symmetry at bottom of chamber
Point(1) = {0, 0, 0, uniform_refine};

// Bottom edge of chamber
Point(2) = {Chamber_Diameter / 2, 0, 0, uniform_refine};

Point(3) = {Chamber_Diameter / 2, Chamber_Height-Ceramic_Metal_Ring_Height, 0, uniform_refine};

// Ceramic Metal Ring Lower Inner Edge
Point(4) = {(Chamber_Diameter / 2) - Ceramic_Metal_Ring_Thickness, Chamber_Height-Ceramic_Metal_Ring_Height, 0, uniform_refine};

// Ceramic Metal Ring Upper Inner Edge
Point(5) = {(Chamber_Diameter / 2) - Ceramic_Metal_Ring_Thickness, Chamber_Height, 0, uniform_refine};

// Ceramic Bottom Outer Edge
Point(6) = {Ceramic_Diameter / 2, Chamber_Height, 0, uniform_refine};

// Ceramic Upper Outer Edge
Point(7) = {Ceramic_Diameter / 2, Ceramic_Thickness + Chamber_Height, 0, uniform_refine};

// Ceramic Upper Inner Edge
Point(8) = {Ceramic_Column_Diameter / 2, Ceramic_Thickness + Chamber_Height, 0, uniform_refine};

// Cermaic Column Top
Point(9) = {Ceramic_Column_Diameter / 2, Ceramic_Thickness + Ceramic_Column_Height + Chamber_Height, 0, uniform_refine};

// Ceramic Column Inner Top
Point(10) = {Ressonator_Pin_Diameter / 2, Ceramic_Thickness + Ceramic_Column_Height + Chamber_Height, 0, uniform_refine};

//Ressonator Pin Upper Edge
Point(11) = {Ressonator_Pin_Diameter / 2, Ceramic_Thickness + Ceramic_Column_Height + Metal_Plate_Thickness + Chamber_Height, 0, uniform_refine};

//Ressontor Pin Upper Axis of Symetry
Point(12) = {0, Ceramic_Thickness + Ceramic_Column_Height + Metal_Plate_Thickness + Chamber_Height, 0, uniform_refine};

//Ressonator Pin Lower Axis of Symetry
Point(13) = {0, Ceramic_Thickness + Ceramic_Column_Height + Metal_Plate_Thickness + Chamber_Height - Ressonator_Pin_Height, 0, uniform_refine};

//Ressonator Pin Lower Edge
Point(14) = {Ressonator_Pin_Diameter / 2, Ceramic_Thickness + Ceramic_Column_Height + Metal_Plate_Thickness + Chamber_Height - Ressonator_Pin_Height, 0, uniform_refine};

// Axis of Symetry at lower edge of Ceramic
Point(15) = {0, Chamber_Height, 0, uniform_refine};


// Lines
//-----------------------------------------------------
// Chamber Axis of symmetry
Line(101) = {15, 1};

// Bottom of Chamber
Line(102) = {1, 2};

// Outer Chamber Edge
Line(103) = {2, 3, 4, 5};

// Ceramic touching metal
Line(104) = {5, 6, 7, 8, 9, 10};

// Pin touching metal
Line(105) = {10, 11, 12};

// Pin Axis of symmetry
Line(106) = {12, 13};

// Pin touching Ceramic
Line(107) = {13, 14, 10};

// Ceramic Axis of symmetry
Line(108) = {13, 15};

// Ceramic to plasma
Line(109) = {15, 5};

//Loops and Surfaces
//-----------------------------------------------------
// Plasma Domain
Line Loop(1) = {101, 102, 103, -109};
Plane Surface(1) = {1};

// Ceramic
Line Loop(2) = {109, 104, -107, 108};
Plane Surface(2) = {2};

// Ressonoator Pin
Line Loop(3) = {105, 106, 107};
Plane Surface(3) = {3};

// Blocks and Boundaries
//-----------------------------------------------------
// Plasma Domain
Physical Surface("Plasma") = {1};

// Ceramic Domain
Physical Surface("Ceramic") = {2};

// Ressonator Pin Domain
Physical Surface("Resonator_Pin") = {3};

// Physical axis of symmetry
Physical Line("axis") = {106, 108, 101};

// Physical metal walls
Physical Line("metal_wall") = {102, 103, 104, 105};
