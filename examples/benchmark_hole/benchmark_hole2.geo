// Gmsh project created on Tue Feb 11 11:38:34 2020
SetFactory("OpenCASCADE");
//+
Point(1) = {-15,0, -15, 1.0};
//+
Point(2) = {15, 0, -15, 1.0};
//+
Point(3) = {15, 0, 15, 1.0};
//+
Point(4) = {-15, 0, 15, 1.0};
//+
Point(5) = {-0.707, 0, -0.707, 1.0};
//+
Point(6) = {0.707, 0, -0.707, 1.0};
//+
Point(7) = {0.707, 0, 0.707, 1.0};
//+
Point(8) = {-0.707, 0, 0.707, 1.0};
//+
Point(9) = {0, -0, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};

//+
Circle(5) = {5, 9, 6};
//+
Circle(6) = {6, 9, 7};
//+
Circle(7) = {7, 9, 8};
//+
Circle(8) = {8, 9, 5};
//+
Line(9) = {1, 5};
//+
Line(10) = {2, 6};
//+
Line(11) = {3, 7};
//+
Line(12) = {4, 8};
//+
Curve Loop(1) = {1, 10, -5, -9};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {10, 6, -11, -2};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {11, 7, -12, -3};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {12, 8, -9, -4};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {5, 6, 7, 8};
//+
Plane Surface(5) = {5};
//+
Extrude {0, 1, 0} {
  Surface{5};Surface{2}; Surface{3}; Surface{4}; Surface{1}; Layers{1};
}
//+
Physical Surface("right", 33) = {13};
//+
Physical Surface("left", 34) = {19};
//+
Physical Surface("top", 35) = {16};
//+
Physical Surface("bottom", 36) = {21};
//+
Physical Surface("back", 37) = {17, 20, 22, 14, 10};
//+
Physical Surface("inner", 38) = {8, 9, 6, 7};
//+
Physical Surface("front", 39) = {4, 3, 2, 1};
//+
Physical Surface("face", 40) = {5};
//+
Physical Volume("rock", 41) = {2, 4, 5, 3};
//+
Physical Volume("tunnel", 42) = {1};
//+
Field[1] = Box;
//+
Field[1].VIn = 0.1;
//+
Field[1].VOut = 2;
//+
Field[1].XMax = 3;
//+
Field[1].XMin = -3;
//+
Field[1].YMax = 1.5;
//+
Field[1].ZMax = 3;
//+
Field[1].ZMin = -3;
//+
Background Field = 1;
