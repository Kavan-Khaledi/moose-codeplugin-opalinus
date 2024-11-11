//+
Point(1) = {-0.5, -0.5, -0.5, 1.0};
//+
Point(2) = {0.5, -0.5, -0.5, 1.0};
//+
Point(3) = {0.5, 0.5, -0.5, 1.0};
//+
Point(4) = {-0.5, 0.5, -0.5, 1.0};

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
Plane Surface(1) = {1};
//+
Transfinite Curve {1, 2, 3, 4} = 2 Using Progression 1;
//+
Extrude {0, 0, 1} {
  Surface{1}; Layers {2};
}
//+
Physical Surface("right", 27) = {17};
//+
Physical Surface("left", 28) = {25};
//+
Physical Surface("top", 29) = {13};
//+
Physical Surface("bottom", 30) = {21};
//+
Physical Surface("front", 31) = {26};
//+
Physical Surface("back", 32) = {1};
//+
Physical Volume("rock", 33) = {1};
