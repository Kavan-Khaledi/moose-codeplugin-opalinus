//+
Point(1) = {-0.015, 0, -0.03, 1.0};
//+
Point(2) = {0, -0.015, -0.03, 1.0};
//+
Point(3) = {0.015, 0, -0.03, 1.0};
//+
Point(4) = {0, 0.015, -0.03, 1.0};
//+
Point(5) = {0, 0, -0.03, 1.0};
//+

//+
Circle(1) = {2, 5, 3};
//+
Circle(2) = {3, 5, 4};
//+
Circle(3) = {4, 5, 1};
//+
Circle(4) = {1, 5, 2};
//+
Curve Loop(1) = {3, 4, 1, 2};
//+
Plane Surface(1) = {1};
//+
Point{5} In Surface{1};

//+
Extrude {0.00, 0.000, 0.06} {
  Surface{1}; 
}
//+
Point{7} In Surface{26};
//+
Physical Surface("top") = {26};
//+
Physical Surface("bottom") = {1};
//+
Physical Surface("confining") = {17, 21, 25, 13};
//+
Physical Volume("rock") = {1};
//+
Physical Point("fixb") = {7};
//+
Physical Point("fixz") = {5};
//+
//+
Field[1] = Box;
//+
Field[1].VIn = 0.005;
//+
Field[1].VOut = 1;
//+
Field[1].XMax = 0.1;
//+
Field[1].XMin = -0.1;
//+
Field[1].YMax = 0.1;
//+
Field[1].YMin = -0.1;
//+
Field[1].ZMax = 0.1;
//+
Field[1].ZMin = -0.1;
//+
Background Field = 1;
