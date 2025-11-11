// Gmsh project created on Tue Oct 10 20:26:15 2023
SetFactory("OpenCASCADE");
//+
D = DefineNumber[ 0.16, Name "Parameters/D" ];
//+
Point(1) = {0, D/2, 0, 1.0};
//+
Point(2) = {D/2, 0, 0, 1.0};
//+
Point(3) = {0, 0, 0, 1.0};
//+
Circle(1) = {2, 3, 1};
//+
Line(2) = {1, 3};
//+
Line(3) = {3, 2};
//+
Curve Loop(1) = {3, 1, 2};
//+
Plane Surface(1) = {1};
//+
Physical Curve("Outer", 4) = {1};
