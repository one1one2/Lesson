Geometry.OldNewReg = 1;


h = 1;

lc = 5e-1;

Point(1) = {0, 0, 0, lc};
Point(2) = {h, 0, 0, lc};
Point(3) = {h, h, 0, lc};
Point(4) = {0, h, 0, lc};
Point(5) = {0, 0 ,h, lc};
Point(6) = {h, 0 ,h, lc};
Point(7) = {h, h ,h, lc};
Point(8) = {0, h ,h, lc};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4}; 
Line(4) = {4,1};
Line(5) = {5,6};
Line(6) = {6,7};
Line(7) = {7,8};
Line(8) = {8,5};
Line(9) = {1,5};
Line(10) = {2,6};
Line(11) = {3,7};
Line(12) = {4,8};

Line Loop(1) = {1,2,3,4};
Line Loop(2) = {5,6,7,8};
Line Loop(3) = {1,10,-5,-9};
Line Loop(4) = {3,12,-7,-11};
Line Loop(5) = {4,9,-8,-12};
Line Loop(6) = {2,11,-6,-10};

Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};
Plane Surface(4) = {4};
Plane Surface(5) = {5};
Plane Surface(6) = {6};

Surface Loop(1) = {1,2,3,4,5,6};

Volume(1) = {1};

Physical Point(1) = {1,2,3,4,5,6,7,8};


Physical Line(1) = {1,2,3,4,5,6,7,8,9,10,11,12};

Physical Surface(1) = {1,2,3,4,5,6};

Recombine Surface{1};

Physical Volume(1) = {1};

//Recombine Volume{1};


/*
Extrude { {0,1,0} , {-0.1,0,0.1} , -Pi/2 } {
Surface{28}; Layers{7}; Recombine;
}

DefineConstant[ angle = {90, Min 0, Max 120, Step 1,
Name "Parameters/Twisting angle"} ];

out[] = Extrude { {-2*h,0,0}, {1,0,0} , {0,0.15,0.25} , angle * Pi / 180 } {
Surface{50}; Layers{10}; Recombine;
};

Physical Volume(101) = {1, 2, out[1]};

Geometry.PointNumbers = 1;
Geometry.Color.Points = Orange;
General.Color.Text = White;
Mesh.Color.Points = {255,0,0};
*/
//Geometry.Color.Surfaces = Geometry.Color.Points;
