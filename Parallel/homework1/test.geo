
lc = 1e-2;

Point(1) = {0, 0, 0, lc};
Point(2) = {.1, 0, 0, lc};
Point(3) = {.1, .1, 0, lc};
Point(4) = {0 , .1, 0 ,lc};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4}; 
Line(4) = {4,1};

Line Loop(1) = {1,2,3,4};

//Plane Surface(1) = {1};

Physical Point(1) = {1,2,3,4};
Physical Line(1) = {1,2,3,4};


out1[] = Extrude { 0 , .1 , 0 } {
	Line{1}; Layers{5};
	Recombine;
};

Physical Surface(1) = {out1[1]};

/*
   out[] = Extrude { 0 , 0, .1 } {
	Surface{out1[1]}; Layers{ 30 };
	Recombine;
	};

Physical Volume(1) = {out[1]};
*/
