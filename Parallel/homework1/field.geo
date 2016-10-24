// Let’s create a simple rectangular geometry
lc = .15;
Point(1) = {0.0,0.0,0,lc}; Point(2) = {1,0.0,0,lc};
Point(3) = {1,1,0,lc};
 Point(4) = {0,1,0,lc};
Point(5) = {0.2,.5,0,lc};
Line(1) = {1,2}; Line(2) = {2,3}; Line(3) = {3,4}; Line(4) = {4,1};
Line Loop(5) = {1,2,3,4}; Plane Surface(6) = {5};


Field[1] = Attractor;
Field[1].NodesList = {5};
Field[1].NNodesByEdge = 100;
Field[1].EdgesList = {2};


Field[2] = Threshold;
Field[2].IField = 1;
Field[2].LcMin = lc / 10;
Field[2].LcMax = lc;
Field[2].DistMin = 0.15;
Field[2].DistMax = 0.5;


Field[3] = MathEval;
Field[3].F = "abs(4*(x-0.5)*(x-0.5)-y)*0.1+0.01";



Field[4] = Attractor;
Field[4].NodesList = {1};
Field[5] = MathEval;
Field[5].F = Sprintf("F4^3 + %g", lc / 100);
Field[6] = Box;
Field[6].VIn = lc / 15;
Field[6].VOut = lc;
Field[6].XMin = 0.3;
Field[6].XMax = 0.6;
Field[6].YMin = 0.3;
Field[6].YMax = 0.6;  
Field[7] = Min;
Field[7].FieldsList = {2, 3, 5, 6};  
Background Field = 2;  

