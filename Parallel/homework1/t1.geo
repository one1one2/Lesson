
lc = 1e-2;

Point(1) = {0, 0, 0, lc};
Point(2) = {1, 0, 0, lc};
Point(3) = {1, 1, 0, lc};
Point(4) = {0 , 1, 0 ,lc};
Point(5) = {0.5,0.5,0,lc};
Point(6) = {0.5,0.75,0,lc};
Point(7) = {0.5,0.25,0,lc};

Line(1) = {1,2};
Line(2) = {3,2};
Line(3) = {3,4}; 
Line(4) = {4,1};
Circle(5) = {6,5,7};
Circle(6) = {7,5,6};

Line Loop(5) = {4,1,-2,3};
Line Loop(6) = {5,6};
Plane Surface(6) = {5,6};

Physical Point(1) = {1,2,3,4};


Physical Line(1) = {1,2,3,4,5,6};

Physical Surface(1) = {6};
