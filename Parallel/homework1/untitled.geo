lc=0.1;

Point(1) = {0, 0, 0, lc};
Point(2) = {0, 0, 0, lc};
Point(3) = {0, 1, 0, lc};
Point(4) = {1, 1, 0, lc};
Point(5) = {1, 0, 0, lc};
Line(1) = {4, 5};
Line(2) = {1, 3};
Line(3) = {3, 4};
Line(4) = {5, 1};
Line Loop(5) = {2, 3, 1, 4};
Plane Surface(6) = {5};
Recombine Surface {6};
Recombine Surface {6};
