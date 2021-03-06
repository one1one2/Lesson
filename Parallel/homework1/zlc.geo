lc = 0.5;
h = 0.4;
//外边界
Point(1) = {0,0,0,lc};
Point(2) = {9,0,0,lc};
Point(3) = {9,5,0,lc};
Point(4) = {0,5,0,lc};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Line Loop(5) = {1,2,3,4};

//Z
Point(101) = {1,1,0,lc};
Point(102) = {3-1.4*h,4-h,0,lc};
Point(103) = {1,4-h,0,lc};
Point(104) = {1,4,0,lc};
Point(105) = {3,4,0,lc};
Point(106) = {1+1.4*h,1+h,0,lc};
Point(107) = {3,1+h,0,lc};
Point(108) = {3,1,0,lc};

Line(101)= {101,102};
Line(102)= {102,103};
Line(103)= {103,104};
Line(104)= {104,105};
Line(105)= {105,106};
Line(106)= {106,107};
Line(107)= {107,108};
Line(108)= {108,101};

Line Loop(109) ={101,102,103,104,105,106,107,108};

//L

Point(201) = {4,1,0,lc};
Point(202) = {4,4,0,lc};
Point(203) = {4+h,4,0,lc};
Point(204) = {4+h,1+h,0,lc};
Point(205) = {6,1+h,0,lc};
Point(206) = {6,1,0,lc};

Line(201)={201,202};
Line(202)={202,203};
Line(203)={203,204};
Line(204)={204,205};
Line(205)={205,206};
Line(206)={206,201};

Line Loop(207) = {201,202,203,204,205,206};
//C

Point(301) = {8,1,0,lc};
Point(302) = {8,4,0,lc};
Point(303) = {8,4-h,0,lc};
Point(304) = {8,1+h,0,lc};
Point(305) = {8,2.5,0,lc};
Point(401) = {6.5,2.5,0,lc};

Circle(301) = {302,305,301};
Line(302) = {302,303};
Circle(303) = {303,305,304};
Line(304) = {304,301};

Line Loop(306) = {-301,302,303,304};

Plane Surface(6) = {5,109,207,306};
