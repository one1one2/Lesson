//定义初始元素，两个点和一条线。
lc = 1e-2;
h = 0.1;
Point(1) = {0, 0, 0, lc};
Point(2) = {h, 0, 0, lc};

Line(1) = {1,2};
//生成一个平面，按照长方形剖分。
out1[] = Extrude { 0 , h , 0 } {
	Line{1}; Layers{5};
	Recombine;
};

//生成长方体
out2[] = Extrude { 0 , 0, h } {
	Surface{out1[1]}; Layers{ 6 };  
	//将已经生成出来带有长方形网格的平面out1[1]来进行extrude，这样
	//生成的三维网格才会随之产生长方体。
	//如果用三角形网格的Surface进行这步操作，会产生三棱柱网格。
	Recombine;
	};
//定义物理实体
Physical Point(1) = {1,2};
Physical Line(1) = {1};
Physical Surface(1) = {out1[1]};
Physical Volume(1) = {out2[1]};

