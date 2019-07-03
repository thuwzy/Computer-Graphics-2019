#pragma once
#include "vec.h"
#include "objects.h"
#include "KDTree.h"
#include "bezier.h"
Vec Cen(50,40.8,-860);
Sphere spheres[] = {//Scene: radius, position, emission, color, material 
 //Sphere(1e5, Vec( 1e5+1-10,40.8,81.6), Vec(),Vec(.75,.25,.25),DIFF),//Left 
 //Sphere(1e5, Vec(-1e5+99,40.8,81.6),Vec(),Vec(.25,.25,.75),DIFF),//Rght 
 //Sphere(1e5, Vec(50,40.8, 1e5),     Vec(),Vec(.75,.75,.75),DIFF),//Back 
 //Sphere(1e5, Vec(50,40.8,-1e5+220), Vec(),Vec(),           DIFF),//Frnt 
 //Sphere(1e5, Vec(50, 1e5, 81.6),    Vec(),Vec(.75,.75,.75),DIFF,"texture/wood.jpg"),//Botm 
 //Sphere(1e5, Vec(50,-1e5+81.6,81.6),Vec(),Vec(.75,.75,.75),DIFF),//Top 
 //Sphere(16.5,Vec(27,16.5,47),       Vec(),Vec(1,1,1)*.999, SPEC),//Mirr 
 //Sphere(16.5,Vec(27,16.5,37),       Vec(),Vec(1,1,1)*.999, SPEC),//Mirr 
 //Sphere(16.5,Vec(73,16.5,78),       Vec(),Vec(1,1,1)*.999, REFR),//Glas
 //Sphere(16.5,Vec(73,16.5,78),       Vec(),Vec(1,1,1)*.999, REFR),//Glas  
 //Sphere(600, Vec(50,681.6-.27,81.6),Vec(12,12,12),  Vec(), DIFF) //Lite 
  //Sphere(10000,Cen+Vec(0,0,-200), Vec(0.00063842, 0.02001478, 0.28923243)*6e-2*8, Vec(.7,.7,1)*.25,  DIFF), // sky

  //Sphere(100000, Vec(50, -100000, 0),  Vec(),Vec(.3,.3,.3),DIFF,"texture/wood.jpg"), // grnd
  //Sphere(110000, Vec(50, -110048.5, 0),  Vec(.9,.5,.05)*4,Vec(),DIFF),// horizon brightener
}; 
KDTree kdtrees[]={
	//KDTree("mesh/dragon/dragon2.txt",REFR, 350,Vec(65,0,25))
	//KDTree("mesh/cup.txt",REFR, 5,Vec(80,15,80))
	//KDTree("mesh/white_oak.txt",DIFF, 0.1,Vec(30,0,40))
	
	//KDTree("mesh/winter_street/model.txt",DIFF, 2000,Vec(35,-45,-150))
	KDTree("mesh/room/Enter a title.txt",DIFF, 50,Vec(35,-45,-150))
};

Bezier beziers[]={
	//Bezier("bezier/bezier.txt",DIFF, 0.21,Vec(10,0,91),"texture/vase.png")
	Bezier("bezier/bezier2.txt",REFR, 0.71,Vec(22,26.5,-10))
	//Bezier("bezier/bezier2.txt",REFR, 0.51,Vec(65,0,40))
};