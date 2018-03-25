#ifndef __WENO_CPP_
#define __WENO_CPP_

#include <cmath>
#include "weno.h"
/**
 * @param v3 cell average u_{i-1}
 * @param v2
 * @param v1 cell average u_{i+1}
 *
 * @return  u point value at i+1/2 within cell i, or left limit at i+1/2
 */


double WENO3CellR(double v3,double v2,double v1){

     double s1,s2; 
     double sum; 
     double w1,w2; 
     double value; 
     double eps=1e-6; 
  
     s1=(v1-v2)*(v1-v2)+eps;  
     s2=(v2-v3)*(v2-v3)+eps; 
  
     s1=2.0/(3.0*s1*s1); 
     s2=1.0/(3.0*s2*s2); 
  
     sum=s1+s2;  
     w1=s1/sum; 
     w2=s2/sum; 
 
     value=w1*(0.5*v1+0.5*v2)+w2*(1.5*v2-0.5*v3); 
  
     return value; 
  
}  
 
 
 
/**  
 *  
 *  
 * @param v1 cell average u_{i-1}   
 * @param v2  
 * @param v3 cell average u_{i+1} 
 *  
 * @return  u point value at i-1/2 within cell i, or right limit at i-1/2
 */ 
  
double WENO3CellL(double v1,double v2,double v3){ 
  
     double s1,s2; 
     double sum; 
     double w1,w2; 
     double value; 
     double eps=1e-6; 
  
     s1=(v1-v2)*(v1-v2)+eps;  
     s2=(v2-v3)*(v2-v3)+eps; 
  
     s1=2.0/(3.0*s1*s1); 
     s2=1.0/(3.0*s2*s2); 
  
     sum=s1+s2;  
     w1=s1/sum; 
     w2=s2/sum; 
 
     value=w1*(0.5*v1+0.5*v2)+w2*(1.5*v2-0.5*v3); 
  
     return value; 
  
}  


/** 
 * 
 * 
 * @param v5 cell average u_{i-2}
 * @param v4 
 * @param v3  cell average u_i , cell i: [i-1/2,i+1/2]
 * @param v2 
 * @param v1 
 * 
 * @return  point value u_{i-1/2}
 */



double WENO5CellL(double v5,double v4,double v3,double v2,double v1){
  double s1,s2,s3;
  double sum;
  double tmp1,tmp2;
  double eps=1e-40;
  double w1,w2,w3;
     
  tmp1=v1-2*v2+v3;
  tmp2=v1-4*v2+3*v3;
  s1=13.0/12*tmp1*tmp1+0.25*tmp2*tmp2;
  s1=0.1/(s1*s1+eps);

  tmp1=v2-2*v3+v4;
  tmp2=v2-v4;
  s2=13.0/12*tmp1*tmp1+0.25*tmp2*tmp2;
  s2=0.6/(s2*s2+eps);

  tmp1=v3-2*v4+v5;
  tmp2=3*v3-4*v4+v5;
  s3=13.0/12*tmp1*tmp1+0.25*tmp2*tmp2;
  s3=0.3/(s3*s3+eps);
  

  sum=s1+s2+s3;
  w1=s1/sum;
  w2=s2/sum;
  w3=s3/sum;

  double v_WENO = w1*(v1/3-7*v2/6+11*v3/6)+w2*(-v2/6+5*v3/6+v4/3)+w3*(v3/3+5*v4/6-v5/6);

  return (v_WENO);
}




/** 
 * 
 * 
 * @param v1 cell average u_{i-2}
 * @param v2 
 * @param v3  cell average u_i , cell i: [i-1/2,i+1/2]
 * @param v4 
 * @param v5 
 * 
 * @return  point value u_{i+1/2}
 */
double WENO5CellR(double v1,double v2,double v3,double v4,double v5){
  double s1,s2,s3;
  double sum;
  double tmp1,tmp2;
  double eps=1e-40;
  double w1,w2,w3;

     
  tmp1=v1-2*v2+v3;
  tmp2=v1-4*v2+3*v3;
  s1=13.0/12*tmp1*tmp1+0.25*tmp2*tmp2;
  s1=0.1/(s1*s1+eps);

  tmp1=v2-2*v3+v4;
  tmp2=v2-v4;
  s2=13.0/12*tmp1*tmp1+0.25*tmp2*tmp2;
  s2=0.6/(s2*s2+eps);

  tmp1=v3-2*v4+v5;
  tmp2=3*v3-4*v4+v5;
  s3=13.0/12*tmp1*tmp1+0.25*tmp2*tmp2;
  s3=0.3/(s3*s3+eps);
     

  sum=s1+s2+s3;
  w1=s1/sum;
  w2=s2/sum;
  w3=s3/sum;     
     

  double v_WENO = w1*(v1/3-7*v2/6+11*v3/6)+w2*(-v2/6+5*v3/6+v4/3)+w3*(v3/3+5*v4/6-v5/6);

  return (v_WENO);

}









/** 
 * 
 * 
 * @param v7 cell average u_{i-3}
 * @param v6 
 * @param v5 
 * @param v4  cell average u_i , cell i: [i-1/2,i+1/2]
 * @param v3 
 * @param v2 
 * @param v1
 * 
 * @return  point value u_{i-1/2}
 */
double WENO7CellL(double v7,double v6,double v5,double v4,double v3,double v2,double v1){
  double s1,s2,s3,s4;
  double sum;
  double tmp1;
  double eps=1e-14;
  double w1,w2,w3,w4;
     
  double D0 = v2-v1;
  double D1 = v3-v2;
  double D2 = v4-v3;
  double D3 = v5-v4;
  double D4 = v6-v5;
  double D5 = v7-v6;

  s1=D0*(547.*D0-2788.*D1+1854.*D2)+D1*(3708.*D1-5188.*D2)
    +2107.*D2*D2;
  s2=D1*(267.*D1-1108.*D2+494.*D3)+D2*(1468.*D2-1428.*D3)
    +547.*D3*D3;
  s3=D2*(547.*D2-1428.*D3+494.*D4)+D3*(1468.*D3-1108.*D4)
    +267.*D4*D4;
  s4=D3*(2107.*D3-5188.*D4+1854.*D5)+D4*(3708.*D4-2788.*D5)
    +547.*D5*D5;

  s1=1.0/(s1*s1+ eps);

  s2 = 12.0/(s2*s2+ eps);

  s3 = 18.0/(s3*s3+ eps);


  s4 = 4.0/(s4*s4+ eps);
     
  sum=s1+s2+s3+s4;
  w1=s1/sum; 
  w2=s2/sum;
  w3=s3/sum;
  w4=s4/sum;

  tmp1 = w1*( -v1*3 + v2*13 - v3*23 + v4*25 )
    +  w2*( v2 - v3*5 + v4*13 + v5*3)
    +  w3*( -v3 + v4*7 + v5*7 - v6 )
    +  w4*( v4*3 + v5*13 - v6*5 + v7 );


  double v_WENO = tmp1/12.0;

  return (v_WENO);

}



/** 
 * 
 * 
 * @param v1 cell average u_{i-3}
 * @param v2 
 * @param v3 
 * @param v4  cell average u_i , cell i: [i-1/2,i+1/2]
 * @param v5 
 * @param v6 
 * @param v7
 * 
 * @return  point value u_{i+1/2}
 */

double WENO7CellR(double v1,double v2,double v3,double v4,double v5,double v6,double v7){
  double s1,s2,s3,s4;
  double sum;
  double tmp1;
  double eps=1e-14;
  double w1,w2,w3,w4;

  double D0 = v2-v1;
  double D1 = v3-v2;
  double D2 = v4-v3;
  double D3 = v5-v4;
  double D4 = v6-v5;
  double D5 = v7-v6;

  s1=D0*(547.*D0-2788.*D1+1854.*D2)+D1*(3708.*D1-5188.*D2)
    +2107.*D2*D2;
  s2=D1*(267.*D1-1108.*D2+494.*D3)+D2*(1468.*D2-1428.*D3)
    +547.*D3*D3;
  s3=D2*(547.*D2-1428.*D3+494.*D4)+D3*(1468.*D3-1108.*D4)
    +267.*D4*D4;
  s4=D3*(2107.*D3-5188.*D4+1854.*D5)+D4*(3708.*D4-2788.*D5)
    +547.*D5*D5;

  s1=1.0/(s1*s1+ eps);

  s2 = 12.0/(s2*s2+ eps);

  s3 = 18.0/(s3*s3+ eps);

  s4 = 4.0/(s4*s4+ eps);

  sum=s1+s2+s3+s4;
  w1=s1/sum; 
  w2=s2/sum;
  w3=s3/sum;
  w4=s4/sum;

  tmp1 = w1*( -v1*3 + v2*13 - v3*23 + v4*25 )
    +  w2*( v2 - v3*5 + v4*13 + v5*3)
    +  w3*( -v3 + v4*7 + v5*7 - v6 )
    +  w4*( v4*3 + v5*13 - v6*5 + v7 );

  double v_WENO = tmp1/12.0;

  return (v_WENO);


}










/** 
 * 
 * 
 * @param v1 cell average u_{i-4}
 * @param v2 
 * @param v3 
 * @param v4 
 * @param v5  cell average u_i , cell i: [i-1/2,i+1/2]
 * @param v6 
 * @param v7
 * @param v8
 * @param v9 
 * 
 * @return  point value u_{i+1/2}
 */

double WENO9CellR(double v1,double v2,double v3,double v4,double v5,double v6,double v7,double v8,double v9){
  double s1,s2,s3,s4,s5;
  double sum;
  double tmp1;
  double eps=1e-40;
  double w1,w2,w3,w4,w5;

  double D0 = v2-v1;
  double D1 = v3-v2;
  double D2 = v4-v3;
  double D3 = v5-v4;
  double D4 = v6-v5;
  double D5 = v7-v6;
  double D6 = v8-v7;
  double D7 = v9-v8;


  s1=D0*(22658.*D0-163185.*D1+201678.*D2-86329.*D3)+
    D1*(297120.*D1-745293.*D2+325158.*D3)+D2*(478980.*D2-
					      433665.*D3)+107918.*D3*D3;
  s2=D1*(6908.*D1-47055*D2+52158.*D3-18079.*D4)+D2*(84600.*
						    D2-196563.*D3+70218.*D4)+D3*(125130.*D3-94935.*D4)+
    22658.*D4*D4;
  s3=D2*(6908.*D2-37185.*D3+30738.*D4-8209.*D5)+D3*(60870.*
						    D3-109413.*D4+30738.*D5)+D4*(60870.*D4-37185.*D5)+
    6908.*D5*D5;
  s4=D3*(22658.*D3-94935.*D4+70218.*D5-18079.*D6)+D4*
    (125130.*D4-196563.*D5+52158.*D6)+D5*(84600.*D5-
					  47055.*D6)+6908.*D6*D6;
  s5=D4*(107918.*D4-433665.*D5+325158.*D6-86329.*D7)+D5*
    (478980.*D5-745293.*D6+201678.*D7)+D6*(297120.*D6-
					   163185.*D7)+22658.*D7*D7;
	
  s1 = fabs(s1) + eps;
  s2 = fabs(s2) + eps;
  s3 = fabs(s3) + eps;
  s4 = fabs(s4) + eps;
  s5 = fabs(s5) + eps;	


  s1 = 1.0/pow(s1,4.0);
  s2 = 20.0/pow(s2,4.0);
  s3 = 60.0/pow(s3,4.0);
  s4 = 40.0/pow(s4,4.0);
  s5 = 5.0/pow(s5,4.0);


  sum=s1+s2+s3+s4+s5;
  w1=s1/sum; 
  w2=s2/sum;
  w3=s3/sum;
  w4=s4/sum;
  w5=s5/sum;

  tmp1 = w1*( v1/5.0 - v2*21.0/20 + (v3+v5)*137.0/60 - v4*163.0/60 )
    +  w2*( -v2/20.0 +v3*17.0/60 - v4*43.0/60 + v5*77.0/60 + v6/5.0)
    +  w3*( v3/30.0 - v4*13.0/60 + v5*47.0/60 + v6*9.0/20 - v7/20.0 )
    +  w4*( -v4/20.0 + v5*9.0/20 + v6*47.0/60 - v7*13.0/60 + v8/30.0 )
    +  w5*( v5/5.0 + v6*77.0/60 - v7*43.0/60 + v8*17.0/60 - v9/20.0 );


  double v_WENO = tmp1;

  return (v_WENO);


}



/** 
 * 
 * 
 * @param v9 cell average u_{i-4}
 * @param v8 
 * @param v7 
 * @param v6 
 * @param v5  cell average u_i , cell i: [i-1/2,i+1/2]
 * @param v4 
 * @param v3
 * @param v2
 * @param v1 
 * 
 * @return  point value u_{i-1/2}
 */

double WENO9CellL(double v9,double v8,double v7,double v6,double v5,double v4,double v3,double v2,double v1){
  double s1,s2,s3,s4,s5;
  double sum;
  double tmp1;
  double eps=1e-40;
  double w1,w2,w3,w4,w5;

  double D0 = v2-v1;
  double D1 = v3-v2;
  double D2 = v4-v3;
  double D3 = v5-v4;
  double D4 = v6-v5;
  double D5 = v7-v6;
  double D6 = v8-v7;
  double D7 = v9-v8;


  s1=D0*(22658.*D0-163185.*D1+201678.*D2-86329.*D3)+
    D1*(297120.*D1-745293.*D2+325158.*D3)+D2*(478980.*D2-
					      433665.*D3)+107918.*D3*D3;
  s2=D1*(6908.*D1-47055*D2+52158.*D3-18079.*D4)+D2*(84600.*
						    D2-196563.*D3+70218.*D4)+D3*(125130.*D3-94935.*D4)+
    22658.*D4*D4;
  s3=D2*(6908.*D2-37185.*D3+30738.*D4-8209.*D5)+D3*(60870.*
						    D3-109413.*D4+30738.*D5)+D4*(60870.*D4-37185.*D5)+
    6908.*D5*D5;
  s4=D3*(22658.*D3-94935.*D4+70218.*D5-18079.*D6)+D4*
    (125130.*D4-196563.*D5+52158.*D6)+D5*(84600.*D5-
					  47055.*D6)+6908.*D6*D6;
  s5=D4*(107918.*D4-433665.*D5+325158.*D6-86329.*D7)+D5*
    (478980.*D5-745293.*D6+201678.*D7)+D6*(297120.*D6-
					   163185.*D7)+22658.*D7*D7;


  s1 = fabs(s1) + eps;
  s2 = fabs(s2) + eps;
  s3 = fabs(s3) + eps;
  s4 = fabs(s4) + eps;
  s5 = fabs(s5) + eps;	




  s1 = 1.0/pow(s1,4.0);
  s2 = 20.0/pow(s2,4.0);
  s3 = 60.0/pow(s3,4.0);
  s4 = 40.0/pow(s4,4.0);
  s5 = 5.0/pow(s5,4.0);
	

  sum=s1+s2+s3+s4+s5;
  w1=s1/sum; 
  w2=s2/sum;
  w3=s3/sum;
  w4=s4/sum;
  w5=s5/sum;

  tmp1 = w1*( v1/5.0 - v2*21.0/20 + (v3+v5)*137.0/60 - v4*163.0/60 )
    +  w2*( -v2/20.0 +v3*17.0/60 - v4*43.0/60 + v5*77.0/60 + v6/5.0)
    +  w3*( v3/30.0 - v4*13.0/60 + v5*47.0/60 + v6*9.0/20 - v7/20.0 )
    +  w4*( -v4/20.0 + v5*9.0/20 + v6*47.0/60 - v7*13.0/60 + v8/30.0 )
    +  w5*( v5/5.0 + v6*77.0/60 - v7*43.0/60 + v8*17.0/60 - v9/20.0 );


  double v_WENO = tmp1;

  return (v_WENO);

}

#endif

