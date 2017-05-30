#ifndef __config_h_
#define __config_h_


#include <iostream>
#include <cmath>
#include <ostream>
#include <fstream>
#include <ctime>
#include <vector>
#include <thread>
#include <functional>
#include <omp.h>
#include <cassert>
#include "Vector.h"


#define gamma 1.4
#define CFL 0.3

enum{ DoubleMachReflection, ForwardStep};

//#define Problem ForwardStep
#define Problem DoubleMachReflection

#define Reconstruct

#define RK


#define limiter minmod



#endif
