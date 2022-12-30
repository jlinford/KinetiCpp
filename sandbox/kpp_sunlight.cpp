#include <cstdio>
#include <cmath>

static double TIME = 0;
static double SUN = 0;

void Update_SUN()
{
double SunRise, SunSet;
double Thour, Tlocal, Ttmp; 
const double PI = 3.14159265358979;  

  SunRise = 4.5;
  SunSet  = 19.5;
  Thour = TIME/3600.0;
  Tlocal = Thour - ((int)Thour/24)*24;

  if ( (Tlocal >= SunRise) && (Tlocal <= SunSet) ) {
    Ttmp = (2.0*Tlocal-SunRise-SunSet)/(SunSet-SunRise);
    if (Ttmp > 0) Ttmp =  Ttmp*Ttmp;
             else Ttmp = -Ttmp*Ttmp;
    SUN = ( 1.0 + cos(PI*Ttmp) )/2.0; 
  } else {
    SUN=0.0;
  }
}

int main(void)
{
    double TSTEP = 60;
    while (TIME < 24*3600) {
        Update_SUN();
        printf("%g,%g\n", TIME, SUN);
        TIME += TSTEP;
    }
    return 0;
}