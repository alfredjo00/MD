#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "t1.h"
#include "t2.h"
#include "t3.h"


int run(int argc, char *argv[])
{
  write_Al_E_pots();
  task2();

  // T (K), deviation, file name label
  temp_pressure_equilibrium_task(700.0 + 273.15, 0.35, (char *) "liq");
  temp_pressure_equilibrium_task(500.0 + 273.15, 0.065, (char *) "solid");
  
	return 0.;
}
