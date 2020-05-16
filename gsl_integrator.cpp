// #include "header_file.h"
// #include "header_file.h"
#include <iostream>
#include <stdio.h>
#include <gsl/gsl_sf_fermi_dirac.h>

// int main(){
//   print_hello_world();
//   std::cout << add(3,5) << std::endl;
//   return 0;
// }

int main (void) {
  double x = 5.0;
  double j = 3;
  double y = gsl_sf_fermi_dirac_2(x);
  printf ("F2(%g) = %.18e\n", x, y);

  double z = gsl_sf_fermi_dirac_int(j, x);
  printf ("F3(%g) = %.18e\n", x, z);
  return 0;
}
