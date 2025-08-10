#include <iostream>

#include "Reaction.h"
#include "constants.h"

using namespace std;

int main() {

  Reaction rxn(1, 6,12, 1,2, 1,1);

  rxn.Run_Kinematics(16.0, 20., 0., 0., 0., 0., 0.);

  rxn.Print();

  return 0;

}
