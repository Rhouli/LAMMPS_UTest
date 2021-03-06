/* Author: Ryan Houlihan
 * BondHarmonicTest()
 *
 *
 */

#ifndef BOND_MORSE_TEST_H
#define BOND_MORSE_TEST_H

#include "mpi.h"
#include <iostream>
#include <math.h>
#include "library.h"
#include "string.h"
#include "stdlib.h"
#include "lammps.h"
#include "atom.h"
#include "atom_vec.h"
#include "comm.h"
#include "input.h"
#include "integrate.h"
#include "domain.h"
#include "force.h"
#include "update.h"
#include "gtest/gtest.h"
#include "bond_morse.h"

class BondMorseTest{
protected:
  BondMorseTest();
  virtual ~BondMorseTest();
  virtual void SetUp();
  virtual void TearDown();
};  

#endif
