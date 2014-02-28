
/* Skylar Rudolph
 * BondFENEExpandTest()
 * Header File for the Fene Bond Tests....
 *
 */

#ifndef BOND_FENE_EXPAND_TEST_H
#define BOND_FENE_EXPAND_TEST_H

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
#include "bond_fene_expand.h"

class BondFENEExpandTest{
protected:
  BondFENEExpandTest();
  virtual ~BondFENEExpandTest();
  virtual void SetUp();
  virtual void TearDown();
};  

#endif
