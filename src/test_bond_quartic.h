/* Author: Ryan Houlihan
 * BondQuarticTest()
 *
 *
 */

#ifndef BOND_QUARTIC_TEST_H
#define BOND_QUARTIC_TEST_H

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
#include "bond_quartic.h"
#include "pair_lj_cut.h"

class BondQuarticTest{
protected:
  BondQuarticTest();
  virtual ~BondQuarticTest();
  virtual void SetUp();
  virtual void TearDown();
};  

#endif
