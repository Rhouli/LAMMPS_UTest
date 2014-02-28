#ifndef BOND_CLASS2_TEST_H
#define BOND_CLASS2_TEST_H

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
#include "CLASS2/bond_class2.cpp"

class BondClass2Test{
protected:
  BondClass2Test();
  virtual ~BondClass2Test();
  virtual void SetUp();
  virtual void TearDown();
};

#endif

