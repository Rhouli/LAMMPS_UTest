/* Ryan Houlihan
 * BondHarmonicTest()
 *
 * The following code tests bond_style bond_harmonic in various ways:
 *
 * Coeff_correctNarg: checks bond->coeff() with correct input parameters
 * 
 * Coeff_highNarg: checks bond->coeff() for right exit status on incorrect input. Narg > expected
 * Coeff_NegativeNarg: checks bond->coeff() for right exit status on incorrect input. Narg = -expected
 * Coeff_lowNarg: checks bond->coeff() for right exit status on incorrect input. Narg < expected
 *
 * Tests the Force, Energy, and Pressure of each of the following simulations 
 * 
 * Bond1Run0_1: Test case with 2 atoms & 1 bond & coeff_argv{1, k = 132.0, r_0 = 2.0}
 * Bond1Run0_2: Test case with 2 atoms & 1 bond & coeff_argv{1, k = 500.0, r_0 = 10.0}
 * 
 * Bond2Run0_1: Test case with 4 atoms & 2 bonds & coeff_argv{1, k = 100.0, r_0 = 1.0}
 * Bond2Run0_2: Test case with 4 atoms & 2 bond & coeff_argv{1, k = 50.0, r_0 = 5.0}
 * 
 * Bond3Run0_1: Test case with 4 atoms & 3 bond & coeff_argv{1, k = 100.0, r_0 = 1.0}
 * Bond3Run0_2: Test case with 4 atoms & 3 bond & coeff_argv{1, k = 50.0, r_0 = 5.0}
 *
 */

#include "test_bond_nonlinear.h"

using namespace std;

using namespace LAMMPS_NS;

namespace {

  // The fixture for testing class Foo.
  class BondNonlinearTest : public ::testing::Test {   
  protected:
    BondNonlinearTest() {
      char *argv[] = {"bond_nonlinear", "-screen", "none", "-log", "none", NULL};
      //char *argv[] = {"bond_morse", "-echo", "screen", "-log", "none", NULL};
      int narg = 5;
      
      //Create a LAMMPS Object
      lmp = new LAMMPS(narg,argv,MPI_COMM_WORLD);
      // Set atom_style
      lmp->input->one("atom_style bond");
      // Modify atom_style parameters
      lmp->input->one("atom_modify map array sort 10000 3.0");
    }

    virtual ~BondNonlinearTest(){
      lmp->destroy();      
    }

    virtual void SetUp() {
      //Initialize simulation domain
      lmp->domain->triclinic = 0;
      lmp->domain->boxlo[0] = -5.0;
      lmp->domain->boxhi[0] =  5.0;
      lmp->domain->boxlo[1] = -5.0;
      lmp->domain->boxhi[1] =  5.0;
      lmp->domain->boxlo[2] = -5.0;
      lmp->domain->boxhi[2] =  5.0;
      lmp->domain->box_exist = 1;
  
      lmp->atom->ntypes = 1;
      lmp->atom->bond_per_atom = 1;
      lmp->atom->angle_per_atom = 0;
      lmp->atom->dihedral_per_atom = 0;
      lmp->atom->improper_per_atom = 0;
      lmp->atom->nbonds = 1;
      lmp->atom->nbondtypes = 1;
      lmp->atom->nangles = 0;
      lmp->atom->ndihedrals = 0;
      lmp->atom->nimpropers = 0;

      // set up two atoms
      lmp->atom->allocate_type_arrays();
      lmp->atom->avec->grow(2);

      lmp->domain->print_box("Created ");
      lmp->domain->set_initial_box();
      lmp->domain->set_global_box();
      lmp->comm->set_procs();
      lmp->domain->set_local_box();
    }

    virtual void TearDown() {
      lmp->input->one("clear"); 
    }
    
    // Objects declared here can be used by all tests in the test case for Foo.
    LAMMPS* lmp;  
  };

  // Check for case with correct value of bcoeff_narg
 TEST_F(BondNonlinearTest, Coeff_CorrectNarg){
    lmp->force->bond = new BondNonlinear(lmp);
    char *bcoeff_argv[] = {"1", "100.0", "2.5", "5.0", NULL};
    int bcoeff_narg; 

    // Check correct input
    bcoeff_narg = 4;
    lmp->force->bond->coeff(bcoeff_narg, bcoeff_argv);
    
    int i;
    for(i = 1; i <=  lmp->atom->nbondtypes; i++){
      EXPECT_DOUBLE_EQ(2.5, lmp->force->bond->equilibrium_distance(i));
      EXPECT_DOUBLE_EQ(1, lmp->force->bond->setflag[i]);
    }
 }

 //////////////////////////////////////////////////////////////////////
 // Multiple check's to make sure bond exits out with incorrect input//
 //////////////////////////////////////////////////////////////////////

  // Check for case with to high a value of bcoeff_narg
 TEST_F(BondNonlinearTest, Coeff_highNarg){
    lmp->force->bond = new BondNonlinear(lmp);
    char *bcoeff_argv[] = {"1", "100.0", "2.5", "5.0", NULL};
    int bcoeff_narg; 

    // Check incorrect input
    bcoeff_narg = 5;
    
    // Check that exit(1) was thrown
    ASSERT_EXIT(lmp->force->bond->coeff(bcoeff_narg, bcoeff_argv), ::testing::ExitedWithCode(1), "");  
 }

  // Check for case with to low a value of bcoeff_narg
 TEST_F(BondNonlinearTest, Coeff_lowNarg){
    lmp->force->bond = new BondNonlinear(lmp);
    char *bcoeff_argv[] = {"1", "100.0", "2.5", "5.0", NULL};
    int bcoeff_narg; 

    // Check incorrect input
    bcoeff_narg = 1;
    
    // Check that exit(1) was thrown
    ASSERT_EXIT(lmp->force->bond->coeff(bcoeff_narg, bcoeff_argv), ::testing::ExitedWithCode(1), "");  
 }
  // Check for case with for a negative value of bcoeff_narg
  TEST_F(BondNonlinearTest,Coeff_NegativeNarg){
    lmp->force->bond = new BondNonlinear(lmp);
    char *bcoeff_argv[] = {"1", "100.0", "2.5", "5.0", NULL};
    int bcoeff_narg; 

    // Check incorrect input
    bcoeff_narg = -4;
    
    // Check that exit(1) was thrown
    ASSERT_EXIT(lmp->force->bond->coeff(bcoeff_narg, bcoeff_argv), ::testing::ExitedWithCode(1), "");  
  }

 //////////////////////////////////////////////////////////////////////
 // Multiple check's to make sure bond compute works properly        //
 //////////////////////////////////////////////////////////////////////

  // Test basic case with 2 atoms & 1 bond & coeff_argv{1, epsilon = 132.0, r0 = 3.0, lamda = 2.0}
  TEST_F(BondNonlinearTest, Bond1Run0_1) {
    // Variables for test
    double epsilon = 132.0;
    double r0 = 3.0;
    double lamda = 2.0;
    double atomNum = 2; 
    double delx = 0.0-0.0;
    double dely = 0.0-0.0; 
    double delz = -1.0 - 1.0; 
    double rsq = delx*delx + dely*dely + delz*delz;
    double r = sqrt(rsq);
    double dr = r - r0;
    double drsq = dr*dr;
    double lamdasq = lamda*lamda;
    double denom = lamdasq - drsq;
    double denomsq = denom*denom;
    double fbond = -epsilon/r * 2.0*dr*lamdasq/denomsq;
    double ebond = epsilon * drsq / denom;

    //*********************************************  
    // Create Bonds and Atoms

    char *lines;
    lines = strdup("1 1 1 0.0 0.0 -1.0\n 2 1 1 0.0 0.0 1.0\n\n");
    // data_atoms(int n, char*) : unpack n lines from Atom section of data file
    lmp->atom->data_atoms(atomNum,lines);
    free(lines);
    lmp->atom->nlocal = atomNum;
    lmp->atom->natoms = atomNum;
  
    if (lmp->atom->map_style) {
      lmp->atom->map_init();
      lmp->atom->map_set();
    }

    // Create bond
    lines = strdup("1 1 1 2\n");
    lmp->atom->data_bonds(1,lines);
    free(lines);

    lmp->input->one("mass 1 1.0");

    // End create bonds and atoms
    //*********************************************  

    lmp->force->bond = new BondNonlinear(lmp);
    char *bcoeff_argv[] = {"1", "132.0", "3.0", "2.0", NULL};
    int bcoeff_narg; 
    
    // Check for correct input
    bcoeff_narg = 4;
    lmp->force->bond->coeff(bcoeff_narg, bcoeff_argv);

    int i;
    for(i = 1; i <=  lmp->atom->nbondtypes; i++){
      EXPECT_DOUBLE_EQ(r0, lmp->force->bond->equilibrium_distance(i));
    }

    // Tests that BondHarmonic::energy works properly 
    lmp->init();
    lmp->update->whichflag=1;
    lmp->update->firststep=0;
    lmp->update->laststep=1;
    lmp->input->one("run 0");

    // Energy for system
    EXPECT_DOUBLE_EQ(2*(ebond*0.5), lmp->force->bond->energy);

    // Forces for atom 1
    EXPECT_DOUBLE_EQ(fbond*delx, lmp->atom[0].f[0][0]); 
    EXPECT_DOUBLE_EQ(fbond*dely, lmp->atom[0].f[0][1]);
    EXPECT_DOUBLE_EQ(fbond*delz, lmp->atom[0].f[0][2]); 
    // Forces for atom 2
    EXPECT_DOUBLE_EQ(-1.0*fbond*delx, lmp->atom[0].f[1][0]); 
    EXPECT_DOUBLE_EQ(-1.0*fbond*dely, lmp->atom[0].f[1][1]);
    EXPECT_DOUBLE_EQ(-1.0*fbond*delz, lmp->atom[0].f[1][2]); 

    // Total computed pressure
    double expecVir[6] = {delx*delx*fbond,
			  dely*dely*fbond,
			  delz*delz*fbond,
			  delx*dely*fbond,
			  delx*delz*fbond,
                          dely*delz*fbond }; 
    for( i = 0; i < 6; i++)
      EXPECT_DOUBLE_EQ(expecVir[i] , lmp->force->bond->virial[i]); 
  }    

  // Test basic case with 2 atoms & 1 bond & coeff_argv{1, d0 = 500.0, alpha = 5.0, r_0 = 10.0}
  TEST_F(BondNonlinearTest, Bond1Run0_2) {
     // Variables for test
    double epsilon = 500.0;
    double r0 = 5.0;
    double lamda = 10.0;
    double atomNum = 2; 
    double delx = 0.0-0.0;
    double dely = 0.0-0.0; 
    double delz = -1.0 - 1.0; 
    double rsq = delx*delx + dely*dely + delz*delz;
    double r = sqrt(rsq);
    double dr = r - r0;
    double drsq = dr*dr;
    double lamdasq = lamda*lamda;
    double denom = lamdasq - drsq;
    double denomsq = denom*denom;
    double fbond = -epsilon/r * 2.0*dr*lamdasq/denomsq;
    double ebond = epsilon * drsq / denom;  

    //*********************************************  
    // Create Bonds and Atoms

    char *lines;
    lines = strdup("1 1 1 0.0 0.0 -1.0\n 2 1 1 0.0 0.0 1.0\n\n");
    // data_atoms(int n, char*) : unpack n lines from Atom section of data file
    lmp->atom->data_atoms(atomNum,lines);
    free(lines);
    lmp->atom->nlocal = atomNum;
    lmp->atom->natoms = atomNum;
  
    if (lmp->atom->map_style) {
      lmp->atom->map_init();
      lmp->atom->map_set();
    }

    // Create bond
    lines = strdup("1 1 1 2\n");
    lmp->atom->data_bonds(1,lines);
    free(lines);

    lmp->input->one("mass 1 1.0");

    // End create bonds and atoms
    //*********************************************  

    lmp->force->bond = new BondNonlinear(lmp);
    char *bcoeff_argv[] = {"1", "500.0", "5.0", "10.0", NULL};
    int bcoeff_narg; 
    
    // Check for correct input
    bcoeff_narg = 4;
    lmp->force->bond->coeff(bcoeff_narg, bcoeff_argv);

    int i;
    for(i = 1; i <=  lmp->atom->nbondtypes; i++){
      EXPECT_DOUBLE_EQ(r0, lmp->force->bond->equilibrium_distance(i));
    }

    // Tests that BondHarmonic::energy works properly 
    lmp->init();
    lmp->update->whichflag=1;
    lmp->update->firststep=0;
    lmp->update->laststep=1;
    lmp->input->one("run 0");

    // Energy for system
    EXPECT_DOUBLE_EQ(2*(ebond*0.5), lmp->force->bond->energy);

    // Forces for atom 1
    EXPECT_DOUBLE_EQ(fbond*delx, lmp->atom[0].f[0][0]); 
    EXPECT_DOUBLE_EQ(fbond*dely, lmp->atom[0].f[0][1]);
    EXPECT_DOUBLE_EQ(fbond*delz, lmp->atom[0].f[0][2]); 
    // Forces for atom 2
    EXPECT_DOUBLE_EQ(-1.0*fbond*delx, lmp->atom[0].f[1][0]); 
    EXPECT_DOUBLE_EQ(-1.0*fbond*dely, lmp->atom[0].f[1][1]);
    EXPECT_DOUBLE_EQ(-1.0*fbond*delz, lmp->atom[0].f[1][2]); 

    // Total computed pressure
    double expecVir[6] = {delx*delx*fbond,
			  dely*dely*fbond,
			  delz*delz*fbond,
			  delx*dely*fbond,
			  delx*delz*fbond,
                          dely*delz*fbond }; 
    for( i = 0; i < 6; i++)
      EXPECT_DOUBLE_EQ(expecVir[i] , lmp->force->bond->virial[i]); 
  }    
  
  // Test basic case with 4 atoms & 2 bond & coeff_argv{1, d0 = 100.0, alpha = 5.0, r_0 = 1.0}
  TEST_F(BondNonlinearTest, Bond2Run0_1) {     
    // Variables for run 
    double epsilon = 100.0;
    double r0 = 5.0;
    double lamda = 1.0;
    double atomNum = 4; 
    double bondNum = 2;

    // Variables for bond between atom 1 and 2
    double delx1 = 0.0-0.0; 
    double dely1 = 0.0-0.0; 
    double delz1 = -1.0 - 1.0; 
    double rsq1 = delx1*delx1 + dely1*dely1 + delz1*delz1;
    double r1 = sqrt(rsq1);
    double dr1 = r1 - r0;
    double drsq1 = dr1*dr1;
    double lamdasq1 = lamda*lamda;
    double denom1 = lamdasq1 - drsq1;
    double denomsq1 = denom1*denom1;
    double fbond1 = -epsilon/r1 * 2.0*dr1*lamdasq1/denomsq1;
    double ebond1 = epsilon * drsq1 / denom1; 

    // Variables for bond between atom 3 and 4
    double delx2 = 0.0-0.0;
    double dely2 = -1.0-1.0; 
    double delz2 = 0.0 - 0.0; 
    double rsq2 = delx2*delx2 + dely2*dely2 + delz2*delz2;
    double r2 = sqrt(rsq2);
    double dr2 = r2 - r0;
    double drsq2 = dr2*dr2;
    double lamdasq2 = lamda*lamda;
    double denom2 = lamdasq2 - drsq2;
    double denomsq2 = denom2*denom2;
    double fbond2 = -epsilon/r2 * 2.0*dr2*lamdasq2/denomsq2;
    double ebond2 = epsilon * drsq2 / denom2;  

    //*********************************************  
    // Create Bonds and Atoms

    char *lines;
    lines = strdup("1 1 1 0.0 0.0 -1.0\n 2 1 1 0.0 0.0 1.0\n 3 1 1 0.0 -1.0 0.0\n 4 1 1 0.0 1.0 0.0\n");
    // data_atoms(int n, char*) : unpack n lines from Atom section of data file
    lmp->atom->data_atoms(atomNum,lines);
    free(lines);
    lmp->atom->nlocal = atomNum;
    lmp->atom->natoms = atomNum;
  
    if (lmp->atom->map_style) {
      lmp->atom->map_init();
      lmp->atom->map_set();
    }

    // Create bond 1
    lines = strdup("1 1 1 2\n");
    lmp->atom->data_bonds(1,lines);
    free(lines);

    // Create bond 2
    lines = strdup("1 1 3 4\n");
    lmp->atom->data_bonds(1,lines);
    free(lines);

    lmp->input->one("mass 1 1.0");

    // End create bonds and atoms
    //*********************************************  

    lmp->force->bond = new BondNonlinear(lmp);
    //  k = 100.0 and r0 = 1
    char *bcoeff_argv[] = {"1", "100.0", "5", "1", NULL};
    int bcoeff_narg; 
    
    // Check for correct input
    bcoeff_narg = 4;
    lmp->force->bond->coeff(bcoeff_narg, bcoeff_argv);

    int i;
    for(i = 1; i <=  lmp->atom->nbondtypes; i++)
      EXPECT_DOUBLE_EQ(r0, lmp->force->bond->equilibrium_distance(i));

    // Tests that BondHarmonic::energy works properly 
    lmp->init();
    lmp->update->whichflag=1;
    lmp->update->firststep=0;
    lmp->update->laststep=1;
    lmp->input->one("run 0");
    EXPECT_DOUBLE_EQ((atomNum/bondNum)*(ebond1*0.5)+(atomNum/bondNum)*(ebond2*0.5), lmp->force->bond->energy); 

    // Test Forces
    // Forces for atom 1
    EXPECT_DOUBLE_EQ(fbond1*delx1, lmp->atom[0].f[0][0]); 
    EXPECT_DOUBLE_EQ(fbond1*dely1, lmp->atom[0].f[0][1]);
    EXPECT_DOUBLE_EQ(fbond1*delz1, lmp->atom[0].f[0][2]); 
    // Forces for atom 2
    EXPECT_DOUBLE_EQ(-1.0*fbond1*delx1, lmp->atom[0].f[1][0]); 
    EXPECT_DOUBLE_EQ(-1.0*fbond1*dely1, lmp->atom[0].f[1][1]);
    EXPECT_DOUBLE_EQ(-1.0*fbond1*delz1, lmp->atom[0].f[1][2]); 

    // Forces for atom 3
    EXPECT_DOUBLE_EQ(fbond2*delx2, lmp->atom[0].f[2][0]); 
    EXPECT_DOUBLE_EQ(fbond2*dely2, lmp->atom[0].f[2][1]);
    EXPECT_DOUBLE_EQ(fbond2*delz2, lmp->atom[0].f[2][2]); 
    // Forces for atom 4
    EXPECT_DOUBLE_EQ(-1.0*fbond2*delx2, lmp->atom[0].f[3][0]); 
    EXPECT_DOUBLE_EQ(-1.0*fbond2*dely2, lmp->atom[0].f[3][1]);
    EXPECT_DOUBLE_EQ(-1.0*fbond2*delz2, lmp->atom[0].f[3][2]); 
    
    // Total computed pressure bond 1 + bond 2
    double expecVir[6] = {delx1*delx1*fbond1 + delx2*delx2*fbond2,
			  dely1*dely1*fbond1 + dely2*dely2*fbond2,
			  delz1*delz1*fbond1 + delz2*delz2*fbond2,
			  delx1*dely1*fbond1 + delx2*dely2*fbond2,
			  delx1*delz1*fbond1 + delx2*delz2*fbond2,
                          dely1*delz1*fbond1 + dely2*delz2*fbond2 }; 
    for( i = 0; i < 6; i++)
      EXPECT_DOUBLE_EQ(expecVir[i] , lmp->force->bond->virial[i]);
  }
  // Test basic case with 4 atoms & 2 bond & coeff_argv{1, d0 = 50.0, alpha = 2.0, r_0 = 5.0}
  TEST_F(BondNonlinearTest, Bond2Run0_2) {     
    // Variables for run 
    double epsilon = 50.0;
    double r0 = 2.0;
    double lamda = 5.0;
    double atomNum = 4; 
    double bondNum = 2;

   // Variables for bond between atom 1 and 2
    double delx1 = 0.0-0.0; 
    double dely1 = 0.0-0.0; 
    double delz1 = -1.0 - 1.0; 
    double rsq1 = delx1*delx1 + dely1*dely1 + delz1*delz1;
    double r1 = sqrt(rsq1);
    double dr1 = r1 - r0;
    double drsq1 = dr1*dr1;
    double lamdasq1 = lamda*lamda;
    double denom1 = lamdasq1 - drsq1;
    double denomsq1 = denom1*denom1;
    double fbond1 = -epsilon/r1 * 2.0*dr1*lamdasq1/denomsq1;
    double ebond1 = epsilon * drsq1 / denom1; 

    // Variables for bond between atom 3 and 4
    double delx2 = 0.0-0.0;
    double dely2 = -1.0-1.0; 
    double delz2 = 0.0 - 0.0; 
    double rsq2 = delx2*delx2 + dely2*dely2 + delz2*delz2;
    double r2 = sqrt(rsq2);
    double dr2 = r2 - r0;
    double drsq2 = dr2*dr2;
    double lamdasq2 = lamda*lamda;
    double denom2 = lamdasq2 - drsq2;
    double denomsq2 = denom2*denom2;
    double fbond2 = -epsilon/r2 * 2.0*dr2*lamdasq2/denomsq2;
    double ebond2 = epsilon * drsq2 / denom2;  

    //*********************************************  
    // Create Bonds and Atoms

    char *lines;
    lines = strdup("1 1 1 0.0 0.0 -1.0\n 2 1 1 0.0 0.0 1.0\n 3 1 1 0.0 -1.0 0.0\n 4 1 1 0.0 1.0 0.0\n");
    // data_atoms(int n, char*) : unpack n lines from Atom section of data file
    lmp->atom->data_atoms(atomNum,lines);
    free(lines);
    lmp->atom->nlocal = atomNum;
    lmp->atom->natoms = atomNum;
  
    if (lmp->atom->map_style) {
      lmp->atom->map_init();
      lmp->atom->map_set();
    }

    // Create bond 1
    lines = strdup("1 1 1 2\n");
    lmp->atom->data_bonds(1,lines);
    free(lines);

    // Create bond 2
    lines = strdup("1 1 3 4\n");
    lmp->atom->data_bonds(1,lines);
    free(lines);

    lmp->input->one("mass 1 1.0");

    // End create bonds and atoms
    //*********************************************  

    lmp->force->bond = new BondNonlinear(lmp);
    //  k = 50.0 and r0 = 5
    char *bcoeff_argv[] = {"1", "50.0", "2", "5", NULL};
    int bcoeff_narg; 
    
    // Check for correct input
    bcoeff_narg = 4;
    lmp->force->bond->coeff(bcoeff_narg, bcoeff_argv);

    int i;
    for(i = 1; i <=  lmp->atom->nbondtypes; i++)
      EXPECT_DOUBLE_EQ(r0, lmp->force->bond->equilibrium_distance(i));

    // Tests that BondHarmonic::energy works properly 
    lmp->init();
    lmp->update->whichflag=1;
    lmp->update->firststep=0;
    lmp->update->laststep=1;
    lmp->input->one("run 0");
    EXPECT_DOUBLE_EQ((atomNum/bondNum)*(ebond1*0.5)+(atomNum/bondNum)*(ebond2*0.5), lmp->force->bond->energy); 

    // Test Forces
    // Forces for atom 1
    EXPECT_DOUBLE_EQ(fbond1*delx1, lmp->atom[0].f[0][0]); 
    EXPECT_DOUBLE_EQ(fbond1*dely1, lmp->atom[0].f[0][1]);
    EXPECT_DOUBLE_EQ(fbond1*delz1, lmp->atom[0].f[0][2]); 
    // Forces for atom 2
    EXPECT_DOUBLE_EQ(-1.0*fbond1*delx1, lmp->atom[0].f[1][0]); 
    EXPECT_DOUBLE_EQ(-1.0*fbond1*dely1, lmp->atom[0].f[1][1]);
    EXPECT_DOUBLE_EQ(-1.0*fbond1*delz1, lmp->atom[0].f[1][2]); 

    // Forces for atom 3
    EXPECT_DOUBLE_EQ(fbond2*delx2, lmp->atom[0].f[2][0]); 
    EXPECT_DOUBLE_EQ(fbond2*dely2, lmp->atom[0].f[2][1]);
    EXPECT_DOUBLE_EQ(fbond2*delz2, lmp->atom[0].f[2][2]); 
    // Forces for atom 4
    EXPECT_DOUBLE_EQ(-1.0*fbond2*delx2, lmp->atom[0].f[3][0]); 
    EXPECT_DOUBLE_EQ(-1.0*fbond2*dely2, lmp->atom[0].f[3][1]);
    EXPECT_DOUBLE_EQ(-1.0*fbond2*delz2, lmp->atom[0].f[3][2]); 
    
    // Total computed pressure bond 1 + bond 2
    double expecVir[6] = {delx1*delx1*fbond1 + delx2*delx2*fbond2,
			  dely1*dely1*fbond1 + dely2*dely2*fbond2,
			  delz1*delz1*fbond1 + delz2*delz2*fbond2,
			  delx1*dely1*fbond1 + delx2*dely2*fbond2,
			  delx1*delz1*fbond1 + delx2*delz2*fbond2,
                          dely1*delz1*fbond1 + dely2*delz2*fbond2 }; 
    for( i = 0; i < 6; i++)
      EXPECT_DOUBLE_EQ(expecVir[i] , lmp->force->bond->virial[i]);
  }

  // Test basic case with 4 atoms & 3 bond & coeff_argv{1, d0 = 100.0, alpha = 3.0, r_0 = 1.0}
  TEST_F(BondNonlinearTest, Bond3Run0_1) {     
    // Variables for run 
    double epsilon= 10.0;
    double r0 = 3.0;
    double lamda = 4.0;
    double atomNum = 4; 
    double bondNum = 3;

    // Variables for bond between atom 1 and 2
    double delx1 = 0.0-0.0; 
    double dely1 = 0.0-0.0; 
    double delz1 = -1.0 - 1.0; 
    double rsq1 = delx1*delx1 + dely1*dely1 + delz1*delz1;
    double r1 = sqrt(rsq1);
    double dr1 = r1 - r0;
    double drsq1 = dr1*dr1;
    double lamdasq1 = lamda*lamda;
    double denom1 = lamdasq1 - drsq1;
    double denomsq1 = denom1*denom1;
    double fbond1 = -epsilon/r1 * 2.0*dr1*lamdasq1/denomsq1;
    double ebond1 = epsilon * drsq1 / denom1; 

    // Variables for bond between atom 3 and 4
    double delx2 = 0.0-0.0;
    double dely2 = -1.0-1.0; 
    double delz2 = 0.0 - 0.0; 
    double rsq2 = delx2*delx2 + dely2*dely2 + delz2*delz2;
    double r2 = sqrt(rsq2);
    double dr2 = r2 - r0;
    double drsq2 = dr2*dr2;
    double lamdasq2 = lamda*lamda;
    double denom2 = lamdasq2 - drsq2;
    double denomsq2 = denom2*denom2;
    double fbond2 = -epsilon/r2 * 2.0*dr2*lamdasq2/denomsq2;
    double ebond2 = epsilon * drsq2 / denom2;  

    // Variables for bond between atom 1 and 4
    double delx3 = 0.0-0.0;
    double dely3 = 0.0-1.0; 
    double delz3 = -1.0 - 0.0; 
    double rsq3 = delx3*delx3 + dely3*dely3 + delz3*delz3;
    double r3 = sqrt(rsq3);
    double dr3 = r3 - r0;
    double drsq3 = dr3*dr3;
    double lamdasq3 = lamda*lamda;
    double denom3 = lamdasq3 - drsq3;
    double denomsq3 = denom3*denom3;
    double fbond3 = -epsilon/r3 * 2.0*dr3*lamdasq3/denomsq3;
    double ebond3 = epsilon * drsq3 / denom3;  

    //*********************************************  
    // Create Bonds and Atoms

    char *lines;
    lines = strdup("1 1 1 0.0 0.0 -1.0\n 2 1 1 0.0 0.0 1.0\n 3 1 1 0.0 -1.0 0.0\n 4 1 1 0.0 1.0 0.0\n");
    // data_atoms(int n, char*) : unpack n lines from Atom section of data file
    lmp->atom->data_atoms(atomNum,lines);
    free(lines);
    lmp->atom->nlocal = atomNum;
    lmp->atom->natoms = atomNum;
  
    if (lmp->atom->map_style) {
      lmp->atom->map_init();
      lmp->atom->map_set();
    }

    // Create bond 1
    lines = strdup("1 1 1 2\n");
    lmp->atom->data_bonds(1,lines);
    free(lines);

    // Create bond 2
    lines = strdup("1 1 3 4\n");
    lmp->atom->data_bonds(1,lines);
    free(lines);

    // Create bond 3
    lines = strdup("1 1 1 4\n");
    lmp->atom->data_bonds(1,lines);
    free(lines);

    lmp->input->one("mass 1 1.0");

    // End create bonds and atoms
    //*********************************************  

    lmp->force->bond = new BondNonlinear(lmp);
    //  k = 100.0 and r0 = 1
    char *bcoeff_argv[] = {"1", "10.0", "3.0", "4.0", NULL};
    int bcoeff_narg; 
    
    // Check for correct input
    bcoeff_narg = 4;
    lmp->force->bond->coeff(bcoeff_narg, bcoeff_argv);

    int i;
    for(i = 1; i <=  lmp->atom->nbondtypes; i++)
      EXPECT_DOUBLE_EQ(r0, lmp->force->bond->equilibrium_distance(i));

    // Tests that BondHarmonic::energy works properly 
    lmp->init();
    lmp->update->whichflag=1;
    lmp->update->firststep=0;
    lmp->update->laststep=1;
    lmp->input->one("run 0");
    EXPECT_DOUBLE_EQ((atomNum/2.0)*(ebond1*0.5)+(atomNum/2.0)*(ebond2*0.5)+(atomNum/2.0)*(ebond3*0.5), lmp->force->bond->energy); 

    // Test Forces
    // Forces for atom 1
    EXPECT_DOUBLE_EQ(fbond1*delx1 + fbond3*delx3, lmp->atom[0].f[0][0]); 
    EXPECT_DOUBLE_EQ(fbond1*dely1 + fbond3*dely3, lmp->atom[0].f[0][1]);
    EXPECT_DOUBLE_EQ(fbond1*delz1 + fbond3*delz3, lmp->atom[0].f[0][2]); 
    // Forces for atom 2
    EXPECT_DOUBLE_EQ(-1.0*fbond1*delx1, lmp->atom[0].f[1][0]); 
    EXPECT_DOUBLE_EQ(-1.0*fbond1*dely1, lmp->atom[0].f[1][1]);
    EXPECT_DOUBLE_EQ(-1.0*fbond1*delz1, lmp->atom[0].f[1][2]); 

    // Forces for atom 3
    EXPECT_DOUBLE_EQ(fbond2*delx2, lmp->atom[0].f[2][0]); 
    EXPECT_DOUBLE_EQ(fbond2*dely2, lmp->atom[0].f[2][1]);
    EXPECT_DOUBLE_EQ(fbond2*delz2, lmp->atom[0].f[2][2]); 
    // Forces for atom 4
    EXPECT_DOUBLE_EQ(-1.0*fbond2*delx2 + -1.0*fbond3*delx3, lmp->atom[0].f[3][0]); 
    EXPECT_DOUBLE_EQ(-1.0*fbond2*dely2 + -1.0*fbond3*dely3, lmp->atom[0].f[3][1]);
    EXPECT_DOUBLE_EQ(-1.0*fbond2*delz2 + -1.0*fbond3*delz3, lmp->atom[0].f[3][2]); 
    
    // Total computed pressure bond 1 + bond 2
    double expecVir[6] = {delx1*delx1*fbond1 + delx2*delx2*fbond2 + delx3*delx3*fbond3,
			  dely1*dely1*fbond1 + dely2*dely2*fbond2 + dely3*dely3*fbond3,
			  delz1*delz1*fbond1 + delz2*delz2*fbond2 + delz3*delz3*fbond3,
			  delx1*dely1*fbond1 + delx2*dely2*fbond2 + delx3*dely3*fbond3,
			  delx1*delz1*fbond1 + delx2*delz2*fbond2 + delx3*delz3*fbond3,
                          dely1*delz1*fbond1 + dely2*delz2*fbond2 + dely3*delz3*fbond3 }; 
    for( i = 0; i < 6; i++)
      EXPECT_DOUBLE_EQ(expecVir[i] , lmp->force->bond->virial[i]);     
  } 
  // Test basic case with 4 atoms & 3 bond & coeff_argv{1, d0 = 50.0, alpha = 9, r_0 = 5.0}
  TEST_F(BondNonlinearTest, Bond3Run0_2) {     
    // Variables for run 
    double epsilon = 50.0;
    double r0 = 9.0;
    double lamda = 5.0;
    double atomNum = 4; 
    double bondNum = 3;

 // Variables for bond between atom 1 and 2
    double delx1 = 0.0-0.0; 
    double dely1 = 0.0-0.0; 
    double delz1 = -1.0 - 1.0; 
    double rsq1 = delx1*delx1 + dely1*dely1 + delz1*delz1;
    double r1 = sqrt(rsq1);
    double dr1 = r1 - r0;
    double drsq1 = dr1*dr1;
    double lamdasq1 = lamda*lamda;
    double denom1 = lamdasq1 - drsq1;
    double denomsq1 = denom1*denom1;
    double fbond1 = -epsilon/r1 * 2.0*dr1*lamdasq1/denomsq1;
    double ebond1 = epsilon * drsq1 / denom1; 

    // Variables for bond between atom 3 and 4
    double delx2 = 0.0-0.0;
    double dely2 = -1.0-1.0; 
    double delz2 = 0.0 - 0.0; 
    double rsq2 = delx2*delx2 + dely2*dely2 + delz2*delz2;
    double r2 = sqrt(rsq2);
    double dr2 = r2 - r0;
    double drsq2 = dr2*dr2;
    double lamdasq2 = lamda*lamda;
    double denom2 = lamdasq2 - drsq2;
    double denomsq2 = denom2*denom2;
    double fbond2 = -epsilon/r2 * 2.0*dr2*lamdasq2/denomsq2;
    double ebond2 = epsilon * drsq2 / denom2;  

    // Variables for bond between atom 1 and 4
    double delx3 = 0.0-0.0;
    double dely3 = 0.0-1.0; 
    double delz3 = -1.0 - 0.0; 
    double rsq3 = delx3*delx3 + dely3*dely3 + delz3*delz3;
    double r3 = sqrt(rsq3);
    double dr3 = r3 - r0;
    double drsq3 = dr3*dr3;
    double lamdasq3 = lamda*lamda;
    double denom3 = lamdasq3 - drsq3;
    double denomsq3 = denom3*denom3;
    double fbond3 = -epsilon/r3 * 2.0*dr3*lamdasq3/denomsq3;
    double ebond3 = epsilon * drsq3 / denom3;  
 
    //*********************************************  
    // Create Bonds and Atoms

    char *lines;
    lines = strdup("1 1 1 0.0 0.0 -1.0\n 2 1 1 0.0 0.0 1.0\n 3 1 1 0.0 -1.0 0.0\n 4 1 1 0.0 1.0 0.0\n");
    // data_atoms(int n, char*) : unpack n lines from Atom section of data file
    lmp->atom->data_atoms(atomNum,lines);
    free(lines);
    lmp->atom->nlocal = atomNum;
    lmp->atom->natoms = atomNum;
  
    if (lmp->atom->map_style) {
      lmp->atom->map_init();
      lmp->atom->map_set();
    }

    // Create bond 1
    lines = strdup("1 1 1 2\n");
    lmp->atom->data_bonds(1,lines);
    free(lines);

    // Create bond 2
    lines = strdup("1 1 3 4\n");
    lmp->atom->data_bonds(1,lines);
    free(lines);

    // Create bond 3
    lines = strdup("1 1 1 4\n");
    lmp->atom->data_bonds(1,lines);
    free(lines);

    lmp->input->one("mass 1 1.0");

    // End create bonds and atoms
    //*********************************************  

    lmp->force->bond = new BondNonlinear(lmp);
    //  k = 50.0 and r0 = 5
    char *bcoeff_argv[] = {"1", "50.0", "9.0", "5.0", NULL};
    int bcoeff_narg; 
    
    // Check for correct input
    bcoeff_narg = 4;
    lmp->force->bond->coeff(bcoeff_narg, bcoeff_argv);

    int i;
    for(i = 1; i <=  lmp->atom->nbondtypes; i++)
      EXPECT_DOUBLE_EQ(r0, lmp->force->bond->equilibrium_distance(i));

    // Tests that BondHarmonic::energy works properly 
    lmp->init();
    lmp->update->whichflag=1;
    lmp->update->firststep=0;
    lmp->update->laststep=1;
    lmp->input->one("run 0");
    EXPECT_DOUBLE_EQ((atomNum/2)*(ebond1*0.5)+(atomNum/2)*(ebond2*0.5) + (atomNum/2)*(ebond3*0.5), lmp->force->bond->energy); 

    // Test Forces
    // Forces for atom 1
    EXPECT_DOUBLE_EQ(fbond1*delx1 + fbond3*delx3, lmp->atom[0].f[0][0]); 
    EXPECT_DOUBLE_EQ(fbond1*dely1 + fbond3*dely3, lmp->atom[0].f[0][1]);
    EXPECT_DOUBLE_EQ(fbond1*delz1 + fbond3*delz3, lmp->atom[0].f[0][2]); 
    // Forces for atom 2
    EXPECT_DOUBLE_EQ(-1.0*fbond1*delx1, lmp->atom[0].f[1][0]); 
    EXPECT_DOUBLE_EQ(-1.0*fbond1*dely1, lmp->atom[0].f[1][1]);
    EXPECT_DOUBLE_EQ(-1.0*fbond1*delz1, lmp->atom[0].f[1][2]); 

   // Forces for atom 3
    EXPECT_DOUBLE_EQ(fbond2*delx2, lmp->atom[0].f[2][0]); 
    EXPECT_DOUBLE_EQ(fbond2*dely2, lmp->atom[0].f[2][1]);
    EXPECT_DOUBLE_EQ(fbond2*delz2, lmp->atom[0].f[2][2]); 
    // Forces for atom 4
    EXPECT_DOUBLE_EQ(-1.0*fbond2*delx2 + -1.0*fbond3*delx3, lmp->atom[0].f[3][0]); 
    EXPECT_DOUBLE_EQ(-1.0*fbond2*dely2 + -1.0*fbond3*dely3, lmp->atom[0].f[3][1]);
    EXPECT_DOUBLE_EQ(-1.0*fbond2*delz2 + -1.0*fbond3*delz3, lmp->atom[0].f[3][2]); 
    
    // Total computed pressure bond 1 + bond 2
    double expecVir[6] = {delx1*delx1*fbond1 + delx2*delx2*fbond2 + delx3*delx3*fbond3,
			  dely1*dely1*fbond1 + dely2*dely2*fbond2 + dely3*dely3*fbond3,
			  delz1*delz1*fbond1 + delz2*delz2*fbond2 + delz3*delz3*fbond3,
			  delx1*dely1*fbond1 + delx2*dely2*fbond2 + delx3*dely3*fbond3,
			  delx1*delz1*fbond1 + delx2*delz2*fbond2 + delx3*delz3*fbond3,
                          dely1*delz1*fbond1 + dely2*delz2*fbond2 + dely3*delz3*fbond3 }; 
    for( i = 0; i < 6; i++)
      EXPECT_DOUBLE_EQ(expecVir[i] , lmp->force->bond->virial[i]);     
 } 
} // namespace
  int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
  }

