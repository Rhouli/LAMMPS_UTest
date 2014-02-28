/* Skylar Rudolph
 * BondFENEExpandTest()
 *
 * The following code tests bond_style bond_fene_expand in various ways:
 *
 * Coeff_correctNarg: checks bond->coeff() with correct input parameters
 * 
 * Coeff_highNarg: checks bond->coeff() for right exit status on incorrect input. Narg > expected
 * Coeff_NegativeNarg: checks bond->coeff() for right exit status on incorrect input. Narg = -expected
 * Coeff_lowNarg: checks bond->coeff() for right exit status on incorrect input. Narg < expected
 *
 * Tests the Force, Energy, and Pressure of each of the following simulations 
 * 
 * Bond1Run0_1: Test case with 2 atoms & 1 bond & coeff= {"1", "20.0", "12.0", "2.0", "5.0", "4.0"};
   Bond1Run0_2: Test case with 2 atoms & 1 bond & coeff= {"1", "40.0", "24.0", "4.0", "10.0", "8.0"};
 * 
 * Bond2Run0_1: Test case with 4 atoms & 2 bonds& coeff = {"1", "20.0", "12.0", "2.0", "5.0", "4.0"};
 * Bond2Run0_2: Test case with 4 atoms & 2 bond & coeff ={"1", "15.0", "13.0", "10.0", "6.0", "5.0"};
 * 
 * Bond3Run0_1: Test case with 4 atoms & 3 bond & coeff ={"1", "15.0", "13.0", "10.0", "6.0", "5.0"};
 * Bond3Run0_2: Test case with 4 atoms & 3 bond & coeff ={"1", "40.0", "25.0", "5.0", "6.0", "5.0"};
 *
 */

#include "test_bond_fene_expand.h"

using namespace std;

using namespace LAMMPS_NS;

//TWO_1_3 = cube root of two.
double TWO_1_3 = pow(2, 1/3);
namespace {

  // The fixture for testing class Foo.
  class BondFENEExpandTest : public ::testing::Test {   
  protected:
    BondFENEExpandTest() {
      char *argv[] = {"bond_fene_expand", "-screen", "none", "-log", "none", NULL};
      int narg = 5;
      
      //Create a LAMMPS Object
      lmp = new LAMMPS(narg,argv,MPI_COMM_WORLD);
      // Set atom_style
      lmp->input->one("atom_style bond");
      // Modify atom_style parameters
      lmp->input->one("atom_modify map array sort 10000 3.0");
    }

    virtual ~BondFENEExpandTest(){
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
 TEST_F(BondFENEExpandTest, Coeff_CorrectNarg){
    lmp->force->bond = new BondFENEExpand(lmp);
	
    char *bcoeff_argv[] = {"1", "30.0", "1.5", "1.0", "1.0", "10.0"};

    double shift = 10.0;
    double sigma = 1.0;
    // Check correct input
    
    int bcoeff_narg = 6;
    lmp->force->bond->coeff(bcoeff_narg, bcoeff_argv);
    
    int i;
    for(i = 1; i <=  lmp->atom->nbondtypes; i++){
      EXPECT_DOUBLE_EQ(0.97*sigma + shift, lmp->force->bond->equilibrium_distance(i));
      EXPECT_DOUBLE_EQ(1, lmp->force->bond->setflag[i]);
    }
 }

 //////////////////////////////////////////////////////////////////////
 // Multiple check's to make sure bond exits out with incorrect input//
 //////////////////////////////////////////////////////////////////////

  // Check for case with to high a value of bcoeff_narg
 TEST_F(BondFENEExpandTest, Coeff_highNarg){
    lmp->force->bond = new BondFENEExpand(lmp);
    char *bcoeff_argv[] = {"1", "100.0", "2.5", "5.0", "10.0", "20.0", NULL};
    int bcoeff_narg; 

    // Check incorrect input
    bcoeff_narg = 7;
    
    // Check that exit(1) was thrown
    ASSERT_EXIT(lmp->force->bond->coeff(bcoeff_narg, bcoeff_argv), ::testing::ExitedWithCode(1), "");  
 }

  // Check for case with to low a value of bcoeff_narg
 TEST_F(BondFENEExpandTest, Coeff_lowNarg){
    lmp->force->bond = new BondFENEExpand(lmp);
    char *bcoeff_argv[] = {NULL};
    int bcoeff_narg; 

    // Check incorrect input
    bcoeff_narg = 1;
    
    // Check that exit(1) was thrown
    ASSERT_EXIT(lmp->force->bond->coeff(bcoeff_narg, bcoeff_argv), ::testing::ExitedWithCode(1), "");  
 }
  // Check for case with for a negative value of bcoeff_narg
  TEST_F(BondFENEExpandTest,Coeff_NegativeNarg){
    lmp->force->bond = new BondFENEExpand(lmp);
    char *bcoeff_argv[] = {"1", "100.0", "2.5", "5.0", "10.0", "11.0"};
    int bcoeff_narg; 

    // Check incorrect input
    bcoeff_narg = -4;
    
    // Check that exit(1) was thrown
    ASSERT_EXIT(lmp->force->bond->coeff(bcoeff_narg, bcoeff_argv), ::testing::ExitedWithCode(1), "");  
  }

 //////////////////////////////////////////////////////////////////////
 // Multiple check's to make sure bond compute works properly        //
 //////////////////////////////////////////////////////////////////////
  // Lowrlogarg: Test special case where rlogarg < .1 but > -3 with 2 atoms & 1 bond & coeff={"1", "30.0", "1.5", "3.0", "2.0", ".5"}; 
 TEST_F(BondFENEExpandTest, Lowrlogarg) {
    // Variables for test
    double sigma   = 2.0;
    double epsilon = 3.0;
    double r0 = 1.5;
    double k = 30.0;
    double shift = 0.5;
    double atomNum = 2.0; 
    double TWO_1_3 = pow(2.0,(1.0/3.0)); 
 
    double delx = 0.0-0.0;
    double dely = 0.0-0.0; 
    double delz = -1.0 - 1.0; 
    
    double rsq = delx*delx + dely*dely + delz*delz;
    double r = sqrt(rsq);
    double rshift = r-shift;
    double rshiftsq = rshift*rshift;
    double r0sq = r0*r0;
    double rlogarg = 0.1;

    double fbond =-k*rshift/rlogarg/r;
    double ebond = -0.5 * k*r0sq*log(rlogarg);
    if (rshiftsq < TWO_1_3*sigma*sigma) {
      double sr2 = (sigma*sigma)/rshiftsq;
      double sr6 = sr2*sr2*sr2;
      ebond += 4.0*epsilon*sr6*(sr6-1.0) + epsilon;    
      fbond += 48.0*epsilon*sr6*(sr6-0.5)/rshift/r;
    }

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

    lmp->force->bond = new BondFENEExpand(lmp);
    char *bcoeff_argv[] = {"1", "30.0", "1.5", "3.0", "2.0", "0.5"};
    int bcoeff_narg; 
    
    // Check for correct input
    bcoeff_narg = 6;
    lmp->force->bond->coeff(bcoeff_narg, bcoeff_argv);

    // Tests that BondHarmonic::energy works properly 
    lmp->init();
    lmp->update->whichflag=1;
    lmp->update->firststep=0;
    lmp->update->laststep=1;
    lmp->input->one("run 0");

    // Energy for system
    EXPECT_DOUBLE_EQ(ebond, lmp->force->bond->energy);

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
    for(int i = 0; i < 6; i++)
      EXPECT_DOUBLE_EQ(expecVir[i] , lmp->force->bond->virial[i]); 
  }
  // VeryLowrlogarg: Test special case where rlogarg < -3 with 2 atoms & 1 bond & coeff={"1", "30.0", ".1", "3.0", "2.0", "0.5"};  
 TEST_F(BondFENEExpandTest, VeryLowrlogarg) {
    // Variables for test
    
    //*********************************************  
    // Create Bonds and Atoms
    double atomNum = 2.0; 
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

    lmp->force->bond = new BondFENEExpand(lmp);
    char *bcoeff_argv[] = {"1", "30.0", "0.1", "3.0", "2.0", "0.5"};
    int bcoeff_narg; 
    
    // Check for correct input
    bcoeff_narg = 6;
    lmp->force->bond->coeff(bcoeff_narg, bcoeff_argv);

    // Tests that BondHarmonic::energy works properly 
    lmp->init();
    lmp->update->whichflag=1;
    lmp->update->firststep=0;
    lmp->update->laststep=1;
    ASSERT_EXIT(lmp->input->one("run 0"), ::testing::ExitedWithCode(1), "");  
  }

  TEST_F(BondFENEExpandTest, Bond1Run0_1) {
    // Variables for test
    double sr2, sr6;
    double atomNum = 2;
    double delx = 0.0-0.0;
    double dely = 0.0-0.0; 
    double delz = -1.0 - 1.0; 

    double k = 20.0; 
    double r0 = 12.0;
    double epsilon = 2.0;
    double sigma = 5.0;
    double shift = 4.0;

    double rsq = delx*delx + dely*dely + delz*delz;
    double r = sqrt(rsq);
    double rshift = r - shift;
    double rshiftsq = rshift*rshift;
    double r0sq = r0 * r0;
    double rlogarg = 1.0 - rshiftsq/r0sq;
   
    double fbond = -1.0 * k * rshift/rlogarg/r;

    if (rshiftsq < TWO_1_3*sigma*sigma) {
      sr2 = sigma*sigma/rshiftsq;
      sr6 = sr2*sr2*sr2;
      fbond += 48.0*epsilon*sr6*(sr6-0.5)/rshift/r;
    }

   double ebond = -0.5 * k*r0sq*log(rlogarg);
      if (rshiftsq < TWO_1_3*sigma*sigma)
        ebond += 4.0*epsilon*sr6*(sr6-1.0) + epsilon;


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

    lmp->force->bond = new BondFENEExpand(lmp);
    char *bcoeff_argv[] = {"1", "20.0", "12.0", "2.0", "5.0", "4.0"};
    int bcoeff_narg; 
    
    // Check for correct input
    bcoeff_narg = 6;
    lmp->force->bond->coeff(bcoeff_narg, bcoeff_argv);

    int i;
    for(i = 1; i <=  lmp->atom->nbondtypes; i++){
      EXPECT_DOUBLE_EQ(0.97*sigma+shift, lmp->force->bond->equilibrium_distance(i));
    }

    // Tests that BondHarmonic::energy works properly 
    lmp->init();
    lmp->update->whichflag=1;
    lmp->update->firststep=0;
    lmp->update->laststep=1;
    lmp->input->one("run 0");

    // Energy for system
    EXPECT_DOUBLE_EQ(ebond, lmp->force->bond->energy);

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

  // Test basic case with 2 atoms
  TEST_F(BondFENEExpandTest, Bond1Run0_2) {
            // Variables for test

    double sr2, sr6;

    double k = 40.0; 
    double r0 = 24.0;
    double epsilon = 4.0;
    double sigma = 10.0;
    double shift = 8.0;


    double atomNum = 2;
    double delx = 0.0-0.0;
    double dely = 0.0-0.0; 
    double delz = -1.0 - 1.0; 
    double rsq = delx*delx + dely*dely + delz*delz;
    double r = sqrt(rsq);
    double rshift = r - shift;
    double rshiftsq = rshift*rshift;
    double r0sq = r0 * r0;
    double rlogarg = 1.0 - rshiftsq/r0sq;
    
    double fbond = -1.0 * k * rshift/rlogarg/r;
    double ebond = -0.5 * k*r0sq*log(rlogarg);
    if (rshiftsq < TWO_1_3*sigma*sigma) {
      sr2 = sigma*sigma/rshiftsq;
      sr6 = sr2*sr2*sr2;
      fbond += 48.0*epsilon*sr6*(sr6-0.5)/rshift/r;
      ebond += 4.0*epsilon*sr6*(sr6-1.0) + epsilon;
    }

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

    lmp->force->bond = new BondFENEExpand(lmp);
    char *bcoeff_argv[] = {"1", "40.0", "24.0", "4.0", "10.0", "8.0"};
    int bcoeff_narg; 
    
    // Check for correct input
    bcoeff_narg = 6;
    lmp->force->bond->coeff(bcoeff_narg, bcoeff_argv);

    int i;
    for(i = 1; i <=  lmp->atom->nbondtypes; i++){
      EXPECT_DOUBLE_EQ(0.97*sigma+shift, lmp->force->bond->equilibrium_distance(i));
    }

    // Tests that BondHarmonic::energy works properly 
    lmp->init();
    lmp->update->whichflag=1;
    lmp->update->firststep=0;
    lmp->update->laststep=1;
    lmp->input->one("run 0");

    // Energy for system
    EXPECT_DOUBLE_EQ(ebond, lmp->force->bond->energy);

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


TEST_F(BondFENEExpandTest, Bond2Run0_1) {     
    // Variables for test

    double k = 20.0; 
    double r0 = 12.0;
    double epsilon = 2.0;
    double sigma = 5.0;
    double shift = 4.0;

    double atomNum = 4.0;
    double bondNum = 2.0;

    // Bond 1
    double delx1 = 0.0-0.0;
    double dely1 = 0.0-0.0; 
    double delz1 = -1.0 - 1.0;

    double rsq1 = delx1*delx1 + dely1*dely1 + delz1*delz1;
    double r1 = sqrt(rsq1);
    double rshift1 = r1 - shift;
    double rshiftsq1 = rshift1*rshift1;
    double r0sq1 = r0 * r0;
    double rlogarg1 = 1.0 - rshiftsq1/r0sq1;

   double fbond1 = -1.0 * k * rshift1/rlogarg1/r1;
   double ebond1 = -0.5 * k*r0sq1*log(rlogarg1);
    if (rshiftsq1 < TWO_1_3*sigma*sigma) {
      double sr2_1 = sigma*sigma/rshiftsq1;
      double sr6_1 = sr2_1*sr2_1*sr2_1;
      fbond1 += 48.0*epsilon*sr6_1*(sr6_1-0.5)/rshift1/r1;
      ebond1 += 4.0*epsilon*sr6_1*(sr6_1-1.0) + epsilon;
    }

    // Bond 2
    double delx2 = 0.0-0.0;
    double dely2 = -1.0-1.0; 
    double delz2 = 0.0 - 0.0; 

    double rsq2 = delx2*delx2 + dely2*dely2 + delz2*delz2;
    double r2 = sqrt(rsq2);
    double rshift2 = r2 - shift;
    double rshiftsq2 = rshift2*rshift2;
    double r0sq2 = r0 * r0;
    double rlogarg2 = 1.0 - rshiftsq2/r0sq2;

   double fbond2 = -1.0 * k * rshift2/rlogarg2/r2;
   double ebond2 = -0.5 * k*r0sq2*log(rlogarg2); 
    if (rshiftsq2 < TWO_1_3*sigma*sigma) {
      double sr2_2 = sigma*sigma/rshiftsq2;
      double sr6_2 = sr2_2*sr2_2*sr2_2;
      fbond2 += 48.0*epsilon*sr6_2*(sr6_2-0.5)/rshift2/r2;
      ebond2 += 4.0*epsilon*sr6_2*(sr6_2-1.0) + epsilon;
    }

    //*********************************************  
    // Create Bonds and Atoms

    char *lines;
    lines = strdup("1 1 1 0.0 0.0 -1.0\n 2 1 1 0.0 0.0 1.0\n 3 1 1 0.0 -1.0 0.0\n 4 1 1 0.0 1.0 0.0\n");
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
    lmp->force->bond = new BondFENEExpand(lmp);

    char *bcoeff_argv[] = {"1", "20.0", "12.0", "2.0", "5.0", "4.0"};
    int bcoeff_narg; 
    
    // Check for correct input
    bcoeff_narg = 6;
    lmp->force->bond->coeff(bcoeff_narg, bcoeff_argv);

    int i;
    for(i = 1; i <=  lmp->atom->nbondtypes; i++){
      EXPECT_DOUBLE_EQ(0.97*sigma+shift, lmp->force->bond->equilibrium_distance(i));
    }

    // Tests that BondHarmonic::energy works properly 
    lmp->init();
    lmp->update->whichflag=1;
    lmp->update->firststep=0;
    lmp->update->laststep=1;
    lmp->input->one("run 0");

    // Energy for system
    EXPECT_DOUBLE_EQ((atomNum/bondNum)*(ebond1*0.5)+(atomNum/bondNum)*(ebond2*0.5), lmp->force->bond->energy); 

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

    // Total computed pressure
        double expecVir[6] = {delx1*delx1*fbond1 + delx2*delx2*fbond2,
			  dely1*dely1*fbond1 + dely2*dely2*fbond2,
			  delz1*delz1*fbond1 + delz2*delz2*fbond2,
			  delx1*dely1*fbond1 + delx2*dely2*fbond2,
			  delx1*delz1*fbond1 + delx2*delz2*fbond2,
                          dely1*delz1*fbond1 + dely2*delz2*fbond2 }; 
    for( i = 0; i < 6; i++)
      EXPECT_DOUBLE_EQ(expecVir[i] , lmp->force->bond->virial[i]);
}  

TEST_F(BondFENEExpandTest, Bond2Run0_2) {     
    // Variables for test
    double k = 15.0; 
    double r0 = 13.0;
    double epsilon = 10.0;
    double sigma = 6.0;
    double shift = 5.0;

    double atomNum = 4.0;
    double bondNum = 2.0;

    // Bond 1
    double delx1 = 0.0-0.0;
    double dely1 = 0.0-0.0; 
    double delz1 = -1.0 - 1.0;

    double rsq1 = delx1*delx1 + dely1*dely1 + delz1*delz1;
    double r1 = sqrt(rsq1);
    double rshift1 = r1 - shift;
    double rshiftsq1 = rshift1*rshift1;
    double r0sq1 = r0 * r0;
    double rlogarg1 = 1.0 - rshiftsq1/r0sq1;

   double fbond1 = -1.0 * k * rshift1/rlogarg1/r1;
   double ebond1 = -0.5 * k*r0sq1*log(rlogarg1);
   if (rshiftsq1 < TWO_1_3*sigma*sigma) {
      double sr2_1 = sigma*sigma/rshiftsq1;
      double sr6_1 = sr2_1*sr2_1*sr2_1;
      fbond1 += 48.0*epsilon*sr6_1*(sr6_1-0.5)/rshift1/r1;
      ebond1 += 4.0*epsilon*sr6_1*(sr6_1-1.0) + epsilon;    
    }
    
    // Bond 2
     double delx2 = 0.0-0.0;
     double dely2 = -1.0-1.0; 
     double delz2 = 0.0 - 0.0; 

     double rsq2 = delx2*delx2 + dely2*dely2 + delz2*delz2;
     double r2 = sqrt(rsq2);
     double rshift2 = r2 - shift;
     double rshiftsq2 = rshift2*rshift2;
     double r0sq2 = r0 * r0;
     double rlogarg2 = 1.0 - rshiftsq2/r0sq2;

     double fbond2 = -1.0 * k * rshift2/rlogarg2/r2;
     double ebond2 = -0.5 * k*r0sq2*log(rlogarg2); 
     if (rshiftsq2 < TWO_1_3*sigma*sigma) {
       double sr2_2 = sigma*sigma/rshiftsq2;
       double sr6_2 = sr2_2*sr2_2*sr2_2;
       fbond2 += 48.0*epsilon*sr6_2*(sr6_2-0.5)/rshift2/r2;
       ebond2 += 4.0*epsilon*sr6_2*(sr6_2-1.0) + epsilon;
     }

     //*********************************************  
     // Create Bonds and Atoms

     char *lines;
     lines = strdup("1 1 1 0.0 0.0 -1.0\n 2 1 1 0.0 0.0 1.0\n 3 1 1 0.0 -1.0 0.0\n 4 1 1 0.0 1.0 0.0\n");
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

     lmp->force->bond = new BondFENEExpand(lmp);
     char *bcoeff_argv[] = {"1", "15.0", "13.0", "10.0", "6.0", "5.0"};
     int bcoeff_narg; 

     // Check for correct input
     bcoeff_narg = 6;
     lmp->force->bond->coeff(bcoeff_narg, bcoeff_argv);

     int i;
     for(i = 1; i <=  lmp->atom->nbondtypes; i++){
       EXPECT_DOUBLE_EQ(0.97*sigma+shift, lmp->force->bond->equilibrium_distance(i));
     }

     // Tests that BondHarmonic::energy works properly 
     lmp->init();
     lmp->update->whichflag=1;
     lmp->update->firststep=0;
     lmp->update->laststep=1;
     lmp->input->one("run 0");

     // Energy for system
     EXPECT_DOUBLE_EQ((atomNum/bondNum)*(ebond1*0.5)+(atomNum/bondNum)*(ebond2*0.5), lmp->force->bond->energy); 

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

     // Total computed pressure
     double expecVir[6] = {delx1*delx1*fbond1 + delx2*delx2*fbond2,
			   dely1*dely1*fbond1 + dely2*dely2*fbond2,
			   delz1*delz1*fbond1 + delz2*delz2*fbond2,
			   delx1*dely1*fbond1 + delx2*dely2*fbond2,
			   delx1*delz1*fbond1 + delx2*delz2*fbond2,
			   dely1*delz1*fbond1 + dely2*delz2*fbond2 }; 
     for( i = 0; i < 6; i++)
       EXPECT_DOUBLE_EQ(expecVir[i] , lmp->force->bond->virial[i]);
 }



   TEST_F(BondFENEExpandTest, Bond3Run0_1) 
   {     
     // Variables for run 
     double atomNum = 4.0; 
     double bondNum = 3.0;

     double k = 15.0; 
     double r0 = 13.0;
     double epsilon = 10.0;
     double sigma = 6.0;
     double shift = 5.0;

     // Bond 1
     double delx1 = 0.0-0.0;
     double dely1 = 0.0-0.0; 
     double delz1 = -1.0 - 1.0;

     double rsq1 = delx1*delx1 + dely1*dely1 + delz1*delz1;
     double r1 = sqrt(rsq1);
     double rshift1 = r1 - shift;
     double rshiftsq1 = rshift1*rshift1;
     double r0sq1 = r0 * r0;
     double rlogarg1 = 1.0 - rshiftsq1/r0sq1;

     double fbond1 = -1.0 * k * rshift1/rlogarg1/r1;
     double ebond1 = -0.5 * k*r0sq1*log(rlogarg1);
     if (rshiftsq1 < TWO_1_3*sigma*sigma) {
       double sr2_1 = sigma*sigma/rshiftsq1;
       double sr6_1 = sr2_1*sr2_1*sr2_1;
       fbond1 += 48.0*epsilon*sr6_1*(sr6_1-0.5)/rshift1/r1;
       ebond1 += 4.0*epsilon*sr6_1*(sr6_1-1.0) + epsilon;
     }

     // Bond 2
     double delx2 = 0.0-0.0; 
     double dely2 = -1.0-1.0;
     double delz2 = 0.0 - 0.0; 

     double rsq2 = delx2*delx2 + dely2*dely2 + delz2*delz2;
     double r2 = sqrt(rsq2);
     double rshift2 = r2 - shift;
     double rshiftsq2 = rshift2*rshift2;
     double r0sq2 = r0 * r0;
     double rlogarg2 = 1.0 - rshiftsq2/r0sq2;

     double fbond2 = -1.0 * k * rshift2/rlogarg2/r2;
     double ebond2 = -0.5 * k*r0sq2*log(rlogarg2);
     if (rshiftsq2 < TWO_1_3*sigma*sigma) {
       double sr2_2 = sigma*sigma/rshiftsq2;
       double sr6_2 = sr2_2*sr2_2*sr2_2;
       fbond2 += 48.0*epsilon*sr6_2*(sr6_2-0.5)/rshift2/r2;
       ebond2 += 4.0*epsilon*sr6_2*(sr6_2-1.0) + epsilon;
     }

     // Bond 3
     double delx3 = 0.0-0.0;
     double delz3 = -1.0 - 0.0; 
     double dely3 = 0.0-1.0;  

     double  rsq3 = delx3*delx3 + dely3*dely3 + delz3*delz3;
     double r3 = sqrt(rsq3);
     double rshift3 = r3 - shift;
     double rshiftsq3 = rshift3*rshift3;
     double r0sq3 = r0 * r0;
     double rlogarg3 = 1.0 - rshiftsq3/r0sq3;

     double fbond3 = -1.0 * k * rshift3/rlogarg3/r3;   
     double ebond3 = -0.5 * k*r0sq3*log(rlogarg3);
     if (rshiftsq3 < TWO_1_3*sigma*sigma) {
       double sr2_3 = sigma*sigma/rshiftsq3;
       double sr6_3 = sr2_3*sr2_3*sr2_3;
       fbond3 += 48.0*epsilon*sr6_3*(sr6_3-0.5)/rshift3/r3;
       ebond3 += 4.0*epsilon*sr6_3*(sr6_3-1.0) + epsilon; 
     }

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

     lmp->force->bond = new BondFENEExpand(lmp);
     //  k = 50.0 and r0 = 5
     char *bcoeff_argv[] = {"1", "15.0", "13.0", "10.0", "6.0", "5.0"};
     int bcoeff_narg; 

     // Check for correct input
     bcoeff_narg = 6;
     lmp->force->bond->coeff(bcoeff_narg, bcoeff_argv);

     int i;
     for(i = 1; i <=  lmp->atom->nbondtypes; i++)
       EXPECT_DOUBLE_EQ(0.97*sigma + shift, lmp->force->bond->equilibrium_distance(i));

     // Tests that BondFENE:Energy is proper.
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


   TEST_F(BondFENEExpandTest, Bond3Run0_2) 
   {     //"40.0", "25.0", "5.0", "6.0", "5.0"
     // Variables for run 
     double atomNum = 4.0; 
     double bondNum = 3.0;

     double k = 40.0; 
     double r0 = 25.0;
     double epsilon = 5.0;
     double sigma = 6.0;
     double shift = 5.0;

    // Bond 1
    double delx1 = 0.0-0.0;
    double dely1 = 0.0-0.0; 
    double delz1 = -1.0 - 1.0;

    double rsq1 = delx1*delx1 + dely1*dely1 + delz1*delz1;
    double r1 = sqrt(rsq1);
    double rshift1 = r1 - shift;
    double rshiftsq1 = rshift1*rshift1;
    double r0sq1 = r0 * r0;
    double rlogarg1 = 1.0 - rshiftsq1/r0sq1;
    
    double fbond1 = -1.0 * k * rshift1/rlogarg1/r1;
    double ebond1 = -0.5 * k*r0sq1*log(rlogarg1);
    if (rshiftsq1 < TWO_1_3*sigma*sigma) {
      double sr2_1 = sigma*sigma/rshiftsq1;
      double sr6_1 = sr2_1*sr2_1*sr2_1;
      fbond1 += 48.0*epsilon*sr6_1*(sr6_1-0.5)/rshift1/r1;
      ebond1 += 4.0*epsilon*sr6_1*(sr6_1-1.0) + epsilon;
    }

    // Bond 2
    double delx2 = 0.0-0.0; 
    double dely2 = -1.0-1.0;
    double delz2 = 0.0 - 0.0; 

    double rsq2 = delx2*delx2 + dely2*dely2 + delz2*delz2;
    double r2 = sqrt(rsq2);
    double rshift2 = r2 - shift;
    double rshiftsq2 = rshift2*rshift2;
    double r0sq2 = r0 * r0;
    double rlogarg2 = 1.0 - rshiftsq2/r0sq2;
    
    double fbond2 = -1.0 * k * rshift2/rlogarg2/r2;
    double ebond2 = -0.5 * k*r0sq2*log(rlogarg2);
    if (rshiftsq2 < TWO_1_3*sigma*sigma) {
      double sr2_2 = sigma*sigma/rshiftsq2;
      double sr6_2 = sr2_2*sr2_2*sr2_2;
      fbond2 += 48.0*epsilon*sr6_2*(sr6_2-0.5)/rshift2/r2;
      ebond2 += 4.0*epsilon*sr6_2*(sr6_2-1.0) + epsilon;
    }
    
    // Bond 3
    double delx3 = 0.0-0.0;
    double delz3 = -1.0 - 0.0; 
    double dely3 = 0.0-1.0;  

    double  rsq3 = delx3*delx3 + dely3*dely3 + delz3*delz3;
    double r3 = sqrt(rsq3);
    double rshift3 = r3 - shift;
    double rshiftsq3 = rshift3*rshift3;
    double r0sq3 = r0 * r0;
    double rlogarg3 = 1.0 - rshiftsq3/r0sq3;

    double fbond3 = -1.0 * k * rshift3/rlogarg3/r3;   
    double ebond3 = -0.5 * k*r0sq3*log(rlogarg3);
    if (rshiftsq3 < TWO_1_3*sigma*sigma) {
      double sr2_3 = sigma*sigma/rshiftsq3;
      double sr6_3 = sr2_3*sr2_3*sr2_3;
      fbond3 += 48.0*epsilon*sr6_3*(sr6_3-0.5)/rshift3/r3;
      ebond3 += 4.0*epsilon*sr6_3*(sr6_3-1.0) + epsilon; 
   }

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

    lmp->force->bond = new BondFENEExpand(lmp);
    //  k = 50.0 and r0 = 5
    char *bcoeff_argv[] = {"1", "40.0", "25.0", "5.0", "6.0", "5.0"};
    int bcoeff_narg; 
    
    // Check for correct input
    bcoeff_narg = 6;
    lmp->force->bond->coeff(bcoeff_narg, bcoeff_argv);

    int i;
    for(i = 1; i <=  lmp->atom->nbondtypes; i++)
      EXPECT_DOUBLE_EQ(0.97*sigma + shift, lmp->force->bond->equilibrium_distance(i));

    // Tests that BondFENE:Energy is proper.
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

