/* Author: Ryan Houlihan
 * BondTableTest():
 *
 * The following code tests bond_style bond_table in various ways:
 * 
 * Settings_CorrectArguments: Checks valid settings input
 * Settings_BadArguments: Checks invalid settings input
 * 
 * Coeff_Correct: Check coeff with corrent input
 * Coeff_RandomTests: Check coeff for various input
 *
 * Tests the Force, Energy, and Pressure of each of the following simulations: 
 *
 * BondLinearRun0_1:
 * Test basic case with 2 atoms & 1 linear table bond
 * input file is as follows
 * N 2 FP 0 0 EQ 0.5                                      
 * 1 0.00 338.0000 1352.0000  
 * 2 0.01 324.6152 1324.9600
 *
 * BondSplineRun0_1:
 * Test basic case with 2 atoms & 1 spline table bond
 * input file is as follows
 * N 2 FP 0 0 EQ 0.5                                      
 * 1 0.00 338.0000 1352.0000  
 * 2 0.01 324.6152 1324.9600
 *
 */

#include "test_bond_table.h"

using namespace std;

using namespace LAMMPS_NS;

namespace {

  // The fixture for testing class BondTableTest
  class BondTableTest : public ::testing::Test {   
  protected:
    BondTableTest() {
      char *argv[] = {"bond_table", "-screen", "none", "-log", "none", NULL};
      //char *argv[] = {"bond_harmonic", "-echo", "screen", "-log", "none", NULL};
      int narg = 5;
      
      //Create a LAMMPS Object
      lmp = new LAMMPS(narg,argv,MPI_COMM_WORLD);
      // Set atom_style
      lmp->input->one("atom_style bond");
      // Modify atom_style parameters
      lmp->input->one("atom_modify map array sort 10000 3.0");
    }

    virtual ~BondTableTest(){
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
      
      // set up atoms
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

  // Check for case with correct value of bsettings_narg
 TEST_F(BondTableTest, Settings_CorrectArguments){
    lmp->force->bond = new BondTable(lmp);
    char *bsettings_argv[] = {"linear", "2", NULL};
    int bsettings_narg; 

    // Check correct input
    bsettings_narg = 2;
    lmp->force->bond->settings(bsettings_narg, bsettings_argv);  

    bsettings_argv = {"spline", "2", NULL};
    bsettings_narg = 2;
    lmp->force->bond->settings(bsettings_narg, bsettings_argv);  
 }

  // Check for various incorrect arguments settings()
  TEST_F(BondTableTest, Settings_BadArguments){
    lmp->force->bond = new BondTable(lmp);
    char *bsettings_argv[] =  {NULL, NULL, NULL};
    int bsettings_narg; 

    // Check for case with incorrect arguments for bsettings_argv
    bsettings_argv = {"harmonic", "1", NULL};
    bsettings_narg = 1;
    ASSERT_EXIT(lmp->force->bond->settings(bsettings_narg, bsettings_argv), ::testing::ExitedWithCode(1), "");  

    // Check for case with incorrect value for bsettings_argv
    bsettings_argv = {"linear", "4", NULL};
    bsettings_narg = 1;
    ASSERT_EXIT(lmp->force->bond->settings(bsettings_narg, bsettings_argv), ::testing::ExitedWithCode(1), "");  

    // Check for case with incorrect value for bsettings_argv
    bsettings_argv = {"linear", "1", NULL};
    bsettings_narg = 1;
    ASSERT_EXIT(lmp->force->bond->settings(bsettings_narg, bsettings_argv), ::testing::ExitedWithCode(1), "");  
       
    // Check for case with incorrect value of bsettings_narg
    bsettings_argv = {"linear", "-2", NULL};
    bsettings_narg = -2;
    ASSERT_EXIT(lmp->force->bond->settings(bsettings_narg, bsettings_argv), ::testing::ExitedWithCode(1), "");  
    // Check for case with incorrect value for bsettings_argv
    bsettings_argv = {"spline", "4", NULL};
    bsettings_narg = 1;
    ASSERT_EXIT(lmp->force->bond->settings(bsettings_narg, bsettings_argv), ::testing::ExitedWithCode(1), "");  

    // Check for case with incorrect value for bsettings_argv
    bsettings_argv = {"spline", "1", NULL};
    bsettings_narg = 1;
    ASSERT_EXIT(lmp->force->bond->settings(bsettings_narg, bsettings_argv), ::testing::ExitedWithCode(1), "");  
       
    // Check for case with incorrect value of bsettings_narg
    bsettings_argv = {"spline", "-2", NULL};
    bsettings_narg = -2;
    ASSERT_EXIT(lmp->force->bond->settings(bsettings_narg, bsettings_argv), ::testing::ExitedWithCode(1), "");  
 }

 //////////////////////////////////////////////////////////////////////
 // Multiple check's to make sure bond exits out with incorrect input//
 //////////////////////////////////////////////////////////////////////

  // Check for case correct coeff arguments
 TEST_F(BondTableTest, Coeff_Correct){
    lmp->force->bond = new BondTable(lmp);

    // Check for linear
    char *bsettings_argv[] = {"linear", "2", NULL};
    int bsettings_narg; 

    bsettings_narg = 2;
    lmp->force->bond->settings(bsettings_narg, bsettings_argv);  

    char *bcoeff_argv[] = {"1", "Input/file.tableTests", "COEFF_CORRECT", NULL};
    int bcoeff_narg; 

    bcoeff_narg = 3;
    lmp->force->bond->coeff(bcoeff_narg, bcoeff_argv);
    EXPECT_EQ(0.5, lmp->force->bond->equilibrium_distance(1));
    
    // Check for spline
    bsettings_argv = {"spline", "2", NULL};
    lmp->force->bond->settings(bsettings_narg, bsettings_argv);  

    lmp->force->bond->coeff(bcoeff_narg, bcoeff_argv);
    EXPECT_EQ(0.5, lmp->force->bond->equilibrium_distance(1));
 }

  // Check for various coeff arguments
 TEST_F(BondTableTest, Coeff_RandomTests){
    lmp->force->bond = new BondTable(lmp);

    // Check for linear
    char *bsettings_argv[] = {"linear", "2", NULL};
    int bsettings_narg; 

    bsettings_narg = 2;
    lmp->force->bond->settings(bsettings_narg, bsettings_argv);  

    char *bcoeff_argv[] = {NULL, NULL, NULL, NULL};
    int bcoeff_narg; 

    // Check incorrect input
    bcoeff_argv = {"1", "Input/file.tableTests", "COEFF_RANDOM_1", NULL};
    bcoeff_narg = 4;    
    ASSERT_EXIT(lmp->force->bond->coeff(bcoeff_narg, bcoeff_argv), ::testing::ExitedWithCode(1), "");    

    // Check incorrect input
    bcoeff_argv = {"1", "Input/file.tableTests", "COEFF_RANDOM_1", NULL};
    bcoeff_narg = 2;    
    ASSERT_EXIT(lmp->force->bond->coeff(bcoeff_narg, bcoeff_argv), ::testing::ExitedWithCode(1), "");    

    // Check incorrect input
    bcoeff_argv = {"1", "Input/file.tableTests", "COEFF_RANDOM_1", NULL};
    bcoeff_narg = -3;    
    ASSERT_EXIT(lmp->force->bond->coeff(bcoeff_narg, bcoeff_argv), ::testing::ExitedWithCode(1), "");    

    // Check incorrect bond table length
    bcoeff_argv = {"1", "Input/file.tableTests", "COEFF_RANDOM_1", NULL};
    bcoeff_narg = 3;    
    ASSERT_EXIT(lmp->force->bond->coeff(bcoeff_narg, bcoeff_argv), ::testing::ExitedWithCode(1), "");

    // Check incorrect bond table input
    bcoeff_argv = {"1", "Input/file.tableTests", "COEFF_RANDOM_2", NULL};
    bcoeff_narg = 3;    
    ASSERT_EXIT(lmp->force->bond->coeff(bcoeff_narg, bcoeff_argv), ::testing::ExitedWithCode(1), "");

    // Check incorrect table
    bcoeff_argv = {"1", "Input/file.tableTests", "COEFF_RANDOM_3", NULL};
    bcoeff_narg = 3;    
    ASSERT_EXIT(lmp->force->bond->coeff(bcoeff_narg, bcoeff_argv), ::testing::ExitedWithCode(1), "");

    // Check incorrect table
    bcoeff_argv = {"1", "Input/file.tableTests", "COEFF_RANDOM_4", NULL};
    bcoeff_narg = 3;    
    ASSERT_EXIT(lmp->force->bond->coeff(bcoeff_narg, bcoeff_argv), ::testing::ExitedWithCode(1), "");

    // Check invalid keyword
    bcoeff_argv = {"1", "Input/file.tableTests", "INCORRECT", NULL};
    bcoeff_narg = 3;    
    ASSERT_EXIT(lmp->force->bond->coeff(bcoeff_narg, bcoeff_argv), ::testing::ExitedWithCode(1), "");
}

 //////////////////////////////////////////////////////////////////////
 // Multiple check's to make sure bond compute works properly        //
 //////////////////////////////////////////////////////////////////////

  /* BondLinearRun0_1:
   * Test basic case with 2 atoms & 1 linear table bond
   * input file is as follows
   * N 2 FP 0 0 EQ 0.5                                      
   * 1 0.00 338.0000 1352.0000  
   * 2 0.01 324.6152 1324.9600
   */
  TEST_F(BondTableTest, BondLinearRun0_1) {
    double delx = 0.0-0.0;
    double dely = 0.0-0.0; 
    double delz = -1.0 - 1.0; 
    double atomNum = 2; 
    int tablelength = 2.0;
    int tlm1 = tablelength - 1.0;
    double rsq = delx*delx + dely*dely + delz*delz;
    double r = sqrt(rsq); 
    double rf[] = {0.00, 0.01};
    double ef[] = {338.0, 324.6152};
    double ff[] = {1352.0, 1324.96};
    double def[tlm1], dff[tlm1];
    for(int x = 0; x < tlm1; x++){
      def[x] = ef[x+1] - ef[x];
      dff[x] = ff[x+1] - ff[x];
    }
    // Input variables
    double r0 = 0.5;
    double ninput = 2;
    double fpflag = 1;
    double rlow = 0; 
    double rhi = .01;   
    double maxN = 1;
    double minN = 2;

    // Force and energy computations
    double delta = (rhi-rlow)/tlm1;
    double invdelta = 1.0/ delta; 
    double deltasq6 = (maxN-minN/6.0);
    double x;
    r > rlow ? x = r : x = rlow;
    if(x > rhi) x = rhi;
    int itable = static_cast<int> ((x-rlow)*invdelta);
    double fraction =  (x-rf[itable])*invdelta;
    double u = ef[itable] + fraction*def[itable];
    double f = ff[itable] + fraction*dff[itable];

    double fbond = f/r;
    double ebond = u;
    lmp->force->bond = new BondTable(lmp);

    // Check for linear
    char *bsettings_argv[] = {"linear", "2", NULL};
    int bsettings_narg; 

    bsettings_narg = 2;
    lmp->force->bond->settings(bsettings_narg, bsettings_argv);  

    char *bcoeff_argv[] = {"1", "Input/file.tableTests", "COEFF_CORRECT", NULL};
    int bcoeff_narg; 

    bcoeff_narg = 3;
    
    lmp->force->bond->coeff(bcoeff_narg, bcoeff_argv);

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

    lmp->input->one("mass * 1.0");

    // End create bonds and atoms
    //*********************************************  

    // Tests that BondTable::energy works properly 
    lmp->init();
    lmp->update->whichflag=1;
    lmp->update->firststep=0;
    lmp->update->laststep=1;
    lmp->input->one("run 0");

    // Energy for system
    EXPECT_EQ(ebond, lmp->force->bond->energy);

    // Forces for atom 1
    EXPECT_EQ(fbond*delx, lmp->atom[0].f[0][0]); 
    EXPECT_EQ(fbond*dely, lmp->atom[0].f[0][1]);
    EXPECT_EQ(fbond*delz, lmp->atom[0].f[0][2]); 

    // Forces for atom 2
    EXPECT_EQ(-1.0*fbond*delx, lmp->atom[0].f[1][0]); 
    EXPECT_EQ(-1.0*fbond*dely, lmp->atom[0].f[1][1]);
    EXPECT_EQ(-1.0*fbond*delz, lmp->atom[0].f[1][2]); 

    // Total computed pressure
    double expecVir[6] = {delx*delx*fbond,
			  dely*dely*fbond,
			  delz*delz*fbond,
			  delx*dely*fbond,
			  delx*delz*fbond,
                          dely*delz*fbond }; 
    for(int i = 0; i < 6; i++)
      EXPECT_EQ(expecVir[i] , lmp->force->bond->virial[i]); 
 
 }    

/* BondSplineRun0_1:
   * Test basic case with 2 atoms & 1 spline table bond
   * input file is as follows
   * N 2 FP 0 0 EQ 0.5                                      
   * 1 0.00 338.0000 1352.0000  
   * 2 0.01 324.6152 1324.9600
   */
  TEST_F(BondTableTest, BondSplineRun0_1) {
    double delx = 0.0-0.0;
    double dely = 0.0-0.0; 
    double delz = -1.0 - 1.0; 
    double atomNum = 2; 
    int tablelength = 2.0;
    int tlm1 = tablelength - 1.0;
    double rsq = delx*delx + dely*dely + delz*delz;
    double r = sqrt(rsq); 
    double rf[] = {0.00, 0.01};
    double ef[] = {338.0, 324.6152};
    double ff[] = {1352.0, 1324.96};
    double def[tlm1], dff[tlm1];
    for(int x = 0; x < tlm1; x++){
      def[x] = ef[x+1] - ef[x];
      dff[x] = ff[x+1] - ff[x];
    }
    // Input variables
    double r0 = 0.5;
    double ninput = 2;
    double fpflag = 1;
    double rlow = 0; 
    double rhi = .01;   
    double maxN = 1;
    double minN = 2;

    // Force and energy computations
    double delta = (rhi-rlow)/tlm1;
    double invdelta = 1.0/ delta; 
    double deltasq6 = (maxN-minN/6.0);
    double x;
    r > rlow ? x = r : x = rlow;
    if(x > rhi) x = rhi;
    int itable = static_cast<int> ((x-rlow)*invdelta);
    double fraction =  (x-rf[itable])*invdelta;
    double b = (x - rf[itable]) *invdelta;
    double a = 1.0-b;
    double u = a * ef[itable] + b * ef[itable+1] + 
      ((a*a*a-a)*pow(ef[itable],2) + (b*b*b-b)*pow(ef[itable+1],2)) * 
      deltasq6;
    
    double f = a * ff[itable] + b * ff[itable+1] + 
      ((a*a*a-a)*pow(ff[itable],2) + (b*b*b-b)*pow(ff[itable+1],2)) * 
      deltasq6;

    double fbond = f/r;
    double ebond = u;
    lmp->force->bond = new BondTable(lmp);

    // Check for linear
    char *bsettings_argv[] = {"spline", "2", NULL};
    int bsettings_narg; 

    bsettings_narg = 2;
    lmp->force->bond->settings(bsettings_narg, bsettings_argv);  

    char *bcoeff_argv[] = {"1", "Input/file.tableTests", "COEFF_CORRECT", NULL};
    int bcoeff_narg; 

    bcoeff_narg = 3;
    
    lmp->force->bond->coeff(bcoeff_narg, bcoeff_argv);

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

    lmp->input->one("mass * 1.0");

    // End create bonds and atoms
    //*********************************************  

    // Tests that BondTable::energy works properly 
    lmp->init();
    lmp->update->whichflag=1;
    lmp->update->firststep=0;
    lmp->update->laststep=1;
    lmp->input->one("run 0");

    // Energy for system
    EXPECT_EQ(ebond, lmp->force->bond->energy);

    // Forces for atom 1
    EXPECT_EQ(fbond*delx, lmp->atom[0].f[0][0]); 
    EXPECT_EQ(fbond*dely, lmp->atom[0].f[0][1]);
    EXPECT_EQ(fbond*delz, lmp->atom[0].f[0][2]); 

    // Forces for atom 2
    EXPECT_EQ(-1.0*fbond*delx, lmp->atom[0].f[1][0]); 
    EXPECT_EQ(-1.0*fbond*dely, lmp->atom[0].f[1][1]);
    EXPECT_EQ(-1.0*fbond*delz, lmp->atom[0].f[1][2]); 

    // Total computed pressure
    double expecVir[6] = {delx*delx*fbond,
			  dely*dely*fbond,
			  delz*delz*fbond,
			  delx*dely*fbond,
			  delx*delz*fbond,
                          dely*delz*fbond }; 
    for(int i = 0; i < 6; i++)
      EXPECT_EQ(expecVir[i] , lmp->force->bond->virial[i]); 
 
 }    
} // namespace
  int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
  }

