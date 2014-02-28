/* Author: Ryan Houlihan
 * BondHybridTest():
 *
 * The following code tests bond_style bond_harmonic in various ways:
 * 
 * Settings_CorrectArguments: Checks valid settings input
 * Settings_BadArguments: Checks invalid settings input
 * 
 * Coeff_Correct: Check coeff with corrent input
 * Coeff_RandomTests: Check coeff for various input
 *
 * Tests the Force, Energy, and Pressure of each of the following simulations: 
 *
 * Style1Bond1Run0_1:
 * Test basic case with 2 atoms & 1 bond and hybrid coeff={"harmonic"};
 * Harmonic coeff_argv{1, k = 132.0, r_0 = 2.0}
 *
 * Bond1Run0_1:
 * Test basic case with 4 atoms & 2 bond and hybrid coeff={"harmonic", "fene"};
 * Harmonic coeff_argv{1, k = 132.0, r_0 = 2.0}
 * Fene coeff={"1", "30.0", "4.0", "3.0", "2.0"};
 *
 * Bond1Run0_2:
 * Test basic case with 4 atoms & 2 bond and hybrid coeff={"nonlinear", "morse"};
 * nonlinear coeff_argv{1, epsilon = 132.0, r0 = 3.0, lamda = 2.0}
 * morse coeff_argv{1, d0 = 132.0, alpha = 3.0, r_0 = 2.0}
 *
 */

#include "test_bond_hybrid.h"

using namespace std;

using namespace LAMMPS_NS;

namespace {

  // The fixture for testing class BondHybridTest
  class BondHybridTest : public ::testing::Test {   
  protected:
    BondHybridTest() {
      char *argv[] = {"bond_hybrid", "-screen", "none", "-log", "none", NULL};
      //char *argv[] = {"bond_harmonic", "-echo", "screen", "-log", "none", NULL};
      int narg = 5;
      
      //Create a LAMMPS Object
      lmp = new LAMMPS(narg,argv,MPI_COMM_WORLD);
      // Set atom_style
      lmp->input->one("atom_style bond");
      // Modify atom_style parameters
      lmp->input->one("atom_modify map array sort 10000 3.0");
    }

    virtual ~BondHybridTest(){
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
 TEST_F(BondHybridTest, Settings_CorrectArguments){
    lmp->force->bond = new BondHybrid(lmp);
    char *bsettings_argv[] = {"harmonic", "morse", NULL};
    int bsettings_narg; 

    // Check correct input
    bsettings_narg = 2;
    lmp->force->bond->settings(bsettings_narg, bsettings_argv);  
 }

  // Check for various incorrect arguments settings()
  TEST_F(BondHybridTest, Settings_BadArguments){
    lmp->force->bond = new BondHybrid(lmp);
    char *bsettings_argv[] =  {NULL, NULL, NULL};
    int bsettings_narg; 

    // Check for case with incorrect "hybrid" arguments for bsettings_argv
    bsettings_argv = {"hybrid", "morse", NULL};
    bsettings_narg = 2;
    ASSERT_EXIT(lmp->force->bond->settings(bsettings_narg, bsettings_argv), ::testing::ExitedWithCode(1), "");  
       
    // Check for case with incorrect value of bsettings_narg
    bsettings_argv = {"harmonic", "morse", NULL};
    bsettings_narg = -2;
    ASSERT_EXIT(lmp->force->bond->settings(bsettings_narg, bsettings_argv), ::testing::ExitedWithCode(1), "");  

    // Check for case with incorrect value of bsettings_narg
    bsettings_argv = {"fene", "harmonic", NULL};
    bsettings_narg = 0;
    ASSERT_EXIT(lmp->force->bond->settings(bsettings_narg, bsettings_argv), ::testing::ExitedWithCode(1), "");  

    // Check for case with incorrect double arguments for bsettings_argv
    bsettings_argv = {"morse", "morse", NULL};
    bsettings_narg = 2;
    ASSERT_EXIT(lmp->force->bond->settings(bsettings_narg, bsettings_argv), ::testing::ExitedWithCode(1), "");  

    // Check for case with "none" incorrect arguments for bsettings_argv
    bsettings_argv = {"none", "nonlinear", NULL};
    bsettings_narg = 2;
    ASSERT_EXIT(lmp->force->bond->settings(bsettings_narg, bsettings_argv), ::testing::ExitedWithCode(1), "");  
 }

  // Still to be done
  // Test of the private char** variable bond->keywords for correct bond styles 

 //////////////////////////////////////////////////////////////////////
 // Multiple check's to make sure bond exits out with incorrect input//
 //////////////////////////////////////////////////////////////////////

  // Check for case correct coeff arguments
 TEST_F(BondHybridTest, Coeff_Correct){
    lmp->force->bond = new BondHybrid(lmp);
    char *bsettings_argv[] = {"harmonic", "morse", NULL};
    int bsettings_narg; 

    // Check correct input
    bsettings_narg = 2;
    lmp->force->bond->settings(bsettings_narg, bsettings_argv);  

    char *bcoeff_argv[] = {"1", "morse", "100.0", "2.5", "5.0", NULL};
    int bcoeff_narg; 

    // Check incorrect input
    bcoeff_narg = 5;
    
    // Check that exit(1) was thrown
    lmp->force->bond->coeff(bcoeff_narg, bcoeff_argv);

    bcoeff_argv = {"1", "harmonic", "20.0", "2.5", NULL};
    bcoeff_narg; 

    // Check incorrect input
    bcoeff_narg = 4;
    
    // Check that exit(1) was thrown
    lmp->force->bond->coeff(bcoeff_narg, bcoeff_argv);
 }

  // Check for various coeff arguments
 TEST_F(BondHybridTest, Coeff_RandomTests){
    lmp->force->bond = new BondHybrid(lmp);
    char *bsettings_argv[] = {"harmonic", "morse", "fene", NULL};
    int bsettings_narg; 

    // Check correct input
    bsettings_narg = 2;
    lmp->force->bond->settings(bsettings_narg, bsettings_argv);  

    char *bcoeff_argv[] = {"1", "morse", "100.0", "2.5", "5.0", NULL};
    int bcoeff_narg; 

    // Check incorrect input number
    bcoeff_argv = {"1", "morse", "100.0", "2.5", "5.0", NULL};
    bcoeff_narg = 4;    
    ASSERT_EXIT(lmp->force->bond->coeff(bcoeff_narg, bcoeff_argv), ::testing::ExitedWithCode(1), "");

    // Check incorrect input number
    bcoeff_argv = {"1", "fene", "100.0", "2.5", "5.0", NULL};
    bcoeff_narg = 5;    
    ASSERT_EXIT(lmp->force->bond->coeff(bcoeff_narg, bcoeff_argv), ::testing::ExitedWithCode(1), "");

   // Check wrong bond style
    bcoeff_argv = {"1", "nonlinear", "100.0", "2.5", "5.0", NULL};
    bcoeff_narg = 5;    
    ASSERT_EXIT(lmp->force->bond->coeff(bcoeff_narg, bcoeff_argv), ::testing::ExitedWithCode(1), "");

    // Check none bond style
    bcoeff_argv = {"1", "none", NULL};
    bcoeff_narg = 2;    
    lmp->force->bond->coeff(bcoeff_narg, bcoeff_argv);

    // Check correct Test
    bcoeff_argv = {"1", "harmonic", "20.0", "2.5", NULL};
    bcoeff_narg = 4;
    lmp->force->bond->coeff(bcoeff_narg, bcoeff_argv);
}

 //////////////////////////////////////////////////////////////////////
 // Multiple check's to make sure bond compute works properly        //
 //////////////////////////////////////////////////////////////////////

  /* Style1Bond1Run0_1:
   * Test basic case with 2 atoms & 1 bond and hybrid coeff={"harmonic"};
   * Harmonic coeff_argv{1, k = 132.0, r_0 = 2.0}
   */
  TEST_F(BondHybridTest, Style1Bond1Run0_1) {
    double delx = 0.0-0.0;
    double dely = 0.0-0.0; 
    double delz = -1.0 - 1.0; 
    double atomNum = 2; 
    double rsq = delx*delx + dely*dely + delz*delz;

    // Harmonic Bond
    double k = 132.0;
    double r0 = 2.0;
    
    double r = sqrt(rsq);
    double dr = r - r0;
    double rk = k*dr; 
    double fbond = -2.0*(rk/r);
    double ebond = rk*dr;
    // End Harmonic Bond

    lmp->force->bond = new BondHybrid(lmp);
    char *bsettings_argv[] = {"harmonic", NULL};
    int bsettings_narg; 

    // Check correct input
    bsettings_narg = 1;

    lmp->force->bond->settings(bsettings_narg, bsettings_argv);  

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

    char *harcoeff_argv[] = {"1", "harmonic", "132.0", "2", NULL};
    int harcoeff_narg; 
    
    harcoeff_narg = 4;
    // Check for correct input
    lmp->force->bond->coeff(harcoeff_narg, harcoeff_argv);

    // Tests that BondHybrid::energy works properly 
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

  /* Bond1Run0_1:
   * Test basic case with 4 atoms & 2 bonds and hybrid coeff={"harmonic", "fene"};
   * Harmonic coeff_argv{1, k = 132.0, r_0 = 2.0}
   * Fene coeff={"1", "30.0", "4.0", "3.0", "2.0"};
   */
  TEST_F(BondHybridTest, Bond1Run0_1) {
      lmp->atom->ntypes = 2;
      lmp->atom->nbonds = 2;
      lmp->atom->nbondtypes = 2;

      double atomNum = 4; 
      double bondNum = 2;

      // Harmonic Bond
      double delx1 = 0.0-0.0;
      double dely1 = 0.0-0.0; 
      double delz1 = -1.0 - 1.0; 
      double H_k = 132.0;
      double H_r0 = 2.0;

      double rsq1 = delx1*delx1 + dely1*dely1 + delz1*delz1;
      double r = sqrt(rsq1);
      double dr = r - H_r0;
      double rk = H_k*dr; 
      double fbond1 = -2.0*(rk/r);
      double ebond1 = rk*dr;
      // End Harmonic Bond

      // Fene Bond
      double delx2 = 0.0-0.0;
      double dely2 = -1.0-1.0; 
      double delz2 = 0.0 - 0.0; 

      double sigma = 2.0;
      double epsilon = 3.0;
      double F_r0 = 4.0;
      double F_k = 30.0;
      double TWO_1_3 = pow(2.0,(1.0/3.0)); 

      double rsq2 = delx2*delx2 + dely2*dely2 + delz2*delz2;
      double r0sq = F_r0 * F_r0;
      double rlogarg = 1.0-rsq2/r0sq;

      double fbond2 = -F_k/rlogarg;
      double ebond2 = -0.5 * F_k*r0sq*log(rlogarg);
      if (rsq2 < TWO_1_3*sigma*sigma) {
	double sr2 = (sigma*sigma)/rsq2;
	double sr6 = sr2*sr2*sr2;
	ebond2 += 4.0*epsilon*sr6*(sr6-1.0) + epsilon;    
	fbond2 += 48.0*epsilon*sr6*(sr6-0.5)/rsq2;
      }
      // End Fene Bond
    
      lmp->force->bond = new BondHybrid(lmp);
      char *bsettings_argv[] = {"harmonic", "fene", NULL};
      int bsettings_narg; 

      // Check correct input
      bsettings_narg = 2;

      lmp->force->bond->settings(bsettings_narg, bsettings_argv);  
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
      lines = strdup("1 2 3 4\n");
      lmp->atom->data_bonds(1,lines);
      free(lines);

    lmp->input->one("mass * 1.0");

    // End create bonds and atoms
    //*********************************************  

    char *harcoeff_argv[] = {"1", "harmonic", "132.0", "2", NULL};
    int harcoeff_narg; 
    
    harcoeff_narg = 4;
    // Check for correct input

    lmp->force->bond->coeff(harcoeff_narg, harcoeff_argv);

    char *fenecoeff_argv[] = {"2", "fene", "30.0", "4.0", "3.0", "2.0", NULL};
    int fenecoeff_narg; 
    
    fenecoeff_narg = 6;
    // Check for correct input
    lmp->force->bond->coeff(fenecoeff_narg, fenecoeff_argv);

   
    // Tests that BondHybrid::energy works properly 
    lmp->init();
    lmp->update->whichflag=1;
    lmp->update->firststep=0;
    lmp->update->laststep=1;
    lmp->input->one("run 0");

    EXPECT_DOUBLE_EQ(ebond1 + ebond2, lmp->force->bond->energy); 

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
    for(int i = 0; i < 6; i++)
      EXPECT_DOUBLE_EQ(expecVir[i] , lmp->force->bond->virial[i]);
  }
  /* Bond1Run0_2:
   * Test basic case with 2 atoms & 1 bond and hybrid coeff={"nonlinear", "morse"};
   * nonlinear coeff_argv{1, epsilon = 132.0, r0 = 3.0, lamda = 2.0}
   * morse coeff_argv{1, d0 = 132.0, alpha = 3.0, r_0 = 2.0}
   */
  TEST_F(BondHybridTest, Bond1Run0_2) {
      lmp->atom->ntypes = 2;
      lmp->atom->nbonds = 2;
      lmp->atom->nbondtypes = 2;

      double atomNum = 4; 
      double bondNum = 2;

    // nonlinear variables for test
    double delx1 = 0.0-0.0;
    double dely1 = 0.0-0.0; 
    double delz1 = -1.0 - 1.0; 

    double epsilon = 132.0;
    double N_r0 = 3.0;
    double lamda = 2.0;

    double rsq1 = delx1*delx1 + dely1*dely1 + delz1*delz1;
    double r1 = sqrt(rsq1);
    double N_dr = r1 - N_r0;
    double drsq = N_dr*N_dr;
    double lamdasq = lamda*lamda;
    double denom = lamdasq - drsq;
    double denomsq = denom*denom;
    double fbond1 = -epsilon/r1* 2.0*N_dr*lamdasq/denomsq;
    double ebond1 = epsilon * drsq / denom;
    // End nonlinear variables

    // morse variables for test
    double delx2 = 0.0-0.0;
    double dely2 = -1.0-1.0; 
    double delz2 = 0.0 - 0.0; 

    double d0 = 132.0;
    double alpha = 3.0;
    double M_r0 = 2.0;
    double rsq2 = delx2*delx2 + dely2*dely2 + delz2*delz2;
    double r2 = sqrt(rsq2);

    double M_dr = r2- M_r0;
    double ralpha = exp(-alpha*M_dr);
    double fbond2 = -2.0*d0*alpha*(1-ralpha)*ralpha/r2; 
    double ebond2 = d0*(1-ralpha)*(1-ralpha); 
    // End morse variables

    lmp->force->bond = new BondHybrid(lmp);
    char *bsettings_argv[] = {"nonlinear", "morse", NULL};
    int bsettings_narg;

    // Check correct input
    bsettings_narg = 2;

    lmp->force->bond->settings(bsettings_narg, bsettings_argv);  
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
    lines = strdup("1 2 3 4\n");
    lmp->atom->data_bonds(1,lines);
    free(lines);

    lmp->input->one("mass * 1.0");
    // End create bonds and atoms
    //*********************************************  
    char *N_coeff_argv[] = {"1","nonlinear", "132.0", "3.0", "2.0", NULL};
    int N_coeff_narg; 
    
    N_coeff_narg = 5;
    // Check for correct input
    lmp->force->bond->coeff(N_coeff_narg, N_coeff_argv);
    char *M_coeff_argv[] = {"2", "morse", "132.0", "3.0", "2.0", NULL};
    int M_coeff_narg; 
    
    M_coeff_narg = 5;
    // Check for correct input
    lmp->force->bond->coeff(M_coeff_narg, M_coeff_argv);
    // Tests that BondHybrid::energy works properly 
    lmp->init();
    lmp->update->whichflag=1;
    lmp->update->firststep=0;
    lmp->update->laststep=1;
    lmp->input->one("run 0");

    EXPECT_DOUBLE_EQ(ebond1 + ebond2, lmp->force->bond->energy); 

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
    for(int i = 0; i < 6; i++)
      EXPECT_DOUBLE_EQ(expecVir[i] , lmp->force->bond->virial[i]);
  }    
} // namespace
  int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
  }

