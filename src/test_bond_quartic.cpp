/* Author: Ryan Houlihan
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
 * BrokenBondRun0
 * Test for case were bond is broken with 2 atoms and 1 bond using pair_style lj_cut
 * Bond Coeff = {"1", "1200.0", "-.55", ".25", "1.3", "34.6878"};
 * Pair Coeff = {"1", "1", "1", "1"};
 * Pair Settings = {"0.5"};
 * 
 * Bond1Run0_1
 * Test for case with 2 atoms and 1 bond using pair_style lj_cut
 * Dont compute pairwise forces by settings rsq > cutoff
 * Bond Coeff = {"1", "1200.0", "-.55", ".25", "3", "34.6878"};
 * Pair Coeff = {"1", "1", "1", "1"};
 * Pair Settings = {"0.05"};
 *
 * Bond2Run0_1
 * Test for case 4 atoms and 2 bond using pair_style lj_cut
 * Compute pairwise forces by settings rsq < cutoff
 * Bond Coeff = {"1", "200.0", "-.55", ".25", "3", "23.6878"};
 * Pair Coeff = {"1", "1", "1", "1"};
 * Pair Settings = {"4"};
 */

#include "test_bond_quartic.h"

using namespace std;

using namespace LAMMPS_NS;

namespace {

  // The fixture for testing class Foo.
  class BondQuarticTest : public ::testing::Test {   
  protected:
    BondQuarticTest() {
      char *argv[] = {"bond_quartic", "-screen", "none", "-log", "none", NULL};
      //char *argv[] = {"bond_quartic", "-echo", "screen", "-log", "none", NULL};
      int narg = 5;
      
      //Create a LAMMPS Object
      lmp = new LAMMPS(narg,argv,MPI_COMM_WORLD);
      // Set atom_style
      lmp->input->one("atom_style bond");
      // Modify atom_style parameters
      lmp->input->one("atom_modify map array sort 10000 3.0");
    }

    virtual ~BondQuarticTest(){
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
 TEST_F(BondQuarticTest, Coeff_CorrectNarg){
    lmp->force->bond = new BondQuartic(lmp);
    char *bcoeff_argv[] = {"1", "100.0", "2.5", "5.0", "3.0", "2.0", NULL};
    int bcoeff_narg; 

    // Check correct input
    bcoeff_narg = 6;
    lmp->force->bond->coeff(bcoeff_narg, bcoeff_argv);
    
    int i;
    for(i = 1; i <=  lmp->atom->nbondtypes; i++){
      EXPECT_DOUBLE_EQ(0.97, lmp->force->bond->equilibrium_distance(i));
      EXPECT_DOUBLE_EQ(1, lmp->force->bond->setflag[i]);
    }
 }

 //////////////////////////////////////////////////////////////////////
 // Multiple check's to make sure bond exits out with incorrect input//
 //////////////////////////////////////////////////////////////////////

  // Check for case with to high a value of bcoeff_narg
 TEST_F(BondQuarticTest, Coeff_highNarg){
    lmp->force->bond = new BondQuartic(lmp);
    char *bcoeff_argv[] = {"1", "100.0", "2.5", "5.0", "3.0", "2.0", NULL};
    int bcoeff_narg; 

    // Check incorrect input
    bcoeff_narg = 7;
    
    // Check that exit(1) was thrown
    ASSERT_EXIT(lmp->force->bond->coeff(bcoeff_narg, bcoeff_argv), ::testing::ExitedWithCode(1), "");  
 }

  // Check for case with to low a value of bcoeff_narg
 TEST_F(BondQuarticTest, Coeff_lowNarg){
    lmp->force->bond = new BondQuartic(lmp);
    char *bcoeff_argv[] = {"1", "100.0", "2.5", "5.0", "3.0", "2.0", NULL};
    int bcoeff_narg; 

    // Check incorrect input
    bcoeff_narg = 1;
    
    // Check that exit(1) was thrown
    ASSERT_EXIT(lmp->force->bond->coeff(bcoeff_narg, bcoeff_argv), ::testing::ExitedWithCode(1), "");  
 }
  // Check for case with for a negative value of bcoeff_narg
  TEST_F(BondQuarticTest,Coeff_NegativeNarg){
    lmp->force->bond = new BondQuartic(lmp);
    char *bcoeff_argv[] = {"1", "100.0", "2.5", "5.0", "3.0", "2.0", NULL};
    int bcoeff_narg; 

    // Check incorrect input
    bcoeff_narg = -6;
    
    // Check that exit(1) was thrown
    ASSERT_EXIT(lmp->force->bond->coeff(bcoeff_narg, bcoeff_argv), ::testing::ExitedWithCode(1), "");  
  }

 //////////////////////////////////////////////////////////////////////
 // Multiple check's to make sure bond compute works properly        //
 //////////////////////////////////////////////////////////////////////

  /* BrokenBondRun0
   * Test for case were bond is broken with 2 atoms and 1 bond using pair_style lj_cut
   * Bond Coeff = {"1", "1200.0", "-.55", ".25", "1.3", "34.6878"};
   * Pair Coeff = {"1", "1", "1", "1"};
   * Pair Settings = {"0.5"};
   */
  TEST_F(BondQuarticTest, BrokenBondRun0) {
    char *lines = strdup("1 1 1 0.0 0.0 -1.0\n 2 1 1 0.0 0.0 1.0\n\n");
    // data_atoms(int n, char*) : unpack n lines from Atom section of data file
    int atomNum = 2; 

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

    // Create bond
    lmp->force->bond = new BondQuartic(lmp);
    char *bcoeff_argv[] = {"1", "1200.0", "-.55", ".25", "1.3", "34.6878", NULL};
    int bcoeff_narg; 
    
    // Check for correct input
    bcoeff_narg = 6;

    lmp->force->bond->coeff(bcoeff_narg, bcoeff_argv);

    // Create pair
    lmp->force->pair = new PairLJCut(lmp);

    char *bsettings_argv[] = {".05", NULL};
    int bsettings_narg = 1;
    bcoeff_argv = {"1", "1", "1", "1", NULL};
    bcoeff_narg = 4;

    lmp->force->pair->settings(bsettings_narg, bsettings_argv);
    lmp->force->pair->coeff(bcoeff_narg, bcoeff_argv);

    lmp->input->one("special_bonds lj 1.0 1.0 1.0");

    // Run test
    lmp->init();
    lmp->update->whichflag=1;
    lmp->update->firststep=0;
    lmp->update->laststep=1;
    lmp->input->one("run 0");

    // Energy for system
    EXPECT_DOUBLE_EQ(0, lmp->force->bond->energy);

    // Forces for atom 1
    EXPECT_DOUBLE_EQ(0, lmp->atom[0].f[0][0]); 
    EXPECT_DOUBLE_EQ(0, lmp->atom[0].f[0][1]);
    EXPECT_DOUBLE_EQ(0, lmp->atom[0].f[0][2]); 

    // Forces for atom 2
    EXPECT_DOUBLE_EQ(0, lmp->atom[0].f[1][0]); 
    EXPECT_DOUBLE_EQ(0, lmp->atom[0].f[1][1]);
    EXPECT_DOUBLE_EQ(0, lmp->atom[0].f[1][2]); 
   
    // Total pressure
    for(int i = 0; i < 6; i++)
      EXPECT_DOUBLE_EQ(0 , lmp->force->bond->virial[i]); 
  }  
  
  /* Bond1Run0_1
   * Test for case with 2 atoms and 1 bond using pair_style lj_cut
   * Dont compute pairwise forces by settings rsq > cutoff
   * Bond Coeff = {"1", "1200.0", "-.55", ".25", "3", "34.6878"};
   * Pair Coeff = {"1", "1", "1", "1"};
   * Pair Settings = {"0.05"};
   */
  TEST_F(BondQuarticTest, Bond1Run0_1) {
     // Variables for test
    double k, b1, b2, rc, u0, cutoff, atomNum, delx, dely, delz;
    double rsq, r, dr, r2, ra, rb, sr2, sr6;
    double fbond, ebond, evdwl, fpair;
    double TWO_1_3 = pow(2.0, (1.0/3.0));

    // Our variables
    k = 1200.0;
    b1 = -0.55;
    b2 = 0.25;
    rc = 3;
    u0 = 34.6878;
    atomNum = 2; 
    cutoff = 0.05;
    delx = 0.0-0.0;
    dely = 0.0-0.0; 
    delz = -1.0 - 1.0; 

    // Energy and force computations
    rsq = delx*delx + dely*dely + delz*delz;
    r = sqrt(rsq);
    dr = r - rc;
    r2 = dr*dr;
    ra = dr-b1;
    rb = dr-b2;
    
    // Compute force and energy for quartic bond
    fbond = -k/r * (r2*(ra+rb) + 2.0*dr*ra*rb);
    if (rsq < TWO_1_3) {
      sr2 = 1.0/rsq;
      sr6 = sr2*sr2*sr2;
      fbond += 48.0*sr6*(sr6-0.5)/rsq;
    }
   
    ebond = k*r2*ra*rb + u0;
    if (rsq < TWO_1_3) ebond += 4.0*sr6*(sr6-1.0) + 1.0;
    
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

    // Create bond
    lmp->force->bond = new BondQuartic(lmp);
    char *bcoeff_argv[] = {"1", "1200.0", "-.55", ".25", "3", "34.6878", NULL};
    int bcoeff_narg; 
    
    // Check for correct input
    bcoeff_narg = 6;

    lmp->force->bond->coeff(bcoeff_narg, bcoeff_argv);

    // Create pair
    lmp->force->pair = new PairLJCut(lmp);

    char *bsettings_argv[] = {"0.05", NULL};
    int bsettings_narg = 1;
    bcoeff_argv = {"1", "1", "1", "1", NULL};
    bcoeff_narg = 4;

    lmp->force->pair->settings(bsettings_narg, bsettings_argv);
    lmp->force->pair->coeff(bcoeff_narg, bcoeff_argv);

    lmp->input->one("special_bonds lj 1.0 1.0 1.0");

    // Run test
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
                          dely*delz*fbond}; 
      
    for(int i = 0; i < 6; i++)
      EXPECT_DOUBLE_EQ(expecVir[i] , lmp->force->bond->virial[i]); 
  }    

  /* Bond1Run0_2
   * Test for case 2 atoms and 1 bond using pair_style lj_cut
   * Compute pairwise forces by settings rsq < cutoff
   * Bond Coeff = {"1", "1200.0", "-.55", ".25", "3", "34.6878"};
   * Pair Coeff = {"1", "1", "1", "1"};
   * Pair Settings = {"5"};
   */
  TEST_F(BondQuarticTest, Bond1Run0_2) {
     // Variables for test
    double k, b1, b2, rc, u0, cutoff, atomNum, delx, dely, delz;
    double rsq, r, dr, r2, ra, rb, sr2, sr6;
    double fbond, ebond, evdwl, fpair;
    double TWO_1_3 = pow(2.0, (1.0/3.0));

    // Our variables
    k = 1200.0;
    b1 = -0.55;
    b2 = 0.25;
    rc = 3;
    u0 = 34.6878;
    cutoff = 5;
    atomNum = 2; 
    delx = 0.0-0.0;
    dely = 0.0-0.0; 
    delz = -1.0 - 1.0; 

    // Energy and force computations
    rsq = delx*delx + dely*dely + delz*delz;
    r = sqrt(rsq);
    dr = r - rc;
    r2 = dr*dr;
    ra = dr-b1;
    rb = dr-b2;
    
    // Compute force and energy for quartic bond
    fbond = -k/r * (r2*(ra+rb) + 2.0*dr*ra*rb);
    if (rsq < TWO_1_3) {
      sr2 = 1.0/rsq;
      sr6 = sr2*sr2*sr2;
      fbond += 48.0*sr6*(sr6-0.5)/rsq;
    }
   
    ebond = k*r2*ra*rb + u0;
    if (rsq < TWO_1_3) ebond += 4.0*sr6*(sr6-1.0) + 1.0;

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

    // Create bond
    lmp->force->bond = new BondQuartic(lmp);
    char *bcoeff_argv[] = {"1", "1200.0", "-.55", ".25", "3", "34.6878", NULL};
    int bcoeff_narg; 
    
    // Check for correct input
    bcoeff_narg = 6;

    lmp->force->bond->coeff(bcoeff_narg, bcoeff_argv);

    // Create pair
    lmp->force->pair = new PairLJCut(lmp);

    char *bsettings_argv[] = {"5", NULL};
    int bsettings_narg = 1;
    bcoeff_argv = {"1", "1", "1", "1", NULL};
    bcoeff_narg = 4;

    lmp->force->pair->settings(bsettings_narg, bsettings_argv);
    lmp->force->pair->coeff(bcoeff_narg, bcoeff_argv);

    lmp->input->one("special_bonds lj 1.0 1.0 1.0");

    // Run test
    lmp->init();
    lmp->update->whichflag=1;
    lmp->update->firststep=0;
    lmp->update->laststep=1;
    lmp->input->one("run 0");
    /*
    // Pair force and  contribution
    fpair = 0;
    if(rsq < pow(cutoff, 2)){
      int i1, i2, itype, jtype;
      i1 = 0;
      i2 = 1;
      itype = 1;
      jtype = 1;
      evdwl = -lmp->force->pair->single(i1,i2,itype,jtype,rsq,1.0,1.0,fpair);
      fpair = -fpair;
    }
    */
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
                          dely*delz*fbond}; 
    for(int i = 0; i < 6; i++)
      EXPECT_DOUBLE_EQ(expecVir[i] , lmp->force->bond->virial[i]); 
  }    

  /* Bond2Run0_1
   * Test for case 4 atoms and 2 bond using pair_style lj_cut
   * Compute pairwise forces by settings rsq < cutoff
   * Bond Coeff = {"1", "200.0", "-.55", ".25", "3", "23.6878"};
   * Pair Coeff = {"1", "1", "1", "1"};
   * Pair Settings = {"4"};
   */
  /*
  TEST_F(BondQuarticTest, Bond2Run0_1) {     
     // Variables for test
    double k, b1, b2, rc, u0, cutoff, atomNum, bondNum;
    double delx1, dely1, delz1, delx2, dely2, delz2;
    double rsq1, r1, dr1, r21, ra1, rb1, sr21, sr61, rsq2, r2, dr2, r22, ra2, rb2, sr22, sr62;
    double fbond1, ebond1, fbond2, ebond2;
    double TWO_1_3 = pow(2.0, (1.0/3.0));

    // Our variables
    k = 200.0;
    b1 = -.55;
    b2 = .25;
    rc = 3;
    u0 = 23.6878;
    cutoff = 4;
    atomNum = 4; 
    bondNum = 2;

    // Bond 1 Computations
    delx1 = 0.0-0.0;
    dely1 = 0.0-0.0; 
    delz1 = -1.0 - 1.0; 

    // Energy and force computations
    rsq1 = delx1*delx1 + dely1*dely1 + delz1*delz1;
    r1 = sqrt(rsq1);
    dr1 = r1 - rc;
    r21 = dr1*dr1;
    ra1 = dr1-b1;
    rb1 = dr1-b2;
    
    // Compute force and energy for quartic bond
    fbond1 = -k/r1 * (r21*(ra1+rb1) + 2.0*dr1*ra1*rb1);
    fprintf(stderr, " fbond1 = %f \n ", fbond1);
    if (rsq1 < TWO_1_3) {
      sr21 = 1.0/rsq1;
      sr61 = sr21*sr21*sr21;
      fbond1 += 48.0*sr61*(sr61-0.5)/rsq1;
    }
   
    ebond1 = k*r21*ra1*rb1 + u0;
    if (rsq1 < TWO_1_3) ebond1 += 4.0*sr61*(sr61-1.0) + 1.0;

    // Bond 2 Computations
    delx2 = 0.0-0.0;
    dely2 = -1.0-1.0; 
    delz2 = 0.0 - 0.0;  

    // Energy and force computations
    rsq2 = delx2*delx2 + dely2*dely2 + delz2*delz2;
    r2 = sqrt(rsq2);
    dr2 = r2 - rc;
    r22 = dr2*dr2;
    ra2 = dr2-b1;
    rb2 = dr2-b2;
    
    // Compute force and energy for quartic bond
    fbond2 = -k/r2 * (r22*(ra2+rb2) + 2.0*dr2*ra2*rb2);
    if (rsq2 < TWO_1_3) {
      sr22 = 1.0/rsq2;
      sr62 = sr22*sr22*sr22;
      fbond2 += 48.0*sr62*(sr62-0.5)/rsq2;
    }
   
    ebond2 = k*r22*ra2*rb2 + u0;
    if (rsq2 < TWO_1_3) ebond2 += 4.0*sr62*(sr62-1.0) + 1.0;
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
    fprintf(stderr,"YEADHHH");
    lmp->input->one("mass 1 1.0");

    // End create bonds and atoms
    //*********************************************  

    // Create bond
    fprintf(stderr,"YEADHHH");
    lmp->force->bond = new BondQuartic(lmp);
    char *bcoeff_argv[] = {"1", "200.0", "-.55", ".25", "3", "23.6878", NULL};
    int bcoeff_narg; 
    
    // Check for correct input
    bcoeff_narg = 6;
    lmp->force->bond->coeff(bcoeff_narg, bcoeff_argv);

    // Create pair
    lmp->force->pair = new PairLJCut(lmp);

    char *bsettings_argv[] = {"4", NULL};
    int bsettings_narg = 1;
    bcoeff_argv = {"1", "1", "1", "1", NULL};
    bcoeff_narg = 4;

    lmp->force->pair->settings(bsettings_narg, bsettings_argv);
    lmp->force->pair->coeff(bcoeff_narg, bcoeff_argv);

    lmp->input->one("special_bonds lj 1.0 1.0 1.0");

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
    for(int i = 0; i < 6; i++)
      EXPECT_DOUBLE_EQ(expecVir[i] , lmp->force->bond->virial[i]);
  } 
*/
} // namespace
  int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
  }

