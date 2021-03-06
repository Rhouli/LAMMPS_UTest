/**
 * Author: Author: Ryan Houlihan
 * This tests the functionality of a bond with a class2 bond_style. 
 */

#include "test_bond_class2.h"

using namespace std;

using namespace LAMMPS_NS;

namespace {

  // The fixture for testing class Foo.
  class BondClass2Test : public ::testing::Test {   
  protected:
    BondClass2Test() {
      char *argv[] = {"bond_class2", "-screen", "none", "-log", "none", NULL};
      int narg = 5;
      
      //Create a LAMMPS Object
      lmp = new LAMMPS(narg,argv,MPI_COMM_WORLD);
      // Set atom_style
      lmp->input->one("atom_style bond");
      // Modify atom_style parameters
      lmp->input->one("atom_modify map array sort 10000 3.0");
    }

    virtual ~BondClass2Test(){
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
 TEST_F(BondClass2Test, Coeff_CorrectNarg){
    lmp->force->bond = new BondClass2(lmp);
    char *bcoeff_argv[] = {"1", "100.0", "2.5", "5.0", "1.0"};
    int bcoeff_narg; 

    double r0 = 100.0;
    // Check correct input
    bcoeff_narg = 5;
    lmp->force->bond->coeff(bcoeff_narg, bcoeff_argv);
    
    int i;
    for(i = 1; i <=  lmp->atom->nbondtypes; i++){
      EXPECT_DOUBLE_EQ(r0, lmp->force->bond->equilibrium_distance(i));
      EXPECT_DOUBLE_EQ(1, lmp->force->bond->setflag[i]);
    }
 }

 //////////////////////////////////////////////////////////////////////
 // Multiple check's to make sure bond exits out with incorrect input//
 //////////////////////////////////////////////////////////////////////

  // Check for case with to high a value of bcoeff_narg
 TEST_F(BondClass2Test, Coeff_highNarg){
    lmp->force->bond = new BondClass2(lmp);
    char *bcoeff_argv[] = {"1", "100.0", "2.5", "5.0", "10", "2", "7"};
    int bcoeff_narg; 

    // Check incorrect input
    bcoeff_narg = 7;
    
    // Check that exit(1) was thrown
    ASSERT_EXIT(lmp->force->bond->coeff(bcoeff_narg, bcoeff_argv), ::testing::ExitedWithCode(1), "");  
 }

  // Check for case with to low a value of bcoeff_narg
 TEST_F(BondClass2Test, Coeff_lowNarg){
    lmp->force->bond = new BondClass2(lmp);
    char *bcoeff_argv[] = {"1", "100.0", "2.5", "5.0", NULL};
    int bcoeff_narg; 

    // Check incorrect input
    bcoeff_narg = 1;
    
    // Check that exit(1) was thrown
    ASSERT_EXIT(lmp->force->bond->coeff(bcoeff_narg, bcoeff_argv), ::testing::ExitedWithCode(1), "");  
 }

  // Check for case with for a negative value of bcoeff_narg
  TEST_F(BondClass2Test,Coeff_NegativeNarg){
    lmp->force->bond = new BondClass2(lmp);
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

 //{"1", R0"10.0", K1"10.0", K2"10.0", K3"10.0"}
  TEST_F(BondClass2Test, Bond1Run0_1) {
    // Variables for test
    double atomNum = 2.0;
    double fbond, ebond =0.0;
    double delx = 0.0-0.0;
    double dely = 0.0-0.0;
    double delz = -1.0-1.0;
    double r0 = 10.0;
    double k2 = 10.0;
    double k3 = 10.0;
    double k4 = 10.0;
    double rsq, r, dr, dr2, dr3, dr4;

    rsq = delx*delx + dely*dely + delz*delz;
    r = sqrt(rsq);
    dr = r - r0;
    dr2 = dr*dr;
    dr3 = dr2*dr;
    dr4 = dr3*dr;
	

	//{"1", "10.0", "10.0", "10.0", "10.0"}
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

    lmp->force->bond = new BondClass2(lmp);
    char *bcoeff_argv[] = {"1", "10.0", "10.0", "10.0", "10.0"};
    int bcoeff_narg; 
    
    // Check for correct input
    bcoeff_narg = 5;
    lmp->force->bond->coeff(bcoeff_narg, bcoeff_argv);

    int i;
    for(i = 1; i <=  lmp->atom->nbondtypes; i++){
      EXPECT_DOUBLE_EQ(r0, lmp->force->bond->equilibrium_distance(i));
    }		


    double de_bond = 2.0*k2*dr + 3.0*k3*dr2 + 4.0*k4*dr3;
    if (r > 0.0) fbond = -de_bond/r;
    else fbond = 0.0;

    ebond = k2*dr2 + k3*dr3 + k4*dr4;


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
  TEST_F(BondClass2Test, Bond1Run0_2) {
        // Variables for test
    double atomNum = 2.0;
    double fbond, ebond =0.0;
    double delx = 0.0-0.0;
    double dely = 0.0-0.0;
    double delz = -1.0-1.0;
    double r0 = 12.0;
    double k2 = 14.0;
    double k3 = 15.0;
    double k4 = 16.0;
    double rsq, r, dr, dr2, dr3, dr4;

    rsq = delx*delx + dely*dely + delz*delz;
    r = sqrt(rsq);
    dr = r - r0;
    dr2 = dr*dr;
    dr3 = dr2*dr;
    dr4 = dr3*dr;
	

	//{"1", "10.0", "10.0", "10.0", "10.0"}
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

    lmp->force->bond = new BondClass2(lmp);
    char *bcoeff_argv[] = {"1", "12.0", "14.0", "15.0", "16.0"};
    int bcoeff_narg; 
    
    // Check for correct input
    bcoeff_narg = 5;
    lmp->force->bond->coeff(bcoeff_narg, bcoeff_argv);

    int i;
    for(i = 1; i <=  lmp->atom->nbondtypes; i++){
      EXPECT_DOUBLE_EQ(r0, lmp->force->bond->equilibrium_distance(i));
    }		


    double de_bond = 2.0*k2*dr + 3.0*k3*dr2 + 4.0*k4*dr3;
    if (r > 0.0) fbond = -de_bond/r;
    else fbond = 0.0;

    ebond = k2*dr2 + k3*dr3 + k4*dr4;


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
  TEST_F(BondClass2Test, Bond2Run0_1) {     
    // Variables for run 
    double atomNum = 4.0;
    double bondNum = 2.0;
    double ebond1, ebond2;
    double fbond1, fbond2;


    double r0 = 12.0;
    double k2 = 14.0;
    double k3 = 15.0;
    double k4 = 16.0;
    double rsq1, r1, dr_1, dr2_1, dr3_1, dr4_1;
    double rsq2, r2, dr_2, dr2_2, dr3_2, dr4_2;

    // Variables for bond between atom 1 and 2
    double delx1 = 0.0-0.0; 
    double dely1 = 0.0-0.0; 
    double delz1 = -1.0 - 1.0; 
    rsq1 = delx1*delx1 + dely1*dely1 + delz1*delz1;
    r1 = sqrt(rsq1);
    dr_1 = r1 - r0;
    dr2_1 = dr_1*dr_1;
    dr3_1 = dr2_1*dr_1;
    dr4_1 = dr3_1*dr_1;

    double de_bond = 2.0*k2*dr_1 + 3.0*k3*dr2_1 + 4.0*k4*dr3_1;
    if (r1 > 0.0) fbond1 = -de_bond/r1;
    else fbond1 = 0.0;

    ebond1 = k2*dr2_1 + k3*dr3_1 + k4*dr4_1;



    // Variables for bond between atom 3 and 4
    double delx2 = 0.0-0.0;
    double dely2 = -1.0-1.0; 
    double delz2 = 0.0 - 0.0; 
    rsq2 = delx2*delx2 + dely2*dely2 + delz2*delz2;
    r2 = sqrt(rsq2);
    dr_2 = r2 - r0;
    dr2_2 = dr_2*dr_2;
    dr3_2 = dr2_2*dr_2;
    dr4_2 = dr3_2*dr_2;

    de_bond = 2.0*k2*dr_2 + 3.0*k3*dr2_2 + 4.0*k4*dr3_2;
    if (r2 > 0.0) fbond2 = -de_bond/r2;
    else fbond2 = 0.0;

    ebond2 = k2*dr2_2 + k3*dr3_2 + k4*dr4_2;
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

    lmp->force->bond = new BondClass2(lmp);
    //  k = 100.0 and r0 = 1
    char *bcoeff_argv[] = {"1", "12.0", "14.0", "15.0", "16.0"};
    int bcoeff_narg; 
    
    // Check for correct input
    bcoeff_narg = 5;
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


  // Test basic case with 4 atoms & 2 bonds
  // Test input == {R0: 14.0, k2: 16.0, k3: 
  TEST_F(BondClass2Test, Bond2Run0_2) 
  {     
  	// Variables for run 
    double atomNum = 4.0;
    double bondNum = 2.0;
    double ebond1, ebond2;
    double fbond1, fbond2;


    double r0 = 14.0;
    double k2 = 16.0;
    double k3 = 18.0;
    double k4 = 20.0;
    double rsq1, r1, dr_1, dr2_1, dr3_1, dr4_1;
    double rsq2, r2, dr_2, dr2_2, dr3_2, dr4_2;

    // Variables for bond between atom 1 and 2
    double delx1 = 0.0-0.0; 
    double dely1 = 0.0-0.0; 
    double delz1 = -1.0 - 1.0; 
    rsq1 = delx1*delx1 + dely1*dely1 + delz1*delz1;
    r1 = sqrt(rsq1);
    dr_1 = r1 - r0;
    dr2_1 = dr_1*dr_1;
    dr3_1 = dr2_1*dr_1;
    dr4_1 = dr3_1*dr_1;

    double de_bond = 2.0*k2*dr_1 + 3.0*k3*dr2_1 + 4.0*k4*dr3_1;
    if (r1 > 0.0) fbond1 = -de_bond/r1;
    else fbond1 = 0.0;

    ebond1 = k2*dr2_1 + k3*dr3_1 + k4*dr4_1;



    // Variables for bond between atom 3 and 4
    double delx2 = 0.0-0.0;
    double dely2 = -1.0-1.0; 
    double delz2 = 0.0 - 0.0; 
    rsq2 = delx2*delx2 + dely2*dely2 + delz2*delz2;
    r2 = sqrt(rsq2);
    dr_2 = r2 - r0;
    dr2_2 = dr_2*dr_2;
    dr3_2 = dr2_2*dr_2;
    dr4_2 = dr3_2*dr_2;

    de_bond = 2.0*k2*dr_2 + 3.0*k3*dr2_2 + 4.0*k4*dr3_2;
    if (r2 > 0.0) fbond2 = -de_bond/r2;
    else fbond2 = 0.0;

    ebond2 = k2*dr2_2 + k3*dr3_2 + k4*dr4_2;
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

    lmp->force->bond = new BondClass2(lmp);
    //  k = 100.0 and r0 = 1
    char *bcoeff_argv[] = {"1", "14.0", "16.0", "18.0", "20.0"};
    int bcoeff_narg; 
    
    // Check for correct input
    bcoeff_narg = 5;
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
  TEST_F(BondClass2Test, Bond3Run0_1) {     
    // Variables for run 
    double atomNum = 4.0;
    double bondNum = 2.0;
    double ebond1, ebond2, ebond3;
    double fbond1, fbond2, fbond3;


    double r0 = 14.0;
    double k2 = 16.0;
    double k3 = 18.0;
    double k4 = 20.0;
    double rsq1, r1, dr_1, dr2_1, dr3_1, dr4_1;
    double rsq2, r2, dr_2, dr2_2, dr3_2, dr4_2;
    double rsq3, r3, dr_3, dr2_3, dr3_3, dr4_3;

    // Variables for bond between atom 1 and 2
    double delx1 = 0.0-0.0; 
    double dely1 = 0.0-0.0; 
    double delz1 = -1.0 - 1.0; 
    rsq1 = delx1*delx1 + dely1*dely1 + delz1*delz1;
    r1 = sqrt(rsq1);
    dr_1 = r1 - r0;
    dr2_1 = dr_1*dr_1;
    dr3_1 = dr2_1*dr_1;
    dr4_1 = dr3_1*dr_1;

    double de_bond = 2.0*k2*dr_1 + 3.0*k3*dr2_1 + 4.0*k4*dr3_1;
    if (r1 > 0.0) fbond1 = -de_bond/r1;
    else fbond1 = 0.0;

    ebond1 = k2*dr2_1 + k3*dr3_1 + k4*dr4_1;

    // Variables for bond between atom 3 and 4
    double delx2 = 0.0-0.0;
    double dely2 = -1.0-1.0; 
    double delz2 = 0.0 - 0.0; 
    rsq2 = delx2*delx2 + dely2*dely2 + delz2*delz2;
    r2 = sqrt(rsq2);
    dr_2 = r2 - r0;
    dr2_2 = dr_2*dr_2;
    dr3_2 = dr2_2*dr_2;
    dr4_2 = dr3_2*dr_2;

    de_bond = 2.0*k2*dr_2 + 3.0*k3*dr2_2 + 4.0*k4*dr3_2;
    if (r2 > 0.0) fbond2 = -de_bond/r2;
    else fbond2 = 0.0;

    ebond2 = k2*dr2_2 + k3*dr3_2 + k4*dr4_2;

    // Variables for bond between atom 1 and 4
    double delx3 = 0.0-0.0;
    double dely3 = 0.0-1.0; 
    double delz3 = -1.0 - 0.0; 
    rsq3 = delx3*delx3 + dely3*dely3 + delz3*delz3;
    r3 = sqrt(rsq3);
    dr_3 = r3 - r0;
    dr2_3 = dr_3*dr_3;
    dr3_3 = dr2_3*dr_3;
    dr4_3 = dr3_3*dr_3;  

    de_bond = 2.0*k2*dr_3 + 3.0*k3*dr2_3 + 4.0*k4*dr3_3;
    if (r3 > 0.0) fbond3 = -de_bond/r3;
    else fbond3 = 0.0;

    ebond3 = k2*dr2_3 + k3*dr3_3 + k4*dr4_3;

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

    lmp->force->bond = new BondClass2(lmp);
    //  k = 100.0 and r0 = 1
    char *bcoeff_argv[] = {"1", "14.0", "16.0", "18.0", "20.0"};
    int bcoeff_narg; 
    
    // Check for correct input
    bcoeff_narg = 5;
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
    EXPECT_DOUBLE_EQ((atomNum/2)*(ebond1*0.5)+(atomNum/2)*(ebond2*0.5)+(atomNum/2)*(ebond3*0.5), lmp->force->bond->energy); 

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


  TEST_F(BondClass2Test, Bond3Run0_2) {     
        // Variables for run 
    double atomNum = 4.0;
    double bondNum = 2.0;
    double ebond1, ebond2, ebond3;
    double fbond1, fbond2, fbond3;


    double r0 = 10.0;
    double k2 = 13.0;
    double k3 = 16.0;
    double k4 = 19.0;
    double rsq1, r1, dr_1, dr2_1, dr3_1, dr4_1;
    double rsq2, r2, dr_2, dr2_2, dr3_2, dr4_2;
    double rsq3, r3, dr_3, dr2_3, dr3_3, dr4_3;

    // Variables for bond between atom 1 and 2
    double delx1 = 0.0-0.0; 
    double dely1 = 0.0-0.0; 
    double delz1 = -1.0 - 1.0; 
    rsq1 = delx1*delx1 + dely1*dely1 + delz1*delz1;
    r1 = sqrt(rsq1);
    dr_1 = r1 - r0;
    dr2_1 = dr_1*dr_1;
    dr3_1 = dr2_1*dr_1;
    dr4_1 = dr3_1*dr_1;

    double de_bond = 2.0*k2*dr_1 + 3.0*k3*dr2_1 + 4.0*k4*dr3_1;
    if (r1 > 0.0) fbond1 = -de_bond/r1;
    else fbond1 = 0.0;

    ebond1 = k2*dr2_1 + k3*dr3_1 + k4*dr4_1;

    // Variables for bond between atom 3 and 4
    double delx2 = 0.0-0.0;
    double dely2 = -1.0-1.0; 
    double delz2 = 0.0 - 0.0; 
    rsq2 = delx2*delx2 + dely2*dely2 + delz2*delz2;
    r2 = sqrt(rsq2);
    dr_2 = r2 - r0;
    dr2_2 = dr_2*dr_2;
    dr3_2 = dr2_2*dr_2;
    dr4_2 = dr3_2*dr_2;

    de_bond = 2.0*k2*dr_2 + 3.0*k3*dr2_2 + 4.0*k4*dr3_2;
    if (r2 > 0.0) fbond2 = -de_bond/r2;
    else fbond2 = 0.0;

    ebond2 = k2*dr2_2 + k3*dr3_2 + k4*dr4_2;

    // Variables for bond between atom 1 and 4
    double delx3 = 0.0-0.0;
    double dely3 = 0.0-1.0; 
    double delz3 = -1.0 - 0.0; 
    rsq3 = delx3*delx3 + dely3*dely3 + delz3*delz3;
    r3 = sqrt(rsq3);
    dr_3 = r3 - r0;
    dr2_3 = dr_3*dr_3;
    dr3_3 = dr2_3*dr_3;
    dr4_3 = dr3_3*dr_3;  

    de_bond = 2.0*k2*dr_3 + 3.0*k3*dr2_3 + 4.0*k4*dr3_3;
    if (r3 > 0.0) fbond3 = -de_bond/r3;
    else fbond3 = 0.0;

    ebond3 = k2*dr2_3 + k3*dr3_3 + k4*dr4_3;

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

    lmp->force->bond = new BondClass2(lmp);
    //  k = 100.0 and r0 = 1
    char *bcoeff_argv[] = {"1", "10.0", "13.0", "16.0", "19.0"};
    int bcoeff_narg; 
    
    // Check for correct input
    bcoeff_narg = 5;
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
    EXPECT_DOUBLE_EQ((atomNum/2)*(ebond1*0.5)+(atomNum/2)*(ebond2*0.5)+(atomNum/2)*(ebond3*0.5), lmp->force->bond->energy); 

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

