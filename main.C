//##############################################################
/**\file main.C

$Source: /home/mann/local/research/1/rshock/RiemannObject/RCS/main.C,v $
$Revision: 1.16 $
$Date: 2010/10/19 22:07:53 $

A main test program for the exact Riemann solver class.
============================================================
*/
#include "iostream"
#include "fstream"
#include "cmath"
#include "unistd.h"
#include "cstdlib"

#include "cloption.h"
#include "riemann.h"

using namespace std;

EqOfState eos;

//==========================================================
int main(const int argc, char* argv[] )
{
  // Check command line

  CommandLineOptions cloption;
  cloption.Get( argc, argv );

  /*
    Read initial data from file
    
    Format: each line has an initial value and then everything
            else on the line is ignored (comment).

    line 1:  gamma
    line 2:  T  total time required
    line 3:  left p
    line 4:  left rho0
    line 5:  left v
    line 6:  right p
    line 7:  right rho0
    line 8:  right vel
    line 9:  leftmost grid position
    line 10:  rightmost grid position
  */

  const int N_LINE_MAX = 256;

  const char* infile = "riemann.dat";

  ifstream s_in;
  s_in.open( infile );
  if( !s_in ){
    cerr << "rshock_riemann: ERROR: unable to open \""
	 << infile << "\"\n";
    exit(1);
  }
 
  double gamma_in;
  s_in >> gamma_in;  s_in.ignore( N_LINE_MAX, '\n' );
  eos.Make(gamma_in);

  double T;
  s_in >> T;  s_in.ignore( N_LINE_MAX, '\n' );

  StateVars lstate;
  s_in >> lstate.p;    s_in.ignore( N_LINE_MAX, '\n' );
  s_in >> lstate.rho0;  s_in.ignore( N_LINE_MAX, '\n' );
  s_in >> lstate.vel;  s_in.ignore( N_LINE_MAX, '\n' );

  StateVars rstate;
  s_in >> rstate.p;    s_in.ignore( N_LINE_MAX, '\n' );
  s_in >> rstate.rho0;  s_in.ignore( N_LINE_MAX, '\n' );
  s_in >> rstate.vel;  s_in.ignore( N_LINE_MAX, '\n' );

  double x_left, x_right;
  s_in >> x_left;   s_in.ignore( N_LINE_MAX, '\n' );
  s_in >> x_right;  s_in.ignore( N_LINE_MAX, '\n' );

  s_in.close();

  // Output some initial indicators

  cout << "Riemann Solver, no tangential velocities\n"
       << "  -Objectified version for Godunov solver\n"
       << "========================================\n";

  cout << "Initial Data:\n"
       << "-------------\n";

  cout << " gamma = " << eos.gamma << '\n'
       << " T     = " << T << '\n';

  /*  
      SPECIFIC INTERNAL ENERGY, SPECIFIC ENTHALPY, SOUND SPEED AND  
      FLOW LORENTZ FACTORS IN THE INITIAL STATES 
  */

  cerr << "main: INFO: make ancillary initial data\n";
  lstate.MakeU( eos);
  lstate.MakeH();
  lstate.MakeCS( eos );
  lstate.MakeW();

  rstate.MakeU( eos );
  rstate.MakeH();
  rstate.MakeCS( eos );
  rstate.MakeW();

#if defined DEBUG
  cout << "Left Initial Data:\n" << lstate << '\n';
  cout << "Right Initial Data:\n" << rstate << '\n';
#endif

  // Useful to print out k (not necessary for the analytic solution)

  const double left_eos_k  = lstate.p/pow(lstate.rho0, 1.0/eos.gamma);
  const double right_eos_k = rstate.p/pow(rstate.rho0, 1.0/eos.gamma);

  cout << "main: INFO: Left equation of state k:  " << left_eos_k << '\n'
       << "            Right equation of state k: " << right_eos_k << '\n';

  // Create Riemann class

  const double x_discontinuity = 0.5 * ( x_left + x_right );
  RiemannExact exact( lstate, rstate, T, x_discontinuity );

  exact.Solve();

  cerr << "main: INFO: completed basic solution\n";

  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  /*
    ---------- 
    SOLUTION ON THE MESH 
    ----------
  */

  const int n = cloption.n;
    
  cerr << "main: INFO: "
       << "compute solution on the mesh: n=" << n << "\n";

  ExactNode* nlist = exact.MakeMeshValues( n, x_left, x_right );

  // dump the vectors

  ofstream s;
  string outfile = "riemann.out";

  s.open( outfile.c_str() );
  if( !s ){
    cerr << "RiemannMartiMuller: ERROR: unable to open \""
	 << outfile << "\"\n";
    exit(1);
  }
  
  s << "############################################################\n"
    << "# Riemann Solver, no tangential velocity\n\n";

  s << "#  n = " << n << '\n'
    << "#  T = " << T << '\n'
    << "#  gamma = " << eos.gamma << '\n'
    << '\n';

  s << "# Left Initial State:\n" << exact.lstate << '\n';
  s << "# Right Initial State:\n" << exact.rstate << '\n';

  s << "# Left equation of state k:  " << left_eos_k << '\n'
    << "# Right equation of state k: " << right_eos_k << '\n';

  s << "# Left wave:\n" << exact.lwave << '\n';
  s << "# Right wave:\n" << exact.rwave << '\n';

  s << "# Positions:\n"
    << "#  x1 = " << exact.x1 << '\n'
    << "#  x2 = " << exact.x2 << '\n'
    << "#  x3 = " << exact.x3 << '\n'
    << "#  x4 = " << exact.x4 << '\n'
    << "#  x5 = " << exact.x5 << '\n';

  s << '\n';

  //--------------------------------------------
  // Integrate D to check for mass conservation

  double mass_future = 0.0;
  double mass_past = 0.0;
  for( int i=0; i<(n-1); ++i ){
    const double deltax = nlist[i+1].x - nlist[i].x;

    const double mean_D = 0.5 * ( nlist[i].D + nlist[i+1].D );
    mass_future += deltax * mean_D;
    
    double D = exact.lstate.rho0 / sqrt( 1.0 - exact.lstate.vel*exact.lstate.vel );
    const double mean_x = 0.5 * ( nlist[i].x + nlist[i+1].x );
    if( mean_x > x_discontinuity ){
      D = exact.rstate.rho0 / sqrt( 1.0 - exact.rstate.vel*exact.rstate.vel );
    }
    mass_past += deltax * D;
  }
  const double mass_diff = mass_future - mass_past;
  const double rel_mass_diff = mass_diff/mass_past;

  s << "# Mass Integrals over whole tube with direct quadrature:\n"
    << "#   mass_past     = " << mass_past << '\n'
    << "#   mass_future   = " << mass_future << '\n'
    << "#   mass_diff     = " << mass_diff << '\n'
    << "#   rel_mass_diff = " << rel_mass_diff << " (" << 100.0*rel_mass_diff << " %)\n";

  s << '\n';

  if( fabs(rel_mass_diff) >0.01 ){
    cerr << "main: WARNING: masses seem to be a bit different!!\n";
    cerr << "  rel_mass_diff = " << rel_mass_diff << " (" << 100.0*rel_mass_diff << " %)\n"
	 << "  mass_diff     = " << mass_diff << '\n'
	 << "  mass_past     = " << mass_past << '\n'
	 << "  mass_future   = " << mass_future << '\n';
  }

  //=============================================================
  // Test integration routines

  IntVars intvars, intvars_past;

  exact.Base_Integrals( n, nlist[0].x, nlist[n-1].x, intvars );
  s << "# Base_Integrals(..) over the tube length\n"
    << "#  Integration with " << n << " quadrature nodes from nlist[0].x="
    << nlist[0].x << " to nlist[n-1].x=" << nlist[n-1].x << '\n'
    << "#  .S = " << intvars.S << '\n'
    << "#  .D = " << intvars.D << " (should be about the same as mass integral above)\n"
    << "#  .E = " << intvars.E << '\n';
  s << endl;

  exact.Integrals( n, intvars );
  exact.Integrals_Past( intvars_past );
  s << "# Integrals(..) from x1 to x5 with " << n << " quadrature nodes\n"
    << "#  .S = " << intvars.S << " (past = " << intvars_past.S << ")\n"
    << "#  .D = " << intvars.D << " (past = " << intvars_past.D << ")\n"
    << "#  .E = " << intvars.E << " (past = " << intvars_past.E << ")\n";
  s << endl;

  exact.Left_Integrals( n, intvars );
  exact.Left_Integrals_Past( intvars_past );
  s << "# Left_Integrals(..) from x1 to x_discontinuity\n"
    << "#  .S = " << intvars.S << " (past = " << intvars_past.S << ")\n"
    << "#  .D = " << intvars.D << " (past = " << intvars_past.D << ")\n"
    << "#  .E = " << intvars.E << " (past = " << intvars_past.E << ")\n";
  s << endl;

  EvolutionVars dleft;
  dleft.S = intvars.S - intvars_past.S;
  dleft.D = intvars.D - intvars_past.D;
  dleft.E = intvars.E - intvars_past.E;

  exact.Right_Integrals( n, intvars );
  exact.Right_Integrals_Past( intvars_past );
  s << "# Right_Integrals(..) from x_discontinuity to x5\n"
    << "#  .S = " << intvars.S << " (past = " << intvars_past.S << ")\n"
    << "#  .D = " << intvars.D << " (past = " << intvars_past.D << ")\n"
    << "#  .E = " << intvars.E << " (past = " << intvars_past.E << ")\n";
  s << endl;

  EvolutionVars dright;
  dright.S = intvars.S - intvars_past.S;
  dright.D = intvars.D - intvars_past.D;
  dright.E = intvars.E - intvars_past.E;

  s << "#======================================================================\n";
  EvolutionVars flux;
  exact.MakeFlux( x_discontinuity, flux );
  s << "# Flux calculations\n"
    << "#   flux.S = " << flux.S << '\n'
    << "#   flux.D = " << flux.D << '\n'
    << "#   flux.E = " << flux.E << '\n';

  double S_transfer = flux.S * T;
  s << "#  S integrals:\n"
    << "#    transferred  = " << S_transfer << " (flux * time)\n"
    << "#    left change  = " << dleft.S << " (difference of left integrals)\n"
    << "#    right change = " << dright.S << " (difference of right integrals)\n";
  double D_transfer = flux.D * T;
  s << "#  Mass (D integrals):\n"
    << "#    transferred  = " << D_transfer << " (flux * time)\n"
    << "#    left change  = " << dleft.D << " (difference of left integrals)\n"
    << "#    right change = " << dright.D << " (difference of right integrals)\n";
  double E_transfer = flux.E * T;
  s << "#  E integrals:\n"
    << "#    transferred  = " << E_transfer << " (flux * time)\n"
    << "#    left change  = " << dleft.E << " (difference of left integrals)\n"
    << "#    right change = " << dright.E << " (difference of right integrals)\n";

  //=============================================================
  // Dump node values for graphing

  s << "#--------------------------------------------------------\n";
  string SPACER = "  ";
  int WIDTH = 15;
  s << "# x(1), p(2), rho0(3), vel(4), u(5), S(6), D(7), E(8)\n";
  for( int i=0; i<n; ++i ){
    s << nlist[i];
  }
  s.close();

  cerr << "main: INFO: output is in \"" << outfile << "\"\n";

  delete[] nlist;
}

