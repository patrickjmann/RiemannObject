//##############################################################
/**\file Cloption.C

Command line options


$Source: /home/mann/local/research/1/rshock/RiemannObject/RCS/Cloption.C,v $
$Revision: 1.5 $
$Date: 2010/08/20 18:21:18 $
*/
//=============================================================================
#include <iostream>
#include <fstream>
#include <cstdlib>

#include "cloption.h"

using namespace std;

//=============================================================================
void CommandLineOptions::Get( const int argc, char* argv[] )
{
  extern int opterr , optind;
  extern char* optarg;
  int opt;
  
  // Defaults

  n = 401;

  // Note: turn off getopt() error checking
  
  opterr = FALSE;
  
  while( ( opt = getopt( argc, argv, "h?dn:" ) ) != EOF ){
    
    switch( opt ){
    case '?':         // getopt() does not recognize the argument
      cerr << "CommandLineOptions: ERROR: invalid command line option\n";
      give_help();
      exit( 0 );
    case 'h':
      give_help();
      exit( 0 );
    case 'd':{
      string datafile = "riemann.dat.test";
      cerr << "Test parameter file generated in \"" << datafile << "\"\n"
	   << "       (should be moved to \"riemann.dat\")\n";
      MakeDataFile( datafile );
      exit( 0 );
    }
    case 'n':{
      n = atoi( optarg );
      break;
    }
    default:
      cerr << "CommandLineOptions: ERROR: invalid command line option\n";
      give_help();
      exit( 0 );
      break;
    }
  }
}
//=============================================================================
/// Create a datafile with default parameters

void CommandLineOptions::MakeDataFile( const string& datafile )
{
  ofstream s( datafile.c_str() );
  if( !s ){
    cerr << "MakeDataFile: ERROR: unable to open \"" << datafile << "\"\n";
    exit(1);
  }

  s << "2.0      eos.gamma\n"
    << "0.25     time\n"
    << "8.0      left p\n"
    << "8.0      left rho0\n"
    << "0.0      left velocity\n"
    << "0.0      right p\n"
    << "2.0      right rho0\n"
    << "0.0      right velocity\n"
    << "-0.3     left grid position\n"
    << " 0.3     right grid position\n";

  s.close();
}
//=============================================================================
void CommandLineOptions::give_help()
{
  cout << "1d Relativistic Riemann Shock with no tangential velocity.\n"
       << "  This version is objectified for inclusion in Godunov method\n";

  cout << '\n'
       << "USAGE: rshock_riemann [-h] [-d] [-n #]\n";

  cout << "  -h:   help and short parameter description\n"
       << "  -d:   generate a test datafile.\n"
       << "  -n #: set the number of nodes in the output grid (default 401)\n";

  cout << '\n'
       << "RESULTS: Various integrated quantities are calculated and presented\n"
       << "         in the output file.  These include initial and final integrals\n"
       << "         which can be compared for conservation.\n"
       << "           -mass and internal energy are conserved.\n"
       << "           -momentum is not conserved due to boundary conditions.\n";

  cout << endl;
}
