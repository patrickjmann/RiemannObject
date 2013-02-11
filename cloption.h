//##############################################################
/**\file cloption.h

$Source: /home/mann/local/research/1/rshock/RiemannObject/RCS/cloption.h,v $
$Revision: 1.2 $
$Date: 2010/05/09 17:40:12 $

\brief Include for command line options.

-No particular options here.  Just used to dump help or data file

\author P. J. Mann
*/
//##############################################################
#if ! defined RIEMANN_CLOPTION_H
#define RIEMANN_CLOPTION_H

#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

#if ! defined FALSE
#define FALSE 0
#define TRUE  1
#endif

//=============================================================
/** Command Line Options
*/

class CommandLineOptions
{
 public:
  int n;

  void give_help();
  CommandLineOptions(){};
  ~CommandLineOptions(){};
  void Get( const int argc, char* argv[] );
  void MakeDataFile( const string& datafile );
};

#endif
