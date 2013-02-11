//##############################################################
/**\file riemann.h

$Source: /home/mann/local/research/1/rshock/RiemannObject/RCS/riemann.h,v $
$Revision: 1.12 $
$Date: 2010/08/20 19:01:03 $

\brief Include file for my C++ version of Marti/Muller Riemann solver

\author P. J. Mann

This is setup as a callable base class which essentially can 
return fluxes on a ray (constant).  Therefore it can be used
in Godunov type codes.

There are various integration routines and vectors of nodes that
are pretty much only used for testing purposes.  In particular
the integration routines return values useful in checking
conservation.
*/
//##############################################################
#if ! defined RIEMANN_H
#define RIEMANN_H

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cmath>

#include "eos.h"

using namespace std;

//==============================================================
/// Identifier for element/node pointers to right or left of current

typedef int LeftRight;
enum { LEFT, RIGHT };

// Field width for output

const int WIDTH = 15;

// Constants

#if ! defined FALSE
#define FALSE 0
#define TRUE  1
#endif

//=========================================================
// The equation of state is used everywhere so just make it external

extern EqOfState eos;

//=========================================================
// Utility functions

inline double max( double a, double b ){
  return ( (a) > (b) ? (a) : (b) );
}

inline double signum( double a ){
  return ( (a) >= 0.0 ? (+1) : (-1) );
}

inline double DSIGN( double a, double b ){
  return fabs(a)*signum(b);
}

/*=============================================================
Original common blocks from the Fortran for reference
-----------------------------------------------------
StateVars: lstate and rstate

      DOUBLE PRECISION RHOL, PL, UL, HL, CSL, VELL, WL, 
     &                 RHOR, PR, UR, HR, CSR, VELR, WR 
      COMMON /STATES/  RHOL, PL, UL, HL, CSL, VELL, WL, 
     &                 RHOR, PR, UR, HR, CSR, VELR, WR

WaveVars: lwave and rwave

      DOUBLE PRECISION RHOLS, ULS, HLS, CSLS, VELLS, VSHOCKL 
      COMMON /LS/      RHOLS, ULS, HLS, CSLS, VELLS, VSHOCKL

      DOUBLE PRECISION RHORS, URS, HRS, CSRS, VELRS, VSHOCKR 
      COMMON /RS/      RHORS, URS, HRS, CSRS, VELRS, VSHOCKR

      DOUBLE PRECISION GAMMA 
      COMMON /ADIND/   GAMMA
--------------------------------------
*/

class StateVars;

//=============================================================
class EvolutionVars
{
 public:
  double S, D, E;
  EvolutionVars(){S=0.0; D=0.0; E=0.0; };
  EvolutionVars( StateVars& );

  EvolutionVars( const double invalue ){
    S=invalue; D=invalue; E=invalue; 
  }
  EvolutionVars( const EvolutionVars& ev ){
    S=ev.S; D=ev.D; E=ev.E; 
  }
  inline void operator = (const EvolutionVars& ev ){
    S=ev.S; D=ev.D; E=ev.E; 
  }
  inline void operator = ( const double invalue ){
    S=invalue; D=invalue; E=invalue;
  }
  inline void operator += ( const EvolutionVars& ev ){
    S+=ev.S; D+=ev.D; E+=ev.E;
  }
  inline void operator += ( const double scalar ){
    S+=scalar; D+=scalar; E+=scalar;
  }
  inline void operator *= ( const double scalar ){
    S*=scalar; D*=scalar; E*=scalar;
  }
};
typedef EvolutionVars IntVars;
//=============================================================
/// Initial States (generally left and right)

class StateVars
{
public:
  double rho0; /// rest density (PJM: rho_0)
  double p;    /// pressure
  double u;    /// specific internal energy (PJM: epsilon)
  double h;    /// specific enthalpy
  double cs;   /// sound speed
  double vel;  /// flow velocity
  double w;    /// Lorentz factor 1/sqrt(1-v^2)

  inline StateVars(){};

  inline StateVars( const StateVars& st ){
    rho0 = st.rho0;
    p = st.p;
    u = st.u;
    h = st.h;
    cs = st.cs;
    vel = st.vel;
    w = st.w;
  }
 
  inline void MakeU( EqOfState& eos ){
    u = eos.EPS( p, rho0 );
  }
  inline void MakeH(){
    h = 1.0 + u + p/rho0;
  }
  inline void MakeCS( EqOfState& eos ){
    cs = eos.CS(p,rho0,h);
  };
  inline void MakeW(){
    w = 1.0 / sqrt( 1.0 - vel*vel );
  }
};
inline ostream& operator << ( ostream& s, StateVars& state )
{
  s << "# rho0 = " << setw(WIDTH) << state.rho0 << '\n'
    << "# p    = " << setw(WIDTH) << state.p << '\n'
    << "# u    = " << setw(WIDTH) << state.u << '\n'
    << "# h    = " << setw(WIDTH) << state.h << '\n'
    << "# cs   = " << setw(WIDTH) << state.cs << '\n'
    << "# vel  = " << setw(WIDTH) << state.vel << '\n'
    << "# w    = " << setw(WIDTH) << state.w << '\n';
  return s;
}
//============================================================
/** Variables in the right and left going waves

These are computed.
 */
class WaveVars
{
 public:
  double rho0;
  double u;
  double h;
  double cs;
  double vel;
  double vshock;

  inline void MakeU( double p ){
    eos.EPS(p,rho0);
  }
};

inline ostream& operator << ( ostream& s, WaveVars& wave )
{
  s << "# rho0   = " << setw(WIDTH) << wave.rho0 << '\n'
    << "# u      = " << setw(WIDTH) << wave.u << '\n'
    << "# h      = " << setw(WIDTH) << wave.h << '\n'
    << "# cs     = " << setw(WIDTH) << wave.cs << '\n'
    << "# vel    = " << setw(WIDTH) << wave.vel << '\n'
    << "# vshock = " << setw(WIDTH) << wave.vshock << '\n';
  return s;
}
//==============================================
/// Node on the grid.  Used for output.

class ExactNode
{
 public:
  double x;
  double rho0; /// rest density
  double p;
  double u;   /// epsilon in my notation
  double vel;
  double h;   /// enthalpy w in my notation
  double w;   /// w=1/sqrt(1-v^2).  a in my notation
  double S, D, E; /// These are the conserved quantities in my notation

  inline void MakeSDEhw(){ Makew(); Makeh(); MakeD(); MakeS(); MakeE(); } // order is important.

 private:
  inline void Makeh(){ h = 1.0 + u + p/rho0; }
  inline void Makew(){ w = 1.0 / sqrt( 1.0 - vel*vel ); }
  inline void MakeS(){ S = rho0*h*w*w*vel; }
  inline void MakeE(){ E = rho0*h*w*w - p - D; }  // D must already be made
  inline void MakeD(){ D = rho0/sqrt(1.0-vel*vel); }
};

inline ostream& operator << ( ostream& s, ExactNode& f )
{
  const char* SPACER = "  ";
  s << setw(WIDTH) << f.x    << SPACER
    << setw(WIDTH) << f.p    << SPACER
    << setw(WIDTH) << f.rho0 << SPACER
    << setw(WIDTH) << f.vel  << SPACER
    << setw(WIDTH) << f.u    << SPACER
    << setw(WIDTH) << f.S    << SPACER
    << setw(WIDTH) << f.D    << SPACER
    << setw(WIDTH) << f.E    << '\n';
  return s;
};
//=========================================================
/* Solution routines.  For details see the source.
 */
class RiemannExact
{
 private:
  RiemannExact(){};

  void base_integration( const int n, const double x_left, const double x_right, IntVars& );
  void base_integration_past( const double x_left, const double x_right, IntVars& );
  double GetDVel( double p );
  double GetP( double pmin, double pmax );

 public:
  StateVars lstate;  /// input
  StateVars rstate;  /// input
  double T;          /// Elapsed time required
  double x_discontinuity;  /// Position of discontinuity

  WaveVars lwave;
  WaveVars rwave;

  double x1, x2, x3, x4, x5;  /// Positions characteristics (left-going or right-going)
  double ps, vels;            /// Shock values

  RiemannExact( StateVars& lstate_in, StateVars& rstate_in, double T_in, double x_discontinuity_in ){
    lstate = lstate_in; rstate = rstate_in; T = T_in; x_discontinuity = x_discontinuity_in;
  }

  /// The solution.  All the real work is done here.

  void Solve();

  /// Utility routines that compute various values from the solution

  void GetExactValues( const double x, ExactNode& );

  void Integrals( const int n, IntVars& );
  void Left_Integrals( const int n, IntVars& );
  void Right_Integrals( const int n, IntVars& );
  inline void Base_Integrals( const int n, const double x_left, const double x_right, IntVars& ivars ){
    base_integration( n, x_left, x_right, ivars );
  }

  inline void Integrals_Past( IntVars& ivars ){
    base_integration_past( x1, x5, ivars );
  }
  inline void Left_Integrals_Past( IntVars& ivars ){
    base_integration_past( x1, x_discontinuity, ivars );
  }
  inline void Right_Integrals_Past( IntVars& ivars ){
    base_integration_past( x_discontinuity, x5, ivars );
  }
  ExactNode* MakeMeshValues( const int n, const double x_left, const double x_right );

  void MakeFlux( const double x, EvolutionVars& );
};
//=========================================================
// Routines used during the solution.

void RareF( const double xi, const double gamma,
	    StateVars& state_in, ExactNode& node, LeftRight LR );

void GetVel( const double gamma, StateVars& state_in, WaveVars& wave,
	     double p, LeftRight LR );

#endif
