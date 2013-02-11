//##############################################################
/**\file Riemann.C

$Source: /home/mann/local/research/1/rshock/RiemannObject/RCS/Riemann.C,v $
$Revision: 1.22 $
$Date: 2012/03/23 15:13:19 $

My translation of the original Fortran Code.

This version is wrapped in a class to be used in Godunov codes.

See <A HREF="http://relativity.livingreviews.org/Articles/lrr-2003-7/">
Numerical Hydrodynamics in Special Relativity</A> by Jose Marti and
Ewald Muller.  The original Fortran was downloaded from the
<A HREF="http://relativity.livingreviews.org/Articles/lrr-2003-7/download/index.html"> download</A> link.

--------------------------------------------------------------
Comments from the Original Fortran code by Marti and Muller

------------
NAME: R I E M A N N 
------------

PURPOSE:
THIS PROGRAM COMPUTES THE SOLUTION OF A 1D   
RELATIVISTIC RIEMANN PROBLEM (FOR CONSTANT-GAMMA IDEAL GASES) WITH  
INITIAL DATA UL IF X<0.5 AND UR IF X>0.5  
IN THE WHOLE SPATIAL DOMAIN [0, 1] 


COMMENTS:
SEE MARTI AND MUELLER, JFM, 1994

WRITTEN BY:     Jose-Maria Marti
Departamento de Astronomia y Astrofisica 
Universidad de Valencia 
46100 Burjassot (Valencia), Spain
jose-maria.marti@uv.es
AND
Ewald Mueller
Max-Planck-Institut fuer Astrophysik
Karl-Schwarzschild-Str. 1
85741 Garching, Germany
emueller@mpa-garching.mpg.de

============================================================
*/
#include <iostream>
#include <fstream>
#include <cmath>
#include <unistd.h>
#include <cstdlib>
#include <limits>

#include "cloption.h"
#include "riemann.h"

using namespace std;

//==========================================================
/** 
Solve for the characteristic positions: x1, x2, x3, x4, x5.

After these positions have been calculated then Riemann:GetExactValues(x,..) can be used to
calculate fluid quantities at a given point x.
*/
void RiemannExact::Solve( )
{
  // No shock if there is no pressure
  // -Pressureless fluid will not allow for a shock.  Just a bunch of particles.
  // -Interesting, what happens if there is no pressure, but there is a density differential?
  // -and what happens if there is no pressure, but velocity differential
  // -Could do with a characteristic pressure to compare to*******************************

  const double EPS = numeric_limits<double>::epsilon();
  //if( fabs(lstate.p - rstate.p) <= EPS )



  // Find a straddle for the GetDVel(..) solver.
  // -that solver uses bisection, inverse quadratic, etc. to find a solution
  //  and needs to start from a straddle.

  //$$cerr << "RiemannExact::Solve: INFO: find pmin and pmax\n";

  const int N_STRADDLE_MAX = 10;

  double pmin = ( lstate.p + rstate.p ) /2.0;
  double pmax = pmin;

  double check = 1.0;

  int istraddle = 0;
  while( check > 0.0 ){   // No straddle, "5" loop in the fortran
    ++istraddle;
    if( istraddle > N_STRADDLE_MAX ){
      cerr << "RiemannExact::Solve: ERROR: unable to find pmin/pmax straddle\n"
	   << " istraddle = " << istraddle << '\n'
	   << " pmin = " << pmin << "  pmax = " << pmax << '\n'
	   << " GetDVel(pmin)=" << GetDVel(pmin) << "  GetDVel(pmax)=" << GetDVel(pmax) << '\n'
	   << " lstate:\n" << lstate << '\n'
	   << " rstate:\n" << rstate << '\n';
      exit(1);
    }

    pmin = 0.5 * max( pmin, 0.0 );
    pmax = 2.0 * pmax;

    double dvel1 = GetDVel( pmin );
    double dvel2 = GetDVel( pmax );

    cerr << "   Iteration: " << istraddle << " pmin=" << pmin << "  pmax=" << pmax
         << " dvel1=" << dvel1 << " dvel2=" << dvel2 << '\n';

    check = dvel1 * dvel2;
  }

  if( pmin >= pmax ){
    cerr << "RiemannExact::Solve: ERROR: straddle has pmin>=pmax\n"
	 << "  pmin=" << pmin << "   pmax=" << pmax << '\n';
    exit(1);
  }
  // Ok, found a straddle and can continue on.

  /*
    --------------------------- 
    PRESSURE AND FLOW VELOCITY IN THE INTERMEDIATE STATES 
    ---------------------------
  */

  /*
  cerr << "RiemannExact::Solve: INFO: "
       << "compute pressure and flow velocity in the intermediate states\n";
  */

  // GetP(..) solves GetDVel(..p..)==0, and starts from the straddle found above
  // -this is  where all the work get's done!

  ps = GetP( pmin, pmax );
  vels = 0.5 * ( lwave.vel + rwave.vel );

  cerr << "RiemannExact::Solve: INFO: after GetP, ps=" << ps << "  vels=" << vels << '\n';

  /*$$
  cerr << "    lwave.vel=" << lwave.vel << " rwave.vel=" << rwave.vel 
       << " ps=" << ps << "  vels=" << vels << '\n';
  */

  /*
    ----------- 
    POSITIONS OF THE WAVES 
    -----------
  */
  // Now that p and flow velocity are known the positions of the waves can be calculated.

  //$$ cerr << "RiemannExact::Solve: INFO: compute positions of the waves\n";

  if( lstate.p >= ps ){    
    x1 = x_discontinuity +
      T * (lstate.vel - lstate.cs )/( 1.0 - lstate.vel*lstate.cs );
    x2 = x_discontinuity +
      T * (vels       - lwave.cs  )/( 1.0 - vels*lwave.cs );
  } else {
    x1 = x_discontinuity + lwave.vshock * T;
    x2 = x1;
  }

  x3 = x_discontinuity + vels * T;

  if( rstate.p >= ps ){
    x4 = x_discontinuity +
      T * ( vels       + rwave.cs  )/( 1.0 + vels*rwave.cs);
    x5 = x_discontinuity +
      T * ( rstate.vel + rstate.cs )/( 1.0 + rstate.vel*rstate.cs );
  } else {
    x4 = x_discontinuity + rwave.vshock*T;
    x5 = x4;
  }

  /*$$
  cerr << "RiemannExact::Solve:INFO: after solution:\n"
       << x1=" << x1 << " x2=" << x2 << " x3=" << x3 << " x4=" << x4 << " x5=" << x5 << '\n';
  cerr << "   T=" << T << " x_discontinuity=" << x_discontinuity << '\n'
       << "   vels=" << vels
       << " lwave.vshock=" << lwave.vshock << " rwave.vshock=" << rwave.vshock << '\n';
  */
}
//=========================================================
/**
  ------- 
  NAME: G E T P 
  -------

  PURPOSE: 
  FIND THE PRESSURE IN THE INTERMEDIATE STATE OF A RIEMANN PROBLEM IN 
  RELATIVISTIC HYDRODYNAMICS 
  
  
  COMMENTS: 
  THIS ROUTINE USES A COMBINATION OF INTERVAL BISECTION AND INVERSE 
  QUADRATIC INTERPOLATION TO FIND THE ROOT IN A SPECIFIED INTERVAL. 
  IT IS ASSUMED THAT DVEL(PMIN) AND DVEL(PMAX) HAVE OPPOSITE SIGNS WITHOUT 
  A CHECK. 
  ADAPTED FROM "COMPUTER METHODS FOR MATHEMATICAL COMPUTATION", 
  BY G. E. FORSYTHE, M. A. MALCOLM, AND C. B. MOLER, 
  PRENTICE-HALL, ENGLEWOOD CLIFFS N.J. 

/retval pressure (PS)
*/

double RiemannExact::GetP( double pmin, double pmax )
{
  const double EPS = numeric_limits<double>::epsilon();
  const double TOL = 10.0*EPS;
  const int N_MAX_ITERATIONS = 100;

  cerr << "GetP: INFO: Entry: machine precision: EPS = " << EPS << "\n"
       << "                           Tolerance: TOL = " << TOL << '\n';
 
  // Test for reasonable input straddle.  Not necessary but useful for development.

  if( pmin >= pmax ){
    cerr << "RiemannExact::GetP: ERROR: pmin >= pmax\n"
	 << "  pmin=" << pmin << "  pmax=" << pmax << '\n';
  }

  /*
    INITIALIZATION: Calculate initial function values for solution iteration.
  */
  
  double a = pmin;
  double b = pmax;
  
  double fa = GetDVel( a );
  double fb = GetDVel( b );

  // Check to see if there is a straddle

  if( fa*fb > 0.0 ){
    cerr << "RiemannExact::GetP: ERROR: the input values do not straddle the zero!\n"
	 << "  a=pmin=" << pmin << "  b=pmax=" << pmax << '\n'
	 << "  fa=" << fa << "  fb=" << fb << '\n';
    exit(1);
  }

  cerr << "GetP: INFO: on entry\n"
       << "  a=" << a << " fa=" << fa << '\n'
       << "  b=" << b << " fb=" << fb << '\n';

  // Initialize the "c" point (for inverse quadratic interpolation)

  double c = a; 
  double fc = fa;
  double diff = b - a;
  double old_diff = diff; 

  // Iterate

  for( int i=1; i<=N_MAX_ITERATIONS; ++i ){

    // Rearrange so |f(b)| <= |f(c)|
    // I don't really understand this.  "a" is used as a temp variable??

    if( fabs(fc) < fabs(fb) ){
      a = b;
      b = c;
      c = a;
      fa = fb;
      fb = fc;
      fc = fa;
    }

    // CONVERGENCE TEST 
    // I don't know why there are factors of 2.0 and 0.5 here?

    double RELATIVE_TOL_PRESSURE = 2.0 * EPS * fabs(b) + 0.5 * TOL;
    double xm = 0.5 * ( c - b );

    //$$ cerr << "xm = " << xm << " ";
    if( fabs(xm) <= RELATIVE_TOL_PRESSURE ){
      cerr << "RiemannExact: CONVERGENCE test passed at i = " << i << ", return\n";
      return b;
    }

    if( fb == 0.0 ){
      cerr << "RiemannExact: CONVERGENCE with fb=0 at i = " << i  << ", return\n";
      return b;
    }

    // IS BISECTION NECESSARY? 

    if( (fabs(old_diff) < RELATIVE_TOL_PRESSURE) || (fabs(fa) <= fabs(fb)) ){
      diff = xm;
      old_diff = diff;
    } else {
  
      // IS QUADRATIC INTERPOLATION POSSIBLE? 
    
      double p, q;
      if( a == c ){
	
	// LINEAR INTERPOLATION 
	
	double fbfa = fb/fa;
	p = 2.0 * xm * fbfa;
	q = 1.0 - fbfa;
      } else {
	
	// INVERSE QUADRATIC INTERPOLATION 

	if( fc == 0 ){
	  return c;
	}
	if( fa == 0 ){
	  return a;
	}
	
	double fafc = fa/fc;
	double fbfc = fb/fc;
	double fbfa = fb/fa;
	
	p = fbfa * ( 2.0*xm*fafc*(fafc-fbfc) - (b-a)*(fbfc-1.0) );
	q = (fafc-1.0)*(fbfc-1.0)*(fbfa-1.0);
      }
      
      // ADJUST SIGNS 
      
      if( p > 0.0 ) q = -q;
      p = fabs(p);
      
      // IS INTERPOLATION ACCEPTABLE? 
      
      if( ( (2.0*p) >= (3.0*xm*q-fabs(RELATIVE_TOL_PRESSURE*q)) ) ||
	  ( p >= fabs(0.5*old_diff*q) ) ){ //no, back to bisection
	// Setup for bisection
	diff = xm;
	old_diff = diff;
      } else {
	// Setup for using the interpolation
	old_diff = diff;
	diff = p/q;
      }
    }
    
    // COMPLETE THE STEP
    // a is replaced by b, and b is given the new (improved) value
    
    a = b;
    fa = fb;

    if( fabs(diff) > RELATIVE_TOL_PRESSURE ){ // Just ensure DIFF>roundoff
      b = b + diff;
    } else {                                  // add a pre-defined small change
      b = b + DSIGN(RELATIVE_TOL_PRESSURE, xm);
    }
    fb = GetDVel( b );

    if( fb*fc > 0.0 ){
      c = a; 
      fc = fa;
      diff = b - a;
      old_diff = diff; 
    }
    cerr << "RiemannExact::GetP: INFO: completed iteration i=" << i
	 << " b=" << setprecision(20) << b << setprecision(6) << " diff=" << diff << '\n'
	 << "   a=" << a << " fa=" << fa << '\n'
	 << "   b=" << b << " fb=" << fb << '\n'
	 << "   c=" << c << " fc=" << fc << '\n';
    cerr << "#--------------------------------------------------------------------------------\n";
  } // Iteration loop

  cerr << "RiemannExact::GetP: ERROR: Did not converge after N_MAX_ITERATIONS=" << N_MAX_ITERATIONS << '\n';
  cerr << "  pmin=" << pmin << " f(pmin)=" << GetDVel(pmin)
       << "  pmax=" << pmax << " f(pmax)=" << GetDVel(pmax) << '\n';
  cerr << "  a=" << a << " fa=" << fa << '\n'
       << "  b=" << b << " fb=" << fb << '\n'
       << "  diff=" << diff << "  old_diff=" << old_diff << '\n';
  exit(1);
}
//=======================================================
/** GetDVel

NAME: G E T D V E L 
----------

PURPOSE: 
COMPUTE THE DIFFERENCE IN FLOW SPEED BETWEEN LEFT AND RIGHT INTERMEDIATE 
STATES FOR GIVEN LEFT AND RIGHT STATES AND PRESSURE 

COMMENTS
NONE

\param[in] p Value of the pressure.  Usually the current iteration estimate.

\retval difference in flow speed between left and right itermediate states (dvel)
*/
double RiemannExact::GetDVel( double p )
{
  GetVel( eos.gamma, lstate, lwave, p, LEFT );
  GetVel( eos.gamma, rstate, rwave, p, RIGHT );

  double dvel = lwave.vel - rwave.vel;
  return dvel;
}
//====================================================================
/** Find solution at a point

This is a general approach for any point.  Might be faster if
there is a priori knowledge of the region the point is in.

Note clever approach for left or right-going shocks. 
  -x1==x2 means right-going
  -x4==x5 means left-going

\param[in] x Position that values are to be calculated.
\param[out] node the values on this node
*/
void RiemannExact::GetExactValues( const double x, ExactNode& node )
{
  /*$$
  cerr << "# RiemannExact::GetExactValues: INFO: On entry x = " << x << '\n'
       << "#    x1=" << x1 << " x2=" << x2 << " x3=" << x3 << " x4=" << x4 << " x5=" << x5 << '\n';
  */

  node.x = x;

  if( x <= x1 ){
    node.p   = lstate.p;
    node.rho0 = lstate.rho0;
    node.vel = lstate.vel;
    node.u   = lstate.u;
    
  } else if( x <= x2 ){

    double xi = ( x - x_discontinuity )/T;
    RareF( xi, eos.gamma, lstate, node, LEFT );
      
  } else if( x <= x3 ){

    node.p   = ps;
    node.rho0 = lwave.rho0;
    node.vel = vels;
    node.u   = lwave.u;

  } else if( x <= x4 ){

    node.p   = ps;
    node.rho0 = rwave.rho0;
    node.vel = vels;
    node.u   = rwave.u;

  } else if( x <= x5 ){

    double xi = ( x - x_discontinuity )/T;
    RareF( xi, eos.gamma, rstate, node, RIGHT );

  } else {

    node.p   = rstate.p;
    node.rho0 = rstate.rho0;
    node.vel = rstate.vel;
    node.u   = rstate.u;
  }
  node.MakeSDEhw();
}
//====================================================================
/*
    --------- 
    NAME: G E T V E L 
    ---------

    PURPOSE: 
    COMPUTE THE FLOW VELOCITY BEHIND A RAREFACTION OR SHOCK IN TERMS OF THE 
    POST-WAVE PRESSURE FOR A GIVEN STATE AHEAD THE WAVE IN A RELATIVISTIC 
    FLOW 
 

    COMMENTS: 
    THIS ROUTINE CLOSELY FOLLOWS THE EXPRESSIONS IN MARTI AND MUELLER, 
    J. FLUID MECH., (1994)


In original code the first set of variables (with A postfix) are
State variables.  The second set (no postfix) are wave vars.

\param[in] gamma Equation of state gamma
\param[in] state Initial state of the region
\param[out] wave variables, including the velocity term
\param[in] p value of the pressure to use in the calculation
\param[in] LR left/right propogating wave
*/
void GetVel( const double gamma, StateVars& state, WaveVars& wave,
	     const double p, const LeftRight LR )
{
  const double gam1 = gamma - 1.0;

  /*$$
  if( p <= 0.0 ){
    cerr << "\nGetVel: ERROR: input p<=0: p=" << p << '\n'
	 << " Input state is:\n" << state << '\n';
    exit(1);
  }
  */

  // Test to see if state.p ~= p
  // -This is from the evolution issues where the solver will not converge
  //  if left and right states are almost the same.

  const double EPS = 10.0 * numeric_limits<double>::epsilon();

  const double p_diff = p - state.p;
  if( fabs(p_diff) < EPS ){
    wave.h = state.h;
    wave.u = state.u;
    wave.rho0 = state.rho0;
    wave.cs = state.cs;
    wave.vshock = 0.0;
    wave.vel = state.vel;
    return;
  }

  /*
    --------------- 
    LEFT OR RIGHT PROPAGATING WAVE 
    ---------------
  */
  double SIGN;

  if( LR == LEFT ){
    SIGN = -1.0;
  }else if( LR == RIGHT ){
    SIGN = 1.0;
  } else {
    cerr << "GetVel: ERROR: Left/Right switch is undefined, LR = " << LR << '\n';
    exit(1);
  }

  if( p > state.p ){

    // SHOCK 

    //$$ cerr << "GetVel: INFO: in shock section \n";    
  
    double a  = 1.0 + gam1 * (state.p - p)/(gamma*p);
    if( a == 0.0 ){
      cerr << "\nGetVel: ERROR: Shock section: a=" << a << '\n'
	   << " Input state is:\n" << state << '\n';
      exit(1);
    }

    double b  = 1.0 - a;
    double c  = state.h * (state.p-p) / state.rho0 - state.h*state.h;
    
    // CHECK FOR UNPHYSICAL ENTHALPIES 
    
    if( c*4.0*a > b*b ){
      cerr << "GetVel: ERROR: in shock section: "
	   << "unphysical specific enthalpy in intermediate state\n";
      exit(1);
    }
    /*
      ----------------------------- 
      SPECIFIC ENTHALPY IN THE POST-WAVE STATE 
      (FROM THE EQUATION OF STATE AND THE TAUB ADIABAT, 
      EQ.(74), MM94) 
      -----------------------------
    */
    wave.h = ( -b + sqrt( b*b - 4.0*a*c ) )/( 2.0*a );
    if( wave.h <= 1.0 ){
      cerr << "GetVel: ERROR: in shock section: "
	   << "wave.h = " << wave.h << " (<1.0) which cannot be used for evaluating p\n";
      exit(1);
    }
       
    // DENSITY IN THE POST-WAVE STATE (FROM EQ.(73), MM94) 

    wave.rho0 = gamma*p / ( gam1*(wave.h-1.0) );
    /*
      ------------------------ 
      SPECIFIC INTERNAL ENERGY IN THE POST-WAVE STATE 
      (FROM THE EQUATION OF STATE) 
      ------------------------
    */
    wave.u = p / ( gam1*wave.rho0 );
    /*
      -------------------------- 
      MASS FLUX ACROSS THE WAVE  
      (FROM THE RANKINE-HUGONIOT RELATIONS, EQ.(71), MM94) 
      --------------------------

      There can be a problem here if p ~= state.p.  Then wave~=state (all variables)
      and the denominator below=0.
    */
    const double denom = state.h/state.rho0 - wave.h/wave.rho0;
    if( denom < 0.0 ){
      cerr << "GetVel: ERROR: in shock section denominator<0\n"
	   << "Usually due to state ~= wave so 0/0 issues.\n";
      exit(1);
    }
    double mflux = SIGN * sqrt( (p-state.p) / denom );
 
    // SHOCK VELOCITY (FROM EQ.(86), MM94 
    
    double mflux2 = mflux * mflux;

    double tmpsa = mflux2 + pow((state.rho0*state.w),2);
    double tmpsb = -state.vel * state.rho0*state.rho0 * state.w*state.w;

    wave.vshock = (-tmpsb + SIGN*mflux2*sqrt( 1.0 + state.rho0*state.rho0/mflux2 ) )
      / tmpsa;
    const double denom_s = 1.0 - wave.vshock*wave.vshock;
    if( denom_s <= 0.0 ){
      cerr << "GetVel: ERROR: in shock section fabs(wave.vshock) >= 1.0\n"
	   << "  wave.vshock=" << wave.vshock << '\n';
    }
    double wshock = 1.0 / sqrt(denom_s);

    // FLOW VELOCITY IN THE POST-SHOCK STATE (FROM EQ.(67), MM94) 

    double tmpvela = wshock * (p-state.p)/mflux + state.h*state.h*state.vel;
    double tmpvelb = state.h*state.w +
      (p-state.p) * ( wshock*state.vel/mflux + 1.0/(state.rho0*state.w) );

    if( tmpvelb == 0.0 ){
      cerr << "GetP: ERROR in shock section.  tmpvelb = " << tmpvelb << '\n';
      exit(1);
    }

    wave.vel = tmpvela / tmpvelb;
    /*
      --------------------- 
      LOCAL SOUND SPEED IN THE POST-SHOCK STATE 
      (FROM THE EQUATION OF STATE) 
      ---------------------
    */
    wave.cs = sqrt( gamma*p/(wave.rho0*wave.h) );
    
  } else {
    // RAREFACTION
 
    /* Compute rho0 in the rarefaction wave.
       This is done by using a polytropic equation of state and computing k from
       the ratio of input (state) values for p and rho0.
       If this state has p=0 (usually from internal energy=0) then the method fails (k=0)
       and I have guessed that the rarefaction wave density is just the input state density.
    */

    if( fabs(state.p) <= 0.0 ){                   // Guess a value
      /*$$$
	 cerr << "GetVel:ERROR: In RAREFACTION section state.p=" << state.p << " (<=0)\n"
	 << "state\n"
	 << "-----\n" << state
	 << "input p=" << p << '\n';
	 exit(1);
      */
      wave.rho0 = state.rho0; // Density behind the rarefaction
    } else {                                      // ok, p>0 so continue on
      /*
	--------------------------- 
	POLYTROPIC CONSTANT OF THE GAS ACROSS THE RAREFACTION 
	---------------------------
	NOTE: if p=0 then this part FAILS.  It assumes p=k*rho0^gamma.
      */

      const double k = state.p / pow(state.rho0,gamma);
      /*
	--------------- 
	DENSITY BEHIND THE RAREFACTION 
	---------------
      */
      wave.rho0 = pow( (p/k), (1.0/gamma) );
      cerr << "GetVel: INFO: Rarefaction section: p=" << p << " state.p=" << state.p
	   << " wave.rho0 = " << wave.rho0 << '\n';
    }
    /*
      ------------------------ 
      SPECIFIC INTERNAL ENERGY BEHIND THE RAREFACTION 
      (FROM THE EQUATION OF STATE) 
      ------------------------
    */
    wave.u = p / ( gam1*wave.rho0 );
    /*
      -------------------- 
      LOCAL SOUND SPEED BEHIND THE RAREFACTION 
      (FROM THE EQUATION OF STATE) 
      --------------------
    */
    wave.cs = sqrt( gamma*p / (wave.rho0 + gamma*p/gam1 ) );
    /*
      ------------------ 
      FLOW VELOCITY BEHIND THE RAREFACTION 
      ------------------
    */
    double sqrt_gam1 = sqrt(gam1);

    double tmp3 = ((sqrt_gam1+state.cs)/(sqrt_gam1-state.cs)) *
      (sqrt_gam1-wave.cs)/(sqrt_gam1+wave.cs);
    cerr << " tmp3=" << tmp3 << " state.cs=" << state.cs << '\n';

    double tmp_exp = -SIGN*2.0/sqrt_gam1;
    cerr << " tmp_exp = " << tmp_exp << '\n';

    double tmp_a = ( ( 1.0 + state.vel )/(1.0 - state.vel ) ) *
      pow( tmp3, tmp_exp);
    cerr << " tmp_a=" << tmp_a << '\n';

    wave.vel = ( tmp_a-1.0)/(tmp_a+1.0);

    cerr << "GetVel: INFO: rarefaction section at end:\n"
	 << "    wave.u=" << wave.u << " wave.vshock=" << wave.vshock << " wave.vel=" << wave.vel
	 << " wave.cs=" << wave.cs << " p=" << p << '\n';
  }
}
//===================================================================
/*
  -------- 
  NAME: R A R E F 
  --------

  PURPOSE: 
  COMPUTE THE FLOW STATE IN A RAREFACTION FOR GIVEN PRE-WAVE STATE 
  
  
  COMMENTS: 
  THIS ROUTINE CLOSELY FOLLOWS THE EXPRESSIONS IN MARTI AND MUELLER, 
  J. FLUID MECH., (1994)

SUBROUTINE RAREF( XI, RHOA, PA, UA, CSA, VELA, S, RHO, P, U, VEL )

RHOA, PA, UA, CSA, VELA are input StateVars (COMMON/STATES/)
S,RHO,P,U,VEL are separate values in the vectors

PJM Note: they get the 'A' ending mixed up!!!  Have to be very
          careful during translation of their code.

/param[in] xi position
/param[in] gamma equation of state
/param[in] state pre-wave state values
/param[out] node values at the position
/param[in] LR left or right propogating wave
*/
void RareF( const double xi, const double gamma,
	    StateVars& state, ExactNode& node, const LeftRight LR )
{ 
  const double gam1 = gamma - 1.0;
  /*
    --------------- 
    LEFT OR RIGHT PROPAGATING WAVE 
    ---------------
  */
  double SIGN;
  if( LR == LEFT ){
    SIGN = 1.0;
  } else if( LR == RIGHT ){
    SIGN = -1.0;
  } else {
    cerr << "RareF: ERROR: LR is neither LEFT nor RIGHT\n";
    exit(1);
  }
  
  double b = sqrt( gamma - 1.0 );
  double c = (b + state.cs)/(b - state.cs );
  double d = -SIGN * b / 2.0;
  double k = (1.0 + xi )/( 1.0 - xi );
  double l = c * pow(k,d);
  double v = pow( (1.0-state.vel)/(1.0+state.vel), d );

  double ocs2 = state.cs;

 label25:;
  double fcs2 = l*v * pow( (1.0+SIGN*ocs2), d ) * (ocs2-b) +
    pow((1.0 - SIGN*ocs2),d) * (ocs2+b);

  double dfdcs2 = l*v * pow( (1.0 + SIGN*ocs2),d ) *
    (1.0 + SIGN*d*(ocs2-b)/(1.0+SIGN*ocs2)) +
    pow( (1.0-SIGN*ocs2), d) *
    (1.0 - SIGN*d*(ocs2+b)/(1.0-SIGN*ocs2));
  
  double cs2 = ocs2 - fcs2/dfdcs2;

  if( fabs( cs2 - ocs2 )/ocs2 > 5.0e-7 ){
    ocs2 = cs2;
    goto label25;
  }

  node.vel = ( xi + SIGN*cs2 )/( 1.0 + SIGN*xi*cs2 );

  double tmp1 = ( cs2*cs2 * (gam1-state.cs*state.cs) )/
    ( state.cs*state.cs*(gam1-cs2*cs2) );

  node.rho0 = state.rho0 * pow( tmp1, (1.0/gam1) );

  node.p = cs2*cs2*gam1 * node.rho0/(( gam1-cs2*cs2 )*gamma);

  node.u = node.p / ( gam1*node.rho0 );
}
//=========================================================
/** Integrals between x_left and x_right

Since it is quite an effort to compute the values all integrals are done in
one sweep.

Uses a simple Trapezoidal rule for now.  Might be good to improve this.

\param[in] n number of integration points (n-1 = number of integration intervals)
\param[in] x_left left boundary
\param[in] x_right right boundary
\param[out] integral Values of the integrals (S,D,E)
*/
void RiemannExact::base_integration( const int n, const double x_left, const double x_right,
				     IntVars& integral )
{
  const double range = x_right - x_left;
  const double deltax = range / double(n-1);

  IntVars work_sum = 0.0;

  // Interior values

  ExactNode node;
  for( int i=1; i<(n-1); ++i ){
    double x = x_left + deltax * double(i);
    GetExactValues( x, node );
    work_sum.S += node.S;
    work_sum.D += node.D;
    work_sum.E += node.E;
  }
  work_sum *= 2.0;

  // Boundary values

  GetExactValues(x_left, node );
  work_sum.S += node.S;
  work_sum.D += node.D;
  work_sum.E += node.E;  

  GetExactValues(x_right, node );
  work_sum.S += node.S;
  work_sum.D += node.D;
  work_sum.E += node.E;
 
  work_sum *= (range/double(2*(n-1)));
  integral = work_sum;
  return;
}
//-----------------------------------------------------------------------
/** Integrals between x_1 and x_discontinuity

\param[in] n number of integration points
\param[out] integral Values of the integrals (S,D,E)
*/

void RiemannExact::Integrals( const int n, IntVars& integral )
{
  const double x_left = x1;
  const double x_right = x5;
  base_integration( n, x_left, x_right, integral );
}

void RiemannExact::Left_Integrals( const int n, IntVars& integral )
{
  const double x_left = x1;
  const double x_right = x_discontinuity;
  base_integration( n, x_left, x_right, integral );
  return;
}

void RiemannExact::Right_Integrals( const int n, IntVars& integral )
{
  const double x_left = x_discontinuity;
  const double x_right = x5;
  base_integration( n, x_left, x_right, integral );
  return;
}
//====================================================================
/// Integration on past (initial) time slice

void RiemannExact::base_integration_past( const double x_left,
					  const double x_right,
					  IntVars& integral )
{
  if( x_left >= x_right ){
    cerr << "RiemannExact::base_integration_past: ERROR: "
	 << " x_left = " << x_left << " >= " << x_right << '\n';
    exit(1);
  }

  if( x_right <= x_discontinuity ){
    const double dx = x_right - x_left;
    EvolutionVars evars( lstate );
    integral.S = evars.S * dx;
    integral.D = evars.D * dx;
    integral.E = evars.E * dx;
    return;
  } else if( x_left >= x_discontinuity ){
    const double dx = x_right - x_left;
    EvolutionVars evars( rstate );
    integral.S = evars.S * dx;
    integral.D = evars.D * dx;
    integral.E = evars.E * dx;
    return;
  } else {
    const double dx_left  = x_discontinuity - x_left;
    const double dx_right = x_right - x_discontinuity;
    EvolutionVars eleft( lstate );
    EvolutionVars eright( rstate );
    integral.S = eleft.S * dx_left + eright.S * dx_right;
    integral.D = eleft.D * dx_left + eright.D * dx_right;
    integral.E = eleft.E * dx_left + eright.E * dx_right;
    return;
  }
}
//====================================================================
/// Make fluxes

void RiemannExact::MakeFlux( const double x, EvolutionVars& flux )
{
  ExactNode exact;
  GetExactValues( x, exact );
  flux.S = exact.S*exact.vel + exact.p;
  flux.D = exact.D*exact.vel;
  flux.E = exact.E*exact.vel + exact.p*exact.vel;
}
//====================================================================
/// Make evolution variables from primitive variables

EvolutionVars::EvolutionVars( StateVars& s )
{
  S = s.rho0*s.h*s.w*s.w*s.vel;
  D = s.rho0/sqrt(1.0-s.vel*s.vel);
  E = s.rho0*s.h*s.w*s.w - s.p - D;  // D must already be made
}
//====================================================================
/** Compute solution on a simple 1d mesh

\retval ExactNode* pointer to array of ExactNodes containing the mesh
*/
ExactNode* RiemannExact::MakeMeshValues( const int n, const double x_left, const double x_right )
{
  if( n <= 5 ){
    cerr << "RiemannExact::MakeMeshValues: ERROR: b too small\n"
	 << "  n=" << n << '\n';
    exit(1);
  }
    
  if( x_right <= x_left ){
    cerr << "RiemannExact::MakeMeshValues: ERROR: x_right <= x_left\n"
	 << "  x_left=" << x_left << "   x_right=" << x_right << '\n';
    exit(1);
  }

  ExactNode* nlist = new ExactNode[n];

  // Values at nodes for graphing purposes

  const double dx = ( x_right - x_left )/double(n-1);

  for( int i=0; i<n; ++i ){
    double x = x_left + double(i) * dx;
    GetExactValues( x, nlist[i] );
  }

  return nlist;
}
