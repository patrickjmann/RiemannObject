//##############################################################
/**\file eqs.h

$Source: /home/mann/local/research/1/rshock/RiemannObject/RCS/eos.h,v $
$Revision: 1.1 $
$Date: 2010/05/09 16:33:58 $

\brief Equation of State Include file for my C++ version of Marti/Muller Riemann solver

\author P. J. Mann
*/
//##############################################################
#if ! defined RIEMANN_EQS_H
#define RIEMANN_EQS_H

//==============================================================
/** Equation of state

  p = pressure
  rho0 = rest density (not including internal energy)
  eps = epsilon = specific internal energy density

For speed there are NO TESTS!
 */
class EqOfState
{
 public:
  double gamma;      /// Polytropic index
  double gam1;       /// gamma-1 (set in constructors)

  inline void Make( double gamma_in ){
    if( gamma_in < 1.0 ){
      cerr << "EqOfState: ERROR: gamma < 1.0\n"; exit(1);
    }
    gamma = gamma_in;
    gam1 = gamma - 1.0;
  };
  inline EqOfState(){
    double default_gamma = 5.0/3.0; Make( default_gamma ); 
  };
  inline EqOfState( double gamma_in ){ Make( gamma_in ); };

  inline double P( double rho0, double eps ){
    return ( gam1*rho0*eps );
  };
  inline double EPS( double p, double rho0 ){
    return ( p/(gam1*rho0) );
  };
  inline double CS( double p, double rho0, double eps ){
    return ( sqrt( gamma * p /( rho0*eps ) ) );
  };
};
#endif
