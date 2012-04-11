# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <ctime>
# include <cmath>
# include <complex>

using namespace std;

# include "normal.hpp"



//****************************************************************************80

float r4_normal ( float a, float b, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    R4_NORMAL returns a scaled pseudonormal R4.
//
//  Discussion:
//
//    The normal probability distribution function (PDF) is sampled,
//    with mean A and standard deviation B.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 June 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, float A, the mean of the PDF.
//
//    Input, float B, the standard deviation of the PDF.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, float R4_NORMAL, a sample of the normal PDF.
//
{
  float value;

  value = a + b * r4_normal_01 ( seed );

  return value;
}
//****************************************************************************80

float r4_normal_01 ( int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    R4_NORMAL_01 returns a unit pseudonormal R4.
//
//  Discussion:
//
//    The standard normal probability distribution function (PDF) has
//    mean 0 and standard deviation 1.
//
//    The Box-Muller method is used, which is efficient, but
//    generates two values at a time.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, float R4_NORMAL_01, a normally distributed random value.
//
{
# define R4_PI 3.1415926

  float r1;
  float r2;
  static int used = -1;
  static int seed1 = 0;
  static int seed2 = 0;
  static int seed3 = 0;
  float x;
  static float y = 0.0;

  if ( used == -1 )
  {
    used = 0;
  }
//
//  If we've used an even number of values so far, generate two more, return one,
//  and save one.
//
  if ( ( used % 2 ) == 0 )
  {
    seed1 = *seed;
    r1 = r4_uniform_01 ( seed );

    if ( r1 == 0.0 )
    {
      cerr << "\n";
      cerr << "R4_NORMAL_01 - Fatal error!\n";
      cerr << "  R4_UNIFORM_01 returned a value of 0.\n";
      exit ( 1 );
    }

    seed2 = *seed;
    r2 = r4_uniform_01 ( seed );
    seed3 = *seed;
    *seed = seed2;

    x = sqrt ( - 2.0 * log ( r1 ) ) * cos ( 2.0 * R4_PI * r2 );
    y = sqrt ( - 2.0 * log ( r1 ) ) * sin ( 2.0 * R4_PI * r2 );
  }
//
//  Otherwise, return the second, saved, value and the corresponding
//  value of SEED.
//
  else
  {
    x = y;
    *seed = seed3;
  }

  used = used + 1;

  return x;
# undef R4_PI
}
//****************************************************************************80

float r4_uniform_01 ( int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    R4_UNIFORM_01 returns a unit pseudorandom R4.
//
//  Discussion:
//
//    This routine implements the recursion
//
//      seed = 16807 * seed mod ( 2**31 - 1 )
//      r4_uniform_01 = seed / ( 2**31 - 1 )
//
//    The integer arithmetic never requires more than 32 bits,
//    including a sign bit.
//
//    If the initial seed is 12345, then the first three computations are
//
//      Input     Output      R4_UNIFORM_01
//      SEED      SEED
//
//         12345   207482415  0.096616
//     207482415  1790989824  0.833995
//    1790989824  2035175616  0.947702
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 November 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Paul Bratley, Bennett Fox, Linus Schrage,
//    A Guide to Simulation,
//    Springer Verlag, pages 201-202, 1983.
//
//    Pierre L'Ecuyer,
//    Random Number Generation,
//    in Handbook of Simulation
//    edited by Jerry Banks,
//    Wiley Interscience, page 95, 1998.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, pages 362-376, 1986.
//
//    Peter Lewis, Allen Goodman, James Miller,
//    A Pseudo-Random Number Generator for the System/360,
//    IBM Systems Journal,
//    Volume 8, pages 136-143, 1969.
//
//  Parameters:
//
//    Input/output, int *SEED, the "seed" value.  Normally, this
//    value should not be 0.  On output, SEED has been updated.
//
//    Output, float R4_UNIFORM_01, a new pseudorandom variate, strictly between
//    0 and 1.
//
{
  int k;
  float r;

  k = *seed / 127773;

  *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

  if ( *seed < 0 )
  {
    *seed = *seed + 2147483647;
  }
//
//  Although SEED can be represented exactly as a 32 bit integer,
//  it generally cannot be represented exactly as a 32 bit real number!
//
  r = ( float ) ( *seed ) * 4.656612875E-10;

  return r;
}
//****************************************************************************80

double r8_normal ( double a, double b, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8_NORMAL returns a scaled pseudonormal R8.
//
//  Discussion:
//
//    The normal probability distribution function (PDF) is sampled,
//    with mean A and standard deviation B.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 November 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, the mean of the PDF.
//
//    Input, double B, the standard deviation of the PDF.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double R8_NORMAL, a sample of the normal PDF.
//
{
  double value;

  value = a + b * r8_normal_01 ( seed );

  return value;
}
//****************************************************************************80

double r8_normal_01 ( int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8_NORMAL_01 returns a unit pseudonormal R8.
//
//  Discussion:
//
//    The standard normal probability distribution function (PDF) has
//    mean 0 and standard deviation 1.
//
//    Because this routine uses the Box Muller method, it requires pairs
//    of uniform random values to generate a pair of normal random values.
//    This means that on every other call, the code can use the second
//    value that it calculated.
//
//    However, if the user has changed the SEED value between calls,
//    the routine automatically resets itself and discards the saved data.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double R8_NORMAL_01, a normally distributed random value.
//
{
# define R8_PI 3.141592653589793

  double r1;
  double r2;
  static int seed1 = 0;
  static int seed2 = 0;
  static int seed3 = 0;
  static int used = 0;
  double v1;
  static double v2 = 0.0;
//
//  If USED is odd, but the input SEED does not match
//  the output SEED on the previous call, then the user has changed
//  the seed.  Wipe out internal memory.
//
  if ( ( used % 2 ) == 1 )
  {
    if ( *seed != seed2 )
    {
      used = 0;
      seed1 = 0;
      seed2 = 0;
      seed3 = 0;
      v2 = 0.0;
    }
  }
//
//  If USED is even, generate two uniforms, create two normals,
//  return the first normal and its corresponding seed.
//
  if ( ( used % 2 ) == 0 )
  {
    seed1 = *seed;

    r1 = r8_uniform_01 ( seed );

    if ( r1 == 0.0 )
    {
      cerr << "\n";
      cerr << "R8_NORMAL_01 - Fatal error!\n";
      cerr << "  R8_UNIFORM_01 returned a value of 0.\n";
      exit ( 1 );
    }

    seed2 = *seed;
    r2 = r8_uniform_01 ( seed );
    seed3 = *seed;
    *seed = seed2;

    v1 = sqrt ( - 2.0 * log ( r1 ) ) * cos ( 2.0 * R8_PI * r2 );
    v2 = sqrt ( - 2.0 * log ( r1 ) ) * sin ( 2.0 * R8_PI * r2 );
  }
//
//  If USED is odd (and the input SEED matched the output value from
//  the previous call), return the second normal and its corresponding seed.
//
  else
  {
    v1 = v2;
    *seed = seed3;
  }

  used = used + 1;

  return v1;
# undef R8_PI
}
//****************************************************************************80

double r8_uniform_01 ( int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8_UNIFORM_01 returns a unit pseudorandom R8.
//
//  Discussion:
//
//    This routine implements the recursion
//
//      seed = 16807 * seed mod ( 2**31 - 1 )
//      r8_uniform_01 = seed / ( 2**31 - 1 )
//
//    The integer arithmetic never requires more than 32 bits,
//    including a sign bit.
//
//    If the initial seed is 12345, then the first three computations are
//
//      Input     Output      R8_UNIFORM_01
//      SEED      SEED
//
//         12345   207482415  0.096616
//     207482415  1790989824  0.833995
//    1790989824  2035175616  0.947702
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Paul Bratley, Bennett Fox, Linus Schrage,
//    A Guide to Simulation,
//    Springer Verlag, pages 201-202, 1983.
//
//    Pierre L'Ecuyer,
//    Random Number Generation,
//    in Handbook of Simulation
//    edited by Jerry Banks,
//    Wiley Interscience, page 95, 1998.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, pages 362-376, 1986.
//
//    Peter Lewis, Allen Goodman, James Miller,
//    A Pseudo-Random Number Generator for the System/360,
//    IBM Systems Journal,
//    Volume 8, pages 136-143, 1969.
//
//  Parameters:
//
//    Input/output, int *SEED, the "seed" value.  Normally, this
//    value should not be 0.  On output, SEED has been updated.
//
//    Output, double R8_UNIFORM_01, a new pseudorandom variate, strictly between
//    0 and 1.
//
{
  int k;
  double r;

  k = *seed / 127773;

  *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

  if ( *seed < 0 )
  {
    *seed = *seed + 2147483647;
  }
//
//  Although SEED can be represented exactly as a 32 bit integer,
//  it generally cannot be represented exactly as a 32 bit real number!
//
  r = ( double ) ( *seed ) * 4.656612875E-10;

  return r;
}
