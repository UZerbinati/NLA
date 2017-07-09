from NLA import *

def r8_uniform_01 ( seed ):

#*****************************************************************************80
#
## R8_UNIFORM_01 returns a unit pseudorandom R8.
#
#  Discussion:
#
#    This routine implements the recursion
#
#      seed = 16807 * seed mod ( 2^31 - 1 )
#      r8_uniform_01 = seed / ( 2^31 - 1 )
#
#    The integer arithmetic never requires more than 32 bits,
#    including a sign bit.
#
#    If the initial seed is 12345, then the first three computations are
#
#      Input     Output      R8_UNIFORM_01
#      SEED      SEED
#
#         12345   207482415  0.096616
#     207482415  1790989824  0.833995
#    1790989824  2035175616  0.947702
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    17 March 2013
#
#  Author:
#
#    John Burkardt
#
#  Reference:
#
#    Paul Bratley, Bennett Fox, Linus Schrage,
#    A Guide to Simulation,
#    Second Edition,
#    Springer, 1987,
#    ISBN: 0387964673,
#    LC: QA76.9.C65.B73.
#
#    Bennett Fox,
#    Algorithm 647:
#    Implementation and Relative Efficiency of Quasirandom
#    Sequence Generators,
#    ACM Transactions on Mathematical Software,
#    Volume 12, Number 4, December 1986, pages 362-376.
#
#    Pierre L'Ecuyer,
#    Random Number Generation,
#    in Handbook of Simulation,
#    edited by Jerry Banks,
#    Wiley, 1998,
#    ISBN: 0471134031,
#    LC: T57.62.H37.
#
#    Peter Lewis, Allen Goodman, James Miller,
#    A Pseudo-Random Number Generator for the System/360,
#    IBM Systems Journal,
#    Volume 8, Number 2, 1969, pages 136-143.
#
#  Parameters:
#
#    Input, integer SEED, the integer "seed" used to generate
#    the output random number.  SEED should not be 0.
#
#    Output, real R, a random value between 0 and 1.
#
#    Output, integer SEED, the updated seed.  This would
#    normally be used as the input seed on the next call.
#
  from sys import exit

  i4_huge = 2147483647

  seed = int ( seed )

  seed = ( seed % i4_huge )

  if ( seed < 0 ):
    seed = seed + i4_huge

  if ( seed == 0 ):
    print ( '' )
    print ( 'R8_UNIFORM_01 - Fatal error!' )
    print ( '  Input SEED = 0!' )
    exit ( 'R8_UNIFORM_01 - Fatal error!' )

  k = ( seed // 127773 )

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ):
    seed = seed + i4_huge

  r = seed * 4.656612875E-10

  return r, seed

def r8_uniform_01_test ( ):

#*****************************************************************************80
#
## R8_UNIFORM_01_TEST tests R8_UNIFORM_01.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    26 July 2014
#
#  Author:
#
#    John Burkardt
#
  import platform

  print ( '' )
  print ( 'R8_UNIFORM_01_TEST' )
  print ( '  Python version: %s' % ( platform.python_version ( ) ) )
  print ( '  R8_UNIFORM_01 produces a sequence of random values.' )

  seed = 123456789

  print ( '' )
  print ( '  Using random seed %d' % ( seed ) )

  print ( '' )
  print ( '  SEED  R8_UNIFORM_01(SEED)' )
  print ( '' )
  for i in range ( 0, 10 ):
    seed_old = seed
    x, seed = r8_uniform_01 ( seed )
    print ( '  %12d  %14f' % ( seed, x ) )

  print ( '' )
  print ( '  Verify that the sequence can be restarted.' )
  print ( '  Set the seed back to its original value, and see that' )
  print ( '  we generate the same sequence.' )

  seed = 123456789
  print ( '' )
  print ( '  SEED  R8_UNIFORM_01(SEED)' )
  print ( '' )

  for i in range ( 0, 10 ):
    seed_old = seed
    x, seed = r8_uniform_01 ( seed )
    print ( '  %12d  %14f' % ( seed, x ) )
#
#  Terminate.
#
  print ( '' )
  print ( 'R8_UNIFORM_01_TEST' )
  print ( '  Normal end of execution.' )
  return

if ( __name__ == '__main__' ):
  from timestamp import timestamp
  timestamp ( )
  r8_uniform_01_test ( )
  timestamp ( )

def wathen_csr ( nx, ny, seed=123456789 ):

#*****************************************************************************80
#
## WATHEN_CSC: Wathen matrix stored in Compressed Sparse Column (CSC) format.
#
#  Discussion:
#
#    When dealing with sparse matrices in MATLAB, it can be much more efficient
#    to work first with a triple of I, J, and X vectors, and only once
#    they are complete, convert to MATLAB's sparse format.
#
#    The Wathen matrix is a finite element matrix which is sparse.
#
#    The entries of the matrix depend in part on a physical quantity
#    related to density.  That density is here assigned random values between
#    0 and 100.
#
#    The matrix order N is determined by the input quantities NX and NY,
#    which would usually be the number of elements in the X and Y directions.
#
#    The value of N is
#
#      N = 3*NX*NY + 2*NX + 2*NY + 1,
#
#    The matrix is the consistent mass matrix for a regular NX by NY grid
#    of 8 node serendipity elements.
#
#    The local element numbering is
#
#      3--2--1
#      |     |
#      4     8
#      |     |
#      5--6--7
#
#    Here is an illustration for NX = 3, NY = 2:
#
#     23-24-25-26-27-28-29
#      |     |     |     |
#     19    20    21    22
#      |     |     |     |
#     12-13-14-15-16-17-18
#      |     |     |     |
#      8     9    10    11
#      |     |     |     |
#      1--2--3--4--5--6--7
#
#    For this example, the total number of nodes is, as expected,
#
#      N = 3 * 3 * 2 + 2 * 2 + 2 * 3 + 1 = 29
#
#    The matrix is symmetric positive definite for any positive values of the
#    density RHO(X,Y).
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    01 September 2014
#
#  Author:
#
#    John Burkardt.
#
#  Reference:
#
#    Nicholas Higham,
#    Algorithm 694: A Collection of Test Matrices in MATLAB,
#    ACM Transactions on Mathematical Software,
#    Volume 17, Number 3, September 1991, pages 289-305.
#
#    Andrew Wathen,
#    Realistic eigenvalue bounds for the Galerkin mass matrix,
#    IMA Journal of Numerical Analysis,
#    Volume 7, Number 4, October 1987, pages 449-457.
#
#  Parameters:
#
#    Input, integer NX, NY, values which determine the size of the matrix.
#
#    Input, integer SEED, the random number seed.
#
#    Output, real A(*,*), the compressed sparxe column version of the matrix.
#
#    Output, integer SEED, the random number seed.
#
  import numpy as np
  from scipy.sparse import coo_matrix

  em = np.array \
  ( \
    ( ( 6.0, -6.0,  2.0, -8.0,  3.0, -8.0,  2.0, -6.0 ), \
      (-6.0, 32.0, -6.0, 20.0, -8.0, 16.0, -8.0, 20.0 ), \
      ( 2.0, -6.0,  6.0, -6.0,  2.0, -8.0,  3.0, -8.0 ), \
      (-8.0, 20.0, -6.0, 32.0, -6.0, 20.0, -8.0, 16.0 ), \
      ( 3.0, -8.0,  2.0, -6.0,  6.0, -6.0,  2.0, -8.0 ), \
      (-8.0, 16.0, -8.0, 20.0, -6.0, 32.0, -6.0, 20.0 ), \
      ( 2.0, -8.0,  3.0, -8.0,  2.0, -6.0,  6.0, -6.0 ), \
      (-6.0, 20.0, -8.0, 16.0, -8.0, 20.0, -6.0, 32.0 ) )\
  )

  node = np.zeros ( 8 )

  st_num = 8 * 8 * nx * ny
  row = np.zeros ( st_num )
  col = np.zeros ( st_num )
  v = np.zeros ( st_num )

  k = 0

  for j in range ( 0, ny ):

    for i in range ( 0, nx ):

      node[0] = ( 3 * ( j + 1 )     ) * nx + 2 * ( j + 1 ) + 2 *  ( i + 1 ) + 1 - 1
      node[1] = node[0] - 1
      node[2] = node[0] - 2

      node[3] = ( 3 * ( j + 1 ) - 1 ) * nx + 2 *  ( j + 1 ) + ( i + 1 ) - 1 - 1
      node[7] = node[3] + 1

      node[4] = ( 3 *  ( j + 1 ) - 3 ) * nx + 2 *  ( j + 1 ) + 2 *  ( i + 1 ) - 3 - 1
      node[5] = node[4] + 1
      node[6] = node[4] + 2

      rho, seed = r8_uniform_01 ( seed )
      rho = 100.0 * rho

      for krow in range ( 0, 8 ):
        for kcol in range ( 0, 8 ):
          row[k] = node[krow]
          col[k] = node[kcol]
          v[k] = rho * em[krow,kcol]
          k = k + 1
#
#  Convert triplet to a Python COO matrix.
#
  a = coo_matrix ( ( v, ( row, col ) ) )
#
#  Convert COO matrix to CSC format.
#
  a = a.tocsr ( )
  A = Matrix(a,Sparse=True)
  return A, seed
