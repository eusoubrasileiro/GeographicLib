/**
 * \file SphericalHarmonic2.hpp
 * \brief Header for GeographicLib::SphericalHarmonic2 class
 *
 * Copyright (c) Charles Karney (2011, 2012) <charles@karney.com> and licensed
 * under the MIT/X11 License.  For more information, see
 * http://geographiclib.sourceforge.net/
 **********************************************************************/

#if !defined(GEOGRAPHICLIB_SPHERICALHARMONIC2_HPP)
#define GEOGRAPHICLIB_SPHERICALHARMONIC2_HPP \
  "$Id: ce4cda614c1966dea65610bc73bc4db562677fa8 $"

#include <vector>
#include <GeographicLib/Constants.hpp>
#include <GeographicLib/SphericalEngine.hpp>
#include <GeographicLib/CircularEngine.hpp>

namespace GeographicLib {

  /**
   * \brief Spherical Harmonic series with two corrections to the coefficients.
   *
   * This classes is similar to SphericalHarmonic, except that the coefficients
   * \e C<sub>\e nm</sub> are replaced by \e C<sub>\e nm</sub> + \e tau'
   * C'<sub>\e nm</sub> + \e tau'' C''<sub>\e nm</sub> (and similarly for \e
   * S<sub>\e nm</sub>).
   *
   * Example of use:
   * \include example-SphericalHarmonic2.cpp
   **********************************************************************/

  class GEOGRAPHIC_EXPORT SphericalHarmonic2 {
  public:
    /**
     * Supported normalizations for associate Legendre polynomials.
     **********************************************************************/
    enum normalization {
      /**
       * Fully normalized associated Legendre polynomials.  See
       * SphericalHarmonic::FULL for documentation.
       *
       * @hideinitializer
       **********************************************************************/
      FULL = SphericalEngine::FULL,
      /**
       * Schmidt semi-normalized associated Legendre polynomials.  See
       * SphericalHarmonic::SCHMIDT for documentation.
       *
       * @hideinitializer
       **********************************************************************/
      SCHMIDT = SphericalEngine::SCHMIDT,
      /// \cond SKIP
      // These are deprecated...
      full = FULL,
      schmidt = SCHMIDT,
      /// \endcond
    };

  private:
    typedef Math::real real;
    SphericalEngine::coeff _c[3];
    real _a;
    unsigned _norm;

  public:
    /**
     * Constructor with a full set of coefficients specified.
     *
     * @param[in] C the coefficients \e C<sub>\e nm</sub>.
     * @param[in] S the coefficients \e S<sub>\e nm</sub>.
     * @param[in] N the maximum degree and order of the sum
     * @param[in] C1 the coefficients \e C'<sub>\e nm</sub>.
     * @param[in] S1 the coefficients \e S'<sub>\e nm</sub>.
     * @param[in] N1 the maximum degree and order of the first correction
     *   coefficients \e C'<sub>\e nm</sub> and \e S'<sub>\e nm</sub>.
     * @param[in] C2 the coefficients \e C''<sub>\e nm</sub>.
     * @param[in] S2 the coefficients \e S''<sub>\e nm</sub>.
     * @param[in] N2 the maximum degree and order of the second correction
     *   coefficients \e C'<sub>\e nm</sub> and \e S'<sub>\e nm</sub>.
     * @param[in] a the reference radius appearing in the definition of the
     *   sum.
     * @param[in] norm the normalization for the associated Legendre
     *   polynomials, either SphericalHarmonic2::FULL (the default) or
     *   SphericalHarmonic2::SCHMIDT.
     *
     * See SphericalHarmonic for the way the coefficients should be stored.  \e
     * N1 and \e N2 should satisfy \e N1 <= \e N and \e N2 <= \e N.
     *
     * The class stores <i>pointers</i> to the first elements of \e C, \e S, \e
     * C', \e S', \e C'', and \e S''.  These arrays should not be altered or
     * destroyed during the lifetime of a SphericalHarmonic object.
     **********************************************************************/
    SphericalHarmonic2(const std::vector<real>& C,
                       const std::vector<real>& S,
                       int N,
                       const std::vector<real>& C1,
                       const std::vector<real>& S1,
                       int N1,
                       const std::vector<real>& C2,
                       const std::vector<real>& S2,
                       int N2,
                       real a, unsigned norm = FULL)
      : _a(a)
      , _norm(norm) {
      if (!(N1 <= N && N2 <= N))
        throw GeographicErr("N1 and N2 cannot be larger that N");
      _c[0] = SphericalEngine::coeff(C, S, N);
      _c[1] = SphericalEngine::coeff(C1, S1, N1);
      _c[2] = SphericalEngine::coeff(C2, S2, N2);
    }

    /**
     * Constructor with a subset of coefficients specified.
     *
     * @param[in] C the coefficients \e C<sub>\e nm</sub>.
     * @param[in] S the coefficients \e S<sub>\e nm</sub>.
     * @param[in] N the degree used to determine the layout of \e C and \e S.
     * @param[in] nmx the maximum degree used in the sum.  The sum over \e n is
     *   from 0 thru \e nmx.
     * @param[in] mmx the maximum order used in the sum.  The sum over \e m is
     *   from 0 thru min(\e n, \e mmx).
     * @param[in] C1 the coefficients \e C'<sub>\e nm</sub>.
     * @param[in] S1 the coefficients \e S'<sub>\e nm</sub>.
     * @param[in] N1 the degree used to determine the layout of \e C' and \e
     *   S'.
     * @param[in] nmx1 the maximum degree used for \e C' and \e S'.
     * @param[in] mmx1 the maximum order used for \e C' and \e S'.
     * @param[in] C2 the coefficients \e C''<sub>\e nm</sub>.
     * @param[in] S2 the coefficients \e S''<sub>\e nm</sub>.
     * @param[in] N2 the degree used to determine the layout of \e C'' and \e
     *   S''.
     * @param[in] nmx2 the maximum degree used for \e C'' and \e S''.
     * @param[in] mmx2 the maximum order used for \e C'' and \e S''.
     * @param[in] a the reference radius appearing in the definition of the
     *   sum.
     * @param[in] norm the normalization for the associated Legendre
     *   polynomials, either SphericalHarmonic2::FULL (the default) or
     *   SphericalHarmonic2::SCHMIDT.
     *
     * The class stores <i>pointers</i> to the first elements of \e C, \e S, \e
     * C', \e S', \e C'', and \e S''.  These arrays should not be altered or
     * destroyed during the lifetime of a SphericalHarmonic object.
     **********************************************************************/
    SphericalHarmonic2(const std::vector<real>& C,
                       const std::vector<real>& S,
                       int N, int nmx, int mmx,
                       const std::vector<real>& C1,
                       const std::vector<real>& S1,
                       int N1, int nmx1, int mmx1,
                       const std::vector<real>& C2,
                       const std::vector<real>& S2,
                       int N2, int nmx2, int mmx2,
                       real a, unsigned norm = FULL)
      : _a(a)
      , _norm(norm) {
      if (!(nmx1 <= nmx && nmx2 <= nmx))
        throw GeographicErr("nmx1 and nmx2 cannot be larger that nmx");
      if (!(mmx1 <= mmx && mmx2 <= mmx))
        throw GeographicErr("mmx1 and mmx2cannot be larger that mmx");
      _c[0] = SphericalEngine::coeff(C, S, N, nmx, mmx);
      _c[1] = SphericalEngine::coeff(C1, S1, N1, nmx1, mmx1);
      _c[2] = SphericalEngine::coeff(C2, S2, N2, nmx2, mmx2);
    }

    /**
     * A default constructor so that the object can be created when the
     * constructor for another object is initialized.  This default object can
     * then be reset with the default copy assignment operator.
     **********************************************************************/
    SphericalHarmonic2() {}

    /**
     * Compute a spherical harmonic sum with two correction terms.
     *
     * @param[in] tau1 multiplier for correction coefficients \e C' and \e S'.
     * @param[in] tau2 multiplier for correction coefficients \e C'' and \e S''.
     * @param[in] x cartesian coordinate.
     * @param[in] y cartesian coordinate.
     * @param[in] z cartesian coordinate.
     * @return \e V the spherical harmonic sum.
     *
     * This routine requires constant memory and thus never throws an
     * exception.
     **********************************************************************/
    Math::real operator()(real tau1, real tau2, real x, real y, real z)
      const throw() {
      real f[] = {1, tau1, tau2};
      real v = 0;
      real dummy;
      switch (_norm) {
      case FULL:
        v = SphericalEngine::Value<false, SphericalEngine::FULL, 3>
          (_c, f, x, y, z, _a, dummy, dummy, dummy);
        break;
      case SCHMIDT:
        v = SphericalEngine::Value<false, SphericalEngine::SCHMIDT, 3>
          (_c, f, x, y, z, _a, dummy, dummy, dummy);
        break;
      }
      return v;
    }

    /**
     * Compute a spherical harmonic sum with two correction terms and its
     * gradient.
     *
     * @param[in] tau1 multiplier for correction coefficients \e C' and \e S'.
     * @param[in] tau2 multiplier for correction coefficients \e C'' and \e S''.
     * @param[in] x cartesian coordinate.
     * @param[in] y cartesian coordinate.
     * @param[in] z cartesian coordinate.
     * @param[out] gradx \e x component of the gradient
     * @param[out] grady \e y component of the gradient
     * @param[out] gradz \e z component of the gradient
     * @return \e V the spherical harmonic sum.
     *
     * This is the same as the previous function, except that the components of
     * the gradients of the sum in the \e x, \e y, and \e z directions are
     * computed.  This routine requires constant memory and thus never throws
     * an exception.
     **********************************************************************/
    Math::real operator()(real tau1, real tau2, real x, real y, real z,
                          real& gradx, real& grady, real& gradz) const throw() {
      real f[] = {1, tau1, tau2};
      real v = 0;
      switch (_norm) {
      case FULL:
        v = SphericalEngine::Value<true, SphericalEngine::FULL, 3>
          (_c, f, x, y, z, _a, gradx, grady, gradz);
        break;
      case SCHMIDT:
        v = SphericalEngine::Value<true, SphericalEngine::SCHMIDT, 3>
          (_c, f, x, y, z, _a, gradx, grady, gradz);
        break;
      }
      return v;
    }

    /**
     * Create a CircularEngine to allow the efficient evaluation of several
     * points on a circle of latitude at fixed values of \e tau1 and \e tau2.
     *
     * @param[in] tau1 multiplier for correction coefficients \e C' and \e S'.
     * @param[in] tau2 multiplier for correction coefficients \e C'' and \e S''.
     * @param[in] p the radius of the circle.
     * @param[in] z the height of the circle above the equatorial plane.
     * @param[in] gradp if true the returned object will be able to compute the
     *   gradient of the sum.
     * @return the CircularEngine object.
     *
     * SphericalHarmonic2::operator()() exchanges the order of the sums in the
     * definition, i.e., sum(n = 0..N)[sum(m = 0..n)[...]] becomes sum(m =
     * 0..N)[sum(n = m..N)[...]].  SphericalHarmonic2::Circle performs the
     * inner sum over degree \e n (which entails about <i>N</i><sup>2</sup>
     * operations).  Calling CircularEngine::operator()() on the returned
     * object performs the outer sum over the order \e m (about \e N
     * operations).  This routine may throw a bad_alloc exception in the
     * CircularEngine constructor.
     *
     * See SphericalHarmonic::Circle for an example of its use.
     **********************************************************************/
    CircularEngine Circle(real tau1, real tau2, real p, real z, bool gradp)
      const {
      real f[] = {1, tau1, tau2};
      switch (_norm) {
      case FULL:
        return gradp ?
          SphericalEngine::Circle<true, SphericalEngine::FULL, 3>
          (_c, f, p, z, _a) :
          SphericalEngine::Circle<false, SphericalEngine::FULL, 3>
          (_c, f, p, z, _a);
        break;
      case SCHMIDT:
      default:                  // To avoid compiler warnings
        return gradp ?
          SphericalEngine::Circle<true, SphericalEngine::SCHMIDT, 3>
          (_c, f, p, z, _a) :
          SphericalEngine::Circle<false, SphericalEngine::SCHMIDT, 3>
          (_c, f, p, z, _a);
        break;
      }
    }

    /**
     * @return the zeroth SphericalEngine::coeff object.
     **********************************************************************/
    const SphericalEngine::coeff& Coefficients() const throw()
    { return _c[0]; }
    /**
     * @return the first SphericalEngine::coeff object.
     **********************************************************************/
    const SphericalEngine::coeff& Coefficients1() const throw()
    { return _c[1]; }
    /**
     * @return the second SphericalEngine::coeff object.
     **********************************************************************/
    const SphericalEngine::coeff& Coefficients2() const throw()
    { return _c[2]; }
  };

} // namespace GeographicLib

#endif  // GEOGRAPHICLIB_SPHERICALHARMONIC2_HPP