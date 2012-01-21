/**
 * \file Accumulator.hpp
 * \brief Header for GeographicLib::Accumulator class
 *
 * Copyright (c) Charles Karney (2010, 2011) <charles@karney.com> and licensed
 * under the MIT/X11 License.  For more information, see
 * http://geographiclib.sourceforge.net/
 **********************************************************************/

#if !defined(GEOGRAPHICLIB_ACCUMULATOR_HPP)
#define GEOGRAPHICLIB_ACCUMULATOR_HPP \
  "$Id: 03b7f4fdb9974c822f98d5f5aab1094b5a9ac8f2 $"

#include <GeographicLib/Constants.hpp>

namespace GeographicLib {

  /**
   * \brief An accumulator for sums.
   *
   * This allow many numbers of floating point type \e T to be added together
   * with twice the normal precision.  Thus if \e T is double, the effective
   * precision of the sum is 106 bits or about 32 decimal places.  The core
   * idea is the error free transformation of a sum, D. E. Knuth, TAOCP, Vol 2,
   * 4.2.2, Theorem B.
   *
   * The implementation follows J. R. Shewchuk,
   * <a href="http://dx.doi.org/10.1007/PL00009321"> Adaptive Precision
   * Floating-Point Arithmetic and Fast Robust Geometric Predicates</a>,
   * Discrete & Computational Geometry 18(3) 305-363 (1997).
   *
   * Approximate timings (summing a vector<double>)
   * - double:               2ns
   * - Accumulator<double>: 23ns
   *
   * In the documentation of the member functions, \e sum stands for the value
   * currently held in the accumulator.
   *
   * Example of use:
   * \include example-Accumulator.cpp
   **********************************************************************/
  template<typename T = Math::real>
  class GEOGRAPHIC_EXPORT Accumulator {
  private:
    // _s + _t accumulates for the sum.
    T _s, _t;
    // Error free transformation of a sum.  Note that t can be the same as one
    // of the first two arguments.
    static inline T sum(T u, T v, T& t) {
      volatile T s = u + v;
      volatile T up = s - v;
      volatile T vpp = s - up;
      up -= u;
      vpp -= v;
      t = -(up + vpp);
      // u + v =       s      + t
      //       = round(u + v) + t
      return s;
    }
    // Same as sum, but requires abs(u) >= abs(v).  This isn't currently used.
    static inline T fastsum(T u, T v, T& t) {
      volatile T s = u + v;
      volatile T vp = s - u;
      t = v - vp;
      return s;
    }
    void Add(T y) throw() {
      // Here's Shewchuk's solution...
      T u;                      // hold exact sum as [s, t, u]
      y  = sum(y, _t,  u);      // Accumulate starting at least significant end
      _s = sum(y, _s, _t);
      // Start is _s, _t decreasing and non-adjacent.  Sum is now (s + t + u)
      // exactly with s, t, u non-adjacent and in decreasing order (except for
      // possible zeros).  The following code tries to normalize the result.
      // Ideally, we want _s = round(s+t+u) and _u = round(s+t+u - _s).  The
      // following does an approximate job (and maintains the decreasing
      // non-adjacent property).  Here are two "failures" using 3-bit floats:
      //
      // Case 1: _s is not equal to round(s+t+u) -- off by 1 ulp
      // [12, -1] - 8 -> [4, 0, -1] -> [4, -1] = 3 should be [3, 0] = 3
      //
      // Case 2: _s+_t is not as close to s+t+u as it shold be
      // [64, 5] + 4 -> [64, 8, 1] -> [64,  8] = 72 (off by 1)
      //                    should be [80, -7] = 73 (exact)
      //
      // "Fixing" these problems is probably not worth the expense.  The
      // representation inevitably leads to small errors in the accumulated
      // values.  The additional errors illustrated here amount to 1 ulp of the
      // less significant word during each addition to the Accumulator and an
      // additional possible error of 1 ulp in the reported sum.
      //
      // Incidentally, the "ideal" representation described above is not
      // canonical, because _s = round(_s + _t) may not be true.  For example,
      // with 3-bit floats:
      //
      // [128, 16] + 1 -> [160, -16] -- 160 = round(145).
      // But [160, 0] - 16 -> [128, 16] -- 128 = round(144).
      //
      if (_s == 0)              // This implies t == 0,
        _s = u;                 // so result is u
      else
        _t += u;                // otherwise just accumulate u to t.
    }
    T Sum(T y) const throw() {
      Accumulator a(*this);
      a.Add(y);
      return a._s;
    }
  public:
    /**
     * Construct from a \e T.  This is not declared explicit, so that you can
     * write <code>Accumulator<double> a = 5;</code>.
     *
     * @param[in] y set \e sum = \e y.
     **********************************************************************/
    Accumulator(T y = T(0)) throw() : _s(y), _t(0) {
      STATIC_ASSERT(!std::numeric_limits<T>::is_integer,
                    "Accumulator type is not floating point");
    }
    /**
     * Set the accumulator to a number.
     *
     * @param[in] y set \e sum = \e y.
     **********************************************************************/
    Accumulator& operator=(T y) throw() { _s = y; _t = 0; return *this; }
    /**
     * Return the value held in the accumulator.
     *
     * @return \e sum.
     **********************************************************************/
    T operator()() const throw() { return _s; }
    /**
     * Return the result of adding a number to \e sum (but don't change \e sum).
     *
     * @param[in] y the number to be added to the sum.
     * @return \e sum + \e y.
     **********************************************************************/
    T operator()(T y) const throw() { return Sum(y); }
    /**
     * Add a number to the accumulator.
     *
     * @param[in] y set \e sum += \e y.
     **********************************************************************/
    Accumulator& operator+=(T y) throw() { Add(y); return *this; }
    /**
     * Subtract a number from the accumulator.
     *
     * @param[in] y set \e sum -= \e y.
     **********************************************************************/
    Accumulator& operator-=(T y) throw() { Add(-y); return *this; }
    /**
     * Multiply accumulator by an integer.  To avoid loss of accuracy, use only
     * integers such that \e n * \e T is exactly representable as a \e T (i.e.,
     * +/- powers of two).  Use \e n = -1 to negate \e sum.
     *
     * @param[in] n set \e sum *= \e n.
     **********************************************************************/
    Accumulator& operator*=(int n) throw() { _s *= n; _t *= n; return *this; }
    /**
     * Test equality of an Accumulator with a number.
     **********************************************************************/
    bool operator==(T y) const throw() { return _s == y; }
    /**
     * Test inequality of an Accumulator with a number.
     **********************************************************************/
    bool operator!=(T y) const throw() { return _s != y; }
    /**
     * Less operator on an Accumulator and a number.
     **********************************************************************/
    bool operator<(T y) const throw() { return _s < y; }
    /**
     * Less or equal operator on an Accumulator and a number.
     **********************************************************************/
    bool operator<=(T y) const throw() { return _s <= y; }
    /**
     * Greater operator on an Accumulator and a number.
     **********************************************************************/
    bool operator>(T y) const throw() { return _s > y; }
    /**
     * Greater or equal operator on an Accumulator and a number.
     **********************************************************************/
    bool operator>=(T y) const throw() { return _s >= y; }
  };

} // namespace GeographicLib

#endif  // GEOGRAPHICLIB_ACCUMULATOR_HPP