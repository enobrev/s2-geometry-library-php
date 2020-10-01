<?php

namespace S2;

class S2 {
    public static function IEEERemainder($dividend, $divisor) {
        return $dividend - ($divisor * round($dividend / $divisor));
    }

    const M_PI = M_PI;
    const M_1_PI = M_1_PI;
    const M_PI_2 = M_PI_2;
    const M_PI_4 = M_PI_4;
    const M_SQRT2 = M_SQRT2;
    const M_E = M_E;

    // Together these flags define a cell orientation. If SWAP_MASK
    // is true, then canonical traversal order is flipped around the
    // diagonal (i.e. i and j are swapped with each other). If
    // INVERT_MASK is true, then the traversal order is rotated by 180
    // degrees (i.e. the bits of i and j are inverted, or equivalently,
    // the axis directions are reversed).
    //  public static final int SWAP_MASK = 0x01;
    //  public static final int INVERT_MASK = 0x02;

    // Number of bits in the mantissa of a double.
    const EXPONENT_SHIFT = 52;
    // Mask to extract the exponent from a double.
    const EXPONENT_MASK = 0x7ff0000000000000;

    /**
     * If v is non-zero, return an integer {@code exp} such that
     * {@code (0.5 <= |v|*2^(-exp) < 1)}. If v is zero, return 0.
     *
     * <p>Note that this arguably a bad definition of exponent because it makes
     * {@code exp(9) == 4}. In decimal this would be like saying that the
     * exponent of 1234 is 4, when in scientific 'exponent' notation 1234 is
     * {@code 1.234 x 10^3}.
     *
     * TODO(dbeaumont): Replace this with "DoubleUtils.getExponent(v) - 1" ?
     */
    public static function exp($v) {
        if ($v == 0) {
            return 0;
        }
        $bits = unpack('N2', strrev(pack('d', $v)));
        $bits = $bits[1] << 32 | $bits[2];
        return ((self::EXPONENT_MASK & $bits) >> self::EXPONENT_SHIFT) - 1022;
    }

    /** Mapping Hilbert traversal order to orientation adjustment mask. */
    private static $POS_TO_ORIENTATION = array(
        S2CellId::SWAP_MASK,
        0,
        0,
        S2CellId::INVERT_MASK + S2CellId::SWAP_MASK
    );

    /**
     * Returns an XOR bit mask indicating how the orientation of a child subcell
     * is related to the orientation of its parent cell. The returned value can
     * be XOR'd with the parent cell's orientation to give the orientation of
     * the child cell.
     *
     * @param int $position position of the subcell in the Hilber traversal,
     *                      in the range [0,3].
     * @return int a bit mask containing some combination of {@link #SWAP_MASK}
     *             and {@link #INVERT_MASK}.
     * @throws \Exception
     */
    public static function posToOrientation($position) {
        if (!(0 <= $position && $position < 4)) throw new \Exception();
        return self::$POS_TO_ORIENTATION[$position];
    }

    /** Mapping from cell orientation + Hilbert traversal to IJ-index. */
    private static $POS_TO_IJ = array(
// 0 1 2 3
        array(0, 1, 3, 2), // canonical order: (0,0), (0,1), (1,1), (1,0)
        array(0, 2, 3, 1), // axes swapped: (0,0), (1,0), (1,1), (0,1)
        array(3, 2, 0, 1), // bits inverted: (1,1), (1,0), (0,0), (0,1)
        array(3, 1, 0, 2), // swapped & inverted: (1,1), (0,1), (0,0), (1,0)
    );

    /**
     * Return the IJ-index of the subcell at the given position in the Hilbert
     * curve traversal with the given orientation. This is the inverse of
     * {@link #ijToPos}.
     *
     * @param int $orientation subcell orientation, in the range [0,3].
     * @param int $position position of the subcell in the Hilbert traversal,
     *                      in the range [0,3].
     * @return int IJ-index where
     * @throws \Exception
     */
    public static function posToIJ($orientation, $position) {
        if (!(0 <= $orientation && $orientation <= 3)) throw new \Exception();
        if (!(0 <= $position && $position <= 3)) throw new \Exception();
        return self::$POS_TO_IJ[$orientation][$position];
    }

    /** Mapping from Hilbert traversal order + cell orientation to IJ-index. */
//  private static final int IJ_TO_POS[][] = {
// (0,0) (0,1) (1,0) (1,1)
//      {0, 1, 3, 2}, // canonical order
//      {0, 3, 1, 2}, // axes swapped
//      {2, 3, 1, 0}, // bits inverted
//      {2, 1, 3, 0}, // swapped & inverted
//  };

    /**
     * Returns the order in which a specified subcell is visited by the Hilbert
     * curve. This is the inverse of {@link #posToIJ}.
     *
     * @param orientation the subcell orientation, in the range [0,3].
     * @param ijIndex the subcell index where
     *     {@code 0->(0,0), 1->(0,1), 2->(1,0), 3->(1,1)}.
     * @return the position of the subcell in the Hilbert traversal, in the range
     *     [0,3].
     * @throws IllegalArgumentException if either parameter is out of bounds.
     */
//  public static final int ijToPos(int orientation, int ijIndex) {
//    Preconditions.checkArgument(0 <= orientation && orientation < 4);
//    Preconditions.checkArgument(0 <= ijIndex && ijIndex < 4);
//    return IJ_TO_POS[orientation][ijIndex];
//  }

// Don't instantiate
    private function __construct() {}

    /**
     * Return a unique "origin" on the sphere for operations that need a fixed
     * reference point. It should *not* be a point that is commonly used in edge
     * tests in order to avoid triggering code to handle degenerate cases. (This
     * rules out the north and south poles.)
     *#/
     * public static S2Point origin() {
     * return new S2Point(0, 1, 0);
     * }
     *
	/**
     * Return true if the given point is approximately unit length (this is mainly
     * useful for assertions).
     */
    public static function isUnitLength(S2Point $p) {
    	return abs($p->norm2() - 1) <= 1e-15;
    }

    /**
     * Return true if edge AB crosses CD at a point that is interior to both
     * edges. Properties:
     *
     *  (1) SimpleCrossing(b,a,c,d) == SimpleCrossing(a,b,c,d) (2)
     * SimpleCrossing(c,d,a,b) == SimpleCrossing(a,b,c,d)
     *#/
     * public static boolean simpleCrossing(S2Point a, S2Point b, S2Point c, S2Point d) {
     * // We compute SimpleCCW() for triangles ACB, CBD, BDA, and DAC. All
     * // of these triangles need to have the same orientation (CW or CCW)
     * // for an intersection to exist. Note that this is slightly more
     * // restrictive than the corresponding definition for planar edges,
     * // since we need to exclude pairs of line segments that would
     * // otherwise "intersect" by crossing two antipodal points.
     *
     * S2Point ab = S2Point.crossProd(a, b);
     * S2Point cd = S2Point.crossProd(c, d);
     * double acb = -ab.dotProd(c);
     * double cbd = -cd.dotProd(b);
     * double bda = ab.dotProd(d);
     * double dac = cd.dotProd(a);
     *
     * return (acb * cbd > 0) && (cbd * bda > 0) && (bda * dac > 0);
     * }
     *
     * /**
     * Return a vector "c" that is orthogonal to the given unit-length vectors "a"
     * and "b". This function is similar to a.CrossProd(b) except that it does a
     * better job of ensuring orthogonality when "a" is nearly parallel to "b",
     * and it returns a non-zero result even when a == b or a == -b.
     *
     *  It satisfies the following properties (RCP == RobustCrossProd):
     *
     *  (1) RCP(a,b) != 0 for all a, b (2) RCP(b,a) == -RCP(a,b) unless a == b or
     * a == -b (3) RCP(-a,b) == -RCP(a,b) unless a == b or a == -b (4) RCP(a,-b)
     * == -RCP(a,b) unless a == b or a == -b
     *#/
     * public static S2Point robustCrossProd(S2Point a, S2Point b) {
     * // The direction of a.CrossProd(b) becomes unstable as (a + b) or (a - b)
     * // approaches zero. This leads to situations where a.CrossProd(b) is not
     * // very orthogonal to "a" and/or "b". We could fix this using Gram-Schmidt,
     * // but we also want b.RobustCrossProd(a) == -b.RobustCrossProd(a).
     * //
     * // The easiest fix is to just compute the cross product of (b+a) and (b-a).
     * // Given that "a" and "b" are unit-length, this has good orthogonality to
     * // "a" and "b" even if they differ only in the lowest bit of one component.
     *
     * // assert (isUnitLength(a) && isUnitLength(b));
     * S2Point x = S2Point.crossProd(S2Point.add(b, a), S2Point.sub(b, a));
     * if (!x.equals(new S2Point(0, 0, 0))) {
     * return x;
     * }
     *
     * // The only result that makes sense mathematically is to return zero, but
     * // we find it more convenient to return an arbitrary orthogonal vector.
     * return ortho(a);
     * }
     *
     * /**
     * Return a unit-length vector that is orthogonal to "a". Satisfies Ortho(-a)
     * = -Ortho(a) for all a.
     */
     public static function ortho(S2Point $a): S2Point {
         // The curr`ent implementation in S2Point has the property we need,
         // i.e. Ort`ho(-a) = -Ortho(a) for all a.
        return $a->ortho();
     }

     /**
     * Return the area of triangle ABC. The method used is about twice as
     * expensive as Girard's formula, but it is numerically stable for both large
     * and very small triangles. The points do not need to be normalized. The area
     * is always positive.
     *
     *  The triangle area is undefined if it contains two antipodal points, and
     * becomes numerically unstable as the length of any edge approaches 180
     * degrees.
     */
     static function area(S2Point $a, S2Point $b, S2Point $c): float {
         // This method is based on l'Huilier's theorem,
         //
         // tan(E/4) = sqrt(tan(s/2) tan((s-a)/2) tan((s-b)/2) tan((s-c)/2))
         //
         // where E is the spherical excess of the triangle (i.e. its area),
         // a, b, c, are the side lengths, and
         // s is the semiperimeter (a + b + c) / 2 .
         //
         // The only significant source of error using l'Huilier's method is the
         // cancellation error of the terms (s-a), (s-b), (s-c). This leads to a
         // *relative* error of about 1e-16 * s / min(s-a, s-b, s-c). This compares
         // to a relative error of about 1e-15 / E using Girard's formula, where E is
         // the true area of the triangle. Girard's formula can be even worse than
         // this for very small triangles, e.g. a triangle with a true area of 1e-30
         // might evaluate to 1e-5.
         //
         // So, we prefer l'Huilier's formula unless dmin < s * (0.1 * E), where
         // dmin = min(s-a, s-b, s-c). This basically includes all triangles
         // except for extremely long and skinny ones.
         //
         // Since we don't know E, we would like a conservative upper bound on
         // the triangle area in terms of s and dmin. It's possible to show that
         // E <= k1 * s * sqrt(s * dmin), where k1 = 2*sqrt(3)/Pi (about 1).
         // Using this, it's easy to show that we should always use l'Huilier's
         // method if dmin >= k2 * s^5, where k2 is about 1e-2. Furthermore,
         // if dmin < k2 * s^5, the triangle area is at most k3 * s^4, where
         // k3 is about 0.1. Since the best case error using Girard's formula
         // is about 1e-15, this means that we shouldn't even consider it unless
         // s >= 3e-4 or so.

         // We use volatile doubles to force the compiler to truncate all of these
         // quantities to 64 bits. Otherwise it may compute a value of dmin > 0
         // simply because it chose to spill one of the intermediate values to
         // memory but not one of the others.
         $sa = $b->angle($c);
         $sb = $c->angle($a);
         $sc = $a->angle($b);
         $s = 0.5 * ($sa + $sb + $sc);
         if ($s >= 3e-4) {
             // Consider whether Girard's formula might be more accurate.
             $s2 = $s * $s;
             $dmin = $s - max($sa, max($sb, $sc));
             if ($dmin < 1e-2 * $s * $s2 * $s2) {
                 // This triangle is skinny enough to consider Girard's formula.
                 $area = self::girardArea($a, $b, $c);
                 if ($dmin < $s * (0.1 * $area)) {
                     return $area;
                 }
             }
         }
         // Use l'Huilier's formula.
         return 4
            * atan(
                sqrt(
                    max(0.0,
                tan(0.5 * $s) * tan(0.5 * ($s - $sa)) * tan(0.5 * ($s - $sb))
                         * tan(0.5 * ($s - $sc)))));
     }

     /**
     * Return the area of the triangle computed using Girard's formula. This is
     * slightly faster than the Area() method above is not accurate for very small
     * triangles.
     */
     public static function girardArea(S2Point $a, S2Point $b, S2Point $c): float {
         // This is equivalent to the usual Girard's formula but is slightly
         // more accurate, faster to compute, and handles a == b == c without
         // a special case.

         $ab = S2Point::crossProd($a, $b);
         $bc = S2Point::crossProd($b, $c);
         $ac = S2Point::crossProd($a, $c);
         return max(0.0, $ab->angle($ac) - $ab->angle($bc) + $bc->angle($ac));
     }

     /**
     * Like Area(), but returns a positive value for counterclockwise triangles
     * and a negative value otherwise.
     */
     public static function signedArea(S2Point $a, S2Point $b, S2Point $c): float {
        return self::area($a, $b, $c) * self::robustCCW($a, $b, $c);
     }

     /* // About centroids:
     * // ----------------
     * //
     * // There are several notions of the "centroid" of a triangle. First, there
     * // // is the planar centroid, which is simply the centroid of the ordinary
     * // (non-spherical) triangle defined by the three vertices. Second, there is
     * // the surface centroid, which is defined as the intersection of the three
     * // medians of the spherical triangle. It is possible to show that this
     * // point is simply the planar centroid projected to the surface of the
     * // sphere. Finally, there is the true centroid (mass centroid), which is
     * // defined as the area integral over the spherical triangle of (x,y,z)
     * // divided by the triangle area. This is the point that the triangle would
     * // rotate around if it was spinning in empty space.
     * //
     * // The best centroid for most purposes is the true centroid. Unlike the
     * // planar and surface centroids, the true centroid behaves linearly as
     * // regions are added or subtracted. That is, if you split a triangle into
     * // pieces and compute the average of their centroids (weighted by triangle
     * // area), the result equals the centroid of the original triangle. This is
     * // not true of the other centroids.
     * //
     * // Also note that the surface centroid may be nowhere near the intuitive
     * // "center" of a spherical triangle. For example, consider the triangle
     * // with vertices A=(1,eps,0), B=(0,0,1), C=(-1,eps,0) (a quarter-sphere).
     * // The surface centroid of this triangle is at S=(0, 2*eps, 1), which is
     * // within a distance of 2*eps of the vertex B. Note that the median from A
     * // (the segment connecting A to the midpoint of BC) passes through S, since
     * // this is the shortest path connecting the two endpoints. On the other
     * // hand, the true centroid is at M=(0, 0.5, 0.5), which when projected onto
     * // the surface is a much more reasonable interpretation of the "center" of
     * // this triangle.
     *
     * /**
     * Return the centroid of the planar triangle ABC. This can be normalized to
     * unit length to obtain the "surface centroid" of the corresponding spherical
     * triangle, i.e. the intersection of the three medians. However, note that
     * for large spherical triangles the surface centroid may be nowhere near the
     * intuitive "center" (see example above).
     *#/
     * public static S2Point planarCentroid(S2Point a, S2Point b, S2Point c) {
     * return new S2Point((a.x + b.x + c.x) / 3.0, (a.y + b.y + c.y) / 3.0, (a.z + b.z + c.z) / 3.0);
     * }
     *
     * /**
     * Returns the true centroid of the spherical triangle ABC multiplied by the
     * signed area of spherical triangle ABC. The reasons for multiplying by the
     * signed area are (1) this is the quantity that needs to be summed to compute
     * the centroid of a union or difference of triangles, and (2) it's actually
     * easier to calculate this way.
     *#/
     * public static S2Point trueCentroid(S2Point a, S2Point b, S2Point c) {
     * // I couldn't find any references for computing the true centroid of a
     * // spherical triangle... I have a truly marvellous demonstration of this
     * // formula which this margin is too narrow to contain :)
     *
     * // assert (isUnitLength(a) && isUnitLength(b) && isUnitLength(c));
     * double sina = S2Point.crossProd(b, c).norm();
     * double sinb = S2Point.crossProd(c, a).norm();
     * double sinc = S2Point.crossProd(a, b).norm();
     * double ra = (sina == 0) ? 1 : (Math.asin(sina) / sina);
     * double rb = (sinb == 0) ? 1 : (Math.asin(sinb) / sinb);
     * double rc = (sinc == 0) ? 1 : (Math.asin(sinc) / sinc);
     *
     * // Now compute a point M such that M.X = rX * det(ABC) / 2 for X in A,B,C.
     * S2Point x = new S2Point(a.x, b.x, c.x);
     * S2Point y = new S2Point(a.y, b.y, c.y);
     * S2Point z = new S2Point(a.z, b.z, c.z);
     * S2Point r = new S2Point(ra, rb, rc);
     * return new S2Point(0.5 * S2Point.crossProd(y, z).dotProd(r),
     * 0.5 * S2Point.crossProd(z, x).dotProd(r), 0.5 * S2Point.crossProd(x, y).dotProd(r));
     * }
     *
     * /**
     * Return true if the points A, B, C are strictly counterclockwise. Return
     * false if the points are clockwise or colinear (i.e. if they are all
     * contained on some great circle).
     *
     *  Due to numerical errors, situations may arise that are mathematically
     * impossible, e.g. ABC may be considered strictly CCW while BCA is not.
     * However, the implementation guarantees the following:
     *
     *  If SimpleCCW(a,b,c), then !SimpleCCW(c,b,a) for all a,b,c.
     *
     * In other words, ABC and CBA are guaranteed not to be both CCW
     *#/
     * public static boolean simpleCCW(S2Point a, S2Point b, S2Point c) {
     * // We compute the signed volume of the parallelepiped ABC. The usual
     * // formula for this is (AxB).C, but we compute it here using (CxA).B
     * // in order to ensure that ABC and CBA are not both CCW. This follows
     * // from the following identities (which are true numerically, not just
     * // mathematically):
     * //
     * // (1) x.CrossProd(y) == -(y.CrossProd(x))
     * // (2) (-x).DotProd(y) == -(x.DotProd(y))
     *
     * return S2Point.crossProd(c, a).dotProd(b) > 0;
     * }
     *
     * /**
     * WARNING! This requires arbitrary precision arithmetic to be truly robust.
     * This means that for nearly colinear AB and AC, this function may return the
     * wrong answer.
     *
     * <p>
     * Like SimpleCCW(), but returns +1 if the points are counterclockwise and -1
     * if the points are clockwise. It satisfies the following conditions:
     *
     *  (1) RobustCCW(a,b,c) == 0 if and only if a == b, b == c, or c == a (2)
     * RobustCCW(b,c,a) == RobustCCW(a,b,c) for all a,b,c (3) RobustCCW(c,b,a)
     * ==-RobustCCW(a,b,c) for all a,b,c
     *
     *  In other words:
     *
     *  (1) The result is zero if and only if two points are the same. (2)
     * Rotating the order of the arguments does not affect the result. (3)
     * Exchanging any two arguments inverts the result.
     *
     *  This function is essentially like taking the sign of the determinant of
     * a,b,c, except that it has additional logic to make sure that the above
     * properties hold even when the three points are coplanar, and to deal with
     * the limitations of floating-point arithmetic.
     *
     *  Note: a, b and c are expected to be of unit length. Otherwise, the results
     * are undefined.
     */
     public static function robustCCW(S2Point $a, S2Point $b, S2Point $c): int {
        return self::robustCCWWithCross($a, $b, $c, S2Point::crossProd($a, $$b));
     }

     /**
     * A more efficient version of RobustCCW that allows the precomputed
     * cross-product of A and B to be specified.
     *
     *  Note: a, b and c are expected to be of unit length. Otherwise, the results
     * are undefined
     */
     public static function robustCCWWithCross(S2Point $a, S2Point $b, S2Point $c, S2Point $aCrossB): int {
        assert (self::isUnitLength($a) && self::isUnitLength($b) && self::isUnitLength($c));

         // There are 14 multiplications and additions to compute the determinant
         // below. Since all three points are normalized, it is possible to show
         // that the average rounding error per operation does not exceed 2**-54,
         // the maximum rounding error for an operation whose result magnitude is in
         // the range [0.5,1). Therefore, if the absolute value of the determinant
         // is greater than 2*14*(2**-54), the determinant will have the same sign
         // even if the arguments are rotated (which produces a mathematically
         // equivalent result but with potentially different rounding errors).
         /** @var float */ $kMinAbsValue = 1.6e-15; // 2 * 14 * 2**-54

         /** @var float */ $det = $aCrossB->dotProd($c);

         // Double-check borderline cases in debug mode.
         // assert ((Math.abs(det) < kMinAbsValue) || (Math.abs(det) > 1000 * kMinAbsValue)
         //    || (det * expensiveCCW(a, b, c) > 0));

         if ($det > $kMinAbsValue) {
            return 1;
         }

         if ($det < -$kMinAbsValue) {
            return -1;
         }

         return self::expensiveCCW($a, $b, $c);
     }

     /**
     * A relatively expensive calculation invoked by RobustCCW() if the sign of
     * the determinant is uncertain.
     */
     private static function expensiveCCW(S2Point $a, S2Point $b, S2Point $c): int {
         // Return zero if and only if two points are the same. This ensures (1).
         if ($a->equals($b) || $b->equals($c) || $c->equals($a)) {
            return 0;
         }

         // Now compute the determinant in a stable way. Since all three points are
         // unit length and we know that the determinant is very close to zero, this
         // means that points are very nearly colinear. Furthermore, the most common
         // situation is where two points are nearly identical or nearly antipodal.
         // To get the best accuracy in this situation, it is important to
         // immediately reduce the magnitude of the arguments by computing either
         // A+B or A-B for each pair of points. Note that even if A and B differ
         // only in their low bits, A-B can be computed very accurately. On the
         // other hand we can't accurately represent an arbitrary linear combination
         // of two vectors as would be required for Gaussian elimination. The code
         // below chooses the vertex opposite the longest edge as the "origin" for
         // the calculation, and computes the different vectors to the other two
         // vertices. This minimizes the sum of the lengths of these vectors.
         //
         // This implementation is very stable numerically, but it still does not
         // return consistent results in all cases. For example, if three points are
         // spaced far apart from each other along a great circle, the sign of the
         // result will basically be random (although it will still satisfy the
         // conditions documented in the header file). The only way to return
         // consistent results in all cases is to compute the result using
         // arbitrary-precision arithmetic. I considered using the Gnu MP library,
         // but this would be very expensive (up to 2000 bits of precision may be
         // needed to store the intermediate results) and seems like overkill for
         // this problem. The MP library is apparently also quite particular about
         // compilers and compilation options and would be a pain to maintain.

         // We want to handle the case of nearby points and nearly antipodal points
         // accurately, so determine whether A+B or A-B is smaller in each case.
         /** @var float */ $sab = ($a->dotProd($b) > 0) ? -1 : 1;
         /** @var float */ $sbc = ($b->dotProd($c) > 0) ? -1 : 1;
         /** @var float */ $sca = ($c->dotProd($a) > 0) ? -1 : 1;
         /** @var S2Point */ $vab = S2Point::add($a, S2Point::mul($b, $sab));
         /** @var S2Point */ $vbc = S2Point::add($b, S2Point::mul($c, $sbc));
         /** @var S2Point */ $vca = S2Point::add($c, S2Point::mul($a, $sca));
         /** @var float */ $dab = $vab->norm2();
         /** @var float */ $dbc = $vbc->norm2();
         /** @var float */ $dca = $vca->norm2();

         // Sort the difference vectors to find the longest edge, and use the
         // opposite vertex as the origin. If two difference vectors are the same
         // length, we break ties deterministically to ensure that the symmetry
         // properties guaranteed in the header file will be true.
         $sign = 0.0;

         if ($dca < $dbc || ($dca == $dbc && $a->lessThan($b))) {
             if ($dab < $dbc || ($dab == $dbc && $a->lessThan($c))) {
                 // The "sab" factor converts A +/- B into B +/- A.
                 $sign = S2Point::crossProd($vab, $vca)->dotProd($a) * $sab; // BC is longest
                 // edge
             } else {
                $sign = S2Point::crossProd($vca, $vbc)->dotProd($c) * $sca; // AB is longest
             // edge
             }
         } else {
             if ($dab < $dca || ($dab == $dca && $b->lessThan($c))) {
                 $sign = S2Point::crossProd($vbc, $vab)->dotProd($b) * $sbc; // CA is longest
                 // edge
             } else {
                 $sign = S2Point::crossProd($vca, $vbc)->dotProd($c) * $sca; // AB is longest
                 // edge
             }
         }
         if ($sign > 0) {
            return 1;
         }
         if ($sign < 0) {
            return -1;

        }

         // The points A, B, and C are numerically indistinguishable from coplanar.
         // This may be due to roundoff error, or the points may in fact be exactly
         // coplanar. We handle this situation by perturbing all of the points by a
         // vector (eps, eps**2, eps**3) where "eps" is an infinitesmally small
         // positive number (e.g. 1 divided by a googolplex). The perturbation is
         // done symbolically, i.e. we compute what would happen if the points were
         // perturbed by this amount. It turns out that this is equivalent to
         // checking whether the points are ordered CCW around the origin first in
         // the Y-Z plane, then in the Z-X plane, and then in the X-Y plane.

         $ccw = self::planarOrderedCCW(
             new R2Vector($a->y, $a->z),
             new R2Vector($b->y, $b->z),
             new R2Vector($c->y, $c->z)
         );

         if ($ccw === 0) {
            $ccw = self::planarOrderedCCW(
                new R2Vector($a->z, $a->x),
                new R2Vector($b->z, $b->x),
                new R2Vector($c->z, $c->x)
            );
            if ($ccw === 0) {
                $ccw = self::planarOrderedCCW(
                    new R2Vector($a->x, $a->y),
                    new R2Vector($b->x, $b->y),
                    new R2Vector($c->x, $c->y)
                );
                assert ($ccw !== 0);
            }
         }
         return $ccw;
     }


     public static function planarCCW(R2Vector $a, R2Vector $b): int {
         // Return +1 if the edge AB is CCW around the origin, etc.
         $sab = ($a->dotProd($b) > 0) ? -1 : 1;
         $vab = R2Vector::add($a, R2Vector::mul($b, $sab));
         $da = $a->norm2();
         $db = $b->norm2();

         if ($da < $db || ($da === $db && $a->lessThan($b))) {
            $sign = $a->crossProd($vab) * $sab;
         } else {
            $sign = $vab->crossProd($b);
         }
         if ($sign > 0) {
            return 1;
         }
         if ($sign < 0) {
            return -1;
         }
         return 0;
     }

     public static function planarOrderedCCW(R2Vector $a, R2Vector $b, R2Vector $c): int {
         $sum = 0;
         $sum += self::planarCCW($a, $b);
         $sum += self::planarCCW($b, $c);
         $sum += self::planarCCW($c, $a);
         if ($sum > 0) {
             return 1;
         }
         if ($sum < 0) {
             return -1;
         }
         return 0;
     }

     /**
     * Return true if the edges OA, OB, and OC are encountered in that order while
     * sweeping CCW around the point O. You can think of this as testing whether
     * A <= B <= C with respect to a continuous CCW ordering around O.
     *
     * Properties:
     * <ol>
     *   <li>If orderedCCW(a,b,c,o) && orderedCCW(b,a,c,o), then a == b</li>
     *   <li>If orderedCCW(a,b,c,o) && orderedCCW(a,c,b,o), then b == c</li>
     *   <li>If orderedCCW(a,b,c,o) && orderedCCW(c,b,a,o), then a == b == c</li>
     *   <li>If a == b or b == c, then orderedCCW(a,b,c,o) is true</li>
     *   <li>Otherwise if a == c, then orderedCCW(a,b,c,o) is false</li>
     * </ol>
     */
     public static function orderedCCW(S2Point $a, S2Point $b, S2Point $c, S2Point $o): bool {
         // The last inequality below is ">" rather than ">=" so that we return true
         // if A == B or B == C, and otherwise false if A == C. Recall that
         // RobustCCW(x,y,z) == -RobustCCW(z,y,x) for all x,y,z.

         $sum = 0;
         if (self::robustCCW($b, $o, $a) >= 0) {
            ++$sum;
         }
         if (self::robustCCW($c, $o, $b) >= 0) {
            ++$sum;
         }
         if (self::robustCCW($a, $o, $c) > 0) {
            ++$sum;
         }
         return $sum >= 2;
     }

     /**
     * Return the angle at the vertex B in the triangle ABC. The return value is
     * always in the range [0, Pi]. The points do not need to be normalized.
     * Ensures that Angle(a,b,c) == Angle(c,b,a) for all a,b,c.
     *
     *  The angle is undefined if A or C is diametrically opposite from B, and
     * becomes numerically unstable as the length of edge AB or BC approaches 180
     * degrees.
     *#/
     * public static double angle(S2Point a, S2Point b, S2Point c) {
     * return S2Point.crossProd(a, b).angle(S2Point.crossProd(c, b));
     * }
     *
     * /**
     * Return the exterior angle at the vertex B in the triangle ABC. The return
     * value is positive if ABC is counterclockwise and negative otherwise. If you
     * imagine an ant walking from A to B to C, this is the angle that the ant
     * turns at vertex B (positive = left, negative = right). Ensures that
     * TurnAngle(a,b,c) == -TurnAngle(c,b,a) for all a,b,c.
     *
     * @param a
     * @param b
     * @param c
     * @return the exterior angle at the vertex B in the triangle ABC
     *#/
     * public static double turnAngle(S2Point a, S2Point b, S2Point c) {
     * // This is a bit less efficient because we compute all 3 cross products, but
     * // it ensures that turnAngle(a,b,c) == -turnAngle(c,b,a) for all a,b,c.
     * double outAngle = S2Point.crossProd(b, a).angle(S2Point.crossProd(c, b));
     * return (robustCCW(a, b, c) > 0) ? outAngle : -outAngle;
     * }
     *
     * /**
     * Return true if two points are within the given distance of each other
     * (mainly useful for testing).
     *#/
     * public static boolean approxEquals(S2Point a, S2Point b, double maxError) {
     * return a.angle(b) <= maxError;
     * }
     *
     * public static boolean approxEquals(S2Point a, S2Point b) {
     * return approxEquals(a, b, 1e-15);
     * }
     *
     * public static boolean approxEquals(double a, double b, double maxError) {
     * return Math.abs(a - b) <= maxError;
     * }
     *
     * public static boolean approxEquals(double a, double b) {
     * return approxEquals(a, b, 1e-15);
     * }
     */
}

/**
 * Defines an area or a length cell metric.
 */

class Metric {
    private $deriv;
    private $dim;

    /**
     * Defines a cell metric of the given dimension (1 == length, 2 == area).
     */
    public function __construct($dim, $deriv) {
        $this->deriv = $deriv;
        $this->dim = $dim;
    }

    /**
     * The "deriv" value of a metric is a derivative, and must be multiplied by
     * a length or area in (s,t)-space to get a useful value.
     */
    public function deriv() {
        return $this->deriv;
    }

    /** Return the value of a metric for cells at the given level. */
     public function getValue($level) {
        return JavaMathHelper::java_Math_ScalB($this->deriv, -$this->dim * $level);
     }

    /**
     * Return the level at which the metric has approximately the given value.
     * For example, S2::kAvgEdge.GetClosestLevel(0.1) returns the level at which
     * the average cell edge length is approximately 0.1. The return value is
     * always a valid level.
     */
    public function getClosestLevel($value) {
        return $this->getMinLevel(M_SQRT2 * $value);
    }

    /**
     * Return the minimum level such that the metric is at most the given value,
     * or S2CellId::kMaxLevel if there is no such level. For example,
     * S2::kMaxDiag.GetMinLevel(0.1) returns the minimum level such that all
     * cell diagonal lengths are 0.1 or smaller. The return value is always a
     * valid level.
     */
    public function getMinLevel($value) {
        if ($value <= 0) {
            return S2CellId::MAX_LEVEL;
        }

// This code is equivalent to computing a floating-point "level"
// value and rounding up.
        $exponent = S2::exp($value / ((1 << $this->dim) * $this->deriv));
        $level = max(
            0,
            min(S2CellId::MAX_LEVEL, -(($exponent - 1) >> ($this->dim - 1)))
        );
// assert (level == S2CellId.MAX_LEVEL || getValue(level) <= value);
// assert (level == 0 || getValue(level - 1) > value);
        return $level;
    }

    /**
     * Return the maximum level such that the metric is at least the given
     * value, or zero if there is no such level. For example,
     * S2.kMinWidth.GetMaxLevel(0.1) returns the maximum level such that all
     * cells have a minimum width of 0.1 or larger. The return value is always a
     * valid level.
     */
    public function getMaxLevel($value) {
        if ($value <= 0) return S2CellId::MAX_LEVEL;

// This code is equivalent to computing a floating-point "level"
// value and rounding down.
        $exponent = S2::exp((1 << $this->dim) * $this->deriv / $value);
        $level = max(
            0,
            min(S2CellId::MAX_LEVEL, (($exponent - 1) >> ($this->dim - 1)))
        );
// assert (level == 0 || getValue(level) >= value);
// assert (level == S2CellId.MAX_LEVEL || getValue(level + 1) < value);
        return $level;
    }
}

class JavaMathHelper {

     /**
      * Maximum exponent a finite {@code double} variable may have.
      * It is equal to the value returned by
      * {@code Math.getExponent(Double.MAX_VALUE)}.
      *
      * @since 1.6
      * @see https://github.com/openjdk/jdk13u/blob/master/src/java.base/share/classes/java/lang/Double.java
      */
     private const DOUBLE_MAX_EXPONENT = 1023;


    /**
     * Minimum exponent a normalized {@code double} variable may
     * have.  It is equal to the value returned by
     * {@code Math.getExponent(Double.MIN_NORMAL)}.
     *
     * @since 1.6
     * @see https://github.com/openjdk/jdk13u/blob/master/src/java.base/share/classes/java/lang/Double.java
     */
     private const DOUBLE_MIN_EXPONENT = -1022;


     /**
      * The number of logical bits in the significand of a
      * {@code double} number, including the implicit bit.
      * @see https://github.com/openjdk/jdk13u/blob/master/src/java.base/share/classes/jdk/internal/math/DoubleConsts.java
      */
     private const DOUBLE_CONSTS_SIGNIFICAND_WIDTH = 53;

    /**
     * Bias used in representing a {@code double} exponent.
     * @see https://github.com/openjdk/jdk13u/blob/master/src/java.base/share/classes/jdk/internal/math/DoubleConsts.java
     */
    private const DOUBLE_CONSTS_EXP_BIAS        = 1023;

    /**
     * Bit mask to isolate the exponent field of a
     * {@code double}.
     * @see https://github.com/openjdk/jdk13u/blob/master/src/java.base/share/classes/jdk/internal/math/DoubleConsts.java
     */
    private const DOUBLE_CONSTS_EXP_BIT_MASK    = 0x7FF0000000000000;

    // Constants used in scalb

    private static function twoToTheDoubleScaleUp(): float {
        return self::powerOfTwoD(512);
    }

    private static function twoToTheDoubleScaleDown(): float {
        return self::powerOfTwoD(-512);
    }

    private static function longBitsToDouble(int $value): float {
        $bin    = pack('q', $value);
        $double = unpack('d', $bin);
        return array_shift($double);
    }

    /**
     * Returns a floating-point power of two in the normal range.
     * @see https://github.com/openjdk/jdk13u/blob/master/src/java.base/share/classes/java/lang/Math.java
     */

    private static function powerOfTwoD(int $n): float {
        return
            self::longBitsToDouble(((($n + self::DOUBLE_CONSTS_EXP_BIAS) <<
                (self::DOUBLE_CONSTS_SIGNIFICAND_WIDTH-1))
               & self::DOUBLE_CONSTS_EXP_BIT_MASK))
        ;
    }

    // https://stackoverflow.com/a/43359819/14651
    public static function uRShift($a, $b)
    {
        if ($b >= 32 || $b < -32) {
            $m = (int)($b/32);
            $b = $b-($m*32);
        }

        if ($b < 0) {
            $b = 32 + $b;
        }

        if ($b == 0) {
            return (($a>>1)&0x7fffffff)*2+(($a>>$b)&1);
        }

        if ($a < 0)
        {
            $a = ($a >> 1);
            $a &= 0x7fffffff;
            $a |= 0x40000000;
            $a = ($a >> ($b - 1));
        } else {
            $a = ($a >> $b);
        }
        return $a;
    }

    /**
     * Returns {@code d} &times;
     * 2<sup>{@code scaleFactor}</sup> rounded as if performed
     * by a single correctly rounded floating-point multiply to a
     * member of the double value set.  See the Java
     * Language Specification for a discussion of floating-point
     * value sets.  If the exponent of the result is between {@link
     * Double#MIN_EXPONENT} and {@link Double#MAX_EXPONENT}, the
     * answer is calculated exactly.  If the exponent of the result
     * would be larger than {@code Double.MAX_EXPONENT}, an
     * infinity is returned.  Note that if the result is subnormal,
     * precision may be lost; that is, when {@code scalb(x, n)}
     * is subnormal, {@code scalb(scalb(x, n), -n)} may not equal
     * <i>x</i>.  When the result is non-NaN, the result has the same
     * sign as {@code d}.
     *
     * <p>Special cases:
     * <ul>
     * <li> If the first argument is NaN, NaN is returned.
     * <li> If the first argument is infinite, then an infinity of the
     * same sign is returned.
     * <li> If the first argument is zero, then a zero of the same
     * sign is returned.
     * </ul>
     *
     * @param d number to be scaled by a power of two.
     * @param scaleFactor power of 2 used to scale {@code d}
     * @return {@code d} &times; 2<sup>{@code scaleFactor}</sup>
     * @since 1.6
     * @see https://github.com/openjdk/jdk13u/blob/master/src/java.base/share/classes/java/lang/Math.java
     */
     public static function java_Math_ScalB($d, int $scaleFactor): float {
         /*
          * This method does not need to be declared strictfp to
          * compute the same correct result on all platforms.  When
          * scaling up, it does not matter what order the
          * multiply-store operations are done; the result will be
          * finite or overflow regardless of the operation ordering.
          * However, to get the correct result when scaling down, a
          * particular ordering must be used.
          *
          * When scaling down, the multiply-store operations are
          * sequenced so that it is not possible for two consecutive
          * multiply-stores to return subnormal results.  If one
          * multiply-store result is subnormal, the next multiply will
          * round it away to zero.  This is done by first multiplying
          * by 2 ^ (scaleFactor % n) and then multiplying several
          * times by 2^n as needed where n is the exponent of number
          * that is a covenient power of two.  In this way, at most one
          * real rounding error occurs.  If the double value set is
          * being used exclusively, the rounding will occur on a
          * multiply.  If the double-extended-exponent value set is
          * being used, the products will (perhaps) be exact but the
          * stores to d are guaranteed to round to the double value
          * set.
          *
          * It is _not_ a valid implementation to first multiply d by
          * 2^MIN_EXPONENT and then by 2 ^ (scaleFactor %
          * MIN_EXPONENT) since even in a strictfp program double
          * rounding on underflow could occur; e.g. if the scaleFactor
          * argument was (MIN_EXPONENT - n) and the exponent of d was a
          * little less than -(MIN_EXPONENT - n), meaning the final
          * result would be subnormal.
          *
          * Since exact reproducibility of this method can be achieved
          * without any undue performance burden, there is no
          * compelling reason to allow double rounding on underflow in
          * scalb.
          */

        // magnitude of a power of two so large that scaling a finite
        // nonzero value by it would be guaranteed to over or
        // underflow; due to rounding, scaling down takes an
        // additional power of two which is reflected here


        $MAX_SCALE = self::DOUBLE_MAX_EXPONENT + -self::DOUBLE_MIN_EXPONENT +
                     self::DOUBLE_CONSTS_SIGNIFICAND_WIDTH + 1;

        // Make sure scaling factor is in a reasonable range

        if($scaleFactor < 0) {
            $scaleFactor = max($scaleFactor, -$MAX_SCALE);
            $scale_increment = -512;
            $exp_delta = self::twoToTheDoubleScaleDown();
        } else {
            $scaleFactor = min($scaleFactor, $MAX_SCALE);
            $scale_increment = 512;
            $exp_delta = self::twoToTheDoubleScaleUp();
        }


        // Calculate (scaleFactor % +/-512), 512 = 2^9, using
        // technique from "Hacker's Delight" section 10-2.
        $t = self::uRShift($scaleFactor >> 9-1, 32 - 9);
        $exp_adjust = (($scaleFactor + $t) & (512 -1)) - $t;

        $d *= self::powerOfTwoD($exp_adjust);
        $scaleFactor -= $exp_adjust;

        while($scaleFactor !== 0) {
            $d *= $exp_delta;
            $scaleFactor -= $scale_increment;
        }

        return $d;
     }
}