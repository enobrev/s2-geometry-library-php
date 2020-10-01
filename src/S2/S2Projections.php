<?php

namespace S2;

use function foo\func;

class S2Projections
{
    /** @var S2Projections  */
    public static $PROJ;

    public static function init() {
        self::$PROJ = S2Projections::S2_QUADRATIC_PROJECTION_OLD();
    }

    public static function S2_LINEAR_PROJECTION() : self {
        return new self(
              4 / (3 * sqrt(3)), // minArea 0.770
              4, // maxArea 4.000
              1.0, // minAngleSpan 1.000
              2, // maxAngleSpan 2.000
              sqrt(2. / 3), // minWidth 0.816
              1.411459345844456965, // avgWidth 1.411
              2 * sqrt(2) / 3, // minEdge 0.943
              1.440034192955603643, // avgEdge 1.440
              2 * sqrt(2) / 3, // minDiag 0.943
              2 * sqrt(2), // maxDiag 2.828
              2.031817866418812674, // avgDiag 2.032
              sqrt(2)  // maxEdgeAspect 1.414
        );
    }
    
    public static function S2_TAN_PROJECTION() : self {
        return new self(
              (S2::M_PI * S2::M_PI) / (4 * sqrt(2)), // minArea 1.745
              S2::M_PI * S2::M_PI / 4, // maxArea 2.467
              S2::M_PI / 2, // minAngleSpan 1.571
              S2::M_PI / 2, // maxAngleSpan 1.571
              S2::M_PI / (2 * sqrt(2)), // minWidth 1.111
              1.437318638925160885, // avgWidth 1.437
              S2::M_PI / (2 * sqrt(2)), // minEdge 1.111
              1.461667032546739266, // avgEdge 1.462
              S2::M_PI * sqrt(2) / 3, // minDiag 1.481
              S2::M_PI * sqrt(2. / 3), // maxDiag 2.565
              2.063623197195635753, // avgDiag 2.064
              sqrt(2)  // maxEdgeAspect 1.414
        );
    }

    public static function S2_QUADRATIC_PROJECTION() : self {
        return new self(
            8 * sqrt(2) / 9, // minArea 1.257
            2.635799256963161491, // maxArea 2.636
            4. / 3, // minAngleSpan 1.333
            1.704897179199218452, // maxAngleSpan 1.705
            2 * sqrt(2) / 3, // minWidth 0.943
            1.434523672886099389, // avgWidth 1.435
            2 * sqrt(2) / 3, // minEdge 0.943
            1.459213746386106062, // avgEdge 1.459
            8 * sqrt(2) / 9, // minDiag 1.257
            2.438654594434021032, // maxDiag 2.439
            2.060422738998471683, // avgDiag 2.060
            1.442615274452682920 // maxEdgeAspect 1.443
        );
    }

    public static function S2_QUADRATIC_PROJECTION_OLD() : self {
        return new self(
            2 * S2::M_SQRT2 / 9, // minArea 0.314
            0.6589498142407903, // maxArea 0.659
            2. / 3, // minAngleSpan 0.667
            0.85244858959960922, // maxAngleSpan 0.852
            S2::M_SQRT2 / 3, // minWidth 0.471
            0.71726183644304969, // avgWidth 0.717
            S2::M_SQRT2 / 3, // minEdge 0.471
            0.72960687319305303, // avgEdge 0.730
            4 * S2::M_SQRT2 / 9, // minDiag 0.629
            1.2193272972170106, // maxDiag 1.219
            1.03021136949923584, // avgDiag 1.030
            1.44261527445268292 // maxEdgeAspect 1.443
        );
    }

    /** @var Metric Minimum area of a cell at level k. */
    public $minArea;

  /** @var Metric Maximum area of a cell at level k. */
  public $maxArea;

  /** @var Metric Average area of a cell at level k. */
  public $avgArea;

  /**
   * Minimum angular separation between opposite edges of a cell at level k. Each cell is bounded by
   * four planes passing through its four edges and the center of the sphere. The angle span metrics
   * relate to the angle between each pair of opposite bounding planes, or equivalently, between the
   * planes corresponding to two different s-values or two different t-values.
   * @var Metric
   */
  public $minAngleSpan;

  /** @var Metric Maximum angular separation between opposite edges of a cell at level k. */
  public $maxAngleSpan;

  /** @var Metric Average angular separation between opposite edges of a cell at level k. */
  public $avgAngleSpan;

  /**
   * Minimum perpendicular angular separation between opposite edges of a cell at level k.
   *
   * <p>The width of a geometric figure is defined as the distance between two parallel bounding
   * lines in a given direction. For cells, the minimum width is always attained between two
   * opposite edges, and the maximum width is attained between two opposite vertices. However, for
   * our purposes we redefine the width of a cell as the perpendicular distance between a pair of
   * opposite edges. A cell therefore has two widths, one in each direction. The minimum width
   * according to this definition agrees with the classic geometric one, but the maximum width is
   * different. (The maximum geometric width corresponds to {@link #maxDiag}.)
   *
   * <p>This is useful for bounding the minimum or maximum distance from a point on one edge of a
   * cell to the closest point on the opposite edge. For example, this is useful when "growing"
   * regions by a fixed distance.
   * @var Metric
   */
  public $minWidth;

  /** @var Metric Maximum perpendicular angular separation between opposite edges of a cell at level k. */
  public $maxWidth;

  /** @var Metric Average perpendicular angular separation between opposite edges of a cell at level k. */
  public $avgWidth;

  /**
   * Minimum angular length of any cell edge at level k. The edge length metrics can also be used to
   * bound the minimum, maximum, or average distance from the center of one cell to the center of
   * one of its edge neighbors. In particular, it can be used to bound the distance between adjacent
   * cell centers along the space-filling Hilbert curve for cells at any given level.
   * @var Metric
   */
  public $minEdge;

  /** @var Metric Maximum angular length of any cell edge at level k. */
  public $maxEdge;

  /** @var Metric Average angular length of any cell edge at level k. */
  public $avgEdge;

  /** @var Metric Minimum diagonal size of cells at level k. */
  public $minDiag;

  /**
   * Maximum diagonal size of cells at level k. The maximum diagonal also happens to be the maximum
   * diameter of any cell, and also the maximum geometric width. So for example, the distance from
   * an arbitrary point to the closest cell center at a given level is at most half the maximum
   * diagonal length.
   * @var Metric
   */
  public $maxDiag;

  /** @var Metric Average diagonal size of cells at level k. */
  public $avgDiag;

  /**
   * Maximum edge aspect ratio over all cells at any level, where the edge aspect ratio of a cell is
   * defined as the ratio of its longest edge length to its shortest edge length.
   * @var float
   */
  public $maxEdgeAspect;

  /**
   * This is the maximum diagonal aspect ratio over all cells at any level, where the diagonal
   * aspect ratio of a cell is defined as the ratio of its longest diagonal length to its shortest
   * diagonal length.
   * @var float
   */
  public $maxDiagAspect; // = sqrt(3); // 1.732

    private function __construct(
        float $minAreaDeriv,
        float $maxAreaDeriv,
        float $minAngleSpanDeriv,
        float $maxAngleSpanDeriv,
        float $minWidthDeriv,
        float $avgWidthDeriv,
        float $minEdgeDeriv,
        float $avgEdgeDeriv,
        float $minDiagDeriv,
        float $maxDiagDeriv,
        float $avgDiagDeriv,
        float $maxEdgeAspect
    ) {
        $this->minArea = new Metric(2, $minAreaDeriv);
        $this->maxArea = new Metric(2, $maxAreaDeriv);
        $this->avgArea = new Metric(2, 4 * S2::M_PI / 6); // ~2.094
        $this->minAngleSpan = new Metric(1, $minAngleSpanDeriv);
        $this->maxAngleSpan = new Metric(1, $maxAngleSpanDeriv);
        $this->avgAngleSpan = new Metric(1, S2::M_PI / 2); // ~1.571
        $this->minWidth = new Metric(1, $minWidthDeriv);
        $this->maxWidth = new Metric(1, $maxAngleSpanDeriv);
        $this->avgWidth = new Metric(1, $avgWidthDeriv);
        $this->minEdge = new Metric(1, $minEdgeDeriv);
        $this->maxEdge = new Metric(1, $maxAngleSpanDeriv);
        $this->avgEdge = new Metric(1, $avgEdgeDeriv);
        $this->minDiag = new Metric(1, $minDiagDeriv);
        $this->maxDiag = new Metric(1, $maxDiagDeriv);
        $this->avgDiag = new Metric(1, $avgDiagDeriv);
        $this->maxEdgeAspect = $maxEdgeAspect;
    }

    const S2_LINEAR_PROJECTION = 1;
    const S2_TAN_PROJECTION = 2;
    const S2_QUADRATIC_PROJECTION = 3;
    const S2_QUADRATIC_PROJECTION_OLD = 4;

    const S2_PROJECTION = self::S2_QUADRATIC_PROJECTION_OLD;

    public static function stToUV($s)
    {
        switch (self::S2_PROJECTION) {
            case self::S2_LINEAR_PROJECTION:
                return 2 * $s - 1;

            case self::S2_TAN_PROJECTION:
              // Unfortunately, tan(M_PI_4) is slightly less than 1.0. This isn't due to a flaw in the
              // implementation of tan(), it's because the derivative of tan(x) at x=pi/4 is 2, and it
              // happens that the two adjacent floating point numbers on either side of the infinite-
              // precision value of pi/4 have tangents that are slightly below and slightly above 1.0 when
              // rounded to the nearest double-precision result.
              $s = tan(S2::M_PI_2 * $s - S2::M_PI_4);
              return $s + (1.0 / (1 << 53)) * $s;

            case self::S2_QUADRATIC_PROJECTION:
              if ($s >= 0.5) {
                return (1 / 3.) * (4 * $s * $s - 1);
              } else {
                return (1 / 3.) * (1 - 4 * (1 - $s) * (1 - $s));
              }

            case self::S2_QUADRATIC_PROJECTION_OLD:
                if ($s >= 0) {
                    return (1 / 3.) * ((1 + $s) * (1 + $s) - 1);
                } else {
                    return (1 / 3.) * (1 - (1 - $s) * (1 - $s));
                }
            default:
                throw new IllegalStateException("Invalid value for S2_PROJECTION");
        }
    }

    public static function uvToST($u)
    {
        switch (self::S2_PROJECTION) {
            case self::S2_LINEAR_PROJECTION:
                return 0.5 * ($u + 1);

            case self::S2_TAN_PROJECTION:
                return (4 * S2::M_1_PI) * atan($u);

            case self::S2_QUADRATIC_PROJECTION:
              if ($u >= 0) {
                return 0.5 * sqrt(1 + 3 * $u);
              } else {
                return 1 - 0.5 * sqrt(1 - 3 * $u);
              }
            case self::S2_QUADRATIC_PROJECTION_OLD:
                if ($u >= 0) {
                    return sqrt(1 + 3 * $u) - 1;
                } else {
                    return 1 - sqrt(1 - 3 * $u);
                }
            default:
                throw new IllegalStateException("Invalid value for S2_PROJECTION");
        }
    }

    /**
     * Convert (face, u, v) coordinates to a direction vector (not necessarily
     * unit length).
     */
    public static function faceUvToXyz($face, $u, $v)
    {
        switch ($face) {
            case 0:
                return new S2Point(1, $u, $v);

            case 1:
                return new S2Point(-$u, 1, $v);

            case 2:
                return new S2Point(-$u, -$v, 1);

            case 3:
                return new S2Point(-1, -$v, -$u);

            case 4:
                return new S2Point($v, -1, -$u);

            default:
                return new S2Point($v, $u, -1);
        }
    }

    public static function validFaceXyzToUv($face, S2Point $p)
    {
        // assert (p.dotProd(faceUvToXyz(face, 0, 0)) > 0);
        switch ($face) {
            case 0:
                $pu = $p->y / $p->x;
                $pv = $p->z / $p->x;
                break;

            case 1:
                $pu = -$p->x / $p->y;
                $pv = $p->z / $p->y;
                break;

            case 2:
                $pu = -$p->x / $p->z;
                $pv = -$p->y / $p->z;
                break;

            case 3:
                $pu = $p->z / $p->x;
                $pv = $p->y / $p->x;
                break;

            case 4:
                $pu = $p->z / $p->y;
                $pv = -$p->x / $p->y;
                break;

            default:
                $pu = -$p->y / $p->z;
                $pv = -$p->x / $p->z;
                break;
        }
        return new R2Vector($pu, $pv);
    }

    public static function xyzToFace(S2Point $p)
    {
        $face = $p->largestAbsComponent();
        if ($p->get($face) < 0) {
            $face += 3;
        }
        return $face;
    }

  public static function faceXyzToUv(int $face, S2Point $p): ?R2Vector {
    if ($face < 3) {
      if ($p->get($face) <= 0) {
        return null;
      }
    } else {
      if ($p->get($face - 3) >= 0) {
        return null;
      }
    }
    return self::validFaceXyzToUv($face, $p);
  }

  public static function getUNorm(int $face, float $u): S2Point {
    switch ($face) {
      case 0:
        return new S2Point($u, -1, 0);
      case 1:
        return new S2Point(1, $u, 0);
      case 2:
        return new S2Point(1, 0, $u);
      case 3:
        return new S2Point(-$u, 0, 1);
      case 4:
        return new S2Point(0, -$u, 1);
      default:
        return new S2Point(0, -1, -$u);
    }
  }

  public static function getVNorm(int $face, float $v): S2Point {
    switch ($face) {
      case 0:
        return new S2Point(-$v, 0, 1);
      case 1:
        return new S2Point(0, -$v, 1);
      case 2:
        return new S2Point(0, -1, -$v);
      case 3:
        return new S2Point($v, -1, 0);
      case 4:
        return new S2Point(1, $v, 0);
      default:
        return new S2Point(1, 0, $v);
    }
  }
/*
  public static S2Point getNorm(int face) {
    return faceUvToXyz(face, 0, 0);
  }
*/
    public static function getUAxis($face)
    {
        switch ($face) {
            case 0:
                return new S2Point(0, 1, 0);

            case 1:
                return new S2Point(-1, 0, 0);

            case 2:
                return new S2Point(-1, 0, 0);

            case 3:
                return new S2Point(0, 0, -1);

            case 4:
                return new S2Point(0, 0, -1);

            default:
                return new S2Point(0, 1, 0);
        }
    }

    public static function getVAxis($face)
    {
        switch ($face) {
            case 0:
                return new S2Point(0, 0, 1);

            case 1:
                return new S2Point(0, 0, 1);

            case 2:
                return new S2Point(0, -1, 0);

            case 3:
                return new S2Point(0, -1, 0);

            case 4:
                return new S2Point(1, 0, 0);

            default:
                return new S2Point(1, 0, 0);
        }
    }
}

S2Projections::init();