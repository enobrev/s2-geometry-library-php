<?php
namespace S2\Tests;

use PHPUnit\Framework\TestCase;

use S2\S2;
use S2\S2Cap;
use S2\S2Cell;
use S2\S2CellId;
use S2\S2CellUnion;
use S2\S2LatLng;
use S2\S2Loop;
use S2\S2Point;
use S2\S2Polygon;
use S2\S2Polyline;
use S2\S2Region;

class GeometryTestCase extends TestCase {

  protected function setUp(): void {
    mt_srand(123456);
  }

    protected static function f_rand(): float {
        return mt_rand(0, 100) / 100;
    }

  public static function assertDoubleNear(float $a, float $b): void {
    self::assertDoubleNearErr($a, $b, 1e-9);
  }

  public static function assertDoubleNearErr(float $a, float $b, float $error): void {
    self::assertTrue($a + $error > $b);
    self::assertTrue($a < $b + $error);
  }

  // maybe these should be put in a special testing util class
  /** Return a random unit-length vector. */
  public function randomPoint(): S2Point {
    return S2Point::normalize(new S2Point(
                                 2 * self::f_rand() - 1,
                                 2 * self::f_rand() - 1,
                                 2 * self::f_rand() - 1));
  }

  /**
   * Return a right-handed coordinate frame (three orthonormal vectors). Returns
   * an array of three points: x,y,z
   */
  public function getRandomFrame(): array {
    $p0 = $this->randomPoint();
    $p1 = S2Point::normalize(S2Point::crossProd($p0, $this->randomPoint()));
    $p2 = S2Point::normalize(S2Point::crossProd($p0, $p1));
    return [$p0, $p1, $p2];
  }

  /**
   * Return a random cell id at the given level or at a randomly chosen level.
   * The distribution is uniform over the space of cell ids, but only
   * approximately uniform over the surface of the sphere.
   */
  public function getRandomCellId(int $level): S2CellId {
    $face = $this->random(S2CellId::NUM_FACES);
    $pos  = mt_rand() & ((1 << (2 * S2CellId::MAX_LEVEL)) - 1);
    return S2CellId::fromFacePosLevel($face, $pos, $level);
  }

  public function getRandomCellIdNoLevel(): S2CellId {
    return $this->getRandomCellId($this->random(S2CellId::MAX_LEVEL + 1));
  }

    protected function random(int $n) {
      if ($n === 0) {
          return 0;
      }

      return mt_rand(0, $n);
  }

  // Pick "base" uniformly from range [0,maxLog] and then return
  // "base" random bits. The effect is to pick a number in the range
  // [0,2^maxLog-1] with bias towards smaller numbers.
  protected function skewed(int $maxLog): int {
    $base = abs(mt_rand()) % ($maxLog + 1);
    // if (!base) return 0; // if 0==base, we & with 0 below.
    //
    // this distribution differs slightly from ACMRandom's Skewed,
    // since 0 occurs approximately 3 times more than 1 here, and
    // ACMRandom's Skewed never outputs 0.
    return mt_rand() & ((1 << $base) - 1);
  }

  /**
   * Checks that "covering" completely covers the given region. If "check_tight"
   * is true, also checks that it does not contain any cells that do not
   * intersect the given region. ("id" is only used internally.)
   */
    protected function checkCovering(S2Region $region, S2CellUnion $covering, bool $checkTight, S2CellId $id): void {
        if (!$id->isValid()) {
            for ($face = 0; $face < 6; ++$face) {
                $this->checkCovering($region, $covering, $checkTight, S2CellId::fromFacePosLevel($face, 0, 0));
            }
          return;
        }

        if (!$region->mayIntersect(new S2Cell($id))) {
            // If region does not intersect id, then neither should the covering.
            if ($checkTight) {
                self::assertTrue(!$covering->intersects($id));
            }

        } else if (!$covering->contains($id)) {
            // The region may intersect id, but we can't assert that the covering
            // intersects id because we may discover that the region does not actually
            // intersect upon further subdivision. (MayIntersect is not exact.)
            self::assertTrue(!$region->contains(new S2Cell($id)));
            self::assertTrue(!$id->isLeaf());
            $end = $id->childEnd();
            for ($child = $id->childBegin(); $child->equals($end); $child = $child->next()) {
                $this->checkCovering($region, $covering, $checkTight, $child);
            }
        }
    }

    protected function getRandomCap(float $minArea, float $maxArea): S2Cap {
        $capArea = $maxArea * pow($minArea / $maxArea, self::f_rand());
        self::assertTrue($capArea >= $minArea && $capArea <= $maxArea);

        // The surface area of a cap is 2*Pi times its height.
        return S2Cap::fromAxisArea($this->randomPoint(), $capArea);
    }

  protected function samplePoint(S2Cap $cap): S2Point {
    // We consider the cap axis to be the "z" axis. We choose two other axes to
    // complete the coordinate frame.

    /* @var S2Point*/ $z = $cap->axis();
    /* @var S2Point*/ $x = $z->ortho();
    /* @var S2Point*/ $y = S2Point::crossProd($z, $x);

    // The surface area of a spherical cap is directly proportional to its
    // height. First we choose a random height, and then we choose a random
    // point along the circle at that height.

    /* @var float */ $h = self::f_rand() * $cap->height();
    /* @var float */ $theta = 2 * S2::M_PI * self::f_rand();
    /* @var float */ $r = sqrt($h * (2 - $h)); // Radius of circle.

    // (cos(theta)*r*x + sin(theta)*r*y + (1-h)*z).Normalize()
    return S2Point::normalize(S2Point::add(
                                 S2Point::add(S2Point::mul($x, cos($theta) * $r), S2Point::mul($y, sin($theta) * $r)),
                                 S2Point::mul($z, (1 - $h))));
  }

    /**
     * @param string $str
     * @param S2Point[] $vertices
     */
  static function parseVertices(string $str, array &$vertices): void {
    if ($str === null) {
        return;
    }

    $tokens = explode(',', $str);
    foreach ($tokens as $token) {
        $colon = strpos($token, ':');
      if ($colon === false) {
          throw new RuntimeException( // IllegalArgumentException(
              "Illegal string:" + token + ". Should look like '35:20'");
      }
      $lat = (float) substr($token, 0, $colon);
      $lng = (float) substr($token, $colon + 1);
      $vertices[] = S2LatLng::fromDegrees($lat, $lng)->toPoint();
    }
  }

  static function makePoint(string $str): S2Point {
    /* @var S2Point[] */  $vertices = [];
    self::parseVertices($str, $vertices);
    return array_shift($vertices);
  }

  static function makeLoop(string $str): S2Loop {
      /* @var S2Point[] */  $vertices = [];
      self::parseVertices($str, $vertices);
    return new S2Loop($vertices);
  }

  static function makePolygon(string $str): S2Polygon {
      /* @var S2Loop[] */  $loops = [];

      $tokens = explode(';', $str);
      $tokens = array_filter($tokens);
      foreach($tokens as $token) {
          $loop = self::makeLoop($token);
          $loop->normalize();
          $loops[] = $loop;
      }

    return new S2Polygon($loops);
  }

  static function makePolyline(string $str): S2Polyline {
    $vertices = [];
    self::parseVertices($str, $vertices);
    return new S2Polyline($vertices);
  }
}
