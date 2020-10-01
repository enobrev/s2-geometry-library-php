<?php
namespace S2\Tests;

require dirname(__DIR__) . '/vendor/autoload.php';
require __DIR__ . '/GeometryTestCase.php';

use S2\S1Angle;
use S2\S2;
use S2\S2Cap;
use S2\S2Cell;
use S2\S2CellId;
use S2\S2CellUnion;
use S2\S2Point;
use S2\S2Projections;
use S2\S2RegionCoverer;

class S2CellUnionTest extends GeometryTestCase {
  //public static Logger logger = Logger.getLogger(S2CellUnionTest.class.getName());

  public function testBasic(): void {
    //logger.info("TestBasic");

    $empty = new S2CellUnion();
    $ids = [];
    $empty->initFromCellIds($ids);
    self::assertEquals(0, $empty->size());

    $face1Id = S2CellId::fromFacePosLevel(1, 0, 0);
    $face1Union = new S2CellUnion();
    $ids[] = $face1Id;
    $face1Union->initFromCellIds($ids);
    self::assertEquals(1, $face1Union->size());
    self::assertEquals($face1Id, $face1Union->cellId(0));

    $face2Id = S2CellId::fromFacePosLevel(2, 0, 0);
    $face2Union = new S2CellUnion();
    $cellids = [];
    $cellids[] = $face2Id->id();
    $face2Union->initFromIds($cellids);
    self::assertEquals(1, $face2Union->size());
    self::assertEquals($face2Id, $face2Union->cellId(0));

    $face1Cell = new S2Cell($face1Id);
    $face2Cell = new S2Cell($face2Id);

    self::assertTrue($face1Union->contains($face1Cell));
    self::assertFalse($face1Union->contains($face2Cell));
  }

  public function testContainsCellUnion(): void {
    // logger.info("TestContainsCellUnion");

    $randomCells = [];
    for ($i = 0; $i < 8; $i++) {
        $randomCells[] = $this->getRandomCellId(S2CellId::MAX_LEVEL);
    }

    $union = new S2CellUnion();
    $union->initFromCellIds($randomCells);

    $randomCells[] = $this->getRandomCellId(S2CellId::MAX_LEVEL);

    $unionPlusOne = new S2CellUnion();
    $unionPlusOne->initFromCellIds($randomCells);

    self::assertTrue($unionPlusOne->contains($union));
    self::assertFalse($union->contains($unionPlusOne));

    // Build the set of parent cells and check containment
    $parents = [];
    foreach ($union->iterator() as $cellId) {
        $parents[] = $cellId->parent();
    }

    $parentUnion = new S2CellUnion();
    $parentUnion->initFromCellIds($parents);

    self::assertTrue($parentUnion->contains($union));
    self::assertFalse($union->contains($parentUnion));
  }

  private function addCells(S2CellId $id, bool $selected, array $input, array $expected): void {
    // Decides whether to add "id" and/or some of its descendants to the
    // test case. If "selected" is true, then the region covered by "id"
    // *must* be added to the test case (either by adding "id" itself, or
    // some combination of its descendants, or both). If cell ids are to
    // the test case "input", then the corresponding expected result after
    // simplification is added to "expected".

    if ($id->equals(S2CellId::none())) {
      // Initial call: decide whether to add cell(s) from each face.
      for ($face = 0; $face < 6; ++$face) {
        $this->addCells(S2CellId::fromFacePosLevel($face, 0, 0), false, $input, $expected);
      }
      return;
    }
    if ($id->isLeaf()) {
      // The rnd.OneIn() call below ensures that the parent of a leaf cell
      // will always be selected (if we make it that far down the hierarchy).
      self::assertTrue($selected);
      $input[] = $id;
      return;
    }
    // The following code ensures that the probability of selecting a cell
    // at each level is approximately the same, i.e. we test normalization
    // of cells at all levels.
    if (!$selected && $this->random(S2CellId::MAX_LEVEL - $id->level()) !== 0) {
      // Once a cell has been selected, the expected output is predetermined.
      // We then make sure that cells are selected that will normalize to
      // the desired output.
      $expected[] = $id;
      $selected = true;
    }

    // With the rnd.OneIn() constants below, this function adds an average
    // of 5/6 * (kMaxLevel - level) cells to "input" where "level" is the
    // level at which the cell was first selected (level 15 on average).
    // Therefore the average number of input cells in a test case is about
    // (5/6 * 15 * 6) = 75. The average number of output cells is about 6.

    // If a cell is selected, we add it to "input" with probability 5/6.
    $added = false;
    if ($selected && $this->random(6) !== 0) {
      $input[] = $id;
      $added = true;
    }
    $numChildren = 0;
    $child = $id->childBegin();
    for ($pos = 0; $pos < 4; ++$pos, $child = $child->next()) {
      // If the cell is selected, on average we recurse on 4/12 = 1/3 child.
      // This intentionally may result in a cell and some of its children
      // being included in the test case.
      //
      // If the cell is not selected, on average we recurse on one child.
      // We also make sure that we do not recurse on all 4 children, since
      // then we might include all 4 children in the input case by accident
      // (in which case the expected output would not be correct).
      if ($this->random($selected ? 12 : 4) == 0 && $numChildren < 3) {
        $this->addCells($child, $selected, $input, $expected);
        ++$numChildren;
      }
      // If this cell was selected but the cell itself was not added, we
      // must ensure that all 4 children (or some combination of their
      // descendents) are added.
      if ($selected && !$added) {
        $this->addCells($child, $selected, $input, $expected);
      }
    }
  }

  public function testNormalize(): void {
    //logger.info("TestNormalize");

    // Try a bunch of random test cases, and keep track of average
    // statistics for normalization (to see if they agree with the
    // analysis above).
    $cellunion = new S2CellUnion();
    $inSum = 0;
    $outSum = 0;
    $kIters = 2000;
    for($i = 0; $i < $kIters; ++$i) {
      $input = [];
      $expected = [];
      $this->addCells(S2CellId::none(), false, $input, $expected);
      $inSum += count($input);
      $outSum += count($expected);
      $cellunion->initFromCellIds($input);
      self::assertEquals($cellunion->size(), count($expected));

      self::assertEquals($expected, $cellunion->cellIds());

      // Test GetCapBound().
      $cap = $cellunion->getCapBound();
      for ($k = 0; $k < $cellunion->size(); ++$k) {
        self::assertTrue($cap->contains(new S2Cell($cellunion->cellId($k))));
      }

      // Test Contains(S2CellId) and Intersects(S2CellId).
      for($j = 0; $j < count($input); ++$j) {
        self::assertTrue($cellunion->contains($input[$j]));
        self::assertTrue($cellunion->intersects($input[$j]));
        if (!$input[$j]->isFace()) {
          self::assertTrue($cellunion->intersects($input[$j]->parent()));
          if ($input[$j]->level() > 1) {
            self::assertTrue($cellunion->intersects($input[$j]->parent()->parent()));
            self::assertTrue($cellunion->intersects($input[$j]->parent(0)));
          }
        }
        if (!$input[$j]->isLeaf()) {
          self::assertTrue($cellunion->contains($input[$j]->childBegin()));
          self::assertTrue($cellunion->intersects($input[$j]->childBegin()));
          self::assertTrue($cellunion->contains($input[$j]->childEnd()->prev()));
          self::assertTrue($cellunion->intersects($input[$j]->childEnd()->prev()));
          self::assertTrue($cellunion->contains($input[$j]->childBegin(S2CellId::MAX_LEVEL)));
          self::assertTrue($cellunion->intersects($input[$j]->childBegin(S2CellId::MAX_LEVEL)));
        }
      }
      for($j = 0; $j < count($expected); ++$j) {
        if (!$expected[$j]->isFace()) {
          self::assertTrue(!$cellunion->contains($expected[$j]->parent()));
          self::assertTrue(!$cellunion->contains($expected[$j]->parent(0)));
        }
      }

      // Test contains(S2CellUnion) and intersects(S2CellUnion)
      $x = [];
      $y = [];
      $xOrY = [];
      $xAndY = [];
      for($j = 0; $j < count($input); ++$j) {
        $inX = $this->random(2) == 0;
        $inY = $this->random(2) == 0;
        if ($inX) {
          $x[] = $input[$j];
        }
        if ($inY) {
            $y[] = $input[$j];
        }
        if ($inX || $inY) {
            $xOrY[] = $input[$j];
        }
      }
      $xCells = new S2CellUnion();
      $yCells = new S2CellUnion();
      $xOrYExpected = new S2CellUnion();
      $xAndYExpected = new S2CellUnion();
      $xCells->initFromCellIds($x);
      $yCells->initFromCellIds($y);
      $xOrYExpected->initFromCellIds($xOrY);

      $xOrYCells = new S2CellUnion();
      $xOrYCells->getUnion($xCells, $yCells);
      self::assertEquals($xOrYExpected, $xOrYCells);

      // Compute the intersection of "x" with each cell of "y",
      // check that this intersection is correct, and append the
      // results to xAndYExpected.
      for($j = 0; $j < $yCells->size(); ++$j) {
        $yId = $yCells->cellId($j);
        $u = new S2CellUnion();
        $u->getIntersection($xCells, $yId);
        for($k = 0; $k < $xCells->size(); ++$k) {
          $xId = $xCells->cellId($k);
          if ($xId->contains($yId)) {
            assertEquals(1, $u->size());
            assertEquals($yId, $u->cellId(0));
          } else if ($yId->contains($xId)) {
            if (!$u->contains($xId)) {
              $u->getIntersection($xCells, $yId);
            }
            assertTrue($u->contains($xId));
          }
        }
        for($k = 0; $k < $u->size(); ++$k) {
          assertTrue($xCells->contains($u->cellId($k)));
          assertTrue($yId->contains($u->cellId($k)));
        }
        $xAndY = array_merge($xAndY, $u->cellIds());
      }
      $xAndYExpected->initFromCellIds($xAndY);

      $xAndYCells = new S2CellUnion();
      $xAndYCells->getIntersection($xCells, $yCells);
      self::assertEquals($xAndYExpected, $xAndYCells);

      $test = [];
      $dummy = [];

      $this->addCells(S2CellId::none(), false, $test, $dummy);
      for($j = 0; $j < count($test); ++$j) {
        $contains = false;
        $intersects = false;
        for($k = 0; $k < count($expected); ++$k) {
          if ($expected[$k]->contains($test[$j])) {
            $contains = true;
          }
          if ($expected[$k]->intersects($test[$j])) {
            $intersects = true;
          }
        }
        self::assertEquals($cellunion->contains($test[$j]), $contains);
        self::assertEquals($cellunion->intersects($test[$j]), $intersects);
      }

    }
  }

  function getMaxAngle(S2CellUnion $covering, S2Point $axis): float {
    $maxAngle = 0;
    for($i = 0; $i < $covering->size(); ++$i) {
      $cell = new S2Cell($covering->cellId($i));
      $cellCap = $cell->getCapBound();
      $angle = $axis->angle($cellCap->axis()) + $cellCap->angle()->radians();
      $maxAngle = max($maxAngle, $angle);
    }
    return $maxAngle;
  }

  public function testExpand(): void {
    //logger.info("TestExpand");

    // This test generates coverings for caps of random sizes, and expands
    // the coverings by a random radius, and then make sure that the new
    // covering covers the expanded cap. It also makes sure that the
    // new covering is not too much larger than expected.

    $coverer = new S2RegionCoverer();
    for($i = 0; $i < 1000; ++$i) {
      $cap = $this->getRandomCap(S2Cell::averageArea(S2CellId::MAX_LEVEL), 4 * S2::M_PI);

      // Expand the cap by a random factor whose log is uniformly distributed
      // between 0 and log(1e2).
      $expandedCap =
          S2Cap::fromAxisHeight($cap->axis(), min(2.0, pow(1e2, $this->f_rand()) * $cap->height()));

      $radius = $expandedCap->angle()->radians() - $cap->angle()->radians();
      $maxLevelDiff = $this->random(8);

      $covering = new S2CellUnion();
      $coverer->setMaxCells(1 + $this->skewed(10));
      $coverer->getCovering($cap, $covering);
      $this->checkCovering($cap, $covering, true, new S2CellId());

      $maxAngle = $this->getMaxAngle($covering, $cap->axis());
      $minLevel = S2CellId::MAX_LEVEL;
      for($j = 0; $j < $covering->size(); ++$j) {
        $minLevel = min($minLevel, $covering->cellId($j)->level());
      }
      $covering->expand(S1Angle::sradians($radius), $maxLevelDiff);
      $this->checkCovering($expandedCap, $covering, false, new S2CellId());

      $expandLevel =
          min($minLevel + $maxLevelDiff, S2Projections::$PROJ->minWidth->getMaxLevel($radius));
      $expandedMaxAngle = $this->getMaxAngle($covering, $cap->axis());

      // If the covering includes a tiny cell along the boundary, in theory the
      // maximum angle of the covering from the cap axis can increase by up to
      // twice the maximum length of a cell diagonal. We allow for an increase
      // of slightly more than this because cell bounding caps are not exact.
      self::assertTrue($expandedMaxAngle - $maxAngle <= 2.01 * S2Projections::$PROJ->maxDiag->getValue($expandLevel));
    }
  }

  public function testLeafCellsCovered(): void {
    $cellUnion = new S2CellUnion();

    // empty union
    self::assertEquals(0, $cellUnion->leafCellsCovered());

    $ids = [];
    $ids[] = S2CellId::fromFacePosLevel(0, (1 << ((S2CellId::MAX_LEVEL << 1) - 1)), S2CellId::MAX_LEVEL);

    // One leaf on face 0.
    $cellUnion->initFromCellIds($ids);
    self::assertEquals(1, $cellUnion->leafCellsCovered());

    // Face 0.
    $ids[] = S2CellId::fromFacePosLevel(0, 0, 0);
    $cellUnion->initFromCellIds($ids);
    self::assertEquals(1 << 60, $cellUnion->leafCellsCovered());

    // Five faces.
    $cellUnion->expand(0);
    self::assertEquals(5 << 60, $cellUnion->leafCellsCovered());

    // Whole world.
    $cellUnion->expand(0);
    self::assertEquals(6 << 60, $cellUnion->leafCellsCovered());

    // Add some disjoint cells.
    $ids[] = S2CellId::fromFacePosLevel(1, 0, 1);
    $ids[] = S2CellId::fromFacePosLevel(2, 0, 2);
    $ids[] = S2CellId::fromFacePosLevel(2, (1 << 60), 2);
    $ids[] = S2CellId::fromFacePosLevel(3, 0, 14);
    $ids[] = S2CellId::fromFacePosLevel(4, (1 << 60), 15);
    $ids[] = S2CellId::fromFacePosLevel(4, 0, 27);
    $ids[] = S2CellId::fromFacePosLevel(5, 0, 30);
    $cellUnion->initFromCellIds($ids);
    $expected = 1 + (1 << 6) + (1 << 30) + (1 << 32) + (2 << 56) + (1 << 58) + (1 << 60);
    self::assertEquals($expected, $cellUnion->leafCellsCovered());
  }

  public function testAverageBasedArea(): void {
    $cellUnion = new S2CellUnion();

    // empty union
    self::assertEquals(0.0, $cellUnion->averageBasedArea());

    $ids = [];
    $ids[] = S2CellId::fromFacePosLevel(1, 0, 1);
    $ids[] = S2CellId::fromFacePosLevel(5, 0, 30);
    $cellUnion->initFromCellIds($ids);

    $expected = S2Cell::averageArea(S2CellId::MAX_LEVEL) * (1 + (1 << 58));
    self::assertEquals($expected, $cellUnion->averageBasedArea());
  }

  public function testApproxArea(): void {
    $cellUnion = new S2CellUnion();

    // empty union
    self::assertEquals(0.0, $cellUnion->approxArea());

    $ids = [];
    $ids[] = S2CellId::fromFacePosLevel(1, 0, 1);
    $ids[] = S2CellId::fromFacePosLevel(5, 0, 30);
    $cellUnion->initFromCellIds($ids);

    $expected = (new S2Cell($ids[0]))->approxArea() + (new S2Cell($ids[1]))->approxArea();
    self::assertEquals($expected, $cellUnion->approxArea());
  }

  public function testExactArea(): void {
    $cellUnion = new S2CellUnion();

    // empty union
    self::assertEquals(0.0, $cellUnion->exactArea());

    $ids = [];
    $ids[] = S2CellId::fromFacePosLevel(1, 0, 1);
    $ids[] = S2CellId::fromFacePosLevel(2, 0, 30);
    $cellUnion->initFromCellIds($ids);

    $expected = (new S2Cell($ids[0]))->exactArea() + (new S2Cell($ids[1]))->exactArea();

    self::assertEquals($expected, $cellUnion->averageBasedArea());
  }
}
