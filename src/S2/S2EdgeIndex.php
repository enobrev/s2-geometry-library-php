<?php

namespace S2;

use S2\utils\Math;

class S2EdgeIndex
{
    /**
     * Thicken the edge in all directions by roughly 1% of the edge length when
     * thickenEdge is true.
     */
    const THICKENING = 0.01;

    /**
     * Threshold for small angles, that help lenientCrossing to determine whether
     * two edges are likely to intersect.
     */
    const MAX_DET_ERROR = 1e-14;

    /**
     * The cell containing each edge, as given in the parallel array
     * <code>edges</code>.
     */
    private $cells;

    /**
     * The edge contained by each cell, as given in the parallel array
     * <code>cells</code>.
     */
    private $edges;

    /**
     * No cell strictly below this level appears in mapping. Initially leaf level,
     * that's the minimum level at which we will ever look for test edges.
     */
    private $minimumS2LevelUsed;

    /**
     * Has the index been computed already?
     * @var bool
     */
    private $indexComputed;

    /**
     * Number of queries so far
     */
    private $queryCount;

    /**
     * Empties the index in case it already contained something.
     */
    public function reset()
    {
        $this->minimumS2LevelUsed = S2CellId::MAX_LEVEL;
        $this->indexComputed = false;
        $this->queryCount = 0;
        $this->cells = null;
        $this->edges = null;
    }

    /**
     * Compares [cell1, edge1] to [cell2, edge2], by cell first and edge second.
     *
     * @param $cell1
     * @param $edge1
     * @param $cell2
     * @param $edge2
     * @return int -1 if [cell1, edge1] is less than [cell2, edge2], 1 if [cell1,
     */
    private static function compare($cell1, $edge1, $cell2, $edge2)
    {
        if ($cell1 < $cell2) {
            return -1;
        } else if ($cell1 > $cell2) {
            return 1;
        } else if ($edge1 < $edge2) {
            return -1;
        } else if ($edge1 > $edge2) {
            return 1;
        } else {
            return 0;
        }
    }

    /** Computes the index (if it has not been previously done). */
     public final function computeIndex() {
         if ($this->indexComputed) {
            return;
         }
         $cellList = [];
         $edgeList = [];
         for ($i = 0; $i < $this->getNumEdges(); ++$i) {
             $from = $this->edgeFrom($i);
             $to = $this->edgeTo($i);
             $cover = [];
             $level = $this->getCovering($from, $to, true, $cover);
             $minimumS2LevelUsed = min($minimumS2LevelUsed, $level);
             foreach ($cover as $cellId) {
                 $cellList[] = $cellId->id();
                 $edgeList[] = $i;
             }
         }
         $this->cells = [];
         $this->edges = [];
         for ($i = 0; $i < count($this->cells); $i++) {
             $this->cells[$i] = $cellList[$i];
             $this->edges[$i] = $edgeList[$i];
         }
         $this->sortIndex();
         $this->indexComputed = true;
     }

     /** Sorts the parallel <code>cells</code> and <code>edges</code> arrays. */
     private function sortIndex() {
         // create an array of indices and sort based on the values in the parallel
         // arrays at each index
         $indices = [];
         for ($i = 0; $i < count($this->cells); $i++) {
            $indices[$i] = $i;
         }
         usort($indices, function(int $index1, int $index2) {
             return S2EdgeIndex::compare(
                 $this->cells[$index1],
                 $this->edges[$index1],
                 $this->cells[$index2],
                 $this->edges[$index2]
             );
         });

         // copy the cells and edges in the order given by the sorted list of indices
         $newCells = [];
         $newEdges = [];
         for ($i = 0; $i < count($indices); $i++) {
            $newCells[$i] = $this->cells[$indices[$i]];
            $newEdges[$i] = $this->edges[$indices[$i]];
         }
         // replace the cells and edges with the sorted arrays
         $this->cells = $newCells;
         $this->edges = $newEdges;
     }

     public function isIndexComputed(): bool {
        return $this->indexComputed;
     }

      /**
     * Tell the index that we just received a new request for candidates. Useful
     * to compute when to switch to quad tree.
     */
     public final function incrementQueryCount() {
        ++$this->queryCount;
     }

     /**
     * If the index hasn't been computed yet, looks at how much work has gone into
     * iterating using the brute force method, and how much more work is planned
     * as defined by 'cost'. If it were to have been cheaper to use a quad tree
     * from the beginning, then compute it now. This guarantees that we will never
     * use more than twice the time we would have used had we known in advance
     * exactly how many edges we would have wanted to test. It is the theoretical
     * best.
     *
     * The value 'n' is the number of iterators we expect to request from this
     * edge index.
     *
     * If we have m data edges and n query edges, then the brute force cost is m
     * * n * testCost where testCost is taken to be the cost of
     * EdgeCrosser.robustCrossing, measured to be about 30ns at the time of this
     * writing.
     *
     * If we compute the index, the cost becomes: m * costInsert + n *
     * costFind(m)
     *
     * - costInsert can be expected to be reasonably stable, and was measured at
     * 1200ns with the BM_QuadEdgeInsertionCost benchmark.
     *
     * - costFind depends on the length of the edge . For m=1000 edges, we got
     * timings ranging from 1ms (edge the length of the polygon) to 40ms. The
     * latter is for very long query edges, and needs to be optimized. We will
     * assume for the rest of the discussion that costFind is roughly 3ms.
     *
     * When doing one additional query, the differential cost is m * testCost -
     * costFind(m) With the numbers above, it is better to use the quad tree (if
     * we have it) if m >= 100.
     *
     * If m = 100, 30 queries will give m*n*testCost = m * costInsert = 100ms,
     * while the marginal cost to find is 3ms. Thus, this is a reasonable thing to
     * do.
     */
     public final function predictAdditionalCalls(int $n) {
         if ($this->indexComputed) {
            return;
         }
         if ($this->getNumEdges() > 100 && ($this->queryCount + $n) > 30) {
            $this->computeIndex();
         }
     }

     /**
     * Overwrite these functions to give access to the underlying data. The
     * function getNumEdges() returns the number of edges in the index, while
     * edgeFrom(index) and edgeTo(index) return the "from" and "to" endpoints of
     * the edge at the given index.
     */
     private $numEdges;

    public function getNumEdges(): int {
         return $this->numEdges;
     }

     public function setNumEdges(int $numEdges): void { // php cannot override methods on the fly, so let's just use a setter/getter
        $this->numEdges = $numEdges;
     }

     /** @var callable */
     private $_edgeFrom;

     protected function edgeFrom(int $index) {
         return ($this->_edgeFrom)($index);
     }

     public function setEdgeFrom(callable $edgeFrom): void { // php cannot override methods on the fly, so let's just use a setter/getter
        $this->_edgeFrom = $edgeFrom;
     }

    /** @var callable */
     private $_edgeTo;

     protected function edgeTo(int $index) {
         return ($this->_edgeTo)($index);
     }

     public function setEdgeTo(callable $edgeTo): void { // php cannot override methods on the fly, so let's just use a setter/getter
        $this->_edgeTo = $edgeTo;
     }
     /*
     * protected abstract function edgeFrom(int index);
     * protected abstract function edgeTo(int index);

     /**
     * Appends to "candidateCrossings" all edge references which may cross the
     * given edge. This is done by covering the edge and then finding all
     * references of edges whose coverings overlap this covering. Parent cells are
     * checked level by level. Child cells are checked all at once by taking
     * advantage of the natural ordering of S2CellIds.
     */
    /**
     * Returns the smallest cell containing all four points, or
     * {@link S2CellId#sentinel()} if they are not all on the same face. The
     * points don't need to be normalized.
     */
     private static function containingCell(S2Point $pa, S2Point $pb, S2Point $pc, S2Point $pd) {
         $a = S2CellId::fromPoint($pa);
         $b = S2CellId::fromPoint($pb);
         $c = S2CellId::fromPoint($pc);
         $d = S2CellId::fromPoint($pd);

         if ($a->face() !== $b->face() || $a->face() !== $c->face() || $a->face() !== $d->face()) {
             return S2CellId::sentinel();
         }

         while (!$a->equals($b) || !$a->equals($c) || !$a->equals($d)) {
             $a = $a->parent();
             $b = $b->parent();
             $c = $c->parent();
             $d = $d->parent();
         }
         return $a;
     }

    public function findCandidateCrossings(S2Point $a, S2Point $b, array &$candidateCrossings) {
        // Preconditions::checkState($this->indexComputed);
        $cover = [];
        $this->getCovering($a, $b, false, $cover);

        // Edge references are inserted into the map once for each covering cell, so
        // absorb duplicates here
        $uniqueSet = new \SplObjectStorage();
        $this->getEdgesInParentCells($cover, $uniqueSet);

        // TODO(user): An important optimization for long query
        // edges (Contains queries): keep a bounding cap and clip the query
        // edge to the cap before starting the descent.
        $this->getEdgesInChildrenCells($a, $b, $cover, $uniqueSet);

        array_filter($candidateCrossings, function() {return false;}); // clearing the referenced array

        foreach($uniqueSet as $value) {
            $candidateCrossings[] = $value;
        }
    }

    /**
     * Returns the smallest cell containing both points, or Sentinel if they are
     * not all on the same face. The points don't need to be normalized.
     */
     private static function containingCell_2(S2Point $pa, S2Point $pb) {
         $a = S2CellId::fromPoint($pa);
         $b = S2CellId::fromPoint($pb);

         if ($a->face() !== $b->face()) {
             return S2CellId::sentinel();
         }

         while (!$a->equals($b) ) {
             $a = $a->parent();
             $b = $b->parent();
         }
         return $a;
     }

     /**
     * Computes a cell covering of an edge. Clears edgeCovering and returns the
     * level of the s2 cells used in the covering (only one level is ever used for
     * each call).
     *
     * If thickenEdge is true, the edge is thickened and extended by 1% of its
     * length.
     *
     * It is guaranteed that no child of a covering cell will fully contain the
     * covered edge.
     */
     private function getCovering(S2Point $a, S2Point $b, bool $thickenEdge, array $edgeCovering) {

         // Selects the ideal s2 level at which to cover the edge, this will be the
         // level whose S2 cells have a width roughly commensurate to the length of
         // the edge. We multiply the edge length by 2*THICKENING to guarantee the
         // thickening is honored (it's not a big deal if we honor it when we don't
         // request it) when doing the covering-by-cap trick.
         $edgeLength = $a->angle($b);
         $idealLevel = S2Projections::$PROJ->minWidth->getMaxLevel($edgeLength * (1 + 2 * self::THICKENING));

         if (!$thickenEdge) {
            $containingCellId = self::containingCell_2($a, $b);
         } else {
             if ($idealLevel === S2CellId::MAX_LEVEL) {
                 // If the edge is tiny, instabilities are more likely, so we
                 // want to limit the number of operations.
                 // We pretend we are in a cell much larger so as to trigger the
                 // 'needs covering' case, so we won't try to thicken the edge.
                 $containingCellId = (new S2CellId(0xFFF0))->parent(3);
             } else {
                 $pq = S2Point::mul(S2Point::minus($b, $a), self::THICKENING);
                 $ortho = S2Point::mul(S2Point::normalize(S2Point::crossProd($pq, $a)), $edgeLength * self::THICKENING);
                 $p = S2Point::minus($a, $pq);
                 $q = S2Point::add($b, $pq);
                 // If p and q were antipodal, the edge wouldn't be lengthened,
                 // and it could even flip! This is not a problem because
                 // idealLevel != 0 here. The farther p and q can be is roughly
                 // a quarter Earth away from each other, so we remain
                 // Theta(THICKENING).
                 $containingCellId = self::containingCell(
                     S2Point::minus($p, $ortho),
                     S2Point::add($p, $ortho),
                     S2Point::minus($q, $ortho),
                     S2Point::add($q, $ortho)
                 );
             }
         }

         // Best case: edge is fully contained in a cell that's not too big.
         if (!$containingCellId->equals(S2CellId::sentinel()) && $containingCellId->level() >= $idealLevel - 2) {
            $edgeCovering[] = $containingCellId;
            return $containingCellId->level();
         }

         if ($idealLevel == 0) {
             // Edge is very long, maybe even longer than a face width, so the
             // trick below doesn't work. For now, we will add the whole S2 sphere.
             // TODO(user): Do something a tad smarter (and beware of the
             // antipodal case).
             for ($cellid = S2CellId::begin(0); !$cellid->equals(S2CellId::end(0)); $cellid = $cellid->next()) {
                $edgeCovering[] = $cellid;
             }
             return 0;
         }
         // TODO(user): Check trick below works even when vertex is at
         // interface
         // between three faces.

         // Use trick as in S2PolygonBuilder.PointIndex.findNearbyPoint:
         // Cover the edge by a cap centered at the edge midpoint, then cover
         // the cap by four big-enough cells around the cell vertex closest to the
         // cap center.
         $middle = S2Point::normalize(S2Point::div(S2Point::add($a, $b), 2));
         $actualLevel = min($idealLevel, S2CellId::MAX_LEVEL - 1);
         S2CellId::fromPoint($middle)->getVertexNeighbors($actualLevel, $edgeCovering);
         return $actualLevel;
     }


    /**
     * A constant holding the minimum value an {@code int} can
     * have, -2<sup>31</sup>.
     */
    private const  INTEGER_MIN_VALUE = 0x80000000;

    /**
     * A constant holding the maximum value an {@code int} can
     * have, 2<sup>31</sup>-1.
     */
    private const INTEGER_MAX_VALUE = 0x7fffffff;

     /**
     * Filters a list of entries down to the inclusive range defined by the given
     * cells, in <code>O(log N)</code> time.
     *
     * @param cell1 One side of the inclusive query range.
     * @param cell2 The other side of the inclusive query range.
     * @return An array of length 2, containing the start/end indices.
     */
     private function getEdges(int $cell1, int $cell2): array {
         // ensure cell1 <= cell2
         if ($cell1 > $cell2) {
             [$cell1, $cell2] = [$cell2, $cell1];
         }
         // The binary search returns -N-1 to indicate an insertion point at index N,
         // if an exact match cannot be found. Since the edge indices queried for are
         // not valid edge indices, we will always get -N-1, so we immediately
         // convert to N.
         return [
             -1 - $this->binarySearch($cell1, self::INTEGER_MIN_VALUE),
             -1 - $this->binarySearch($cell2, self::INTEGER_MAX_VALUE)
         ];
    }

     private function binarySearch(int $cell, int $edge) {
         $low = 0;
         $high = count($this->cells) - 1;
         while ($low <= $high) {
             $mid = JavaMathHelper::uRShift(($low + $high), 1);
             $cmp = self::compare($this->cells[$mid], $this->edges[$mid], $cell, $edge);
             if ($cmp < 0) {
                $low = $mid + 1;
             } else if ($cmp > 0) {
                $high = $mid - 1;
             } else {
                return $mid;
             }
         }
         return -($low + 1);
     }

     /**
     * Adds to candidateCrossings all the edges present in any ancestor of any
     * cell of cover, down to minimumS2LevelUsed. The cell->edge map is in the
     * variable mapping.
     */
     private function getEdgesInParentCells(array $cover, \SplObjectStorage $candidateCrossings) {
         // Find all parent cells of covering cells.
         $parentCells = new \SplObjectStorage();
         foreach ($cover as $coverCell) {
             for ($parentLevel = $coverCell->level() - 1; $parentLevel >= $this->minimumS2LevelUsed; --$parentLevel) {
                 $addingCell = $coverCell->parent($parentLevel);
                 if ($parentCells->contains($addingCell)) {
                     break;
                 }

                 $parentCells->attach($addingCell);
             }
         }

         // Put parent cell edge references into result.
         foreach ($parentCells as $parentCell) {
             $bounds = $this->getEdges($parentCell->id(), $parentCell->id());
             for ($i = $bounds[0]; $i < $bounds[1]; $i++) {
                 $candidateCrossings->attach($this->edges[$i]);
             }
         }
     }

     /**
     * Returns true if ab possibly crosses cd, by clipping tiny angles to zero.
     */
     private static function lenientCrossing(S2Point $a, S2Point $b, S2Point $c, S2Point $d) {
         // assert (S2.isUnitLength(a));
         // assert (S2.isUnitLength(b));
         // assert (S2.isUnitLength(c));

         $acb = S2Point::crossProd($a, $c)->dotProd($b);
         $bda = S2Point::crossProd($b, $d)->dotProd($a);
         if (abs($acb) < self::MAX_DET_ERROR || abs($bda) < self::MAX_DET_ERROR) {
            return true;
         }
         if ($acb * $bda < 0) {
            return false;
         }
         $cbd = S2Point::crossProd($c, $b)->dotProd($d);
         $dac = S2Point::crossProd($c, $a)->dotProd($c);
         if (abs($cbd) < self::MAX_DET_ERROR || abs($dac) < self::MAX_DET_ERROR) {
             return true;
         }
         return ($acb * $cbd >= 0) && ($acb * $dac >= 0);
     }

     /**
     * Returns true if the edge and the cell (including boundary) intersect.
     */
     private static function edgeIntersectsCellBoundary(S2Point $a, S2Point $b, S2Cell $cell) {
         $vertices = [];
         for ($i = 0; $i < 4; ++$i) {
            $vertices[$i] = $cell->getVertex($i);
         }
         for ($i = 0; $i < 4; ++$i) {
             $fromPoint = $vertices[$i];
             $toPoint = $vertices[($i + 1) % 4];
             if (self::lenientCrossing($a, $b, $fromPoint, $toPoint)) {
                return true;
             }
         }
         return false;
     }

     /**
     * Appends to candidateCrossings the edges that are fully contained in an S2
     * covering of edge. The covering of edge used is initially cover, but is
     * refined to eliminate quickly subcells that contain many edges but do not
     * intersect with edge.
     */
     private function getEdgesInChildrenCells(S2Point $a, S2Point $b, array $cover, \SplObjectStorage $candidateCrossings) {
         // Put all edge references of (covering cells + descendant cells) into
         // result.
         // This relies on the natural ordering of S2CellIds.
         $children = null;
         while (count($cover) > 0) {
             $cell = array_pop($cover);
             $bounds = $this->getEdges($cell->rangeMin()->id(), $cell->rangeMax()->id());
             if ($bounds[1] - $bounds[0] <= 16) {
                 for ($i = $bounds[0]; $i < $bounds[1]; $i++) {
                     $candidateCrossings->attach($this->edges[$i]);
                 }
             } else {
                 // Add cells at this level
                 $bounds = $this->getEdges($cell->id(), $cell->id());
                 for ($i = $bounds[0]; $i < $bounds[1]; $i++) {
                    $candidateCrossings->attach($this->edges[$i]);
                 }
                 // Recurse on the children -- hopefully some will be empty.
                 if ($children === null) {
                     $children = [];
                     for ($i = 0; $i < 4; ++$i) {
                        $children[$i] = new S2Cell();
                     }
                 }
                 (new S2Cell($cell))->subdivide($children);
                 foreach ($children as $child) {
                 // TODO(user): Do the check for the four cells at once,
                 // as it is enough to check the four edges between the cells. At
                 // this time, we are checking 16 edges, 4 times too many.
                 //
                 // Note that given the guarantee of AppendCovering, it is enough
                 // to check that the edge intersect with the cell boundary as it
                 // cannot be fully contained in a cell.
                     if (self::edgeIntersectsCellBoundary($a, $b, $child)) {
                         $cover[] = $child->id();
                     }
                 }
             }
         }
     }

}

/*
        * An iterator on data edges that may cross a query edge (a,b). Create the
        * iterator, call getCandidates(), then hasNext()/next() repeatedly.
        *
        * The current edge in the iteration has index index(), goes between from()
        * and to().
        */

class DataEdgeIterator
{
    /**
     * The structure containing the data edges.
     * @var S2EdgeIndex
     */
    private $edgeIndex;

    /**
     * Tells whether getCandidates() obtained the candidates through brute force
     * iteration or using the quad tree structure.
     */
    private $isBruteForce;

    /**
     * Index of the current edge and of the edge before the last next() call.
     */
    private $currentIndex;

    /**
     * Cache of edgeIndex.getNumEdges() so that hasNext() doesn't make an extra
     * call
     */
    private $numEdges;

    /**
     * All the candidates obtained by getCandidates() when we are using a
     * quad-tree (i.e. isBruteForce = false).
     *
     * @var int[]
     */
     public $candidates;

     /**
     * Index within array above. We have: currentIndex =
     * candidates.get(currentIndexInCandidates).
     * @var int
     */
     private $currentIndexInCandidates;

     public function __construct(S2EdgeIndex $edgeIndex) {
         $this->edgeIndex = $edgeIndex;
         $this->candidates = [];
     }

     /**
     * Initializes the iterator to iterate over a set of candidates that may
     * cross the edge (a,b).
     */
     public function getCandidates(S2Point $a, S2Point $b): void {
         $this->edgeIndex->predictAdditionalCalls(1);
         $isBruteForce = !$this->edgeIndex->isIndexComputed();
         if ($isBruteForce) {
             $this->edgeIndex->incrementQueryCount();
             $this->currentIndex = 0;
             $this->numEdges = $this->edgeIndex->getNumEdges();
         } else {
             $this->edgeIndex->findCandidateCrossings($a, $b, $this->candidates);
             $this->currentIndexInCandidates = 0;
             if (count($this->candidates) > 0) {
                $this->currentIndex = $this->candidates[0];
             }
         }
     }

     /**
     * Index of the current edge in the iteration.
     */
     public function index(): int {
        //Preconditions.checkState(hasNext());
        return $this->currentIndex;
     }
     /**
     * False if there are no more candidates; true otherwise.
     */
     public function hasNext(): bool {
        if ($this->isBruteForce) {
            return ($this->currentIndex < $this->numEdges);
        } else {
            return $this->currentIndexInCandidates < count($this->candidates);
        }
     }

     /**
     * Iterate to the next available candidate.
     */
     public function next(): void {
         //Preconditions.checkState(hasNext());
         if ($this->isBruteForce) {
            ++$this->currentIndex;
         } else {
            ++$this->currentIndexInCandidates;
            if ($this->currentIndexInCandidates < count($this->candidates)) {
                $this->currentIndex = $this->candidates[$this->currentIndexInCandidates];
            }
         }
     }

}
