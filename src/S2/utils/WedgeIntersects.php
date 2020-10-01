<?php
    namespace S2\utils;

    use S2\S2;
    use S2\S2Point;

    class WedgeIntersects implements WedgeRelation {
        /**
         * Given two edge chains (see WedgeRelation above), this function returns -1
         * if the region to the left of A intersects the region to the left of B,
         * and 0 otherwise. Note that regions are defined such that points along a
         * boundary are contained by one side or the other, not both. So for
         * example, if A,B,C are distinct points ordered CCW around a vertex O, then
         * the wedges BOA, AOC, and COB do not intersect.
         */
         public function test(S2Point $a0, S2Point $ab1, S2Point $a2, S2Point $b0, S2Point $b2): int {
         // For A not to intersect B (where each loop interior is defined to be
         // its left side), the CCW edge order around ab1 must be a0 b2 b0 a2.
         // Note that it's important to write these conditions as negatives
         // (!OrderedCCW(a,b,c,o) rather than Ordered(c,b,a,o)) to get correct
         // results when two vertices are the same.
             return (S2::orderedCCW($a0, $b2, $b0, $ab1) && S2::orderedCCW($b0, $a2, $a0, $ab1) ? 0 : -1);
         }
     }