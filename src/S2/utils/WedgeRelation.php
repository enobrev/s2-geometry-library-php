<?php
    namespace S2\utils;

    use S2\S2Point;

    interface WedgeRelation {
        function test(S2Point $a0, S2Point $ab1, S2Point $a2, S2Point $b0, S2Point $b2): int;
    }