module Math.Beta.Incomplete (incompleteBeta) where

import Math.ContinuedFraction
import Math.Sequence.Converge
import Math.Gamma

absEps, relEps :: Fractional a => a
useQ :: Num a => a
absEps  = 0
relEps  = 1e-15
useQ    = 3000

lnBeta z w = lnGamma z + lnGamma w - lnGamma (z+w)

checkBounds fName f x a b
    | x <  0    = error (fName ++ ": x < 0")
    | x >  1    = error (fName ++ ": x > 1")
    | a <  0    = error (fName ++ ": a < 0")
    | b <  0    = error (fName ++ ": b < 0")
    | otherwise = f x a b

incompleteBeta :: (Gamma a, Ord a) => a -> a -> a -> a
incompleteBeta = checkBounds "incompleteBeta" uncheckedIncompleteBeta

uncheckedIncompleteBeta x a b
    | x == 0            = 0
    | x == 1            = 1
    | all (>useQ) [a,b] = incompleteBetaByQ        x a b
    | x < (a+1)/(a+b+2) = incompleteBetaByCF       x a b
    | otherwise         = iSymm incompleteBetaByCF x a b

iSymm i x a b = 1 - i (1-x) b a

incompleteBetaByQ    x a b = error "incompleteBeta: implement quadrature version"
incompleteBetaByCF x a b = c * z / a
    where
        c   = exp (a*log x + b*log(1-x) - lnBeta a b)
        z   = convergeTo absEps relEps (lentz zCF)
        zCF = betaCF x a b

-- This computes a value Z such that
-- B_x(a,b) = (x**a * (1-x)**b / a) * Z
-- I_x(a,b) = B_x(a,b)/B(a,b) = (x**a * (1-x)**b / (a * beta a b)) * Z
-- TODO: this needs major cleanup
betaCF x a b = gcf 0 (map (\d -> (d,1)) ds)
    where
        ap1 = a+1
        ap2 = a+2
        apb = a+b
        bm1 = b-1
        ds = 1 : concat
            [ [ negate (m+a)*(m+apb)*x / ((m2+a)*(m2+ap1))
              , (m+1)*(bm1-m)*x / ((m2+ap1)*(m2+ap2))
              ]
            | m <- iterate (1+) 0
            , let m2 = 2*m
            ]

interleave []     ys = ys
interleave (x:xs) ys = x : interleave ys xs
