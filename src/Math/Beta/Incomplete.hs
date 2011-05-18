{-# LANGUAGE FlexibleContexts #-}
module Math.Beta.Incomplete 
    ( incompleteBeta
    , checkBounds
    , uncheckedIncompleteBeta
    ) where

import Math.ContinuedFraction
import Math.Sequence.Converge
import Math.Gamma
import Math.Quadrature.Gaussian
import qualified Data.Vector.Generic as V

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

incompleteBeta :: (RealFloat a, Gamma a, Ord a, V.Vector v (a,a)) => QRule v a -> a -> a -> a -> a
incompleteBeta = checkBounds "incompleteBeta" . uncheckedIncompleteBeta

uncheckedIncompleteBeta qRule x a b
    | x == 0            = 0
    | x == 1            = 1
    | all (>useQ) [a,b] = incompleteBetaByQ qRule  x a b
    | x < (a+1)/(a+b+2) = incompleteBetaByCF       x a b
    | otherwise         = iSymm incompleteBetaByCF x a b

iSymm i x a b = 1 - i (1-x) b a

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

incompleteBetaByQ qRule x a b
    | x < muA   = 1 - betaQ x xu
    | otherwise =     betaQ xl x
    where
        -- xl and xh are points above and below x, respectively, which
        -- are far enough that they contain all the significant mass of the
        -- integrand, but (hopefully) close enough that they have decent
        -- resolution.  Only one or the other will be used; the one on the
        -- same side of the integrand's maximum ('muA') as 'x'.
        xl = max 0 (min (muA - 10*t) (x - 5*t))
        xu = min 1 (max (muA + 10*t) (x + 5*t))
        
        muA = a / (a + b)
        muB = b / (a + b)
        t = sqrt (muA*muB/(1+a+b))
        
        -- compute the actual integral using the specified endpoints.  In
        -- order to avoid underflow (which is pretty much guaranteed 
        -- otherwise), the terms in the integrand are divided by their
        -- maxima so that the maximum value of 'f' will be 1.  Corresponding
        -- terms are introduced under 'exp' to cancel that out.
        betaQ x0 x1
            = exp ((a-1) * log muA - lnGamma a + (b-1) * log muB - lnGamma b + lnGamma (a + b))
            * estimate (integrateRange qRule f x0 x1)
            where
                am1 = a-1; lnMuA = log muA
                bm1 = b-1; lnMuB = log muB
                
                f t = exp (am1 * (log t - lnMuA) + bm1 * (log (1-t) - lnMuB))
