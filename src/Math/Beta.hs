module Math.Beta
    ( Beta(..)
    ) where

import Math.Beta.Incomplete
import Math.Gamma hiding (beta)
import GHC.Float (float2Double, double2Float)
import qualified Data.Vector.Unboxed as U
import Math.Quadrature.Gaussian

-- |Minimum implementation: {lnBeta or beta} and {i_ or b_}
class Floating a => Beta a where
    lnBeta :: a -> a -> a
    lnBeta z w = log (beta z w)
    
    beta :: a -> a -> a
    beta z w = exp (lnBeta z w)
    
    i_ :: a -> a -> a -> a
    i_ x a b = b_ x a b / beta a b
    
    b_ :: a -> a -> a -> a
    b_ x a b = i_ x a b * beta a b

instance Beta Float where
    lnBeta z w = double2Float (lnBeta (float2Double z) (float2Double w))
    beta z w = double2Float (beta (float2Double z) (float2Double w))
    i_ x a b = double2Float (i_ (float2Double x) (float2Double a) (float2Double b))
    b_ x a b = double2Float (b_ (float2Double x) (float2Double a) (float2Double b))
instance Beta Double where
    lnBeta z w = lnGamma z + lnGamma w - lnGamma (z+w)
    i_ = checkBounds "i_" (incompleteBeta doubleQRule)

doubleQRule :: QRule U.Vector Double
doubleQRule = gaussLegendre 0 1 18 1e-16
