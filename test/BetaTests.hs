module BetaTests where

import Control.Applicative
import Data.List
import Math.Beta
import System.Random (Random)
import Test.QuickCheck

newtype X a = X a deriving Show
instance (Arbitrary a, Random a, Num a) => Arbitrary (X a) where
    arbitrary = X <$> choose (0,1)

prop_within eps a b = abs (a-b) <= eps

prop_withinR eps a b = abs (a-b) <= eps * s
    where s = max (abs a) (abs b)

prop_complement_i eps (X x) (Positive a) (Positive b)
    = prop_within eps 1 (i_ x a b + i_ (1-x) b a)

prop_shift_i eps (X x) (Positive a) (Positive b)
    = prop_withinR eps
        (i_ x (a+1) b + delta)
        (i_ x a b)
    where
        delta = exp lnDelta
        lnDelta = a * log x + b * log (1-x) - log a - lnBeta a b

prop_shift_b eps (X x) (Positive a) (Positive b)
    = prop_withinR eps
        ((a+b) * b_ x (a+1) b + delta)
        (a * b_ x a b)
    where
        mu = a / (a + b)
        delta = exp lnDelta
        lnDelta = a * log x + b * log (1-x)
