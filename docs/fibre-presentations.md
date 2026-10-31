# Selecting exact split-torsion presentations

For a fully split cubic and a class in the rational 2-torsion Kummer image,
the auxiliary cover has a known rational point over that torsion point. Its
component at the corresponding cubic root is zero; the other two components
are the checked rational square roots of `(x(T)-theta)/delta`. Using the CRT
coordinate direction at that root and one other CRT direction to parametrize
the conic gives a biquadratic branch quartic. The matrix is an invertible change
of the existing conic parameters, and all arithmetic labels and signed maps
are verified by the same presentation checker.

The fibre assembler chooses between these verified presentations and the default
one, switching only for strictly fewer primes in the **full recomputed
conservative support**. Ties retain the reference presentation: smaller
coefficients can have much worse local reduction. A regression with rational
4-torsion demonstrates this failure of coefficient size as a tie-breaker.
The selector does not delete primes from a support set or infer local solubility
from a small model. For `y^2=x^3-289x`, this avoids artificial
leading-coefficient primes such as 13 and 89. Every relative candidate and local
condition for the selected presentation is still computed. This is a limited
exact presentation optimization, not general minimization/reduction.

See [the global fibre contract](four-selmer-fibres.md) for equivalence and evidence conventions.
