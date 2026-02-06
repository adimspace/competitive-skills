# Competitive Mathematics — International Level (IMO & Beyond)

A comprehensive breakdown of skills and topics required to compete at the level of the **International Mathematical Olympiad (IMO)**, Putnam, and similar elite competitions.

---

## 1. Algebra

### 1.1 Inequalities
- **Classical inequalities**: AM-GM, Cauchy-Schwarz, Power Mean, Rearrangement, Chebyshev's Sum
- **Advanced inequalities**: Schur's inequality, Muirhead's inequality, Karamata (majorization)
- **Techniques**: SOS (Sum of Squares) decomposition, tangent line trick, smoothing/convexity (Jensen's inequality), normalization, homogenization
- **Substitution methods**: trigonometric substitutions, symmetric substitutions, the uvw method (Schur + Newton)

### 1.2 Polynomials
- Roots and factorization over Q, R, C
- **Vieta's formulas** (Newton's identities, power sums)
- Irreducibility criteria: Eisenstein's criterion, rational root theorem, reduction modulo p
- Symmetric polynomials and the Fundamental Theorem of Symmetric Polynomials
- Polynomials with integer values, interpolation (Lagrange)
- Resultants and discriminants (advanced)

### 1.3 Functional Equations
- Cauchy's functional equation and its variants
- **Injective / surjective / bijective** deductions
- Substitution strategies (setting x=0, x=y, etc.)
- Monotonicity and continuity arguments
- Functional equations over integers vs. rationals vs. reals
- Involutions and periodic functions
- Pointwise vs. global approaches

### 1.4 Sequences and Series
- Linear recurrences (characteristic equation, generating functions)
- Nonlinear recurrences and closed-form solutions
- Convergence arguments, bounding techniques
- Telescoping sums and products
- Periodicity detection in sequences

### 1.5 Abstract Algebra (advanced / Putnam level)
- Groups, rings, fields — basic definitions and properties
- Finite fields and their applications
- Permutation groups, group actions
- Polynomial rings, quotient rings

---

## 2. Combinatorics

### 2.1 Counting Techniques
- **Bijections** and constructive counting
- **Double counting** (combinatorial identities proved two ways)
- **Principle of Inclusion-Exclusion (PIE)** — and Mobius inversion
- **Stars and bars**, multiset counting
- Generating functions (ordinary and exponential)
- Catalan numbers, Stirling numbers, Bell numbers, Fibonacci identities
- Lattice path counting, ballot problems

### 2.2 Pigeonhole Principle
- Basic and generalized pigeonhole
- Infinite pigeonhole (Ramsey-flavored)
- Applications to number theory (e.g., sums modulo n)

### 2.3 Graph Theory (Combinatorial)
- **Ramsey theory**: R(3,3), general bounds, Ramsey numbers on graphs
- **Extremal graph theory**: Turan's theorem, forbidden subgraphs
- **Graph coloring**: chromatic number, edge coloring, list coloring
- Trees and spanning trees, Cayley's formula
- Planarity, Euler's formula (V - E + F = 2)
- Matchings and Hall's theorem
- Hamiltonian and Eulerian paths/circuits
- Degree sequences, Erdos-Gallai theorem

### 2.4 Combinatorial Geometry
- Convex hull arguments and Helly's theorem
- Picking's theorem and lattice point problems
- Coloring of the plane, tiling problems
- Extremal problems in discrete geometry

### 2.5 Extremal & Probabilistic Combinatorics
- Extremal set theory: Sperner's theorem, Dilworth's theorem, sunflower lemma
- The **probabilistic method** (Lovasz Local Lemma, first/second moment methods)
- Turán-type problems
- Anti-chains, chains, and partially ordered sets

### 2.6 Invariants and Monovariants
- Parity arguments
- Coloring invariants (checkerboard, modular coloring)
- **Monovariants**: quantities that only increase (or decrease)
- Weight function / potential function arguments
- Extremal principle (choosing the min/max element)

### 2.7 Combinatorial Games & Strategies
- Nim and Sprague-Grundy theory
- Strategy stealing arguments
- Pairing strategies
- Games on graphs and boards

### 2.8 Algorithms and Processes
- Greedy algorithms in combinatorial proofs
- Inductive constructions
- Iterative processes and convergence arguments

---

## 3. Geometry

### 3.1 Euclidean Geometry — Triangle Theory
- **Cevians**: medians, altitudes, angle bisectors, symmedians
- Notable points: centroid (G), orthocenter (H), circumcenter (O), incenter (I), excenters
- Nine-point circle, Euler line, Simson line
- Stewart's theorem, angle bisector theorem, mass point geometry
- **Ceva's theorem** and **Menelaus' theorem** (trigonometric forms too)

### 3.2 Circle Geometry
- Power of a point, radical axes, radical center
- **Cyclic quadrilaterals**: Ptolemy's theorem, properties of inscribed angles
- Tangent lines, tangent-tangent and secant-tangent relations
- Mixtilinear incircles
- Apollonius circle, coaxial circles

### 3.3 Projective Geometry
- Cross-ratio and harmonic conjugates
- **Poles and polars** with respect to a circle/conic
- Harmonic division, harmonic bundles
- Projective transformations, perspectivities
- Desargues' theorem, Pascal's theorem, Brianchon's theorem

### 3.4 Inversive Geometry
- Circle inversion: definition, properties, invariants
- Mapping lines/circles to lines/circles
- Applications to tangency problems (Ptolemy via inversion)
- Inversion distance formula

### 3.5 Transformations
- **Isometries**: translations, rotations, reflections, glide reflections
- **Similarities**: spiral similarities (rotation + homothety)
- **Homothety**: homothetic center, ratio
- Composition of transformations
- Applications to proving concyclicity and collinearity

### 3.6 Trigonometric Methods
- Law of sines, law of cosines, extended law of sines
- Area formulas (Heron's formula, trigonometric area)
- Trigonometric identities in geometry
- Trigonometric form of Ceva and Menelaus

### 3.7 Coordinate & Analytic Methods
- **Cartesian coordinates**: strategic coordinate placement, computation
- **Complex numbers**: rotation, spiral similarity, collinearity/concyclicity conditions
- **Barycentric coordinates**: displacement vectors, distance formula, equation of a line, circle equations
- **Vectors**: dot product, cross product, area computations

### 3.8 Geometric Inequalities
- Triangle inequality and its generalizations
- Erdos-Mordell inequality
- Isoperimetric inequality
- Geometric optimization (Fermat point, etc.)

---

## 4. Number Theory

### 4.1 Divisibility & Primes
- GCD, LCM, and the Euclidean algorithm
- **Fundamental Theorem of Arithmetic** (unique prime factorization)
- Bezout's identity and linear Diophantine equations
- Sieve of Eratosthenes, properties of primes
- Lifting the exponent lemma (LTE)

### 4.2 Modular Arithmetic
- Congruences, residue classes
- **Fermat's Little Theorem**, **Euler's theorem**
- **Chinese Remainder Theorem (CRT)**
- Wilson's theorem
- Modular inverses, solving linear congruences
- **Order of an element** modulo n, **primitive roots**

### 4.3 p-adic Valuations (v_p)
- Definition and properties of the p-adic valuation
- Legendre's formula (v_p(n!))
- Applications to divisibility problems
- Lifting the exponent lemma (detailed use)

### 4.4 Quadratic Residues
- **Legendre symbol**, Euler's criterion
- **Quadratic reciprocity** (Gauss's law)
- Jacobi symbol
- Quadratic residues modulo primes and prime powers

### 4.5 Diophantine Equations
- **Pythagorean triples** (parametric solutions)
- **Pell's equation** (x^2 - Dy^2 = 1) and continued fractions
- Fermat's method of infinite descent
- **Vieta jumping** (a.k.a. root flipping) — IMO 1988 P6 technique
- Zsygmondy's theorem
- Sophie Germain identity
- Sum of two/four squares theorems

### 4.6 Arithmetic Functions
- Euler's totient function phi(n)
- Divisor function d(n), sum of divisors sigma(n)
- Mobius function mu(n), Mobius inversion formula
- Multiplicativity and Dirichlet convolution

### 4.7 Advanced Topics
- Hensel's lemma (lifting solutions mod p to mod p^k)
- Cyclotomic polynomials
- Algebraic integers and norms (Gaussian integers Z[i], Eisenstein integers)
- Bertrand's postulate and prime gaps
- Legendre's conjecture applications

---

## 5. Cross-Cutting Skills

### 5.1 Problem-Solving Strategies
- Working backwards, extreme cases
- Small case analysis and pattern recognition
- Generalization and specialization
- Symmetry exploitation
- Wishful thinking and construction

### 5.2 Proof Techniques
- Direct proof, proof by contradiction, proof by contrapositive
- **Mathematical induction** (weak, strong, structural, transfinite)
- Infinite descent
- Constructive vs. existential proofs
- Well-ordering principle

### 5.3 Writing & Communication
- Clear, rigorous mathematical writing
- Proper logical flow and justification of every step
- Handling edge cases explicitly
- Knowing what level of rigor is expected

---

## Recommended Progression

| Stage | Focus |
|-------|-------|
| **Foundation** | Basic inequalities, Euclidean geometry, modular arithmetic, counting |
| **Intermediate** | Functional equations, projective geometry, generating functions, p-adic valuations |
| **Advanced** | Probabilistic method, Vieta jumping, barycentric/complex coordinates, advanced NT |
| **IMO-ready** | Full synthesis, speed, elegance, multi-topic problems |
