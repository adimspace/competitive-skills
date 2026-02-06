# Professional Applications of Competitive Math & Programming Skills

For each skill from the competitive math and competitive programming skill lists, this document describes where it is used in modern professional work. Skills with no meaningful professional application outside of pure research or competition are marked **None (competition-only)**.

---
---

# PART I — Competitive Mathematics Skills

---

## 1. Algebra

### 1.1 Inequalities

| Skill | Professional Application |
|-------|------------------------|
| Classical inequalities (AM-GM, Cauchy-Schwarz, Power Mean, etc.) | **Optimization & ML** — bounding loss functions, proving convergence rates for gradient descent, deriving PAC-learning bounds. **Finance** — risk bound estimation (e.g., portfolio variance bounds). **Operations Research** — proving optimality of LP relaxations. |
| Advanced inequalities (Schur, Muirhead, Karamata) | **None (competition-only)** — too specialized; the general mindset of bounding quantities transfers to optimization research but these specific inequalities are almost never invoked. |
| Techniques (SOS, tangent line trick, Jensen's, homogenization) | **Convex optimization** — Jensen's inequality is foundational in ML (proof of EM algorithm convergence, KL-divergence properties, variational inference). SOS decomposition is used in **polynomial optimization and control theory** (Lyapunov stability analysis via SOS programming). |
| Substitution methods (trig, symmetric, uvw) | **None (competition-only)** — contest-specific manipulation skills. |

### 1.2 Polynomials

| Skill | Professional Application |
|-------|------------------------|
| Roots and factorization | **Control systems engineering** — pole-zero analysis of transfer functions. **Signal processing** — filter design (placing roots of transfer function). **Computer algebra systems** — symbolic computation engines (Mathematica, Maple, SymPy). |
| Vieta's formulas, Newton's identities | **Signal processing** — relationship between filter coefficients and poles/zeros. **Algebraic geometry** — computational algebraic geometry tools. Niche but real. |
| Irreducibility criteria (Eisenstein, etc.) | **Cryptography** — constructing irreducible polynomials over finite fields for AES (Rijndael), CRC codes, and elliptic curve cryptography. **Coding theory** — BCH and Reed-Solomon error-correcting codes. |
| Symmetric polynomials | **Physics** — quantum mechanics (symmetric functions of eigenvalues). **Statistics** — theory of U-statistics and symmetric estimators. Niche. |
| Interpolation (Lagrange) | **Cryptography** — Shamir's Secret Sharing scheme directly uses Lagrange interpolation. **Numerical methods** — curve fitting, FEM (finite element method). **Computer graphics** — spline interpolation for animation curves. |
| Resultants and discriminants | **Robotics** — kinematics equations (eliminating variables from polynomial systems). **Computational algebraic geometry** — automated theorem proving, CAD (cylindrical algebraic decomposition). Niche. |

### 1.3 Functional Equations

| Skill | Professional Application |
|-------|------------------------|
| Cauchy's equation and variants | **None (competition-only)** — the classification of additive/multiplicative functions is a foundational result but not directly applied in industry. |
| Injectivity/surjectivity deductions | **Software engineering** — the concept of bijections is fundamental in database design (unique keys), serialization/deserialization, and hash function analysis, but the olympiad-style deduction technique itself is competition-only. |
| Substitution strategies, pointwise approaches | **None (competition-only)** |
| Functional equations over various domains | **Economics (very niche)** — characterizing utility functions, but this is academic micro-economics, not industry work. |
| Involutions and periodic functions | **None (competition-only)** |

### 1.4 Sequences and Series

| Skill | Professional Application |
|-------|------------------------|
| Linear recurrences, characteristic equation | **Algorithm analysis** — analyzing runtime of divide-and-conquer algorithms (Master theorem is a recurrence). **Finance** — modeling compound growth, amortization schedules. **Signal processing** — IIR (infinite impulse response) filter design. |
| Generating functions | **Probability theory** — probability generating functions and moment generating functions in statistics & ML. **Algorithm analysis** — analyzing expected complexity of algorithms. **Combinatorics in CS** — analytic combinatorics for average-case analysis (Flajolet-Sedgewick framework). |
| Convergence, bounding techniques | **Numerical computing** — convergence analysis of iterative solvers (Newton's method, conjugate gradient). **ML** — proving convergence of training algorithms. |
| Telescoping sums and products | **Finance** — compound return decomposition. Minor utility as a computational trick in various fields. |
| Periodicity detection | **Signal processing** — pitch detection in audio. **Time series analysis** — seasonality detection in data science. **Cryptanalysis** — detecting periods in pseudo-random sequences. |

### 1.5 Abstract Algebra

| Skill | Professional Application |
|-------|------------------------|
| Groups, rings, fields | **Cryptography** — foundational for RSA (group theory of Z/nZ), ECC (elliptic curve groups), AES (field GF(2^8)). **Coding theory** — all linear codes live in vector spaces over finite fields. **Quantum computing** — quantum error correction uses group theory extensively. |
| Finite fields | **Cryptography** — AES, elliptic curve cryptography, Diffie-Hellman over finite fields. **Networking** — CRC checksums, Reed-Solomon in QR codes, CD/DVD/Blu-ray error correction. **Database systems** — consistent hashing. |
| Permutation groups, group actions | **Chemistry** — molecular symmetry analysis (Burnside/Polya for isomer counting). **Computer graphics** — symmetry groups for procedural generation. **Quantum computing** — representation theory for quantum algorithms. |
| Polynomial rings, quotient rings | **Cryptography** — lattice-based (post-quantum) cryptography uses polynomial rings R_q = Z_q[x]/(x^n + 1). **Coding theory** — cyclic codes are ideals in polynomial quotient rings. |

---

## 2. Combinatorics

### 2.1 Counting Techniques

| Skill | Professional Application |
|-------|------------------------|
| Bijections, constructive counting | **Software engineering** — encoding/decoding schemes, database normalization. **Combinatorial testing** — enumerating test configurations. |
| Double counting | **None (competition-only)** — a proof technique, not directly applied. |
| Inclusion-Exclusion (PIE) | **Database query optimization** — computing cardinality of JOIN/UNION operations. **Probability** — calculating P(A ∪ B ∪ C) in reliability engineering. **Data science** — deduplication logic, Venn diagram-based analytics. **Security** — estimating collision probabilities. |
| Stars and bars, multiset counting | **Statistics** — multinomial distributions, sampling with replacement. **Resource allocation** — distributing workloads, inventory problems in operations research. |
| Generating functions | (Same as §1.4 — probability, algorithm analysis, analytic combinatorics.) |
| Catalan, Stirling, Bell, Fibonacci | **Computer science** — Catalan numbers appear in parsing (number of valid parenthesizations, BST shapes). Stirling numbers in database partition optimization. Bell numbers in set partition algorithms. Fibonacci in data structures (Fibonacci heaps). Mostly the concepts matter, not the numbers themselves. |
| Lattice path counting, ballot problems | **Finance** — ballot problems relate to the reflection principle used in pricing barrier options (Black-Scholes extensions). **Queueing theory** — analyzing queue behavior. Niche. |

### 2.2 Pigeonhole Principle

| Skill | Professional Application |
|-------|------------------------|
| Basic and generalized pigeonhole | **Distributed systems** — birthday paradox for collision estimation (hash table sizing, UUID collision probability). **Networking** — proving impossibility results (CAP theorem intuitions). **Compression** — proving limits of lossless compression. |
| Infinite pigeonhole | **None (competition-only)** — theoretical tool. |
| Applications to number theory | **None (competition-only)** |

### 2.3 Graph Theory (Combinatorial)

| Skill | Professional Application |
|-------|------------------------|
| Ramsey theory | **None (competition-only)** — mostly theoretical. Some distant application in communication complexity and theoretical CS. |
| Extremal graph theory (Turan's theorem) | **None (competition-only)** — research-level combinatorics and theoretical computer science. |
| Graph coloring | **Compiler design** — register allocation is a graph coloring problem. **Scheduling** — exam/meeting scheduling (conflict graphs). **Map rendering** — geographic map coloring. **Frequency assignment** — wireless channel allocation. |
| Trees, spanning trees, Cayley's formula | **Network design** — minimum cost network topologies. **Electrical engineering** — Kirchhoff's matrix tree theorem for circuit analysis. |
| Planarity, Euler's formula | **VLSI/PCB design** — checking if a circuit can be laid out without wire crossings. **Graph visualization** — planar layout algorithms. Niche. |
| Matchings, Hall's theorem | **Job/resource allocation** — assigning workers to tasks (the mathematical foundation). **Kidney exchange programs** — matching donors to recipients. **Online advertising** — ad-slot assignment. |
| Hamiltonian / Eulerian paths | **Logistics** — route optimization (TSP variants). **DNA sequencing** — de Bruijn sequences for genome assembly use Eulerian paths. **Snow plow routing** — Chinese postman problem (Eulerian). |
| Degree sequences | **Network science** — characterizing real-world networks (social, biological). Configuration model for random graph generation. Niche. |

### 2.4 Combinatorial Geometry

| Skill | Professional Application |
|-------|------------------------|
| Convex hull arguments, Helly's theorem | **Machine learning** — SVM (support vector machines) rely on convex hull separation. **Computational geometry** — collision detection in games/robotics. |
| Pick's theorem, lattice point problems | **Computer graphics** — rasterization (deciding which pixels a polygon covers). Niche. |
| Coloring of the plane, tiling | **None (competition-only)** — beautiful math, but no direct industrial application. |
| Extremal problems in discrete geometry | **None (competition-only)** |

### 2.5 Extremal & Probabilistic Combinatorics

| Skill | Professional Application |
|-------|------------------------|
| Sperner's theorem, Dilworth's theorem | **Database theory** — antichain/chain decomposition concepts appear in query optimization and concurrency control. Dilworth's theorem relates to minimum path cover in DAGs (scheduling). Niche. |
| Probabilistic method (LLL, moment methods) | **Randomized algorithm design** — proving existence of efficient hash families, expander graphs for networking, derandomization. **Theoretical CS research.** |
| Turán-type problems | **None (competition-only)** |
| Partially ordered sets (posets) | **Version control** — Git's commit DAG is a poset. **Distributed systems** — vector clocks, causal ordering of events (Lamport clocks). **Task scheduling** — dependency resolution. |

### 2.6 Invariants and Monovariants

| Skill | Professional Application |
|-------|------------------------|
| Parity arguments | **Error detection** — parity bits in data transmission. **Formal verification** — proving program properties. Simple but ubiquitous concept. |
| Coloring invariants | **None (competition-only)** |
| Monovariants | **Formal verification** — loop variants (proving termination of programs). **Distributed systems** — proving convergence of self-stabilizing algorithms. **Control theory** — Lyapunov functions are continuous-world monovariants. |
| Weight/potential function arguments | **Amortized analysis** — analyzing data structure operations (e.g., splay tree amortized cost). **Algorithm design** — potential method is a standard technique in algorithm analysis textbooks and real system performance analysis. |
| Extremal principle | **None (competition-only)** — a proof/problem-solving strategy. |

### 2.7 Combinatorial Games & Strategies

| Skill | Professional Application |
|-------|------------------------|
| Nim, Sprague-Grundy theory | **Game AI development** — directly applicable in designing AI for combinatorial games. **Mechanism design** — auction theory and strategic interaction modeling. Niche. |
| Strategy stealing | **None (competition-only)** — elegant proof technique but not applied industrially. |
| Pairing strategies | **None (competition-only)** |
| Games on graphs | **AI/Reinforcement learning** — game-tree search (Alpha-Beta, MCTS). **Economics** — modeling strategic behavior on networks. |

### 2.8 Algorithms and Processes

| Skill | Professional Application |
|-------|------------------------|
| Greedy algorithms in proofs | **Everywhere in software engineering** — the greedy mindset is one of the most practical algorithmic paradigms (scheduling, caching, Huffman coding, network routing). |
| Inductive constructions | **Formal verification** — structural induction for proving correctness of recursive programs. **Compiler design** — inductive definitions of language grammars. |
| Iterative processes and convergence | **ML/AI** — convergence analysis of gradient descent, EM algorithm, iterative solvers. **Scientific computing** — convergence of PDE solvers. |

---

## 3. Geometry

### 3.1 Euclidean Geometry — Triangle Theory

| Skill | Professional Application |
|-------|------------------------|
| Cevians, medians, altitudes, angle bisectors | **None (competition-only)** — classic Euclidean geometry has very little direct professional application. |
| Notable points (centroid, orthocenter, circumcenter, incenter) | **Centroid** is used in **computational geometry, physics simulations** (center of mass), and **clustering** (k-means centroid). The rest are competition-only. |
| Nine-point circle, Euler line, Simson line | **None (competition-only)** |
| Stewart's theorem, angle bisector theorem, mass points | **None (competition-only)** — mass point geometry has a cute connection to barycentric coordinates used in graphics, but the technique itself is not used. |
| Ceva's theorem, Menelaus' theorem | **None (competition-only)** |

### 3.2 Circle Geometry

| Skill | Professional Application |
|-------|------------------------|
| Power of a point, radical axes | **None (competition-only)** |
| Cyclic quadrilaterals, Ptolemy's theorem | **None (competition-only)** |
| Tangent lines and relations | **CAD software (minor)** — tangent calculations in curve-fitting and tool-path generation. Very basic tangent math is used, not the olympiad theory. |
| Mixtilinear incircles | **None (competition-only)** |
| Apollonius circle, coaxial circles | **None (competition-only)** |

### 3.3 Projective Geometry

| Skill | Professional Application |
|-------|------------------------|
| Cross-ratio, harmonic conjugates | **Computer vision** — cross-ratio is invariant under perspective projection; used in camera calibration and 3D reconstruction from 2D images. |
| Poles and polars | **None (competition-only)** — some use in classical optics but not modern industry. |
| Harmonic division, bundles | **None (competition-only)** |
| Projective transformations | **Computer vision & AR/VR** — homography estimation (mapping one plane to another), image stitching (panoramas), augmented reality marker tracking. This is heavily used. |
| Desargues', Pascal's, Brianchon's theorems | **None (competition-only)** |

### 3.4 Inversive Geometry

| Skill | Professional Application |
|-------|------------------------|
| Circle inversion and all related | **None (competition-only)** — some application in conformal mapping (2D fluid dynamics, electrostatics), but this is niche physics/engineering research, not mainstream professional work. |

### 3.5 Transformations

| Skill | Professional Application |
|-------|------------------------|
| Isometries (translations, rotations, reflections) | **Computer graphics & game engines** — every frame in every video game applies these. **Robotics** — rigid body transformations for robot arm kinematics. **CAD** — part placement and assembly modeling. Core skill. |
| Similarities (spiral similarities) | **Image processing** — image registration, template matching with scale+rotation invariance. |
| Homothety | **Computer graphics** — scaling operations. **Cartography** — map scaling. Basic but universal. |
| Composition of transformations | **Computer graphics** — transformation matrices, scene graph traversal (parent-child transforms). **Robotics** — forward/inverse kinematics chains. Critical skill. |
| Applications to concyclicity/collinearity | **None (competition-only)** |

### 3.6 Trigonometric Methods

| Skill | Professional Application |
|-------|------------------------|
| Law of sines, law of cosines | **Surveying & GPS** — triangulation and trilateration for positioning. **Game development** — distance and angle calculations. **Structural engineering** — force decomposition in trusses. |
| Area formulas (Heron's, trig area) | **GIS (Geographic Information Systems)** — land area computation. **Computer graphics** — triangle rasterization. |
| Trigonometric identities | **Signal processing** — Fourier analysis, modulation/demodulation in telecommunications. **Control systems** — Bode plots, frequency response analysis. **Physics simulations** — wave mechanics, harmonic motion. Core skill in many engineering fields. |
| Trig forms of Ceva and Menelaus | **None (competition-only)** |

### 3.7 Coordinate & Analytic Methods

| Skill | Professional Application |
|-------|------------------------|
| Cartesian coordinates | **Universal** — every field that involves spatial computation: game dev, robotics, GIS, CAD, physics simulations, data visualization. |
| Complex numbers (geometric) | **Electrical engineering** — AC circuit analysis (impedance as complex numbers). **Signal processing** — FFT, frequency domain analysis. **Control theory** — root locus, Nyquist plots. **Quantum computing** — quantum states are complex vectors. |
| Barycentric coordinates | **Computer graphics** — triangle rasterization in GPU pipelines (every pixel on screen is computed using barycentric interpolation). **Finite element method** — basis functions on triangular meshes. Important and directly used. |
| Vectors (dot product, cross product) | **Universal in engineering** — physics simulations, game engines, robotics, computer graphics, structural analysis. One of the most broadly applicable math skills. |

### 3.8 Geometric Inequalities

| Skill | Professional Application |
|-------|------------------------|
| Triangle inequality | **Metric spaces** — foundational axiom in database indexing (metric trees, VP-trees for nearest neighbor search). **Networking** — triangle inequality violations in network coordinates (Vivaldi system). |
| Erdos-Mordell inequality | **None (competition-only)** |
| Isoperimetric inequality | **Materials science / architecture (niche)** — optimizing surface-area-to-volume ratios. **ML theory** — isoperimetric inequalities on graphs for understanding graph partitioning (Cheeger's inequality). |
| Geometric optimization (Fermat point) | **Facility location** — the Fermat-Weber problem is a classic operations research problem (where to place a warehouse to minimize total transport cost). |

---

## 4. Number Theory

### 4.1 Divisibility & Primes

| Skill | Professional Application |
|-------|------------------------|
| GCD, LCM, Euclidean algorithm | **Cryptography** — RSA key generation uses extended GCD. **Systems programming** — scheduling (LCM of periods), fraction simplification. **Music technology** — rhythm LCM for polyrhythm computation. |
| Fundamental Theorem of Arithmetic | **Cryptography** — the entire RSA cryptosystem is built on the difficulty of factoring. **Database design** — normalization theory has conceptual parallels. |
| Bezout's identity, linear Diophantine equations | **Cryptography** — computing modular inverses for RSA and ECC. **Scheduling** — solving integer constraints in resource allocation. |
| Sieve of Eratosthenes | **Cryptography** — generating large primes. **Number-theoretic libraries** — fundamental building block. |
| Lifting the exponent lemma | **None (competition-only)** |

### 4.2 Modular Arithmetic

| Skill | Professional Application |
|-------|------------------------|
| Congruences, residue classes | **Cryptography** — the entire field. **Hashing** — hash functions are modular arithmetic. **Checksums** — ISBN, credit card validation (Luhn), CRC. **Databases** — hash partitioning, consistent hashing. |
| Fermat's Little Theorem, Euler's theorem | **Cryptography** — RSA encryption/decryption directly uses Euler's theorem. **Performance** — efficient modular exponentiation. |
| Chinese Remainder Theorem | **Cryptography** — RSA-CRT (4x speedup of RSA decryption in practice). **Distributed systems** — secret sharing. **Signal processing** — residue number systems for parallel computation. |
| Wilson's theorem | **None (competition-only)** — occasionally used in primality proofs but not in practice. |
| Modular inverses | **Cryptography** — key computation in RSA, ECC, Diffie-Hellman. **Competitive programming** — common, but also used in any system doing modular arithmetic. |
| Order, primitive roots | **Cryptography** — Diffie-Hellman key exchange requires generators (primitive roots) of cyclic groups. **Random number generation** — full-period linear congruential generators use primitive roots. |

### 4.3 p-adic Valuations

| Skill | Professional Application |
|-------|------------------------|
| All p-adic topics | **None (competition-only)** — p-adic analysis has applications in theoretical physics (p-adic string theory) and pure mathematics research, but no mainstream professional application. |

### 4.4 Quadratic Residues

| Skill | Professional Application |
|-------|------------------------|
| Legendre symbol, Euler's criterion | **Cryptography** — Solovay-Strassen primality test, Tonelli-Shanks square root algorithm (used in ECC point decompression). |
| Quadratic reciprocity | **Cryptography research** — understanding structure of quadratic residues for cryptanalysis. Niche. |
| Jacobi symbol | **Cryptography** — Solovay-Strassen test, some factoring algorithms. |
| QR modulo prime powers | **None (competition-only)** — too specialized. |

### 4.5 Diophantine Equations

| Skill | Professional Application |
|-------|------------------------|
| Pythagorean triples | **Computer graphics (minor)** — generating rational points on circles. Very niche. |
| Pell's equation, continued fractions | **Continued fractions** are used in **numerical methods** (best rational approximations — used in calendar systems, gear ratio design). Pell's equation itself is competition-only. |
| Infinite descent | **None (competition-only)** — a proof technique. |
| Vieta jumping | **None (competition-only)** |
| Zsygmondy's theorem | **None (competition-only)** |
| Sophie Germain identity | **None (competition-only)** |
| Sum of squares theorems | **None (competition-only)** — beautiful number theory with no direct application. |

### 4.6 Arithmetic Functions

| Skill | Professional Application |
|-------|------------------------|
| Euler's totient phi(n) | **Cryptography** — RSA directly: d = e^(-1) mod phi(n). Computing phi(n) quickly is equivalent to factoring n. |
| Divisor functions d(n), sigma(n) | **None in mainstream** — some use in algorithm complexity analysis and number-theoretic computation libraries. |
| Mobius function, Mobius inversion | **Combinatorics in CS** — inclusion-exclusion over lattices, used in some counting algorithms. **Number-theoretic transforms** — computing multiplicative functions efficiently. Niche. |
| Dirichlet convolution | **None (competition-only)** — theoretical number theory. |

### 4.7 Advanced Topics

| Skill | Professional Application |
|-------|------------------------|
| Hensel's lemma | **Cryptography (niche)** — lifting solutions in p-adic analysis for some factoring algorithms. |
| Cyclotomic polynomials | **Cryptography** — cyclotomic fields are used in some lattice-based (post-quantum) cryptographic schemes. |
| Algebraic integers (Gaussian, Eisenstein) | **Cryptography research** — some advanced number-theoretic algorithms use algebraic integers. Niche. |
| Bertrand's postulate, prime gaps | **None (competition-only)** |
| Legendre's conjecture | **None (competition-only)** |

---

## 5. Cross-Cutting Skills

### 5.1 Problem-Solving Strategies

| Skill | Professional Application |
|-------|------------------------|
| Working backwards, extreme cases | **Debugging** — root cause analysis, working backwards from a failure. **Systems design** — designing from the desired output backwards. **Universal professional skill.** |
| Small case analysis, pattern recognition | **Data science** — exploratory data analysis. **Machine learning** — feature engineering through pattern recognition. **Software engineering** — testing with small inputs first. Universal. |
| Generalization and specialization | **Software architecture** — abstraction, designing general interfaces from specific requirements. **Mathematical modeling** — building models by specializing general frameworks. Universal. |
| Symmetry exploitation | **Algorithm optimization** — eliminating redundant computation. **Physics simulations** — exploiting symmetry to reduce problem size. **ML** — weight sharing in CNNs exploits translational symmetry. |
| Wishful thinking and construction | **Product engineering** — prototype-first development. **API design** — "design the interface you wish you had." |

### 5.2 Proof Techniques

| Skill | Professional Application |
|-------|------------------------|
| Direct proof, contradiction, contrapositive | **Formal verification** — proving correctness of critical systems (aerospace, medical devices, smart contracts). **Protocol design** — proving security properties. |
| Mathematical induction | **Formal verification** — proving properties of recursive/iterative programs, especially in theorem provers (Coq, Isabelle, Lean). **Algorithm correctness** — the standard technique for proving loop invariants. |
| Infinite descent | **None (competition-only)** as a standalone technique. |
| Constructive vs. existential proofs | **Software engineering** — constructive proofs correspond to algorithms (Curry-Howard correspondence). Important in type theory and functional programming (Haskell, Rust type system). |
| Well-ordering principle | **Formal verification** — proving termination of programs. |

### 5.3 Writing & Communication

| Skill | Professional Application |
|-------|------------------------|
| Rigorous mathematical writing | **Technical documentation** — writing specs, RFCs, design docs. **Research** — academic papers, patent applications. **ML research** — theorem-proof style is standard in ML theory papers. |
| Logical flow, edge case handling | **Code review** — reviewing logic for correctness. **QA engineering** — systematic test case design. Universal engineering skill. |

---
---

# PART II — Competitive Programming Skills

---

## 1. Data Structures

### 1.1 Foundational Structures

| Skill | Professional Application |
|-------|------------------------|
| Arrays, dynamic arrays | **Universal** — the most fundamental data structure in all of software engineering. |
| Linked lists | **Operating systems** — kernel data structures (Linux kernel uses linked lists extensively). **Memory allocators** — free-list management. **Undo systems** — editor undo/redo chains. Used less in application-level code due to cache inefficiency. |
| Stacks (monotonic stack) | **Compilers** — parsing expressions, call stack management. **Web browsers** — navigation history. **Monotonic stack** specifically is niche (mainly competition), but stack-based algorithms appear in histogram problems in image processing. |
| Queues, deques, monotonic deque | **Operating systems** — process scheduling, BFS in network packet processing. **Streaming systems** — message queues (Kafka, RabbitMQ use queue abstractions). **Sliding window max/min** — real-time signal processing, financial data (rolling max/min of stock prices). |
| Priority queues / heaps | **Operating systems** — process scheduling (priority schedulers). **Networking** — packet prioritization, Dijkstra in routing protocols (OSPF). **Databases** — top-K queries, merge step in external sort. **Event simulation** — discrete event simulation engines. Very common. |
| Hash tables | **Universal** — dictionaries/maps are the most-used non-trivial data structure in software. Databases (hash indices), caches (memcached, Redis), compilers (symbol tables), networking (routing tables), deduplication. |

### 1.2 Tree-Based Structures

| Skill | Professional Application |
|-------|------------------------|
| Binary Search Trees | **Databases** — the concept underlies all ordered index structures. **Standard libraries** — Java TreeMap, C++ std::map. |
| Self-balancing BSTs (AVL, Red-Black) | **Language runtimes** — Red-Black trees power Java's TreeMap, C++ std::map/set, Linux CFS scheduler. **Databases** — in-memory indices. **Operating systems** — virtual memory management (Linux uses RB-trees for VMAs). |
| Treaps (implicit treap) | **Niche** — used in some specialized systems (e.g., some persistent data structure libraries). Implicit treap concepts appear in rope data structures for **text editors** (e.g., Xi editor). |
| Splay trees | **Caching systems** — splay trees naturally bring frequently accessed items to the root (self-adjusting). **Networking** — some implementations of link-cut trees for dynamic graph algorithms. Niche in practice. |
| Order-statistic trees | **Databases** — rank/select queries (finding the k-th smallest element, PERCENTILE functions). **Finance** — order-book implementations in trading systems. |

### 1.3 Segment Trees

| Skill | Professional Application |
|-------|------------------------|
| Basic segment tree | **Databases** — range aggregate queries (SUM, MIN, MAX over ranges). The concept is embedded in database query engines, though typically as B-tree variants rather than literal segment trees. **Gaming** — range queries for collision detection. |
| Lazy propagation | **Niche** — not commonly implemented directly in production, but the concept of deferred computation appears in **lazy evaluation** (Haskell, Spark), **write-back caches**, and **batched updates** in databases. |
| Persistent segment trees | **Version control systems** — conceptually related to immutable/persistent data (Git's content-addressable storage). **Functional programming** — persistent data structures in Clojure, Haskell. **Databases** — MVCC (Multi-Version Concurrency Control) uses similar ideas. |
| 2D segment trees, segment tree beats, merge sort tree | **None (competition-only)** — too specialized for production use. |
| Dynamic/implicit segment trees | **Niche** — sparse range-query scenarios in some analytics systems, but in practice people use different approaches (e.g., skip lists, B-trees). |

### 1.4 Binary Indexed Trees (Fenwick Trees)

| Skill | Professional Application |
|-------|------------------------|
| Point update + prefix sum | **Statistics** — computing cumulative frequency distributions efficiently. **Inversion counting** — used in some ranking and recommendation algorithms. Simpler to implement than segment trees for prefix-sum-type problems. |
| Range update + range query | **Niche** — same domain as segment trees but less commonly seen in production. |
| 2D Fenwick trees | **Niche** — 2D cumulative statistics in some scientific computing or image processing applications. |

### 1.5 Other Advanced Structures

| Skill | Professional Application |
|-------|------------------------|
| Tries | **Autocomplete & search engines** — prefix-based suggestion systems (Google search, IDE autocomplete). **Networking** — IP routing tables (longest prefix match). **Spell checkers** — dictionary lookup. **NLP** — morphological analysis. Very practical. |
| Disjoint Set Union (DSU) | **Image processing** — connected component labeling. **Network connectivity** — determining if nodes are connected (Kruskal's in network design). **Social networks** — community detection. **Compilers** — type unification in type inference. Widely used. |
| Sparse tables | **Database engines** — fast RMQ (range minimum query) for certain index structures. **Bioinformatics** — LCA queries on suffix trees for genome analysis. |
| Sqrt decomposition, Mo's algorithm | **None (competition-only)** — too specialized; real systems use different approaches for range queries. |
| Euler tour, HLD, centroid decomposition | **Niche** — Euler tour flattening is used in some database systems for XML/JSON querying (nested set model). HLD and centroid decomposition are competition-only in practice. |
| Link-Cut Trees | **Research** — dynamic connectivity algorithms. Not used in mainstream production code. |
| Li Chao trees | **None (competition-only)** |

---

## 2. Graph Algorithms

### 2.1 Graph Representations & Traversals

| Skill | Professional Application |
|-------|------------------------|
| Adjacency list/matrix/edge list | **Universal in graph-based systems** — social networks (Facebook/LinkedIn graph), knowledge graphs (Google), dependency graphs (package managers like npm, Maven). |
| BFS | **Web crawlers** — Googlebot uses BFS-like strategies. **Social networks** — friend-of-friend suggestions, six-degrees-of-separation queries. **Networking** — broadcast protocols. **Game AI** — pathfinding in unweighted grids. |
| DFS | **Compilers** — control flow analysis, dead code detection. **File systems** — recursive directory traversal. **Web crawlers** — deep crawling. **Garbage collectors** — reachability analysis (mark phase of mark-and-sweep). |
| DFS timestamps | **Compilers** — dominance analysis, loop detection in CFGs. **Databases** — nested set model for hierarchical data. |

### 2.2 Shortest Paths

| Skill | Professional Application |
|-------|------------------------|
| Dijkstra's algorithm | **GPS / Maps** — Google Maps, Waze, Apple Maps (with heuristic enhancements like A*). **Network routing** — OSPF protocol uses Dijkstra. **Robotics** — path planning. **Logistics** — vehicle routing. One of the most practically important algorithms. |
| Bellman-Ford / SPFA | **Networking** — RIP (Routing Information Protocol) uses Bellman-Ford. **Finance** — detecting negative-weight cycles corresponds to arbitrage detection in currency exchange markets. |
| Floyd-Warshall | **Network analysis** — all-pairs shortest paths for small-to-medium networks (pre-computation for routing). **Transitive closure** — computing reachability in dependency graphs. |
| Johnson's algorithm | **Niche** — used in sparse all-pairs scenarios; most practical systems use Dijkstra from each source or specialized techniques. |
| 0-1 BFS | **Niche** — useful in some grid-based game AI and robotics path planning with free/costly movements. |

### 2.3 Minimum Spanning Trees

| Skill | Professional Application |
|-------|------------------------|
| Kruskal's / Prim's / Boruvka's | **Network design** — designing minimum-cost telecommunication networks, power grids, road networks. **Clustering** — single-linkage clustering (cutting MST edges). **Image segmentation** — Felzenszwalb's algorithm uses MST. **Approximation algorithms** — MST-based 2-approximation for metric TSP. |
| MST properties (cycle, cut) | **Network engineering** — understanding redundancy and backup paths. |
| Minimum spanning arborescence | **Niche** — dependency resolution in directed networks, some NLP (Eisner's algorithm for dependency parsing uses similar ideas). |
| Second-best MST, bottleneck ST | **Network reliability** — finding backup topologies. **Bottleneck path** problems in logistics. |

### 2.4 Tree Algorithms

| Skill | Professional Application |
|-------|------------------------|
| LCA (Lowest Common Ancestor) | **Version control** — Git merge-base computation (finding common ancestor of branches). **Databases** — XML/JSON path queries in hierarchical data. **Bioinformatics** — phylogenetic tree analysis (finding common ancestor of species). **File systems** — computing relative paths. |
| Tree diameter | **Network design** — determining the worst-case latency in a network tree topology. |
| Rerooting technique | **Niche** — some network optimization and tree-based analytics. Mostly competition-only. |
| HLD, centroid decomposition | **None (competition-only)** in mainstream. Some research use in algorithms for dynamic graphs. |
| Virtual tree / auxiliary tree | **None (competition-only)** |
| Tree isomorphism, hashing | **Cheminformatics** — molecular graph comparison. **Compiler optimization** — CSE (common subexpression elimination) compares expression trees. |
| Prufer sequence | **None (competition-only)** — a theoretical bijection. |

### 2.5 Connectivity & Components

| Skill | Professional Application |
|-------|------------------------|
| Strongly Connected Components | **Web analysis** — structure of the web graph (bow-tie model). **Compiler optimization** — finding loops in control flow graphs. **Dependency analysis** — detecting circular dependencies in software modules. |
| Bridges and articulation points | **Network reliability** — identifying single points of failure in networks. **Infrastructure engineering** — critical link/node analysis. **Social networks** — identifying broker nodes. Directly applicable. |
| Biconnected components | **Network design** — ensuring 2-connectivity (fault tolerance). Similar to bridges analysis. |
| Euler path/circuit | **DNA sequencing** — de Bruijn graph assembly. **Circuit design** — PCB trace routing. **Snow plow / postal routing** — Chinese postman problem. |
| Topological sorting | **Build systems** — Make, Gradle, Bazel (dependency-ordered compilation). **Package managers** — npm, pip (install order). **Spreadsheets** — cell dependency evaluation order. **CI/CD pipelines** — job ordering. Very common and practical. |
| Condensation (DAG of SCCs) | **Compiler optimization** — simplifying control flow for optimization passes. **Dependency analysis** — collapsing circular dependencies into single units. |

### 2.6 Network Flow

| Skill | Professional Application |
|-------|------------------------|
| Max-flow / Min-cut | **Image segmentation** — graph-cut-based image segmentation (GrabCut in OpenCV). **Network capacity planning** — determining maximum throughput. **Supply chain** — maximum goods flow through a distribution network. |
| Flow algorithms (Dinic's, Push-relabel) | **Operations research** — transportation and assignment problems in logistics. **Telecommunications** — bandwidth allocation. **Airlines** — crew scheduling. Dinic's is the go-to in practice for max-flow. |
| Min-cost max-flow | **Supply chain optimization** — minimizing transportation cost while meeting demand. **Assignment problems** — optimal task-worker assignment with costs. **Energy markets** — optimal power flow. |
| Hungarian algorithm | **Computer vision** — object tracking (matching detections across frames, e.g., SORT tracker). **Resource allocation** — optimal assignment of tasks to machines. **Ride-sharing** — matching riders to drivers. Very practical. |
| Circulation with demands | **Supply chain** — modeling networks with minimum throughput requirements. Niche but important in OR. |

### 2.7 Matching

| Skill | Professional Application |
|-------|------------------------|
| Bipartite matching (Hopcroft-Karp, Kuhn's) | **Job platforms** — matching candidates to positions. **Kidney exchange** — donor-recipient matching (Nobel Prize-winning application). **Advertising** — ad-to-slot matching. **Ride-sharing** — rider-driver matching. |
| General matching (Blossom) | **Chemistry** — Kekulé structure enumeration. **Scheduling** — some non-bipartite assignment problems. Niche. |
| Konig's theorem | **Operations research** — minimum vertex cover problems in bipartite networks (e.g., minimum guards to cover all corridors). |
| Stable matching (Gale-Shapley) | **Medical residency** — NRMP (National Resident Matching Program) uses Gale-Shapley. **School choice** — student-school assignment. **Kidney exchange**. Nobel Prize in Economics (2012) for this work. |
| Weighted bipartite matching | **Logistics** — optimal assignment with costs. **Computer vision** — point set registration. |

### 2.8 Special Graph Problems

| Skill | Professional Application |
|-------|------------------------|
| 2-SAT | **Configuration management** — checking satisfiability of boolean constraints (package dependency resolution with conflicts). **EDA (Electronic Design Automation)** — VLSI placement constraints. |
| Graph coloring | **Compiler design** — register allocation (Chaitin's algorithm). **Scheduling** — exam timetabling, frequency assignment in wireless networks. |
| Dominator tree | **Compilers** — SSA (Static Single Assignment) form construction uses dominance frontiers. Critical in modern compiler infrastructure (LLVM, GCC). |
| Functional graphs, cycle detection | **Cryptography** — Pollard's rho algorithm for factoring. **Distributed systems** — deadlock detection. **Linked list problems** — detecting cycles in data structures (Floyd's). |
| Cactus graphs, chordal graphs, interval graphs | **Niche** — interval graphs appear in scheduling theory and bioinformatics (interval scheduling). Chordal graphs in sparse matrix computation (fill-reducing orderings for Cholesky decomposition). |

---

## 3. Dynamic Programming

### 3.1 Classical DP

| Skill | Professional Application |
|-------|------------------------|
| 1D DP (Fibonacci, coin change, etc.) | **Finance** — optimal coin/denomination systems, pricing strategies. **Resource allocation** — budget distribution. Foundation for all DP-based professional applications. |
| 2D DP (LCS, edit distance) | **Bioinformatics** — DNA/protein sequence alignment (BLAST, Smith-Waterman). **NLP** — spell checkers, diff tools (Unix diff). **Version control** — Git diff uses LCS variants. Very widely used. |
| Knapsack | **Finance** — portfolio optimization (selecting assets with budget constraint). **Logistics** — container loading, cutting stock problem. **Cloud computing** — VM placement (bin packing variants). **Advertising** — budget allocation across campaigns. |
| LIS (Longest Increasing Subsequence) | **Statistics** — Patience sorting, longest monotone subsequence in data analysis. **Bioinformatics** — gene order comparison. Niche. |
| Interval DP | **Compiler optimization** — optimal expression evaluation order. **Operations research** — optimal BST construction for search engines (given access frequencies). |
| Counting DP | **Combinatorics in production** — computing number of valid configurations, test coverage counting. **Probability** — computing exact probabilities via counting. |

### 3.2 Bitmask DP

| Skill | Professional Application |
|-------|------------------------|
| TSP via bitmask DP | **Logistics** — exact TSP solutions for small instances (< 20 cities). **Robotics** — optimal tour planning for warehouse robots. **PCB design** — drill path optimization. In practice, used for small n and as a subroutine in branch-and-bound. |
| Assignment problems | **Operations research** — job-shop scheduling for small instances. **Hardware** — FPGA configuration optimization. |
| Subset enumeration, broken profile | **None (competition-only)** for the most part — too specialized. |

### 3.3 DP on Trees

| Skill | Professional Application |
|-------|------------------------|
| Subtree DP | **Network design** — optimal facility placement on tree networks. **Compiler optimization** — optimal code generation for expression trees. **Organizational analytics** — metrics aggregation on org hierarchy. |
| Rerooting DP | **Niche** — network centroid finding, some tree-based analytics. |
| Independent set / vertex cover on trees | **Network security** — optimal sensor/monitor placement on tree-structured networks. |

### 3.4 Digit DP

| Skill | Professional Application |
|-------|------------------------|
| All digit DP topics | **Niche** — counting numbers with specific properties (useful in some analytics and number-theoretic programming). Not a standard professional tool. |

### 3.5 DP on DAGs

| Skill | Professional Application |
|-------|------------------------|
| Longest/shortest path on DAGs | **Project management** — critical path method (CPM) for finding project duration = longest path in DAG. **Build systems** — determining critical build path. Very practical. |
| Counting paths | **Reliability engineering** — counting failure paths in fault trees. **Web analytics** — counting user journey paths. |

### 3.6 SOS DP

| Skill | Professional Application |
|-------|------------------------|
| Subset sum transforms | **None (competition-only)** in mainstream. Some use in theoretical CS research and specialized combinatorial algorithms. |

### 3.7 DP Optimizations

| Skill | Professional Application |
|-------|------------------------|
| Convex Hull Trick | **Finance** — optimizing piecewise-linear cost functions. **Operations research** — convex cost flow. Niche but real. |
| Knuth's, D&C, SMAWK optimizations | **None (competition-only)** — too specialized. The underlying ideas (exploiting monotonicity, convexity) are general optimization principles. |
| Aliens trick / WQS binary search | **None (competition-only)** |
| Slope trick | **None (competition-only)** |
| Hirschberg's algorithm | **Bioinformatics** — space-efficient sequence alignment (important when aligning very long DNA sequences). Also used in diff tools. |

### 3.8 Other DP Types

| Skill | Professional Application |
|-------|------------------------|
| Probability / Expected value DP | **Finance** — option pricing (binomial tree model). **Reinforcement learning** — value iteration, policy iteration (Bellman equations are DP). **Queueing theory** — analyzing system performance. |
| DP with matrix exponentiation | **Population modeling** — Leslie matrix for age-structured populations. **Markov chains** — computing n-step transition probabilities. **Finance** — credit rating migration. |
| Broken profile / plug DP | **None (competition-only)** |
| DP on permutations | **None (competition-only)** |
| Connection profile DP (Steiner tree) | **Network design (niche)** — Steiner tree problem appears in VLSI routing and telecommunication network design. Exact DP is for small instances; heuristics are used at scale. |

---

## 4. String Algorithms

### 4.1 Hashing

| Skill | Professional Application |
|-------|------------------------|
| Polynomial hashing (Rabin-Karp) | **Plagiarism detection** — document fingerprinting (Turnitin-like systems use rolling hashes). **Data deduplication** — content-defined chunking (rsync, cloud storage dedup). **Malware detection** — signature matching. |
| Double hashing, anti-hash techniques | **Security** — hash collision resistance analysis, designing robust hash functions. |
| Palindrome/substring hashing | **Bioinformatics** — finding repeated motifs in DNA. **NLP** — text analysis. |

### 4.2 Pattern Matching

| Skill | Professional Application |
|-------|------------------------|
| KMP | **Text editors / IDEs** — find/replace functionality. **Log analysis** — searching patterns in log files (grep internals). **Network security** — intrusion detection systems (Snort, Suricata). |
| Z-function | **Same domains as KMP** — text search, though KMP is more commonly implemented in production. |
| Aho-Corasick | **Network security** — multi-pattern matching in IDS/IPS (scanning network traffic against thousands of malware signatures simultaneously). **Content filtering** — web filtering, spam detection. **Bioinformatics** — searching for multiple DNA motifs. Very practical. |
| Rabin-Karp | **Plagiarism detection** — document similarity (winnowing algorithm). **Data dedup** — rsync, Docker layer dedup. |

### 4.3 Suffix Structures

| Skill | Professional Application |
|-------|------------------------|
| Suffix array + LCP array | **Bioinformatics** — genome indexing (BWA aligner for DNA sequencing uses suffix arrays). **Search engines** — full-text indexing. **Data compression** — BWT-based compression (bzip2). |
| Suffix tree | **Bioinformatics** — genome assembly, repeat finding, multiple sequence alignment. **Text analytics** — document clustering. Less common than suffix arrays due to memory overhead. |
| Suffix automaton | **Niche** — substring statistics in some text analytics and bioinformatics tools. Mostly competition-only in practice. |

### 4.4 Palindromes

| Skill | Professional Application |
|-------|------------------------|
| Manacher's algorithm | **Bioinformatics** — DNA palindromes (restriction enzyme sites are palindromic). **NLP** — niche text analysis. |
| Palindromic tree (Eertree) | **None (competition-only)** |

### 4.5 Other String Techniques

| Skill | Professional Application |
|-------|------------------------|
| Trie-based algorithms | **Search engines** — autocomplete, prefix-based search. **Networking** — IP routing (longest prefix match in routers). **NLP** — dictionary operations, morphological analysis. Highly practical. |
| Burrows-Wheeler Transform, FM-index | **Data compression** — bzip2 uses BWT. **Bioinformatics** — BWA, Bowtie (DNA read aligners) use FM-index. One of the most impactful data structures in genomics. |
| Lyndon factorization, Booth's algorithm | **Niche** — lexicographically smallest rotation has application in data normalization (canonical forms for cyclic strings in databases). |
| Run enumeration | **None (competition-only)** |
| Bitap algorithm | **Text search** — Unix agrep tool for approximate matching. **NLP** — fuzzy string matching in search engines. |

---

## 5. Mathematics for CP

### 5.1 Number Theory

| Skill | Professional Application |
|-------|------------------------|
| Sieve of Eratosthenes | **Cryptography** — generating prime tables for key generation. **Number theory libraries** — foundational. |
| Modular arithmetic, exponentiation | **Cryptography** — RSA, Diffie-Hellman, ECC. **Blockchain** — hashing and digital signatures. **Checksums** — CRC, hash functions. Universal in security. |
| Extended Euclidean algorithm | **Cryptography** — computing modular inverses for RSA key generation. |
| CRT | **Cryptography** — RSA-CRT optimization. **Parallel computing** — residue number systems for hardware arithmetic. |
| Euler's totient, multiplicative functions | **Cryptography** — RSA. **Algorithm analysis** — counting coprime pairs. |
| Discrete logarithm (BSGS, Pohlig-Hellman) | **Cryptography** — breaking Diffie-Hellman (the hardness of discrete log is the security assumption). **Cryptanalysis** — security analysis. |
| Miller-Rabin primality test | **Cryptography** — standard primality test used in all major crypto libraries (OpenSSL, libsodium, Java BigInteger). Directly used in production. |
| Pollard's rho factorization | **Cryptography** — factoring algorithms. **Security auditing** — testing RSA key strength. |

### 5.2 Combinatorics

| Skill | Professional Application |
|-------|------------------------|
| Binomial coefficients, Lucas' theorem | **Statistics** — computing exact probabilities for hypothesis testing. **ML** — combinatorial probability calculations. |
| PIE | **Database analytics** — computing union cardinalities, deduplication. **A/B testing** — overlap analysis. **Reliability** — system failure probability. |
| Burnside / Polya enumeration | **Chemistry** — counting distinct molecular structures (isomers). **Computer graphics** — counting distinct tilings/patterns under symmetry. **Puzzle design** — enumerating distinct configurations. |
| Catalan, Stirling, Bell numbers | **CS theory** — analyzing number of BST shapes, parsing trees, set partitions. Mostly theoretical/educational. |
| Generating functions | **Probability** — PGFs and MGFs in statistics. **Algorithm analysis** — average-case analysis (Knuth-style). **Physics** — partition functions in statistical mechanics. |

### 5.3 Linear Algebra

| Skill | Professional Application |
|-------|------------------------|
| Gaussian elimination | **Scientific computing** — solving systems of equations (FEM, circuit analysis, structural mechanics). **ML** — least squares regression, PCA (via SVD). **Computer graphics** — ray-tracing (ray-plane intersections). One of the most practically important algorithms in all of computing. |
| Matrix exponentiation | **Markov chains** — computing steady-state probabilities (PageRank). **Population dynamics** — Leslie matrix models. **Finance** — credit rating transitions over time. |
| XOR linear basis | **Networking** — network coding, RAID parity calculations. **Cryptography** — analyzing linear properties of ciphers. |
| Systems over GF(2) | **Coding theory** — LDPC codes, turbo codes. **Cryptography** — cryptanalysis of stream ciphers. **Networking** — erasure coding in distributed storage (Reed-Solomon). |
| Kirchhoff's matrix tree theorem | **Network analysis** — counting spanning trees (network reliability). **Chemistry** — graph-theoretic molecular descriptors. Niche. |

### 5.4 Probability & Expected Value

| Skill | Professional Application |
|-------|------------------------|
| Linearity of expectation | **Universal in data science / ML** — expected value calculations, A/B test analysis. One of the most useful probability tools. |
| Geometric distribution, coupon collector | **Marketing** — "how many ads until a user converts?" **Gaming** — loot box probability analysis. **Distributed systems** — expected time to collect all responses. |
| Markov chains | **NLP** — language models (n-gram models), text generation. **Finance** — credit rating migration, stock price models. **Queueing theory** — system performance analysis. **Web** — PageRank is a Markov chain. |
| DP for expected values | **Reinforcement learning** — Bellman equations, value iteration, Q-learning. **Finance** — decision trees for investment analysis. Core to modern AI. |

### 5.5 Game Theory

| Skill | Professional Application |
|-------|------------------------|
| Sprague-Grundy, Nim | **Game AI** — solving combinatorial games exactly. Niche professional application. |
| Games on graphs/DAGs | **AI** — game tree search (chess, Go engines). **Economics** — mechanism design, auction theory (Google ad auctions). **Multi-agent systems** — strategic interaction modeling. |

### 5.6 Polynomials & Transforms

| Skill | Professional Application |
|-------|------------------------|
| FFT | **Signal processing** — audio processing (equalizers, compression, noise cancellation). **Telecommunications** — OFDM modulation (WiFi, 4G/5G). **Image processing** — frequency domain filtering. **Scientific computing** — solving PDEs. One of the top 10 most impactful algorithms ever. |
| NTT | **Cryptography** — polynomial multiplication in lattice-based (post-quantum) cryptography (CRYSTALS-Kyber uses NTT). **Big integer arithmetic** — used internally in arbitrary-precision math libraries. |
| Karatsuba multiplication | **Big number libraries** — GMP (GNU Multiple Precision) uses Karatsuba for medium-sized multiplications. **Cryptography** — big integer operations in RSA. |
| Polynomial division/interpolation | **Cryptography** — Shamir's secret sharing, Reed-Solomon error correction. **Signal processing** — filter design. **Finance** — yield curve interpolation. |
| Formal power series | **Symbolic computation** — computer algebra systems (Mathematica, Maple). **Combinatorics research** — analytic combinatorics. Niche. |
| Walsh-Hadamard Transform | **Coding theory** — decoding Reed-Muller codes. **ML** — fast kernel methods. **Quantum computing** — Hadamard gates. |
| Subset convolution | **None (competition-only)** |

---

## 6. Computational Geometry

### 6.1 Primitives

| Skill | Professional Application |
|-------|------------------------|
| Point, vector, line, segment representations | **Universal in spatial computing** — game engines (Unity, Unreal), CAD (AutoCAD, SolidWorks), GIS (ArcGIS, QGIS), robotics, AR/VR. |
| Cross product (orientation test) | **Computer graphics** — back-face culling, polygon winding order. **Robotics** — determining turn direction. **GIS** — spatial relationship queries. |
| Dot product | **Computer graphics** — lighting calculations (Lambert's cosine law), projection. **ML** — similarity measures, attention mechanisms. **Physics engines** — collision response. |
| Line/segment intersection | **CAD** — interference checking. **GIS** — road intersection detection, boundary analysis. **Game engines** — ray casting, collision detection. |
| Point-to-line/segment distance | **GIS** — proximity queries ("how far is this point from the road?"). **Robotics** — obstacle avoidance. |

### 6.2 Convex Hull

| Skill | Professional Application |
|-------|------------------------|
| Graham scan, Andrew's monotone chain | **Computer vision** — object bounding, shape analysis. **GIS** — territory delineation, minimum enclosing region. **Robotics** — obstacle representation. **Statistics** — data outlier detection (peeling convex hulls). |
| Dynamic convex hull | **Finance** — maintaining efficient frontier as assets change. **Computational geometry libraries** — CGAL. |
| Convex hull trick (DP) | **Finance (niche)** — optimizing over piecewise-linear functions. Mostly competition-only. |

### 6.3 Sweep Line & Plane Algorithms

| Skill | Professional Application |
|-------|------------------------|
| Line sweep (closest pair, rectangle union, segment intersection) | **VLSI/EDA** — design rule checking (detecting overlapping wires). **GIS** — map overlay operations. **Game engines** — broad-phase collision detection (sweep and prune). |
| Event-based processing | **Event-driven systems** — the general paradigm is foundational in UI frameworks, game loops, simulation engines. |
| Voronoi via sweep (Fortune's) | **GIS** — nearest facility analysis, territory planning. **Robotics** — path planning (Voronoi diagrams give maximum-clearance paths). **Materials science** — modeling crystal grain boundaries. |
| Angular sweep | **Robotics** — radar/lidar processing, visibility computation. |

### 6.4 Polygon Operations

| Skill | Professional Application |
|-------|------------------------|
| Shoelace formula | **GIS** — computing land parcel areas. **Surveying** — area computation from GPS coordinates. **Computer graphics** — polygon area for physics (mass/inertia). |
| Point-in-polygon | **GIS** — "which district is this address in?" (geo-fencing). **Game engines** — click detection, spatial queries. **Logistics** — delivery zone checking. Very common. |
| Minkowski sum | **Robotics** — configuration space computation (robot path planning around obstacles). **CAD** — offset curves and surfaces. **3D printing** — collision checking. |
| Polygon clipping (Sutherland-Hodgman) | **Computer graphics** — clipping polygons against view frustum (every GPU does this). **GIS** — map clipping operations. **CAD** — boolean operations on shapes. |
| Rotating calipers | **Computer vision** — minimum bounding box for object detection. **GIS** — optimal enclosing rectangle for spatial data. |

### 6.5 Advanced Geometry

| Skill | Professional Application |
|-------|------------------------|
| Half-plane intersection | **Robotics** — visibility regions, sensor coverage computation. **Linear programming** — feasible region visualization (intersection of half-planes = LP feasible region). |
| Voronoi diagram, Delaunay triangulation | **Mesh generation** — FEM (finite element method) in engineering simulation (Delaunay is the standard). **GIS** — spatial interpolation, territory planning. **Game development** — procedural terrain generation. **Meteorology** — weather station interpolation. Widely used. |
| Smallest enclosing circle (Welzl's) | **GIS** — finding minimum enclosing region for a set of locations. **Wireless networks** — optimal base station placement. |
| k-d trees | **ML** — k-nearest neighbor search. **Game engines** — spatial indexing, ray tracing acceleration. **Databases** — multidimensional range queries. **Robotics** — point cloud processing (PCL library). Very practical. |
| Pick's theorem | **Computer graphics** — counting lattice points in rasterization. Niche. |
| Great circle distance | **Aviation / shipping** — flight and shipping route optimization. **GIS** — accurate distance computation on Earth. **Ride-sharing** — distance calculations (Haversine formula). Directly used in production. |

---

## 7. General Algorithmic Techniques

### 7.1 Searching & Sorting

| Skill | Professional Application |
|-------|------------------------|
| Binary search | **Universal** — used everywhere: database indices (B-tree search), version control (git bisect), debugging (binary search on commits), API rate limiting, capacity planning ("what's the max load?"). Arguably the single most versatile algorithmic technique. |
| Ternary search | **ML** — hyperparameter optimization on unimodal functions. **A/B testing** — finding optimal threshold. Niche. |
| Sorting algorithms (merge sort, quickselect) | **Universal** — database ORDER BY, file systems, UI rendering (z-ordering). **Statistics** — quickselect for median finding. **Data processing** — every data pipeline sorts. |
| Coordinate compression | **Databases** — dictionary encoding in columnar stores (Parquet, ORC). **Data visualization** — mapping continuous data to discrete bins. |
| Parallel binary search | **None (competition-only)** as a specific technique, but the concept of parallelizing search appears in distributed systems. |

### 7.2 Greedy Algorithms

| Skill | Professional Application |
|-------|------------------------|
| Activity/interval scheduling | **Operating systems** — CPU scheduling (Shortest Job First). **Calendar apps** — meeting room booking. **Cloud computing** — VM scheduling. |
| Huffman coding | **Data compression** — DEFLATE (used in ZIP, gzip, PNG), JPEG (Huffman + DCT). One of the most used algorithms in practice (every compressed file uses it). |
| Fractional knapsack | **Finance** — portfolio allocation, resource budgeting. **Cloud** — fractional resource allocation. |
| Exchange argument | **Algorithm design** — proving greedy correctness. The mindset is valuable in any algorithm/system design role. |
| Matroid intersection | **Operations research (niche)** — theoretical foundation for greedy optimization. Academic. |

### 7.3 Divide and Conquer

| Skill | Professional Application |
|-------|------------------------|
| Merge sort applications | **Databases** — external sort (merge-sort based). **Analytics** — inversion counting for rank correlation (Kendall's tau). |
| D&C on segment tree / merge sort tree | **None (competition-only)** |
| CDQ divide and conquer | **None (competition-only)** |
| Centroid decomposition | **None (competition-only)** in production. |

### 7.4 Bit Manipulation

| Skill | Professional Application |
|-------|------------------------|
| Bitwise operations | **Systems programming** — flags, permissions (Unix file permissions), protocol headers. **Embedded systems** — register manipulation, hardware control. **Networking** — subnet masking, IP operations. **Databases** — bitmap indices. |
| XOR properties, XOR basis | **Networking** — RAID parity (XOR-based). **Cryptography** — stream cipher design. **Error detection** — simple parity checks. |
| Popcount, lowest set bit | **Databases** — bitmap index cardinality. **Chess engines** — bitboard move generation (Stockfish). **Compilers** — bit manipulation optimizations. |
| Iterating over submasks | **None (competition-only)** |
| Gray code | **Rotary encoders** — hardware position sensing. **Error reduction** — minimizing bit errors in ADC. **Genetic algorithms** — Gray code encoding for smoother fitness landscapes. |

### 7.5 Randomized Algorithms

| Skill | Professional Application |
|-------|------------------------|
| Randomized hashing | **Security** — hash DoS protection (SipHash in Python, Rust). **Databases** — randomized hash functions for load balancing. |
| Randomized pivot (quickselect/quicksort) | **Standard libraries** — virtually all production sorting algorithms use randomized pivoting (or its equivalent, introsort). **Databases** — approximate median/percentile. |
| Randomized geometry | **Computer graphics** — Monte Carlo ray tracing (every modern movie uses this). **Robotics** — RRT (Rapidly-exploring Random Trees) for path planning. |
| Birthday paradox | **Security** — collision resistance requirements for hash functions (explains why we need 256-bit hashes, not 128-bit). **Distributed systems** — UUID collision probability analysis. |
| Schwartz-Zippel lemma | **Cryptography** — zero-knowledge proofs (used in blockchain/ZK-SNARKs). **Verification** — polynomial identity testing in symbolic computation. |
| Treaps | **Persistent data structures** — some functional language runtimes. Niche in production vs. Red-Black trees. |

### 7.6 Offline & Online Techniques

| Skill | Professional Application |
|-------|------------------------|
| Mo's algorithm | **None (competition-only)** |
| Sqrt decomposition | **None (competition-only)** — the blocking idea appears in database systems (block-based I/O) but as a different concept. |
| Persistent data structures | **Functional programming** — Clojure's persistent vectors, Haskell's immutable data. **Databases** — MVCC (PostgreSQL, Oracle). **Version control** — Git's content-addressed storage. **Blockchain** — immutable state history. |
| Fractional cascading | **Computational geometry libraries** — speeding up multi-level queries. Niche. |
| CDQ, rollback | **Databases** — transaction rollback (UNDO logs). The general rollback concept is universal; CDQ is competition-only. |

### 7.7 Miscellaneous Techniques

| Skill | Professional Application |
|-------|------------------------|
| Meet in the middle | **Cryptanalysis** — MITM attacks on double encryption (reason why Double-DES is insecure). **Security** — brute-force optimization for key recovery. |
| Two pointers / sliding window | **Streaming analytics** — real-time aggregation (rolling average, moving sum). **Networking** — TCP sliding window protocol. **Data science** — time-series analysis (rolling statistics). Very common. |
| Small-to-large merging | **Databases** — merge strategies in LSM trees (LevelDB, RocksDB, Cassandra). The concept of merging smaller into larger is directly used. |
| Virtual tree construction | **None (competition-only)** |
| Interactive problems | **API design** — request-response protocol design, adaptive systems. The mindset of "querying an oracle" is relevant to API integration and testing. |
| Constructive algorithms | **Software engineering** — building valid configurations (e.g., generating valid test data, constructing schedules that satisfy constraints). |
| Output-sensitive algorithms | **Databases** — query planning considers output size. **Graphics** — rendering only visible objects (frustum culling). |

---

## 8. Implementation & Practical Skills

### 8.1 Language Mastery

| Skill | Professional Application |
|-------|------------------------|
| C++ STL mastery | **Systems programming** — high-performance servers, game engines, embedded systems. **Finance** — HFT (high-frequency trading) systems. **ML frameworks** — TensorFlow, PyTorch internals are C++. |
| Performance optimization (fast I/O, cache-friendly code) | **Systems engineering** — database engine development, game engine optimization, embedded systems. **HFT** — every nanosecond matters. **ML infrastructure** — training pipeline optimization. |
| Precision handling | **Scientific computing** — numerical stability in simulations. **Finance** — decimal precision for monetary calculations. **Graphics** — floating-point precision in coordinate systems. |

### 8.2 Debugging & Testing

| Skill | Professional Application |
|-------|------------------------|
| Stress testing, random test generation | **Quality assurance** — fuzz testing (AFL, libFuzzer), property-based testing (QuickCheck, Hypothesis). **Security** — vulnerability discovery. Industry standard practice. |
| Edge case identification | **Software engineering** — universal skill. Boundary conditions, null handling, overflow detection. Critical for writing robust production code. |
| Assert-based validation | **Production systems** — assertion-based monitoring, contract programming (Eiffel, Rust debug_assert). **Testing** — unit test assertions. |
| Binary search on test cases | **Debugging** — git bisect, binary search debugging. Very practical. |

### 8.3 Problem-Solving Strategy

| Skill | Professional Application |
|-------|------------------------|
| Identifying problem type | **Software architecture** — choosing the right tool/pattern for the job. **System design** — matching problems to known solutions (caching, sharding, queuing, etc.). |
| Time complexity estimation | **System design interviews & production** — capacity planning, choosing algorithms for given constraints, predicting system behavior under load. Critical for backend engineering. |
| Managing time under pressure | **On-call / incident response** — debugging production issues under time pressure. **Project management** — prioritizing tasks with deadlines. |

### 8.4 Common Contest Patterns

| Skill | Professional Application |
|-------|------------------------|
| "Binary search on answer" pattern | **System design** — capacity planning ("what's the maximum QPS we can handle?"), auto-scaling thresholds. **ML** — threshold tuning for classifiers. |
| "Range update/query" pattern | **Analytics** — time-range aggregation queries in dashboards. **Finance** — computing portfolio values over date ranges. |
| Pattern recognition in general | **Software engineering** — recognizing when a problem maps to a known algorithm is the core competitive programming skill that transfers most directly to professional work. It is the difference between reinventing a slow solution and applying an efficient known one. |
