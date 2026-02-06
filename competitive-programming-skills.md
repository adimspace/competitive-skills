# Competitive Programming — International Level (IOI, ICPC & Beyond)

A comprehensive breakdown of skills and topics required to compete at the level of the **International Olympiad in Informatics (IOI)**, **ICPC World Finals**, and top-tier platforms (Codeforces, AtCoder, TopCoder).

---

## 1. Data Structures

### 1.1 Foundational Structures
- Arrays, dynamic arrays (vectors)
- Linked lists (singly, doubly)
- **Stacks**: monotonic stack, stack-based parsing
- **Queues**: standard, deque, monotonic deque (sliding window max/min)
- **Priority queues** (binary heap, d-ary heap)
- Hash tables / hash maps (collision handling, custom hashing to avoid anti-hash attacks)

### 1.2 Tree-Based Structures
- **Binary Search Trees (BST)**: balancing concepts
- **Self-balancing BSTs**: AVL trees, Red-Black trees
- **Treaps** (tree + heap, implicit treap for array operations)
- **Splay trees** (amortized O(log n), link-cut applications)
- **Order-statistic trees** (k-th smallest, rank queries)

### 1.3 Segment Trees
- **Basic segment tree**: range queries (sum, min, max, GCD)
- **Lazy propagation**: range updates + range queries
- **Persistent segment trees**: versioned queries, functional updates
- **Segment tree with coordinate compression**
- **2D segment trees** (segment tree of segment trees)
- **Segment tree beats** (Ji driver segment tree — range chmin/chmax)
- **Dynamic / implicit segment trees** (sparse indices)
- **Merge sort tree** (segment tree with sorted lists in nodes)

### 1.4 Binary Indexed Trees (Fenwick Trees)
- Point update, prefix sum queries
- Range update + point query (difference array trick)
- Range update + range query (two BITs)
- 2D Fenwick trees
- Order-statistic via BIT (compressed coordinates)

### 1.5 Other Advanced Structures
- **Tries**: standard, compressed (Patricia trie), binary trie (XOR queries)
- **Disjoint Set Union (DSU / Union-Find)**: path compression, union by rank/size, rollback DSU
- **Sparse tables**: O(1) range minimum queries (idempotent operations)
- **Sqrt decomposition**: block decomposition for queries
- **Mo's algorithm**: offline range queries in O((N+Q)√N)
- **Mo's algorithm on trees**
- **Euler tour** (flatten tree into array for subtree/path queries)
- **Heavy-Light Decomposition (HLD)**: path queries/updates on trees
- **Centroid decomposition**: distance queries, path-counting on trees
- **Link-Cut Trees**: dynamic tree connectivity, access/cut/link in amortized O(log n)
- **Li Chao trees**: dynamic convex hull trick, line container

---

## 2. Graph Algorithms

### 2.1 Graph Representations & Traversals
- Adjacency list, adjacency matrix, edge list
- **BFS**: shortest path in unweighted graphs, multi-source BFS, 0-1 BFS
- **DFS**: connected components, cycle detection, DFS tree, back/forward/cross edges
- DFS timestamps (tin/tout), subtree queries

### 2.2 Shortest Paths
- **Dijkstra's algorithm**: O((V+E) log V) with priority queue, edge cases with negative weights
- **Bellman-Ford**: negative weight detection, SPFA (Shortest Path Faster Algorithm)
- **Floyd-Warshall**: all-pairs shortest paths in O(V^3)
- **Johnson's algorithm**: all-pairs in sparse graphs (reweighting + Dijkstra)
- **0-1 BFS**: deque-based, for edge weights 0 or 1
- Dial's algorithm, shortest path in DAGs (topological order relaxation)

### 2.3 Minimum Spanning Trees
- **Kruskal's algorithm** (greedy + DSU)
- **Prim's algorithm** (greedy + priority queue)
- **Boruvka's algorithm** (useful in parallel/MST-related problems)
- MST properties: cycle property, cut property
- **Minimum spanning arborescence**: Edmonds' algorithm (directed MST)
- Second-best MST, bottleneck spanning tree

### 2.4 Tree Algorithms
- **Lowest Common Ancestor (LCA)**: binary lifting O(N log N + Q log N), sparse table on Euler tour O(N log N + Q), Farach-Colton-Bender O(N + Q)
- **Tree diameter**: two BFS / DFS approach
- **Rerooting technique** (re-rooting DP)
- **Heavy-Light Decomposition**: path queries on trees
- **Centroid decomposition**: distance-related queries
- **Virtual tree (auxiliary tree)**: subset of important nodes
- Tree isomorphism, hashing trees
- Prufer sequence (tree encoding)

### 2.5 Connectivity & Components
- **Strongly Connected Components (SCC)**: Tarjan's algorithm, Kosaraju's algorithm
- **Bridges and articulation points**: Tarjan's bridge-finding
- **Biconnected components** (block-cut tree)
- **2-edge-connected / 2-vertex-connected components**
- **Euler path/circuit**: Hierholzer's algorithm (directed and undirected)
- Topological sorting (Kahn's BFS, DFS-based)
- Condensation graph (DAG of SCCs)

### 2.6 Network Flow
- **Max-flow / Min-cut theorem**
- **Ford-Fulkerson method**: DFS-based augmenting paths
- **Edmonds-Karp**: BFS-based, O(VE^2)
- **Dinic's algorithm**: O(V^2 E), O(E√V) for unit-capacity graphs
- **Push-relabel** (FIFO, highest-label): O(V^2 √E)
- **Min-cost max-flow**: SPFA-based / Bellman-Ford augmentation, successive shortest paths
- **Hungarian algorithm** (assignment problem)
- Circulation with demands (lower bounds on edges)
- Applications: project selection, closure problems, baseball elimination

### 2.7 Matching
- **Bipartite matching**: Hopcroft-Karp O(E√V), Kuhn's algorithm O(VE)
- **Maximum matching in general graphs**: Edmond's blossom algorithm
- **Konig's theorem**: min vertex cover = max matching (bipartite)
- Dilworth's theorem (antichain decomposition)
- Stable matching (Gale-Shapley algorithm)
- Weighted bipartite matching (Hungarian algorithm / min-cost flow)

### 2.8 Special Graph Problems
- **2-SAT**: implication graph + SCC
- **Graph coloring**: bipartiteness check, chromatic polynomial (small graphs)
- **Dominator tree**: immediate dominators in directed graphs
- **Functional graphs**: successor graphs, cycle detection (Floyd's / Brent's)
- **Cactus graphs**: special structure, DP on cacti
- **Chordal graphs**, **interval graphs** — recognition and exploitation

---

## 3. Dynamic Programming

### 3.1 Classical DP
- **1D DP**: Fibonacci, staircase, coin change, rod cutting
- **2D DP**: grid paths, LCS (Longest Common Subsequence), edit distance
- **Knapsack**: 0/1, unbounded, bounded, multi-dimensional
- **LIS** (Longest Increasing Subsequence): O(n log n) with patience sorting
- **Interval DP**: matrix chain multiplication, optimal BST, stone merging
- **Counting DP**: number of ways, paths, partitions

### 3.2 Bitmask DP
- Subset enumeration: iterating over subsets of a bitmask in O(3^n)
- Traveling Salesman Problem (TSP): O(2^n * n^2)
- Hamiltonian path/cycle via bitmask
- Assignment problems
- DP on subsets of a set (broken profile, profile DP)

### 3.3 DP on Trees
- Subtree DP (rooted tree problems)
- **Rerooting DP** (compute answer for every root efficiently)
- Tree DP with knapsack (merge children)
- Independent set on trees, vertex cover on trees
- Tree diameter, tree center via DP

### 3.4 Digit DP
- Counting numbers in a range [L, R] with specific digit properties
- Sum of digits, digit frequency, divisibility constraints
- Leading zeros handling, tight/free states

### 3.5 DP on DAGs
- Longest/shortest path on DAGs (topological order)
- Counting paths on DAGs
- DP combined with SCC condensation

### 3.6 SOS DP (Sum over Subsets)
- Subset sum transform in O(n * 2^n)
- Superset sum, Zeta/Mobius transforms
- Applications: counting pairs with AND/OR properties

### 3.7 DP Optimizations
- **Knuth's optimization**: O(n^2) for certain interval DPs (quadrangle inequality)
- **Divide and Conquer optimization**: O(n m log n) when the optimal split point is monotone
- **Convex Hull Trick (CHT)**: linear DP with convex/concave cost, Li Chao tree
- **SMAWK algorithm**: totally monotone matrix optimization
- **Aliens trick (Lambda optimization / WQS binary search)**: converting constrained DP to unconstrained
- **Slope trick**: DP on convex piecewise-linear functions
- **Hirschberg's algorithm**: space-efficient LCS/edit distance

### 3.8 Other DP Types
- **Probability / Expected value DP**
- **DP with matrix exponentiation** (linear recurrences in O(k^3 log n))
- **Broken profile DP** (plug DP / tiling problems)
- **DP on permutations**
- **Connection profile DP** (Steiner tree, Hamiltonian path on grid)
- **DP on convex hull** (inserting points, maintaining hull properties)

---

## 4. String Algorithms

### 4.1 Hashing
- **Polynomial hashing** (Rabin-Karp): rolling hash for substring matching
- Double hashing to reduce collision probability
- Anti-hash attacks: randomized bases and moduli
- Hashing for palindrome detection, substring comparison

### 4.2 Pattern Matching
- **KMP (Knuth-Morris-Pratt)**: failure function, O(n+m) pattern search
- **Z-function**: Z-array, pattern matching, string period
- **Aho-Corasick automaton**: multi-pattern matching, failure links, dictionary matching
- **Rabin-Karp** (hash-based multi-pattern)

### 4.3 Suffix Structures
- **Suffix array**: O(n log n) or O(n) construction, LCP array (Kasai's algorithm)
- **Suffix tree**: Ukkonen's algorithm, applications (longest repeated substring, etc.)
- **Suffix automaton (SAM)**: online construction, substring counting, shortest non-occurring string
- Generalized suffix structures (multiple strings)

### 4.4 Palindromes
- **Manacher's algorithm**: all maximal palindromic substrings in O(n)
- **Palindromic tree (Eertree)**: all distinct palindromic substrings
- Palindrome hashing techniques

### 4.5 Other String Techniques
- **Trie-based algorithms**: autocomplete, XOR maximization
- **Burrows-Wheeler Transform** (BWT) and FM-index
- Lyndon factorization, Booth's algorithm (lexicographically smallest rotation)
- Duval's algorithm, Lyndon words
- Minimum / maximum suffix, lexicographic comparisons via suffix array
- Run enumeration (Runs theorem)
- Bitap algorithm (for approximate matching)

---

## 5. Mathematics for CP

### 5.1 Number Theory
- **Sieve of Eratosthenes**: standard, linear sieve, segmented sieve
- **Modular arithmetic**: modular inverse (Fermat / extended GCD), modular exponentiation
- **Extended Euclidean algorithm**: solving ax + by = gcd(a,b)
- **Chinese Remainder Theorem** (CRT): system of congruences
- **Euler's totient function** phi(n), totient sieve
- Multiplicative functions, Mobius function, Mobius inversion
- **Discrete logarithm**: Baby-step Giant-step (BSGS), Pohlig-Hellman
- **Primitive roots**: existence, computation
- **Miller-Rabin** primality test (deterministic for n < 3.3 * 10^24 with fixed bases)
- **Pollard's rho** factorization

### 5.2 Combinatorics
- **Binomial coefficients**: Pascal's triangle, Lucas' theorem, precomputation with modular inverse
- **Principle of Inclusion-Exclusion** (PIE)
- **Burnside's lemma / Polya enumeration theorem** (counting under symmetry)
- Catalan numbers, Stirling numbers (first/second kind), Bell numbers, partition numbers
- Derangements, subfactorials
- **Generating functions** (ordinary, exponential)
- Twelvefold way (balls into boxes framework)
- Stars and bars, compositions, partitions of integers

### 5.3 Linear Algebra
- **Gaussian elimination**: solving systems, finding rank, determinant, matrix inverse
- **Matrix exponentiation**: computing n-th term of linear recurrences in O(k^3 log n)
- Linear basis (XOR basis for maximum XOR subset)
- Systems of linear equations over GF(2) (binary / XOR systems)
- Kirchhoff's matrix tree theorem (counting spanning trees)

### 5.4 Probability & Expected Value
- Linearity of expectation
- Geometric distribution, coupon collector problem
- Markov chains and absorbing states
- DP for expected values, conditional expectation

### 5.5 Game Theory
- **Sprague-Grundy theorem**: nim-values, game sums
- **Nim** and its variants (Wythoff's game, Moore's Nim)
- Games on graphs, games on DAGs
- Multi-pile games, impartial vs. partisan games

### 5.6 Polynomials & Transforms
- **FFT (Fast Fourier Transform)**: polynomial multiplication in O(n log n)
- **NTT (Number Theoretic Transform)**: FFT over modular arithmetic
- **Karatsuba multiplication** (O(n^1.585))
- Polynomial division, interpolation (Lagrange)
- Formal power series: inverse, logarithm, exponentiation, composition
- Walsh-Hadamard Transform (XOR/AND/OR convolution)
- Subset convolution

---

## 6. Computational Geometry

### 6.1 Primitives
- Point, vector, line, segment, ray representations
- **Cross product**: orientation test (CCW / CW / collinear)
- **Dot product**: angle computation, projection
- Line-line intersection, segment-segment intersection
- Point-segment distance, point-line distance

### 6.2 Convex Hull
- **Graham scan** O(n log n)
- **Andrew's monotone chain** O(n log n)
- Jarvis march (gift wrapping) O(nh)
- Dynamic convex hull (online insertion)
- Convex hull trick (for DP optimization)

### 6.3 Sweep Line & Plane Algorithms
- **Line sweep**: closest pair of points, rectangle union area, segment intersection
- **Event-based processing**: opening/closing events
- Plane sweep for Voronoi diagrams (Fortune's algorithm)
- Angular sweep

### 6.4 Polygon Operations
- **Shoelace formula** (area of simple polygon)
- Point-in-polygon tests (ray casting, winding number)
- Polygon convexity check
- **Minkowski sum** of convex polygons
- Polygon clipping (Sutherland-Hodgman)
- Rotating calipers (diameter, width, closest/farthest pair on convex hull)

### 6.5 Advanced Geometry
- **Half-plane intersection** (O(n log n))
- **Voronoi diagram** and **Delaunay triangulation**
- Smallest enclosing circle (Welzl's algorithm, randomized O(n))
- k-d trees (nearest neighbor queries in 2D/3D)
- Pick's theorem (lattice point counting)
- Great circle distance (spherical geometry)

---

## 7. General Algorithmic Techniques

### 7.1 Searching & Sorting
- **Binary search**: on arrays, on answer (parametric search), real-valued binary search
- **Ternary search**: unimodal functions
- **Sorting algorithms**: merge sort (inversions count), quick select (k-th element O(n))
- **Coordinate compression**
- **Parallel binary search**: answering multiple binary-search queries simultaneously

### 7.2 Greedy Algorithms
- Activity selection, interval scheduling, interval partitioning
- Huffman coding
- Fractional knapsack
- Exchange argument proofs
- Matroid intersection (advanced)

### 7.3 Divide and Conquer
- Merge sort applications (inversion counting, closest pair)
- Divide and conquer on segment tree / merge sort tree
- CDQ divide and conquer (offline multidimensional queries)
- Centroid decomposition as divide and conquer on trees

### 7.4 Bit Manipulation
- Bitwise operations (AND, OR, XOR, NOT, shifts)
- **XOR properties**: linear independence, XOR basis / Gaussian elimination over GF(2)
- Popcount (Hamming weight), lowest set bit tricks
- Iterating over submasks
- Gray code

### 7.5 Randomized Algorithms
- Randomized hashing (anti-hack)
- Randomized pivot (quickselect, quicksort)
- Randomized algorithms for geometry (smallest enclosing circle)
- Birthday paradox applications
- Schwartz-Zippel lemma (polynomial identity testing)
- Treaps (randomized BSTs)

### 7.6 Offline & Online Techniques
- **Mo's algorithm**: offline range queries O((N+Q)√N), with updates O(N^(5/3))
- **Square root decomposition**: block-based query answering
- **Persistent data structures**: immutable updates, versioned queries
- **Fractional cascading**: multi-level binary search
- **Offline processing with divide and conquer (CDQ)**
- **Rollback** technique (undo operations in DSU, etc.)

### 7.7 Miscellaneous Techniques
- **Meet in the middle**: splitting search space in half, O(2^(n/2)) instead of O(2^n)
- **Two pointers / sliding window**: subarray problems, pair-finding
- **Small-to-large merging** (DSU on tree / Euler tour merging)
- **Virtual tree construction** (auxiliary tree with key nodes)
- **Interactive problems**: binary search with queries, adaptive strategies
- **Constructive algorithms**: building solutions that meet constraints
- **Output-sensitive algorithms**: performance depends on output size

---

## 8. Implementation & Practical Skills

### 8.1 Language Mastery (C++ / Python / Java)
- **C++ STL**: vector, set, map, unordered_map, priority_queue, bitset, __builtin_* functions
- **Pragmas and optimizations**: #pragma GCC optimize, fast I/O (ios::sync_with_stdio)
- Template-based coding (pre-written libraries)
- Avoiding TLE: constant factor optimization, cache-friendly access
- Precision handling (long long vs int, __int128, floating point pitfalls)

### 8.2 Debugging & Testing
- Stress testing: brute force vs. optimized, random test generation
- Edge case identification (n=0, n=1, maximum constraints)
- Assert-based validation
- Binary search on test cases (finding minimal failing input)

### 8.3 Problem-Solving Strategy
- Reading and parsing problem statements carefully
- Identifying problem type (greedy, DP, graph, math, data structure)
- Estimating time complexity from constraints (n <= 10^5 suggests O(n log n), etc.)
- Choosing the right algorithm for the constraint bounds
- Managing time in contests (easy first, skip and return)

### 8.4 Common Contest Patterns
- "Answer is small" — brute force or meet in the middle
- "Queries offline" — Mo's algorithm, CDQ, persistent structures
- "Range update/query" — segment tree, BIT, difference arrays
- "Minimize/maximize with constraint" — binary search on answer + greedy check
- "Count pairs with property" — two pointers, sorting, transforms
- "Shortest path with state" — Dijkstra on layered graph / state-space BFS

---

## Recommended Progression

| Stage | Focus |
|-------|-------|
| **Foundation** | Sorting, BFS/DFS, basic DP, greedy, binary search, number theory basics |
| **Intermediate** | Segment trees, graph algorithms (Dijkstra, MST, SCC), string hashing, combinatorics |
| **Advanced** | Flows, advanced DP optimizations, suffix structures, computational geometry, FFT/NTT |
| **World Finals / IGM** | Persistent structures, link-cut trees, DP on subsets, advanced flows, contest strategy |
