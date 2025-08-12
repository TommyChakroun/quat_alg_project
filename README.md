# Quaternion Algebra Project

A SageMath implementation for working with quaternion algebras, including algorithms for isomorphisms, maximal orders, and various algebraic computations. See the internship report in the files for a presentation of the mathematical problems.

## Project Structure

src
├── iso_isotropicsubspcae
│   └── quadratic_equivalence.m
├── iso_solving_quadratics
│   ├── solving_quadratics.m
│   └── solving_quadratics.sage
├── iso_splitting_algebra
│   ├── explicit_iso_matrix_ring.sage
│   ├── explicit_iso_quat_alg.sage
│   ├── maximal_orders
│   │   ├── find_maximal_orders.sage
│   │   └── maximal_orders_utilities.sage
│   ├── minimal_ideals
│   │   ├── idempotent_to_ideals.sage
│   │   ├── minimal_ideals_from_magma.sage
│   │   ├── minimal_ideals_from_sage.sage
│   │   └── minimal_ideals_manually.sage
│   └── rank_one
│       ├── rank_one_MnFp.sage
│       ├── rank_one_MnQ.sage
│       └── rank_one_real.sage
└── quaternion_recognition
    ├── identify_quaternion_algebra.sage
    ├── identify_standard_involution.sage
    └── normalize_quadratic.sage

## Requirements

- **SageMath** (version 10.6, Release Date: 2025-03-31)
- All dependencies are included within SageMath

## Getting Started

To run this project:

1. Open your terminal and navigate to the `quat_alg_project` directory
2. Start SageMath by typing `sage` and pressing Enter
3. Load all project files at the `sage:` prompt by typing `load("load_project.sage")` and pressing Enter

This loading process may take a few seconds. Upon successful completion, all project functions and variables will be available in your SageMath environment.

### Example Session

```
[tommy@tommy-hplaptop15db0xxx quat_alg_project]$ sage
┌────────────────────────────────────────────────────────────────────┐
│ SageMath version 10.6, Release Date: 2025-03-31                    │
│ Using Python 3.13.3. Type "help()" for help.                       │
└────────────────────────────────────────────────────────────────────┘
sage: load("load_project.sage")
============================================
---   Start loading quat_alg_project   ---
============================================

--> Loading Source Code...

--> Loading Test Functions...

============================================
---  All project and test files loaded.  ---
---   Environment is ready to use.       ---
============================================
sage: 
```

## Usage Example

Here is a basic example demonstrating maximal order computation:

```sage
sage: A = QuaternionAlgebra(QQ, -1, -1)
sage: B = QuaternionAlgebra(QQ, -1, -2)
sage: C = tensor(A, opposite(B))
sage: %time Zbasis_L = max_order(C)
CPU times: user 4.28 s, sys: 70.5 ms, total: 4.35 s
Wall time: 4.06 s
sage: discriminant(C, Zbasis_L)
1
sage: is_order(C, Zbasis_L)
True
sage: len(Zbasis_L)
16


sage: A = QuaternionAlgebra(QQ,-69565111111111111116,-465656533333333336635)
sage: A.ramified_places()
([3, 151, 6389, 134917, 1097729, 238505243, 5280987012817],
 [Ring morphism:
    From: Rational Field
    To:   Real Field with 53 bits of precision
    Defn: 1 |--> 1.00000000000000])
sage: B = QuaternionAlgebra(QQ,[3, 151, 6389, 134917, 1097729, 238505243, 5280987012817],[1/2])
sage: B
Quaternion Algebra (-97, -539890808015802473001758328641975308911) with base ring Rational Field
sage: %time iso_quat_alg(A,B)
CPU times: user 2.97 s, sys: 2.74 ms, total: 2.97 s
Wall time: 3 s
(True,
 [1,
  3051505944662292589324023606525881581992/6109654102414479896365283261459*i - 14460173029262124110/62986124767159586560466837747*j + 109794820220326674708/6109654102414479896365283261459*k,
  -52374902494206141125821144016205699232545333865409583914094504623068977/43844931822703108744970588984500902943970259414258312936269794*i + 28352834537936478837330486790463326599054296301417/226004803209809838891600974146911870845207522753908829568401*j + 81414424032369471418810813441/1043561314078444167871600909789049455618*k,
  -238669359081341878828146709616831594349687739020095350981201942644203269418757063/21922465911351554372485294492250451471985129707129156468134897*i - 1324828360187091476869435344670327065963284786643000885454794/226004803209809838891600974146911870845207522753908829568401*j - 110399940432773669468529148179629866821/521780657039222083935800454894524727809*k])
sage: 
```

For more examples, see the files in the `examples/` directory.

## Project Organization

- **`src/`**: Main mathematical algorithms and core functionality
- **`utilities/`**: Helper functions used across multiple modules
- **`database/`**: Functions for generating and managing algebraic databases
- **`tests/`**: Comprehensive test suites for validating functionality
- **`examples/`**: Demonstration scripts showing library usage
- **`benchmarks/`**: This directory contains time comparisons of the different algorithms.
- **`docs/`**: This directory holds the mathematical report associated with the project and notes on some SageMath bugs encountered during the implementation

## Mathematical and Implementation Conventions

- **Morphisms**: A morphism between two algebras A → B (with canonical bases) is represented as a list of dim(A) elements of B, representing f(a₁), ..., f(aₙ) where a₁, ..., aₙ form a basis of A
  
- **Tensor products**: A morphism from A ⊗ B → C is represented as a dictionary with keys (i,j) and values in C, where i ranges from 1 to dim(A) and j ranges from 1 to dim(B)

- **Lattices**: If B is a Q-algebra, a Z-lattice I is represented by a Z-basis `Zbasis_I`, which is a Q-basis of B such that:
  ```
  I = Z·e₁ ⊕ ... ⊕ Z·eₙ
  ```

- **Orders**: Orders are represented as lattices with the assumption that they are stable under multiplication and contain the identity element


## Author

**Tommy Chakroun**

Developed as part of an M1 internship project (2025) at the Mathematics Department of Virginia Tech, supervised by Travis Morrison.