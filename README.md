# Quaternion Algebra Project

A SageMath implementation for working with quaternion algebras, including algorithms for isomorphisms, maximal orders, and various algebraic computations. See the internship report in the files for a presentation of the mathematical problems.

## Project Structure

```
.
├── benchmarks
│   ├── iso_quat_alg
│   │   ├── iso_quat_alg_0.txt
│   │   ├── iso_quat_alg_1.txt
│   │   ├── iso_quat_alg_2.txt
│   │   ├── iso_quat_alg_3.txt
│   │   ├── iso_quat_alg_4.txt
│   │   ├── iso_quat_alg_5.txt
│   │   ├── iso_quat_alg_6.txt
│   │   ├── iso_quat_alg_7.txt
│   │   ├── iso_quat_alg_8.txt
│   │   └── iso_quat_alg_9.txt
│   └── maximal_order
│       ├── maximal_order_tensor_Bp_inf.txt
│       ├── max_order_from_inseparable.txt
│       ├── sage_vs_magma_may_29.txt
│       ├── sage_vs_magma_may_31.txt
│       ├── standard_test_M4Q_june_1.txt
│       ├── standard_test_M4Q_june_2.txt
│       ├── standard_test_M4Q_june_3.txt
│       ├── standard_test_M4Q_may_28.txt
│       ├── standard_test_M4Q_may_30.txt
│       ├── standard_test_M4Q_may_31.txt
│       └── starting_order_M4Q.txt
├── database
│   ├── db_iso_matrix_ring.txt
│   ├── db_iso_quat_alg.txt
│   ├── generate_data_base.sage
│   └── read_data_base.txt
├── docs
│   ├── report.pdf
│   └── sage_bugs
│       └── bug_central_idp_sag.txt
├── examples
│   └── examples.txt
├── global_variables.sage
├── load_project.sage
├── README.md
├── src
│   ├── count_matrices_MnFp
│   │   └── nb_irreducible_matrices_MnFp.sage
│   ├── isomorphism
│   │   ├── explicit_iso_matrix_ring.sage
│   │   ├── explicit_iso_quat_alg_equations.sage
│   │   └── explicit_iso_quat_alg.sage
│   ├── maximal_orders
│   │   ├── find_maximal_orders.sage
│   │   └── maximal_orders_utilities.sage
│   ├── minimal_ideals
│   │   ├── idempotent_to_ideals.sage
│   │   ├── minimal_ideals_from_magma.sage
│   │   ├── minimal_ideals_from_sage.sage
│   │   └── minimal_ideals_manually.sage
│   ├── quaternion_recognition
│   │   ├── identify_quaternion_algebra.sage
│   │   ├── identify_standard_involution.sage
│   │   └── normalize_quadratic.sage
│   └── rank_one
│       ├── rank_one_MnFp.sage
│       ├── rank_one_MnQ.sage
│       ├── rank_one_norm_equation.sage
│       └── rank_one_real.sage
├── tests
│   ├── test_formal
│   │   └── formal_expression.sage
│   ├── tests_iso
│   │   ├── test_id_quat_alg.sage
│   │   ├── test_iso_matrix_ring.sage
│   │   ├── test_iso_quat_alg.sage
│   │   ├── test_normalize_quadratic.sage
│   │   └── test_randomization.sage
│   ├── tests_magma
│   │   └── norm_equation.m
│   ├── tests_max_order
│   │   ├── test_central_idempotents.sage
│   │   ├── test_diverse_max_orders.sage
│   │   ├── test_given_max_order_in_quat_alg.sage
│   │   └── test_sage_vs_magma.sage
│   └── tests_rank_one
│       ├── test_rank_one_MnZ.sage
│       └── variables.sage
└── utilities
    ├── algebra_type.sage
    └── utilities.sage
```

## Requirements

- **SageMath** (version 10.6, Release Date: 2025-03-31)
- All dependencies are included within SageMath

## Getting Started

To run this project:

1. Open your terminal and navigate to the `quat_alg_project` directory
2. Start SageMath by typing `sage` and pressing Enter
3. Load all project files at the `sage:` prompt by typing `load("load_all_files.sage")` and pressing Enter

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
```

For more examples, see the files in the `examples/` directory.

## Project Organization

- **`core/`**: Main mathematical algorithms and core functionality
- **`utilities/`**: Helper functions used across multiple modules
- **`database/`**: Functions for generating and managing algebraic databases
- **`maximal_orders/`**: Specialized algorithms for computing maximal orders
- **`minimal_ideals/`**: Algorithms for working with minimal ideals
- **`tests/`**: Comprehensive test suites for validating functionality
- **`examples/`**: Demonstration scripts showing library usage

## Mathematical and Implementation Conventions

- **Morphisms**: A morphism between two algebras A → B (with canonical bases) is represented as a list of dim(A) elements of B, representing f(a₁), ..., f(aₙ) where a₁, ..., aₙ form a basis of A
  
- **Tensor products**: A morphism from A ⊗ B → C is represented as a dictionary with keys (i,j) and values in C, where i ranges from 1 to dim(A) and j ranges from 1 to dim(B)

- **Lattices**: If B is a Q-algebra, a Z-lattice I is represented by a Z-basis `Zbasis_I`, which is a Q-basis of B such that:
  ```
  I = Z·e₁ ⊕ ... ⊕ Z·eₙ
  ```

- **Orders**: Orders are represented as lattices with the assumption that they are stable under multiplication and contain the identity element

## Available Functions

After loading the project, use these helper functions:
- `list_function()`: Display all available project functions
- `list_global_variable()`: Display all global variables

## Author

**Tommy Chakroun**

Developed as part of an M1 internship project (2025) at the Mathematics Department of Virginia Tech, supervised by Travis Morrison.