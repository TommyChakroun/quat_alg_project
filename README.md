# Quaternion Algebra Project

A SageMath implementation for working with quaternion algebras, including algorithms for isomorphisms, maximal orders, and various algebraic computations.

## Project Structure

```
.
├── core
│   ├── explicit_iso_matrix_ring.sage
│   ├── explicit_iso_quat_alg.sage
│   ├── identify_quaternion_algebra.sage
│   ├── identify_standard_involution.sage
│   ├── normalize_quadratic.sage
│   └── rank_one.sage
├── database
│   ├── database_utilities.sage
│   ├── db_iso_matrix_ring.txt
│   ├── db_iso_quat_alg_from_tensor_iso.txt
│   ├── db_iso_quat_alg.txt
│   ├── generate_iso_matrix_ring.sage
│   ├── generate_iso_quat_alg.sage
│   ├── max_order_M4Q.txt
│   └── read_data_base.txt
├── examples
│   ├── bug_central_idp_sag.txt
│   ├── estimate_randomization.txt
│   ├── ex_1.txt
│   ├── maximal_order_tensor_Bp_inf.txt
│   ├── sage_vs_magma_may_29.txt
│   ├── sage_vs_magma_may_31.txt
│   ├── standard_test_M4Q_may_28.txt
│   ├── standard_test_M4Q_may_30.txt
│   ├── standard_test_M4Q_may_31.txt
│   └── starting_order_M4Q.txt
├── load_all_files.sage
├── maximal_orders
│   ├── find_maximal_orders.sage
│   ├── maximal_orders_utilities.sage
│   ├── minimal_ideals_from_magma.sage
│   ├── minimal_ideals_from_sage.sage
│   └── minimal_ideals_manually.sage
├── README.md
├── tests
│   ├── test_central_idempotents.sage
│   ├── test_id_quat_alg.sage
│   ├── test_iso_matrix_ring.sage
│   ├── test_iso_quat_alg.sage
│   ├── test_normalize_quadratic.sage
│   └── test_randomization.sage
├── tests_max_order
│   ├── comments.pdf
│   ├── test_diverse_max_orders.sage
│   ├── test_given_max_order_in_quat_alg.sage
│   ├── test_max_order_structure_constants.sage
│   └── test_sage_vs_magma.sage
└── utilities
    ├── algebra_type.sage
    └── utilities.sage
```

## Requirements

- **SageMath** (version 10.6, Release Date: 2025-03-3)
- All dependencies are included within SageMath

## How to Run

### Important Setup

**Always launch SageMath from the project root directory** for all file paths to work correctly:

```bash
cd quat_alg_project
sage
```

### Quick Start - Load All Functions

To load all project functions at once:

```sage
sage: load("load_all_files.sage")
```

This will load all core functions, utilities, and algorithms, making them available in your SageMath session.

### Running Individual Components

You can also load specific modules:

```sage
# Load core algorithms
sage: load("core/identify_quaternion_algebra.sage")
sage: load("core/explicit_iso_quat_alg.sage")

# Load utilities
sage: load("utilities/utilities.sage")

# Load database functions
sage: load("database/database_utilities.sage")

# Load maximal order algorithms
sage: load("maximal_orders/find_maximal_order.sage")
```

### Running Tests

Execute test files to verify functionality:

```sage
sage: load("tests/test_id_quat_alg.sage")
sage: load("tests/test_iso_quat_alg.sage")
sage: load("tests_max_order/test_max_order_in_tensor_quat_alg.sage")
```

### Running Examples

```sage
sage: load("examples/example_script.sage")  # Replace with actual example files
```

## Usage

After loading the functions (via `load_all_files.sage` or individual modules), you can use the various quaternion algebra functions directly in your SageMath session.

## Project Organization

- **`core/`**: Contains the main mathematical algorithms and core functionality
- **`utilities/`**: Helper functions used across multiple modules
- **`database/`**: Functions for generating and managing algebraic databases
- **`maximal_orders/`**: Specialized algorithms for computing maximal orders
- **`tests/`**: Comprehensive test suites for validating functionality
- **`examples/`**: Demonstration scripts showing how to use the library


## Author
Tommy Chakroun

Developed as part of an M1 internship project (2025) at the math department of Virginia Tech supervised by Travis Morrison.