# Quaternion Algebra Project

A SageMath implementation for working with quaternion algebras, including algorithms for isomorphisms, maximal orders, and various algebraic computations.

## Project Structure

```
quat_alg_project/
├── README.md
├── load_all_files.sage          # Loads all project functions
├── core/                        # Core algorithms and functions
│   ├── explicit_iso_matrix_ring.sage
│   ├── explicit_iso_quat_alg.sage
│   ├── identify_quaternion_algebra.sage
│   ├── identify_standard_involution.sage
│   ├── normalize_quadratic.sage
│   └── rank_one.sage
├── utilities/                   # Helper functions and utilities
│   ├── algebra_type.sage
│   ├── utilities.sage
│   └── load_all_utilities.sage
├── database/                    # Database generation and utilities
│   ├── database_utilities.sage
│   ├── generate_iso_matrix_ring.sage
│   ├── generate_iso_quat_alg.sage
│   ├── db_iso_matrix_ring.txt
│   ├── db_iso_quat_alg_from_tensor_iso.txt
│   └── db_iso_quat_alg.txt
├── maximal_orders/             # Maximal order computations
│   ├── find_maximal_orders.sage
│   ├── maximal_orders_utilities.sage
│   ├── minimal_ideals_from_magma.sage
│   ├── minimal_ideals_from_sage.sage
│   └── minimal_ideals_manually.sage
├── tests/                      # Test files
│   ├── test_central_idempotents.sage
│   ├── test_id_quat_alg.sage
│   ├── test_iso_matrix_ring.sage
│   ├── test_iso_quat_alg.sage
│   ├── test_normalize_quadratic.sage
│   └── test_randomization.sage
├── tests_max_order/           # Specific tests for maximal orders
│   ├── test_diverse_max_orders.sage
│   ├── test_max_order_in_tensor_quat_alg.sage
│   ├── test_max_order_in_tensor_of_quat_alg.sage
│   └── test_max_order_structure_constants.sage
├── examples/                  # Usage examples and demonstrations
│   ├── estimate_randomization.txt
│   ├── ex_1.txt
│   ├── max_order_in_tensor_quat_alg.txt
│   ├── standard_test_M4Q.txt
│   └── starting_order_M4Q.txt
└── comments.pdf              # Additional documentation
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