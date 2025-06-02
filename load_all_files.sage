# load_all_files.sage
# Loads all functions from the quaternion algebra project
print("Loading all quaternion algebra project files...")

# Load utilities first (dependencies for other modules)
print("- Loading utilities...")
load("utilities/algebra_type.sage")
load("utilities/utilities.sage")

# Load core functions
print("- Loading core algorithms...")
load("core/explicit_iso_matrix_ring.sage")
load("core/explicit_iso_quat_alg.sage")
load("core/identify_quaternion_algebra.sage")
load("core/identify_standard_involution.sage")
load("core/normalize_quadratic.sage")
load("core/rank_one.sage")

# Load database functions
print("- Loading database functions...")
load("database/database_utilities.sage")
load("database/generate_iso_matrix_ring.sage")
load("database/generate_iso_quat_alg.sage")

# Load maximal order functions
print("- Loading maximal order algorithms...")
load("maximal_orders/maximal_orders_utilities.sage")
load("maximal_orders/minimal_ideals_manually.sage")
load("maximal_orders/minimal_ideals_from_sage.sage")
load("maximal_orders/minimal_ideals_from_magma.sage")
load("maximal_orders/find_maximal_orders.sage")

print("- Loading general test scripts (except those with global constants)...")
load("tests/test_central_idempotents.sage")
load("tests/test_id_quat_alg.sage")
load("tests/test_iso_matrix_ring.sage")
load("tests/test_iso_quat_alg.sage")

print("-Loading maximal order test scripts (except those with global constants)...")
load("tests_max_order/test_diverse_max_orders.sage")
load("tests_max_order/test_sage_vs_magma.sage")


print("✓ All files loaded successfully!")
print("You can now use all quaternion algebra functions.")
print("See examples in the examples/ folder.")