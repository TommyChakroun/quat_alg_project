# File: load_project.sage

print("============================================")
print("---   Start loading quat_alg_project   ---")
print("============================================")

# --- Loading Source Code ---
print("\n--> Loading Source Code...")
load("global_variables.sage")
load("utilities/utilities.sage")
load("utilities/algebra_type.sage")
load("src/quaternion_recognition/identify_standard_involution.sage")
load("src/quaternion_recognition/normalize_quadratic.sage")
load("src/quaternion_recognition/identify_quaternion_algebra.sage")
load("src/iso_solving_quadratics/solving_quadratics.sage")
load("src/iso_splitting_algebra/maximal_orders/find_maximal_orders.sage")
load("src/iso_splitting_algebra/maximal_orders/maximal_orders_utilities.sage")
load("src/iso_splitting_algebra/minimal_ideals/idempotent_to_ideals.sage")
load("src/iso_splitting_algebra/minimal_ideals/minimal_ideals_from_magma.sage")
load("src/iso_splitting_algebra/minimal_ideals/minimal_ideals_from_sage.sage")
load("src/iso_splitting_algebra/minimal_ideals/minimal_ideals_manually.sage")
load("src/iso_splitting_algebra/rank_one/rank_one_MnFp.sage")
load("src/iso_splitting_algebra/rank_one/rank_one_MnQ.sage")
load("src/iso_splitting_algebra/rank_one/rank_one_real.sage")
load("src/iso_splitting_algebra/explicit_iso_matrix_ring.sage")
load("src/iso_splitting_algebra/explicit_iso_quat_alg.sage")


print("\n============================================")
print("---  All project and test files loaded.  ---")
print("---   Environment is ready to use.       ---")
print("============================================")
