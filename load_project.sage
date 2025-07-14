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
load("src/minimal_ideals/idempotent_to_ideals.sage")
load("src/minimal_ideals/minimal_ideals_from_sage.sage")
load("src/minimal_ideals/minimal_ideals_manually.sage")
load("src/minimal_ideals/minimal_ideals_from_magma.sage")
load("src/maximal_orders/maximal_orders_utilities.sage")
load("src/maximal_orders/find_maximal_orders.sage")
load("src/rank_one/rank_one_MnFp.sage")
load("src/rank_one/rank_one_MnQ.sage")
load("src/rank_one/rank_one_real.sage")
load("src/rank_one/rank_one_norm_equation.sage")
load("src/isomorphism/explicit_iso_matrix_ring.sage")
load("src/isomorphism/explicit_iso_quat_alg.sage")
load("src/isomorphism/explicit_iso_quat_alg_equations.sage")

# --- Loading Test Function Definitions ---
print("\n--> Loading Test Functions...")
# Isomorphism Tests
load("tests/tests_iso/test_id_quat_alg.sage")
load("tests/tests_iso/test_iso_matrix_ring.sage")
load("tests/tests_iso/test_iso_quat_alg.sage")
load("tests/tests_iso/test_normalize_quadratic.sage")
load("tests/tests_iso/test_randomization.sage")
# Maximal Order Tests
load("tests/tests_max_order/test_diverse_max_orders.sage")
load("tests/tests_max_order/test_given_max_order_in_quat_alg.sage")
load("tests/tests_max_order/test_sage_vs_magma.sage")
load("tests/tests_max_order/test_central_idempotents.sage")
## Rank One Tests
load("tests/tests_rank_one/test_rank_one_MnZ.sage")


print("\n============================================")
print("---  All project and test files loaded.  ---")
print("---   Environment is ready to use.       ---")
print("============================================")


def list_function():
    print("Functions available in the project:\n")
    functions = [
        "get_coefficients",
        "is_in_field",
        "dimension",
        "change_matrix",
        "coordinate",
        "evaluate",
        "is_linear_iso",
        "is_algebra_homomorphism",
        "random_invertible_matrix",
        "matrix_in_basis",
        "",
        "structure_constants",
        "finite_dimensional_algebra_format",
        "tensor",
        "opposite",
        "quat_alg_mixed_with_table",
        "",
        "has_standard_involution",
        "reduced_trace_si",
        "standard_involution",
        "reduced_norm",
        "",
        "is_isomorphic_to_quaternion_algebra",
        "quaternion_structure",
        "",
        "quat_alg_iso_from_matrix_ring_iso",
        "explicit_iso_quat_alg",
        "",
        "matrix_ring_iso_from_algebra_iso",
        "basis_of_translated",
        "matrix_ring_iso_from_rank_one",
        "",
        "normalize_quadratic_form_over_Z_loc_2",
        "normalize_symetric_matrix_over_Z_loc_2",
        "",
        "right_rank",
        "min_rank",
        "rank_one",
        "",
        "read_line",
        "read_one_block",
        "read_file",
        "random_iso_quat_alg",
        "generate_quaternion_iso_db",
        "generate_matrix_ring_iso_db",
        "regenerate_data_base",
        "",
        "strictly_bigger_order_local",
        "strictly_bigger_order_local_printers",
        "strictly_bigger_order",
        "is_maximal_order",
        "max_order_containing_order_local",
        "max_order_containing_order",
        "max_order",
        "max_order_tensor_quat_alg",
        "",
        "Z_mod_span_by_rows",
        "intersection_Z_mod",
        "kernel_mod_non_split",
        "kernel_mod",
        "lattice_LLL",
        "is_in_lattice",
        "is_sub_lattice",
        "is_order",
        "are_equals_lattices",
        "left_order",
        "reduced_trace",
        "discriminant",
        "finite_algebra_from_order",
        "quotient_algebra_ideal",
        "kernel_Z_mod_map",
        "",
        "center_perso"
        "idempotents_to_ideals",
        "",
        "central_idempotents_magma",
        "minimal_ideals_magma",
        "",
        "structure_constants_subspace",
        "subalgebra_from_subspace",
        "orthogonal_decomposition_sage",
        "idempotents_commutative_sage",
        "central_idempotents_sage",
        "minimal_ideals_sage",
        "",
        "minimal_polynomial",
        "split_from_idemp",
        "idempotents_commutative_perso_split",
        "idempotents_commutative_perso_direct",
        "central_idempotents_perso",
        "minimal_ideals_perso",
        "minimal_ideals_commutative",
        "minimal_ideals_manually",
        "",
        "test_1",
        "test_2",
        "test_3",
        "test_4",
        "test_5",
        "test_6",
        "test_7",
        "test_8",
        "test_9",
        "test_10",
        "test_11",
        "test_12",
        "test_13",
        "test_14",
        "test_15",
        "test_16"
    ]
    
    for name in functions:
        if name == "":
            print()
        else :
            print("-", name)
    return None

def list_global_variable():
    print("ListIsoQuatAlgDb")
    print("ListIsoM4QDb")
    print()
    print("MatriceRingM4Q")
    print()
    print("RadomInvertibleMatrix_16x16_example_1")
    print("RadomInvertibleMatrix_16x16_example_2")
    print("RadomInvertibleMatrix_16x16_example_3")
    print()
    print("RandomBasisM4Q_example_1")
    print("RandomBasisM4Q_example_2")
    print("RandomBasisM4Q_example_3")
    print()
    print("BasisRandomOrder_M4Q_example_1")
    print("BasisRandomOrder_M4Q_example_2")
    print("BasisRandomOrder_M4Q_example_3")
    print()
    print("MixedM4Q_example_1")
    print("MixedM4Q_example_2")
    print("MixedM4Q_example_3")
    return None