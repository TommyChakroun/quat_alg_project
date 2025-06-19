load("utilities/utilities.sage")
load("utilities/algebra_type.sage")
load("src/maximal_orders/maximal_orders_utilities.sage")
load("src/maximal_orders/find_maximal_orders.sage")

import time
from datetime import timedelta

def format_time(seconds):
    """Format time duration in a readable way"""
    if seconds < 1:
        return f"{seconds*1000:.2f} ms"
    elif seconds < 60:
        return f"{seconds:.3f} s"
    else:
        return str(timedelta(seconds=int(seconds)))

def print_section(title, width=60):
    """Print a nicely formatted section header"""
    print("=" * width)
    print(f" {title.center(width-2)} ")
    print("=" * width)
    print()

def print_subsection(title, width=50):
    """Print a nicely formatted subsection header"""
    print("-" * width)
    print(f" {title} ")
    print("-" * width)

def test_13(a, b, c, d,parallel = False):
    print_section("QUATERNION ALGEBRA MAXIMAL ORDER COMPUTATION")
    
    # ==================== PREPARATION ====================
    print_section("PREPARATION", 50)
    
    prep_start = time.time()
    
    print(f"• Creating Quaternion Algebra A: ({a}, {b} | Q)")
    A = QuaternionAlgebra(QQ, a, b)
    
    print(f"• Creating Quaternion Algebra B: ({c}, {d} | Q)")
    B = QuaternionAlgebra(QQ, c, d)
    
    print(f"• Computing opposite algebra B^op")
    Bop = opposite(B)
    
    print(f"• Computing tensor product C = A ⊗ B^op")
    C = tensor(A, Bop)
    print(f"  → Dimension of C: {C.dimension()}")
    
    print(f"• Computing structure constants...")
    T = C.table()
    structure_constants = [T[j][i][k] for i in range(16) for j in range(16) for k in range(16)]
    print(f"  → Structure constants computed ({len(structure_constants)} elements)")
    
    print(f"• Getting canonical basis of C")
    basis_C = C.basis()
    
    print(f"• Computing canonical order O (left order of canonical basis)...")
    Zbasis_O = left_order(C, basis_C)
    
    prep_time = time.time() - prep_start
    
    print(f"\n📋 CANONICAL ORDER O:")
    print(f"   Z-basis of O:")
    print(Zbasis_O)
    print(f"   Z-basis elements: {len(Zbasis_O)}")
    print(f"   Discriminant: {discriminant(C, Zbasis_O)}")
    
    print(f"\n⏱️  Preparation completed in {format_time(prep_time)}")
    
    # ==================== GOAL ====================
    print_section("GOAL: COMPUTE A MAXIMAL ORDER L IN C", 60)
    
    # ==================== COMPETITION ====================
    print_section("SAGE vs MAGMA COMPARISON", 60)
    
    # ==================== MAGMA ====================
    print_subsection("🔥 MAGMA COMPUTATION")
    
    magma_start = time.time()
    
    print("• Preparing MAGMA code...")
    magma_code = (
        f"Q := Rationals();\n"
        f"structure_constants := [ {', '.join(map(str, structure_constants))} ];\n"
        f"C := AssociativeAlgebra<Q, 16| structure_constants>;\n"
        f"L:= MaximalOrder(C);\n"
        f"basis:= Basis(L);\n"
        f"[Eltseq(e) : e in basis];\n"
    )
    
    print("• Executing MAGMA computation...")
    magma_output = magma.eval(magma_code)
    
    print("• Converting MAGMA output to Sage format...")
    convert_output = sage_eval(magma_output)
    Zbasis_L_Magma = [sum(convert_output[i][k]*basis_C[k] for k in range(16)) for i in range(16)]
    
    magma_time = time.time() - magma_start
    
    print(f"\n📊 MAGMA RESULTS:")
    print(f"   Z-basis of L_MAGMA:")
    print(Zbasis_L_Magma)
    print(f"   Z-basis elements: {len(Zbasis_L_Magma)}")
    print(f"   Is Order: {is_order(C,Zbasis_L_Magma)}")
    print(f"   Discriminant: {discriminant(C, Zbasis_L_Magma)}")
    print(f"   O ⊆ L_MAGMA: {is_sub_lattice(C, Zbasis_O, Zbasis_L_Magma)}")
    print(f"   ⏱️  MAGMA time: {format_time(magma_time)}")
    
    # ==================== SAGE ====================
    print_subsection("🐍 SAGE COMPUTATION")
    
    sage_start = time.time()
    
    print("• Computing maximal order using Sage...")
    Zbasis_L_SAGE = max_order(C,parallel = parallel)
    
    sage_time = time.time() - sage_start
    
    print(f"\n📊 SAGE RESULTS:")
    print(f"   Z-basis of L_SAGE:")
    print(Zbasis_L_SAGE)
    print(f"   Z-basis elements: {len(Zbasis_L_SAGE)}")
    print(f"   Is Order: {is_order(C,Zbasis_L_SAGE)}")
    print(f"   Discriminant: {discriminant(C, Zbasis_L_SAGE)}")
    print(f"   O ⊆ L_SAGE: {is_sub_lattice(C, Zbasis_O, Zbasis_L_SAGE)}")
    print(f"   ⏱️  SAGE time: {format_time(sage_time)}")
    
    # ==================== SUMMARY ====================
    print_section("PERFORMANCE SUMMARY", 60)
    
    total_time = prep_time + magma_time + sage_time
    
    print(f"📈 TIMING BREAKDOWN:")
    print(f"   • Preparation:     {format_time(prep_time):>12}")
    print(f"   • MAGMA compute:   {format_time(magma_time):>12}")
    print(f"   • SAGE compute:    {format_time(sage_time):>12}")
    print(f"   • Total runtime:   {format_time(total_time):>12}")
    
    # Speed comparison
    if magma_time > 0 and sage_time > 0:
        if magma_time < sage_time:
            speedup = sage_time / magma_time
            print(f"\n🏆 MAGMA is {speedup:.2f}x faster than SAGE")
        else:
            speedup = magma_time / sage_time
            print(f"\n🏆 SAGE is {speedup:.2f}x faster than MAGMA")
    
    # Verify results match
    disc_magma = discriminant(C, Zbasis_L_Magma)
    disc_sage = discriminant(C, Zbasis_L_SAGE)
    lattices_equal = are_equals_lattices(C, Zbasis_L_Magma, Zbasis_L_SAGE)
    
    print(f"\n🔍 VERIFICATION:")
    print(f"   • Same discriminants: {disc_magma == disc_sage}")
    print(f"   • Equal lattices: {lattices_equal}")
    
    if disc_magma == disc_sage and lattices_equal:
        print(f"   ✅ Both methods found the SAME maximal order with discriminant {disc_magma}")
    elif disc_magma == disc_sage and not lattices_equal:
        print(f"   ⚠️  Same discriminant ({disc_magma}) but DIFFERENT lattices!")
        print(f"       This suggests both are maximal but not identical")
    elif lattices_equal and disc_magma != disc_sage:
        print(f"   🤔 Same lattices but different discriminants - this shouldn't happen!")
    else:
        print(f"   ❌ DIFFERENT results: MAGMA discriminant={disc_magma}, SAGE discriminant={disc_sage}")
        print(f"       Lattices are also different")
    
    print("\n" + "=" * 60)
    print(" COMPUTATION COMPLETED ")
    print("=" * 60)