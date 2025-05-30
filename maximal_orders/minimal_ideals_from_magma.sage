load("maximal_orders/maximal_orders_utilities.sage")
load("utilities/algebra_type.sage")
load("utilities/utilities.sage")


#------------------------------------------------------------------------------------------
#
#            MINIMAL IDEALS FROM INDEMPOTENTS
#           using the magma function CentralIdempotents(A)
#-----------------------------------------------------------------------------------------


def central_orthogonal_idempotents_magma(A):
    """
    INPUT:
        - A: a finite-dimensional algebra over a finite field F_p.
    
    OUTPUT:
        - idempotents: list of primitive central orthogonal idempotents of A.
    """
    F = A.base_ring()
    dim_A = A.dimension()
    p = F.order()
    basis_A = list(A.basis())

    # Build structure constants as flat list L
    L = []
    for e in basis_A:
        for f in basis_A:
            x = e * f
            coords = get_coefficients(x, A)  # This function must return coordinates in the basis
            L.extend([int(c) for c in coords])

    # MAGMA code must be flat (no indentation!)
    magma_code = (
        f"F := GF({p});\n"
        f"L := [ {', '.join(map(str, L))} ];\n"
        f"A := AssociativeAlgebra<F, {dim_A} | L>;\n"
        f"Ids := CentralIdempotents(A);\n"
        f"[Eltseq(e) : e in Ids];\n"
    )


    # Call MAGMA (assumes `magma()` function is defined in your environment)
    magma_output = magma.eval(magma_code)

    # Parse MAGMA's output into a list of lists of integers
    import ast
    try:
        convert_output = ast.literal_eval(magma_output.strip())
    except Exception as e:
        raise ValueError(f"Failed to parse MAGMA output: {magma_output}") from e

    # Convert back to elements of A using the algebra's basis
    idempotents = []
    for coords in convert_output:
        e = sum(F(coords[i]) * basis_A[i] for i in range(dim_A))
        idempotents.append(e)

    return idempotents


def minimal_ideals_magma(A):
    basis_A = list(A.basis())
    F = A.base_ring()
    dim_A = dimension(A)

    CentralIdempotents = central_orthogonal_idempotents_magma(A)
    
    MinIdealsList =  []
    for e in CentralIdempotents :
        SetGeneratorsIdeal = [a*e for a in basis_A]

        ## Write vector of SetGeneratorsIdeal in each row of a matrix
        rows = []
        for x in SetGeneratorsIdeal:
            rows.append(get_coefficients(x,A))
        M = Matrix(F,rows)
        ## Extract free family of SetGeneratorsIdeal

        basis_Ai = []
        independent_vectors = M.row_space().basis_matrix()
        for v in independent_vectors:
            basis_Ai.append(sum(v[i] * basis_A[i] for i in range(dim_A)))
    
        MinIdealsList.append(basis_Ai)
    return MinIdealsList
