// Corrected Benchmarking Script: Interleaved Trials

num_trials := 10;
factor_times := [];
iso_times := [];

printf "Running Interleaved Benchmarks...\n";

for i in [1..num_trials] do
    // --- Factorization Trial ---
    // Generate a unique number for this trial (using offset 2*i)
    p_factor := NextPrime(2^78 + 2*i*1000);
    d_factor := -(24^2 + p_factor*53^2 + p_factor*27^2);
    N_factor := p_factor^2 * d_factor;

    t_factor := Cputime();
    Factorization(N_factor);
    time_factor := Cputime(t_factor);
    Append(~factor_times, time_factor);
    printf "Trial %o (Factorization): %o seconds\n", i, time_factor;

    // --- Isotropic Subspace Trial ---
    // Generate a DIFFERENT unique number (using offset 2*i+1)
    p_iso := NextPrime(2^78 + (2*i+1)*1000);
    d_iso := -(24^2 + p_iso*53^2 + p_iso*27^2);
    M_iso := DiagonalMatrix(Rationals(), [1, p_iso, p_iso, d_iso]);

    t_iso := Cputime();
    IsotropicSubspace(M_iso);
    time_iso := Cputime(t_iso);
    Append(~iso_times, time_iso);
    printf "Trial %o (Isotropic Subspace): %o seconds\n", i, time_iso;
    
end for;


// --- Analyze the Results ---

avg_factor_time := &+factor_times / num_trials;
avg_iso_time := &+iso_times / num_trials;

printf "\n--- Results ---\n";
print "Average Factorization Time:", avg_factor_time;
print "Average Isotropic Subspace Time:", avg_iso_time;
print "Average Solver Overhead (IsoTime - FactorTime):", avg_iso_time - avg_factor_time;