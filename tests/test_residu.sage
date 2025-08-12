# SageMath code using matplotlib for the final plot

# 1. Define the function using Sage's native capabilities
def f(N):
    """
    Returns the smallest prime p such that N is not a square mod p.
    """
    p = 3
    while legendre_symbol(N, p) != -1:
        p = next_prime(p)
    return p

# 2. Generate the data points cleanly
# For f(N)
M = 10000
points_f = [(N, f(N)) for N in range(2, M+1) if not is_square(N)]

# For the comparison curve
# Using .n() to ensure numerical values for plotting
points_curve = [(N, (3 * (log(N))^2).n()) for N in sxrange(2, M, 0.5)]

# 3. Unpack the data for matplotlib
# This is the standard way to separate coordinates for plotting
N_vals_f, y_vals_f = zip(*points_f)
N_vals_curve, y_vals_curve = zip(*points_curve)

# 4. Use matplotlib for plotting
import matplotlib.pyplot as plt

plt.figure(figsize=(12, 8))  # Set the figure size

# Plot the scatter points for f(N)
plt.scatter(N_vals_f, y_vals_f, color='blue', label='$f(N)$', s=20)

# Plot the smooth curve
plt.plot(N_vals_curve, y_vals_curve, color='red', label='$3(\\log N)^2$', linestyle='--')

# Add titles and labels for clarity
plt.title('Plot of $f(N)$ vs $3(\\log N)^2$ (Using Matplotlib)')
plt.xlabel('$N$')
plt.ylabel('$y$')
plt.legend()
plt.grid(True)

# Display the plot in the Sage worksheet
plt.show()