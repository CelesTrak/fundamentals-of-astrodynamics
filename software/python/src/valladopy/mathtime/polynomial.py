import numpy as np


def quadric(a, b, c, option="include_imaginary"):
    """Solves for the roots of a quadratic equation ax^2 + bx + c = 0.

    Args:
        a (float): Coefficient of x^2
        b (float): Coefficient of x
        c (float): Constant term
        option (str): Determines the type of roots to return:
            'include_imaginary' - returns all roots including complex
            'real_only' - returns only real roots
            'unique_real' - returns unique real roots, filtering out duplicates

    Returns:
        tuple: Contains real and imaginary parts of two roots:
            - r1r (float): Real part of the first root
            - r1i (float): Imaginary part of the first root
            - r2r (float): Real part of the second root
            - r2i (float): Imaginary part of the second root
    """
    roots = np.roots([a, b, c])

    # Get roots based on option
    if option == "real_only":  # Only real roots
        roots = roots[np.isreal(roots)]
    elif option == "unique_real":  # Only unique real roots
        real_roots = roots[np.isreal(roots)]
        roots = np.unique(real_roots)

    # Ensure there are two roots for consistent output
    roots = np.pad(
        roots, (0, max(0, 2 - len(roots))), "constant", constant_values=np.nan
    )

    return roots[0].real, roots[0].imag, roots[1].real, roots[1].imag
