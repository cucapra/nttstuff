from sympy.ntheory import isprime, primitive_root
import math


def check_eq(a, b):
    """Assert that two arrays are equal everywhere."""
    for i, (x, y) in enumerate(zip(a, b)):
        assert x == y, f'difference at element {i}: {x}, {y}'
    print('ok!')


def gen_omegas(n, q):
    # Generate an omega: g^k (mod q) for a generator of the field, g.
    g = primitive_root(q)
    k = (q - 1) // n
    print(q, k)
    omega = (g ** k) % q
    assert 0 <= omega and omega < q

    # Generate pre-computed omega array (also from Shunning).
    omegas = [1]
    for i in range(n):
        omegas.append((omegas[i] * omega) % q)
    for i in range(n):
        assert (omegas[n - i] * omegas[i]) % q == 1
    omegas = omegas[:n]  # Drop the last, needless value.

    return omegas


def inversed(omegas, q):
    def multiplicative_inverse(a, q):
        # We want to find `i` such that `(x * i) mod q == 1`,
        # which which is guaranteed if `x` is a unit of the
        # multiplicative group under `q`.
        #
        # Reference:
        # https://en.wikipedia.org/wiki/Extended_Euclidean_algorithm
        prime = q
        Y = 0
        X = 1
        if a == 1:
            return 0
        while a > 1:
            quotient = a // q
            t = q
            q = a % q
            a = t
            t = Y
            Y = X - quotient * Y
            X = t
        return X if X >= 0 else X + prime

    return [multiplicative_inverse(x, q) for x in omegas]


def reverse_bits(number, bit_length):
    # Reverses the bits of `number` up to `bit_length`.
    reversed = 0
    for i in range(0, bit_length):
        if (number >> i) & 1: reversed |= 1 << (bit_length - 1 - i)
    return reversed


def find_prime(n):
    """Pick a prime number k*n+1 for some k."""
    for k in range(1, 100):
        p = k * n + 1
        if isprime(p):
            return p
    assert False, 'prime not found!!!1111one'


def _data(arr, bitwidth=32):
    """Make a FuTIL-ready JSON data dict."""
    return {'data': arr, 'bitwidth': bitwidth}
