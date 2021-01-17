import random
import math
from sympy.ntheory import isprime
from ntt_utils import check_eq, reverse_bits, gen_omegas, get_omegas_inversed


def cooley_tukey_ntt_opt(inp, n, q, omegas):
    """Cooley-Tukey DIT algorithm with an extra optimization.
    We can avoid computing bit reversed order with each call by
    pre-computing the omegas in bit-reversed order.

    Requires:
     `omegas` are provied in bit-reversed order.
     `n` is a power of two.
     `q` is equivalent to `1 mod 2n`.

    Reference:
       https://www.microsoft.com/en-us/research/wp-content/
       uploads/2016/05/RLWE-1.pdf
    """
    a = inp.copy()

    assert q % (2 * n) == 1, f'{q} is not equivalent to 1 mod {2 * n}'
    assert (n & (n - 1) == 0) and n > 0, f'n: {n} is not a power of 2.'

    t = n
    m = 1
    while m < n:
        t >>= 1
        for i in range(0, m):
            j1 = i * (t << 1)
            j2 = j1 + t - 1
            S = omegas[m + i]
            for j in range(j1, j2 + 1):
                U = a[j]
                V = a[j + t] * S
                a[j] = (U + V) % q
                a[j + t] = (U - V) % q
        m <<= 1
    return a


def gentleman_sande_intt_opt(inp, n, q, inv_omegas):
    """Gentleman-Sande INTT butterfly algorithm.
    Assumes that inverse omegas are stored in bit-reversed order.
    Reference:
       https://www.microsoft.com/en-us/research/wp-content/
       uploads/2016/05/RLWE-1.pdf
    """
    a = inp.copy()
    t = 1
    m = n
    while (m > 1):
        j1 = 0
        h = m // 2
        for i in range(h):
            j2 = j1 + t - 1
            S = inv_omegas[h + i]
            for j in range(j1, j2 + 1):
                U = a[j]
                V = a[j + t]
                a[j] = (U + V) % q
                a[j + t] = ((U - V) * S) % q
            j1 = j1 + t * 2
        t *= 2
        m //= 2

    return [(i // n) % q for i in a]


def gen_bit_reversed_omegas(n, q):
    # Generate powers of omega in bit-reversed order.
    omegas = gen_omegas(n, q)

    for i in range(n):
        rev_i = reverse_bits(i, n.bit_length() - 1)
        if rev_i > i:
            omegas[i], omegas[rev_i] = omegas[rev_i], omegas[i]

    return omegas


def find_modulus(n, min_modulus):
    """
    We need to find a number `q` that
    satisfies the following properties:
    (1) q = 2 * n * k + 1, for some k.
    (2) q is a prime number.
    (3) q >= min_modulus > n
    (4) q is equivalent to 1 mod 2n.

    By Dirichlet's theorem, this is guaranteed to exist.
    """
    # Pre-condition for minimum modulus.
    assert min_modulus > n

    for k in range(1, n ** 2):
        q = k * 2 * n + 1  # (1)
        if not isprime(q) or q < min_modulus:
            # (2), (3)
            continue
        if q % (2 * n) != 1:
            # (4)
            continue
        return q


def round_trip_ntt(n, input_upper_bound):
    """Round trips through ntt and inverse ntt,
    verifying intt(ntt(a)) = a.

    `n` is the length of our input,
    and `input_upper_bound` is the maximum value in our input.
    """

    # Pick a sufficiently large minimum modulus to avoid overflow.
    min_modulus = input_upper_bound * n + 1

    a = [random.randint(0, input_upper_bound) for _ in range(n)]
    q = find_modulus(n, min_modulus)

    br_omegas = gen_bit_reversed_omegas(n, q)
    ntt_res = cooley_tukey_ntt_opt(a, n, q, br_omegas)

    inv_br_omegas = get_omegas_inversed(br_omegas, q)
    intt_res = gentleman_sande_intt_opt(ntt_res, n, q, inv_br_omegas)

    check_eq(a, intt_res)


if __name__ == '__main__':
    round_trip_ntt(n=1024, input_upper_bound=1000)
