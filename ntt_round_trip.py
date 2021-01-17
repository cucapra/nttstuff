import random
import math
from sympy.ntheory import isprime
from ntt_utils import check_eq, reverse_bits, gen_omegas, inversed


def cooley_tukey_ntt_opt(inp, n, q, phis):
    """Cooley-Tukey DIT algorithm with an extra optimization.
    We can avoid computing bit reversed order with each call by
    pre-computing the phis in bit-reversed order.

    Requires:
     `phis` are provided in bit-reversed order.
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
            S = phis[m + i]
            for j in range(j1, j2 + 1):
                U = a[j]
                V = a[j + t] * S
                a[j] = (U + V) % q
                a[j + t] = (U - V) % q
        m <<= 1
    return a


def gentleman_sande_intt_opt(inp, n, q, inv_phis):
    """Gentleman-Sande INTT butterfly algorithm.
    Assumes that inverse phis are stored in bit-reversed order.
    Reference:
       https://www.microsoft.com/en-us/research/wp-content/
       uploads/2016/05/RLWE-1.pdf
    """
    a = inp.copy()
    t = 1
    m = n
    while (m > 1):
        j1 = 0
        h = m >> 1
        for i in range(h):
            j2 = j1 + t - 1
            S = inv_phis[h + i]
            for j in range(j1, j2 + 1):
                U = a[j]
                V = a[j + t]
                a[j] = (U + V) % q
                a[j + t] = ((U - V) * S) % q
            j1 += (t << 1)
        t <<= 1
        m >>= 1

    return [(i // n) % q for i in a]


def get_bit_reversed(c, n, q):
    cc = c.copy()
    for i in range(n):
        rev_i = reverse_bits(i, n.bit_length() - 1)
        if rev_i > i:
            cc[i], cc[rev_i] = cc[rev_i], cc[i]

    return cc


def gen_phis(omegas, q):
    def legendre(x, q):
        return pow(x, (q - 1) // 2, q)

    def tonelli_shanks(x, q):
        # Finds the square-root of `x mod q`.
        # Source: https://rosettacode.org/wiki/Tonelli-Shanks_algorithm
        Q = q - 1
        s = 0
        while Q % 2 == 0:
            Q //= 2
            s += 1
        if s == 1:
            return pow(x, (q + 1) // 4, q)
        for z in range(2, q):
            if q - 1 == legendre(z, q):
                break
        c = pow(z, Q, q)
        r = pow(x, (Q + 1) // 2, q)
        t = pow(x, Q, q)
        m = s
        t2 = 0
        while (t - 1) % q != 0:
            t2 = (t * t) % q
            for i in range(1, m):
                if (t2 - 1) % q == 0:
                    break
                t2 = (t2 * t2) % q
            b = pow(c, 1 << (m - i - 1), q)
            r = (r * b) % q
            c = (b * b) % q
            t = (t * c) % q
            m = i
        return r

    return [tonelli_shanks(x, q) for x in omegas]


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

    for k in range(1, n ** n):
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
    min_modulus = input_upper_bound * n
    q = find_modulus(n, min_modulus)

    a = [random.randint(0, input_upper_bound) for _ in range(n)]

    omegas = gen_omegas(n, q)
    phis = gen_phis(omegas, q)
    br_phis = get_bit_reversed(phis, n, q)
    ntt_res = cooley_tukey_ntt_opt(a, n, q, br_phis)

    inv_phis = inversed(phis, q)
    inv_phis = get_bit_reversed(inv_phis, n, q)
    intt_res = gentleman_sande_intt_opt(ntt_res, n, q, inv_phis)

    check_eq(a, intt_res)


if __name__ == '__main__':
    round_trip_ntt(n=1024, input_upper_bound=1000)
