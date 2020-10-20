from sympy.discrete.transforms import ntt
from sympy.ntheory import isprime, primitive_root
import random
import json
import argparse


def find_prime(n):
    """Pick a prime number k*n+1 for some k."""
    for k in range(1, 100):
        p = k * n + 1
        if isprime(p):
            return p
    assert False, 'prime not found!!!1111one'


def naive_ntt(inp, P, omegas):
    """Simple NTT from Shunning."""
    N = len(inp)
    ret = [0] * N
    for i in range(N):
        for j in range(N):
            ret[i] = (ret[i] + inp[j] * omegas[(i * j) % N]) % P
    return ret


def _data(arr, bitwidth=32):
    """Make a FuTIL-ready JSON data dict."""
    return {'data': arr, 'bitwidth': bitwidth}


def check_eq(a, b):
    """Assert that two arrays are equal everywhere."""
    for i, (x, y) in enumerate(zip(a, b)):
        assert x == y, 'difference at element {}'.format(i)
    print('ok!')


def gen_omegas(n, p):
    # Generate an omega: g^k (mod p) for a generator of the field, g.
    g = primitive_root(p)
    k = (p - 1) // n
    omega = (g ** k) % p

    # Generate pre-computed omega array (also from Shunning).
    omegas = [1]
    for i in range(n):
        omegas.append(omegas[i] * omega % p)
    for i in range(n):
        assert omegas[n - i] * omegas[i] % p == 1
    omegas = omegas[:n]  # Drop the last, needless value.

    return omegas


def run_ntt(n, dump_input, dump_output, indata=None):
    if indata:
        p = indata['prime0']['data'][0]
        a = indata['inp0']['data']
    else:
        p = find_prime(n)
        a = [random.randint(0, 1000) for _ in range(n)]

    # SymPy's NTT generates omega stuff internally.
    sympy_res = ntt(a, prime=p)

    if indata:
        omegas = indata['omegas0']['data']
    else:
        omegas = gen_omegas(n, p)

    if dump_input:
        print(json.dumps({
            'inp0': _data(a),
            'prime0': _data([p]),
            'omegas0': _data(omegas),
            'ret0': sympy_res if dump_output else _data([0] * n),
        }, indent=2, sort_keys=True))
        return

    naive_res = naive_ntt(a, p, omegas)

    check_eq(naive_res, sympy_res)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('num', metavar='N', type=int, nargs='?',
                        default=2048, help='input size')
    parser.add_argument('-d', dest='dump', action='store_true', default=False,
                        help='dump inputs')
    parser.add_argument('-i', dest='input', default=None, help='input data')
    parser.add_argument('-o', dest='output', action='store_true',
                        default=False, help='dump results (requires -d)')
    args = parser.parse_args()

    if args.input:
        with open(args.input) as f:
            indata = json.load(f)
    else:
        indata = None

    run_ntt(args.num, args.dump, args.output, indata)


if __name__ == '__main__':
    main()
