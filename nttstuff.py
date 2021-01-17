from sympy.discrete.transforms import ntt
import random
import math
import json
import argparse
from ntt_utils import check_eq, _data, find_prime, reverse_bits, gen_omegas


def naive_ntt(inp, P, omegas):
    """Simple NTT from Shunning."""
    N = len(inp)
    ret = [0] * N
    for i in range(N):
        for j in range(N):
            ret[i] = (ret[i] + inp[j] * omegas[(i * j) % N]) % P
    return ret


def cooley_tukey_ntt(inp, P, omegas):
    """Cooley-Tukey NTT algorithm."""
    ret = inp
    N = len(ret)
    bit_length = N.bit_length() - 1

    for i in range(N):
        rev_i = reverse_bits(i, bit_length)
        if rev_i > i:
            ret[i] ^= ret[rev_i]
            ret[rev_i] ^= ret[i]
            ret[i] ^= ret[rev_i]

    M = 2
    iters = int(math.log2(N))
    for _ in range(iters):
        for i in range(0, N, M):
            g = 0
            for j in range(0, M // 2):
                k = i + j + (M // 2)
                U = ret[i + j]
                V = ret[k] * omegas[g]
                ret[i + j] = (U + V) % P
                ret[k] = (U - V) % P
                g = g + N // M
        M = M * 2

    return ret


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
    cooley_tukey_res = cooley_tukey_ntt(a, p, omegas)

    check_eq(sympy_res, cooley_tukey_res)
    check_eq(sympy_res, naive_res)


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
