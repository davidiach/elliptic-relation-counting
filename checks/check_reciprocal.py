"""Finite checks of a quotient-free divisibility predicate for C38's curve equations."""
from math import gcd
from random import Random

import arithmetic_alternating_chord as chord


def cyclic_parameters(p, radix, length):
    power = radix**length
    sign = 1 if gcd(p, power - 1) == 1 else -1
    wrap = power - sign
    assert gcd(p, wrap) == 1 and wrap > 16*p*p
    return sign, wrap, pow(p, -1, wrap)


def closure(coefficients, radix, sign, initial_carry, p):
    """Check the wraparound carry and range, inferring each private output digit."""
    carry = initial_carry
    digits = []
    for coefficient in coefficients:
        carry, digit = divmod(coefficient + carry, radix)
        digits.append(digit)
    return initial_carry == sign*carry and chord.evaluate(digits, radix) < 16*p


def wrap_product(left, right, sign):
    length = len(left)
    folded = [0] * length
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            degree = i+j
            folded[degree % length] += a*b*(sign if degree >= length else 1)
    return folded


def scalar_examples():
    count = 0
    for p in (13, 17, 29, 31):
        radix, length = 4, 8
        sign, wrap, multiplier = cyclic_parameters(p, radix, length)
        digits = chord.encode(multiplier, length, radix)
        for value in range(16*p*p):
            coefficients = wrap_product(chord.encode(value, length, radix), digits, sign)
            raw = chord.evaluate(coefficients, radix)
            initial = sign*(raw//wrap)
            accepted = closure(coefficients, radix, sign, initial, p)
            assert accepted == (value % p == 0)
            for alternative in (initial-1, initial+1):
                assert not closure(coefficients, radix, sign, alternative, p)
            count += 1
    print(f"All {count:,} bounded scalar inputs agree with divisibility across both wrap signs.")


def coefficient_ring(radix, length):
    prime = 2*length+1
    while not chord.prime(prime):
        prime += 2*length
    nonresidue = next(a for a in range(2, prime) if pow(a, (prime-1)//2, prime) == prime-1)
    modulus = prime
    while modulus <= 256*length*length*radix**3:
        modulus *= prime
    phi = modulus//prime*(prime-1)
    root = pow(nonresidue, phi//(2*length), modulus)
    assert pow(root, length, modulus) == modulus-1
    return modulus, root


def transformed_equations(data, inputs, sign, multiplier):
    """Run the two Fourier stages on the source limbs and return signed closure coefficients."""
    par = dict(data["par"])
    length, radix, rows, columns = (par[key] for key in ("L", "B", "L2", "L1"))
    modulus, root = coefficient_ring(radix, length)
    omega, twist = root**2 % modulus, root if sign == -1 else 1
    par.update(modulus=modulus, omega=omega)
    row_root = pow(omega, columns, modulus)
    twiddles = [pow(twist, j, modulus) for j in range(length)]

    def boundary(digits):
        twisted = [digit*factor % modulus for digit, factor in zip(digits, twiddles)]
        return chord.H(twisted, par)

    sources = {name: boundary(chord.limbs(x, par)) for name, x in inputs.items()}
    public = {name: boundary(chord.limbs(x, par)) for name, x in data["constants"].items()}
    multiplier_boundary = boundary(chord.encode(multiplier, length, radix))
    equations = chord.arithmetic_count(data["m"])
    outputs = [[[0]*columns for _ in range(rows)] for _ in range(equations)]
    for column in range(columns):
        def second_stage(table):
            return chord.dft([row[column] for row in table], row_root, modulus)
        vals = {name: second_stage(table) for name, table in sources.items()}
        const = {name: second_stage(table) for name, table in public.items()}
        multiplier_values = second_stage(multiplier_boundary)
        frequencies = []
        for j in range(rows):
            at_source = {name: values[j] for name, values in vals.items()}
            at_public = {name: values[j] for name, values in const.items()}
            raw = chord.operations(at_source, at_public, data["m"], mode="special")[:equations]
            frequencies.append([(value+8*at_public["p"]**2)*multiplier_values[j] % modulus
                                for value in raw])
        for eq, table in enumerate(outputs):
            values = chord.idft([row[eq] for row in frequencies], row_root, modulus)
            for row, value in enumerate(values):
                table[row][column] = value
    coefficients = []
    for table in outputs:
        raw = [0]*length
        for row, boundary_words in enumerate(table):
            for j, value in enumerate(chord.H_inverse_row(boundary_words, row, par)):
                index = rows*j+row
                signed = value*pow(twiddles[index], -1, modulus) % modulus
                raw[index] = signed if 2*signed < modulus else signed-modulus
        coefficients.append(raw)
    return coefficients


def check_curve_equations():
    rng = Random(20261120)
    checked = positive = negative = 0
    for bits, length, offset, arity in ((4, 8, 3, 3), (5, 8, 15, 4), (7, 16, 27, 5)):
        data, valid, _ = chord.make_instance(bits, length, offset, arity, rng, "general")
        par, constants = data["par"], data["constants"]
        radix, p = par["B"], constants["p"]
        sign, wrap, multiplier = cyclic_parameters(p, radix, length)
        full_digits = chord.encode(multiplier, length, radix)
        source = {name: valid[name] for name in chord.source_names(arity, "special")}
        for trial in range(30):
            inputs = source if trial == 0 else {name: rng.randrange(p) for name in source}
            got = transformed_equations(data, inputs, sign, multiplier)
            unshifted = chord.operations(inputs, constants, arity, mode="special")
            in_polys = {name: chord.limbs(x, par) for name, x in inputs.items()}
            public_polys = {name: chord.limbs(x, par) for name, x in constants.items()}
            direct = chord.operations(in_polys, public_polys, arity, vector=True, mode="special")
            shift = chord.scale(8, chord.conv(public_polys["p"], public_polys["p"]))
            for eq, coefficients in enumerate(got):
                expected = wrap_product(chord.add(direct[eq], shift), full_digits, sign)
                assert coefficients == expected
                carry = sign*(chord.evaluate(coefficients, radix)//wrap)
                accepts = closure(coefficients, radix, sign, carry, p)
                assert accepts == (unshifted[eq] % p == 0)
                assert trial or accepts
                positive += int(accepts)
                negative += int(not accepts)
                checked += 1
        print(f"Field {p}, {arity} slots, radix {radix}, wrap sign {sign:+}: direct and Fourier coefficients agree.")
    assert positive and negative
    print(f"{checked:,} curve-equation checks: {positive} zero and {negative} nonzero residuals.")


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    scalar_examples()
    check_curve_equations()
