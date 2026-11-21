"""Finite C38/simplex interface checks, separate from asymptotic claims."""
from collections import Counter
from copy import deepcopy
from itertools import product
import random
import common
import c38_simplex as bridge


def main():
    common.verify_sources()
    rng = random.Random(2026103101)
    counts = Counter()
    records = []
    for b, L, c in ((4, 16, 3), (7, 16, 27), (8, 32, 5)):
        for m in (2, 3, 4, 5, 6):
            data, values, factors = bridge.arithmetic.make_instance(b, L, c, m, rng, 'general')
            original = bridge.arithmetic.witness(data, values)
            assert bridge.arithmetic.accept(data, original)
            for t in (1, 2, 3, 4):
                codec = bridge.Codec(t, data['par']['modulus'])
                encoded = bridge.encode_assignment(data, original, codec)
                decoded = bridge.decode_assignment(data, encoded, codec)
                assert decoded == original
                assert bridge.accept(data, encoded, codec)
                assert bridge.arithmetic.recover_tuple(data, decoded) == factors
                ledger = bridge.ledger(data, codec)
                assert ledger['words'] == ledger['exact_array_count']*data['par']['L']
                records.append(dict(p=data['constants']['p'], m=m, t=t, **ledger))
                counts['actual_C38_simplex_instances'] += 1
                counts['Fourier_words_recoded'] += ledger['words']
                counts['local_first_predicates'] += data['par']['L2']
                counts['local_middle_predicates'] += data['par']['L1']
                counts['local_closure_predicates'] += data['par']['L2']//data['par']['R']
                for kind in ('out_of_domain', 'noncanonical', 'wrong_length', 'changed_output', 'bad_state'):
                    bad = deepcopy(encoded)
                    if kind == 'out_of_domain':
                        bad['hs']['x1'][0][0][0] = (t+1,)*t
                    elif kind == 'noncanonical':
                        bad['hs']['x1'][0][0] = codec.encode(codec.modulus)
                    elif kind == 'wrong_length':
                        bad['hs']['x1'][0][0].append(codec.alphabet.unrank(0))
                    elif kind == 'changed_output':
                        value = codec.decode(bad['outs'][0][0][0])
                        bad['outs'][0][0][0] = codec.encode((value+1) % codec.modulus)
                    elif codec.Q > 2:
                        boundary = min(bad['carry_ports'][0])
                        bad['carry_ports'][0][boundary][0] = codec.alphabet.unrank(2)
                    else:
                        continue
                    assert not bridge.accept(data, bad, codec), kind
                    counts['rejected_'+kind] += 1

    # Enumerate every encoded word for small domains, including invalid tails.
    for t in (1, 2, 3):
        for modulus in (3, 7, 17):
            codec = bridge.Codec(t, modulus)
            for digits in product(range(codec.Q), repeat=codec.width):
                word = [codec.alphabet.unrank(digit) for digit in digits]
                expected = sum(digit*codec.Q**i for i, digit in enumerate(digits))
                try:
                    got = codec.decode(word)
                except ValueError:
                    assert expected >= modulus
                    counts['exhaustive_noncanonical_rejections'] += 1
                else:
                    assert expected < modulus and got == expected
                    assert codec.encode(got) == word
                    counts['exhaustive_word_bijections'] += 1

    common.save('c38_simplex.json', dict(status='passed', counts=dict(counts), cases=records,
        scope='Actual C38 source/middle/closure predicates, canonical joint-simplex words, unchanged binary side states and recovered curve tuples. No complete C38 network contraction or asymptotic packing is executed.'))


if __name__ == '__main__':
    main()
