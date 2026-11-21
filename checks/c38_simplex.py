"""Canonical simplex symbols on actual C38 Fourier boundaries.

This adapter evaluates the inherited local predicates after decoding their
incident words. It adds no conversion edges. It is not a relation searcher.
"""
from copy import deepcopy
from math import comb
import arithmetic_alternating_chord as arithmetic
from simplex_contraction import Alphabet


class Codec:
    def __init__(self, t, modulus):
        self.alphabet = Alphabet(t, t)
        self.Q = comb(2*t, t)
        self.modulus = modulus
        self.width = arithmetic.digits_needed(modulus, self.Q)

    def encode(self, value):
        if not 0 <= value < self.Q**self.width:
            raise ValueError('word does not fit')
        output = []
        for _ in range(self.width):
            value, digit = divmod(value, self.Q)
            output.append(self.alphabet.unrank(digit))
        return output

    def rank(self, symbol):
        if not isinstance(symbol, (tuple, list)) or any(type(v) is not int for v in symbol):
            raise ValueError('not an integer simplex symbol')
        return self.alphabet.rank(symbol)

    def decode(self, word):
        if len(word) != self.width:
            raise ValueError('wrong word width')
        value = 0
        for symbol in reversed(word):
            value = value*self.Q + self.rank(symbol)
        if value >= self.modulus:
            raise ValueError('noncanonical Fourier word')
        return value

    def encode_bit(self, bit):
        if type(bit) is not int or bit not in (0, 1):
            raise ValueError('nonbinary state')
        return self.alphabet.unrank(bit)

    def decode_bit(self, symbol):
        bit = self.rank(symbol)
        if bit not in (0, 1):
            raise ValueError('nonbinary state symbol')
        return bit


def state_map(value, transform):
    if isinstance(value, dict):
        return {key: state_map(item, transform) for key, item in value.items()}
    if isinstance(value, list):
        return [state_map(item, transform) for item in value]
    return transform(value)


def encode_assignment(data, assignment, codec):
    output = deepcopy(assignment)
    for key, rows in assignment['hs'].items():
        output['hs'][key] = [[codec.encode(arithmetic.read_word(word, data['par']))
                              for word in row] for row in rows]
    output['outs'] = [[[codec.encode(arithmetic.read_word(word, data['par']))
                       for word in row] for row in rows] for rows in assignment['outs']]
    for field in ('quotients', 'carry_ports', 'range', 'nonzero'):
        output[field] = state_map(assignment[field], codec.encode_bit)
    return output


def decode_states(value, codec):
    # State containers use lists/dicts; every leaf joint symbol is a tuple.
    if isinstance(value, tuple):
        return codec.decode_bit(value)
    if isinstance(value, dict):
        return {key: decode_states(item, codec) for key, item in value.items()}
    if isinstance(value, list):
        return [decode_states(item, codec) for item in value]
    raise ValueError('malformed state container')


def decode_assignment(data, assignment, codec):
    output = deepcopy(assignment)
    for key, rows in assignment['hs'].items():
        output['hs'][key] = [[arithmetic.word(codec.decode(word), data['par'])
                              for word in row] for row in rows]
    output['outs'] = [[[arithmetic.word(codec.decode(word), data['par'])
                       for word in row] for row in rows] for rows in assignment['outs']]
    for field in ('quotients', 'carry_ports', 'range', 'nonzero'):
        output[field] = decode_states(assignment[field], codec)
    return output


def accept(data, assignment, codec):
    try:
        return arithmetic.accept(data, decode_assignment(data, assignment, codec))
    except (AssertionError, ValueError, TypeError, KeyError, IndexError):
        return False


def ledger(data, codec):
    """Enumerate actual endpoints; replace each full word, retain every side edge."""
    old = arithmetic.graph_incidence(data)
    old_fourier = old['categories']['five_symbol_fourier']
    assert old_fourier % data['par']['W5'] == 0
    words = old_fourier // data['par']['W5']
    auxiliary = old['digit_edges'] - old_fourier
    return dict(words=words, symbols_per_word=codec.width, alphabet_size=codec.Q,
                Fourier_symbol_edges=words*codec.width, auxiliary_symbol_edges=auxiliary,
                total_symbol_edges=words*codec.width+auxiliary,
                old_five_symbol_edges=old['digit_edges'],
                exact_array_count=8*data['m']-7+3*((data['m']-1)//2),
                endpoint_rule='hs[name][j2][t1]: F[j2]--M[t1]; outs[eq][j2][t1]: M[t1]--E[j2//R]')
