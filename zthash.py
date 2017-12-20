import copy
from bitstring import Bits

_BITS = (Bits('0b0'), Bits('0b1'))
_MITM_DEFAULT_MAX_LENGTH = 32


class ZTHashParams(object):
    def __init__(self, generators):
        assert len(generators) == 2, 'Weird amount of generators given'
        gen0 = copy.copy(generators[0])
        gen0.set_immutable()
        gen1 = copy.copy(generators[1])
        gen1.set_immutable()
        self.generators = (gen0, gen1)

        inverse_gen0 = gen0.inverse()
        inverse_gen0.set_immutable()
        inverse_gen1 = gen1.inverse()
        inverse_gen1.set_immutable()
        self.inverse_generators = (inverse_gen0, inverse_gen1)

        assert gen0.parent() == gen1.parent(), 'Generators don\'t belong to the same MatrixSpace'
        assert gen0.parent().dims() == (2, 2), 'Weird MatrixSpace dimensions'
        self.matrix_space = gen0.parent()
        self.initial_value = self.matrix_space.one()
        self.initial_value.set_immutable()

        ring = gen0.base_ring()
        assert ring.is_field(), 'Base ring isn\'t a field'
        self.field = ring


class ZTHash(object):
    def __init__(self, zthash_params):
        self._params = zthash_params
        self._state = self._params.initial_value

    def update(self, data):
        for bit in data:
            self._state *= self._params.generators[bit]

    def digest(self):
        state_copy = copy.copy(self._state)
        state_copy.set_immutable()
        return state_copy

    def reset(self):
        self._state = self._params.initial_value

    def compute(self, data):
        self.reset()
        self.update(data)
        return self.digest()


def mitm(zthash_params, state_to_suffix, except_bitstrings=(), max_length=_MITM_DEFAULT_MAX_LENGTH):
    preimage = state_to_suffix.get(zthash_params.initial_value)
    if preimage is not None:
        return preimage
    
    # Starting conditions
    state_to_prefix = {zthash_params.initial_value: Bits()}
    advance_forward = True
    current_length = 0    

    while current_length <= max_length:

        if advance_forward:
            # Advance the prefixes (forwards)
            new_state_to_prefix = {}
            for state, prefix in state_to_prefix.iteritems():
                for bit in _BITS:
                    new_state = state * zthash_params.generators[bit[0]]
                    new_prefix = prefix + bit
                    new_state_to_prefix[new_state] = new_prefix
                    
                    # Check whether we have an intersection
                    suffix = state_to_suffix.get(new_state)
                    if suffix is None:
                        continue
                    preimage = new_prefix + suffix
                    if preimage not in except_bitstrings:
                        return preimage
            state_to_prefix = new_state_to_prefix
            
        else:
            # Advance the suffixes (backwards)
            new_state_to_suffix = {}
            for state, suffix in state_to_suffix.iteritems():
                for bit in _BITS:
                    new_state = zthash_params.inverse_generators[bit[0]] * state
                    new_suffix = bit + suffix
                    new_state_to_suffix[new_state] = new_suffix

                    # Check whether we have an intersection
                    prefix = state_to_prefix.get(new_state)
                    if prefix is None:
                        continue
                    preimage = prefix + new_suffix
                    if preimage not in except_bitstrings:
                        return preimage
            state_to_suffix = new_state_to_suffix

        advance_forward = not advance_forward
        current_length += 1


def mitm_preimage(zthash_params, digest, except_bitstrings=(), max_length=_MITM_DEFAULT_MAX_LENGTH):
    state_to_suffix = {digest: Bits()}
    return mitm(zthash_params=zthash_params,
                state_to_suffix=state_to_suffix,
                except_bitstrings=except_bitstrings,
                max_length=max_length)


def mitm_second_preimage(zthash_params, bitstring, except_bitstrings=(), max_length=_MITM_DEFAULT_MAX_LENGTH):
    zthash = ZTHash(zthash_params)
    digest = zthash.compute(bitstring)
    return mitm_preimage(zthash_params=zthash_params,
                         digest=digest,
                         except_bitstrings=(bitstring,) + except_bitstrings,
                         max_length=max_length)
