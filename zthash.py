import copy
from bitstring import Bits

_BITS = (Bits('0b0'), Bits('0b1'))


class SL2HashParams(object):
    def __init__(self, generators):
        assert len(generators) == 2, 'Weird amount of generators given'
        gen0 = copy.copy(generators[0])
        gen0.set_immutable()
        gen1 = copy.copy(generators[1])
        gen1.set_immutable()
        self.generators = (gen0, gen1)

        assert gen0.parent() == gen1.parent(), 'Generators don\'t belong to the same MatrixSpace'
        assert gen0.parent().dims() == (2, 2), 'Weird MatrixSpace dimensions'
        self.matrix_space = gen0.parent()
        self.initial_value = self.matrix_space.one()
        self.initial_value.set_immutable()

        ring = gen0.base_ring()
        assert ring.is_field(), 'Base ring isn\'t a field'
        self.field = ring


class SL2Hash(object):
    def __init__(self, sl2hash_params):
        self._params = sl2hash_params
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


def mitm():
    pass

def mitm_preimage(sl2hash_params, digest):
    empty = Bits()
    if digest == sl2hash_params.initial_value:
        return empty

    # Starting conditions
    state_to_preimage_previous = {sl2hash_params.initial_value: empty}
    state_to_preimage_current = {sl2hash_params.generators[i]: _BITS[i] for i in (0, 1)}

    while True:
        state_to_preimage_next = {}
        for state, preimage in state_to_preimage_current.iteritems():
            other_half = state.inverse() * digest
            other_half.set_immutable()

            other_half_preimage = state_to_preimage_previous.get(other_half) or state_to_preimage_current.get(
                other_half)
            if other_half_preimage is not None:
                return preimage + other_half_preimage

            for bit in _BITS:
                new_state = state * sl2hash_params.generators[bit[0]]
                new_state.set_immutable()
                new_preimage = preimage + bit
                state_to_preimage_next[new_state] = new_preimage

        state_to_preimage_previous = state_to_preimage_current
        state_to_preimage_current = state_to_preimage_next


def mitm_second_preimage():
    pass