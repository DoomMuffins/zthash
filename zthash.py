# coding=utf-8
"""
This package provides a research-oriented pure Python implementation of
ZT-hash (Zémor-Tillich), including utility functions and algorithms to
examine the properties of particular instantiations of ZT-hash.

Classes:

ZTHashParams -- The parameters (keys) of a particular function from the ZT-hash family.
ZTHash -- A stateful object providing a stream-based hash interface.
"""

import copy
from bitstring import Bits

_BITS = (Bits('0b0'), Bits('0b1'))
_MITM_DEFAULT_MAX_LENGTH = 32


class ZTHashParams(object):
    """
    Encapsulates the keys of a particular function from the ZT-hash family.
    """
    def __init__(self, generators):
        assert len(generators) == 2, 'Weird amount of generators given'
        # Copying possibly-mutable objects that belong to the caller
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
    """
    Provides a both a streaming and a single-call interface to ZT-hash calculations.
    """
    def __init__(self, zthash_params):
        self._params = zthash_params
        self._state = self._params.initial_value

    def update(self, data):
        """ Feed bits into the hash. """
        for bit in data:
            self._state *= self._params.generators[bit]

    def digest(self):
        """ Get the current digest. """
        state_copy = copy.copy(self._state)
        state_copy.set_immutable()
        return state_copy

    def reset(self):
        """ Reset the state for a new hash computation. """
        self._state = self._params.initial_value

    def compute(self, data):
        """ Compute the hash of a whole bitstring. """
        self.reset()
        self.update(data)
        return self.digest()


def mitm(zthash_params, end_states, except_bitstrings=(), max_length=_MITM_DEFAULT_MAX_LENGTH):
    """
    Perform a meet-in-the-middle attack trying to arrive at one of the desired end-states.
    :param tshash_params: The TSHash parameters.
    :param end_states: A dictionary of desired end states at which we want to arrive.
    :param except_bitstrings: Preimages to exclude.
    :param max_length: Maximum depth of BFS traversal across the graph.
    :return: A preimage, or None if no preimage was found.
    """

    empty = Bits()
    if zthash_params.initial_value in end_states and empty not in except_bitstrings:
        return empty
    
    # Starting conditions
    state_to_prefix = {zthash_params.initial_value: empty}
    state_to_suffix = {end_state: empty for end_state in end_states}
    advance_forward = True
    current_length = 0    

    while current_length <= max_length:

        if advance_forward:
            # Advance the prefixes (forwards)
            new_state_to_prefix = {}
            for state, prefix in state_to_prefix.iteritems():
                for bit in _BITS:
                    new_state = state * zthash_params.generators[bit[0]]
                    new_state.set_immutable()
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
                    new_state = state * zthash_params.inverse_generators[bit[0]]
                    new_state.set_immutable()
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

    return None


def mitm_preimage(zthash_params, digest, except_bitstrings=(), max_length=_MITM_DEFAULT_MAX_LENGTH):
    """ Performs a meet-in-the-middle attack to find a preimage for a digest. """
    return mitm(zthash_params=zthash_params,
                end_states=(digest, ),
                except_bitstrings=except_bitstrings,
                max_length=max_length)


def mitm_second_preimage(zthash_params, bitstring, except_bitstrings=(), max_length=_MITM_DEFAULT_MAX_LENGTH):
    """ Performs a meet-in-the-middle attack to find a second preimage to the image of a given bitstring. """
    digest = ZTHash(zthash_params).compute(bitstring)
    return mitm_preimage(zthash_params=zthash_params,
                         digest=digest,
                         except_bitstrings=(bitstring,) + except_bitstrings,
                         max_length=max_length)
