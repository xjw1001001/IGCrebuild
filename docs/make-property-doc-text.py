"""
A helper script for writing the docs.

"""
from __future__ import print_function, division, absolute_import

import sys

from jsonctmctree.common_unpacking_ex import gen_valid_extended_properties


def main():

    print('.. automatically generated using the python script')
    print('..', ' '.join(sys.argv))
    print()
    print()

    reduction_names = (
            'observation_reduction',
            'edge_reduction',
            'state_reduction')

    response_names = (
            'observation index',
            'edge index',
            "'unraveled' state index",
            )

    s = '    '

    for extended_property in gen_valid_extended_properties():
        prefix = extended_property[:3]
        core_property = extended_property[-4:]
        observation_code, edge_code, state_code = prefix

        # define the lower case name
        # to be used for the extended property label and heading
        name = extended_property.lower()

        # write the link and a line break
        print('.. _' + name + ':')
        print()

        # write the section header
        print(name)
        print('^' * 7)
        print()

        # define the list of applicable reductions
        reductions = []
        for code, reduction in zip(prefix, reduction_names):
            if code == 'w':
                reductions.append(reduction)
        if core_property.lower() == 'tran':
            reductions.append('transition_reduction')

        # write the list of required custom reductions if any
        print('user-specified reductions:')
        print()
        if reductions:
            for reduction in reductions:
                print('* :ref:`%s`' % reduction)
        else:
            print('This extended property')
            print('is not associated with any user-specified reduction.')
        print()

        # define the shape of the response array
        responses = []
        for code, response in zip(prefix, response_names):
            if code == 'd':
                responses.append(response)
        if core_property.lower() == 'node':
            responses.append('node index')

        # say something about the interpretation of the output
        print('response shape:')
        print()
        if responses:
            for i, response in enumerate(responses):
                print('* axis %d of the response is the %s' % (i, response))
        else:
            print('The response associated with this extended property')
            print('is a single number rather than an array.')
        print()

        # extra line break for good measure
        print()

if __name__ == '__main__':
    main()
