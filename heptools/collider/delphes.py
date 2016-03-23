# -*- coding: utf-8 -*-
''' Utilities for using Delphes.

This module require for libDelphes.so to be in LD_LIBRARY_PATH

Examples:
    A code to load the first 1000 of the photons in the specified root file::

        from heptools.collider import delphes

        @delphes.activate
        def photon(i, photons):
            if i > 1000:
                delphes.deactivate(photon)

            for p1, p2 in itertools.combinations(photons, 2):
                do_something()

        delphes.load_root_file(SOME_ROOT_FILE)
'''

from ROOT import gSystem
gSystem.Load('libDelphes')

from ROOT import TChain, ExRootTreeReader

import inspect

DELPHES_BRANCH = 'Delphes'

delphes_name = ['Event', 'Particle', 'Track', 'Tower', 'EFlowTrack', 'EFlowPhoton',
                'EFlowNeutralHadron', 'GenJet', 'Jet', 'Electron', 'Photon', 'Muon',
                'MissingET', 'ScalarHT']

callback_functions = {}

def standardize(snake):
    '''This function returns the "standardized" name.

    Args: 
        snake (str): a name to standardize

    Returns:
        str: The standardized name

    Examples:
        >>> standardize('gen_particles')
        'GenParticle'
        >>> standardize('muons')
        'Muon'
    '''
    return ''.join(x.title() for x in snake.split('_'))[:-1]

def load_root_file(root_file):
    '''This function starts to load the specified root file.
    
    Each callbacks are called once this function is called.
    If the number of the callbacks becomes zero while loading,
    this function stops to load the root file anymore.

    Args:
        root_file (str): a path to the root file to load
    '''
    chain = TChain(DELPHES_BRANCH)
    chain.Add(root_file)
    tree_reader = ExRootTreeReader(chain)
    branches = {n: tree_reader.UseBranch(n) for n in delphes_name}

    for i in xrange(tree_reader.GetEntries()):
        if not callback_functions:
            return

        tree_reader.ReadEntry(i)

        callback_copy = {f: callback_functions[f] for f in callback_functions}

        for callback in callback_copy:
            branch_name = callback_functions[callback]
            kwargs = {b: branches[standardize(b)] for b in branch_name}
            callback(i, **kwargs)

def activate(f):
    '''A decorator to register a callback function.

    The first argument of each callbacks are expected to be an integer
    corresponding to the index of the events.
    The others are for the particle arrays. heptools determines the correspondance
    from the name of the arguments. See the example above.

    Args:
        f (callable): a function to wrap
    '''
    args = inspect.getargspec(f).args
    args = [a for a in args if standardize(a) in delphes_name]
    callback_functions[f] = args
    return f

def deactivate(f):
    '''Remove the function from the callback list.

    Args:
        f (callable): a function to remove
    '''
    del callback_functions[f]
