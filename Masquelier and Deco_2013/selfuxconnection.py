#from brian.connection.base import *
#from brian.connection import *
from brian import *

__all__ = [
         'SelfUXConnection'
         ]

class SelfUXConnection(Connection):
    '''
    
    Initialised with arguments:
    
    ``source``, NeuronGroup.
    ``depression``, boolean indicating if x should be updated
    ``facilitation``, boolean indicating if u should be updated
    ``delay``
        Only homogeneous delays are allowed.
    
    The benefit of this class is that it has no storage requirements and is optimised for
    this special case.
    '''
    @check_units(delay=second)
    def __init__(self, source, facilitation, depression, delay=0 * msecond, U=0.5):
        self.source = source # pointer to source group
        self.target = source # pointer to target group
        self.depression = depression
        self.facilitation = facilitation        
        source.set_max_delay(delay)
        self.delay = int(delay / source.clock.dt) # Synaptic delay in time bins
        if facilitation:
            self.n_u_f = self.target.get_var_index('u_f') # index of u_f
        if depression:
            self.n_x_d = self.target.get_var_index('x_d') # index of u_d
        self.U = U
    def propagate(self, spikes):
        '''
        Propagates the spikes to the target.
        '''
        if self.facilitation:
            self.target._S[self.n_u_f, spikes] += self.U * (1-self.target._S[self.n_u_f, spikes])
        if self.depression:
            if self.facilitation:
                self.target._S[self.n_x_d, spikes] *= (1-self.target._S[self.n_u_f, spikes])
            else:  
                self.target._S[self.n_x_d, spikes] *= (1-self.U)

    def compress(self):
        pass







