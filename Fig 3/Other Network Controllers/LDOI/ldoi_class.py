'''
used by ldoi_aug.py to help organize

TODO (2/2/2024)
- add_input_driven_row() should consider mutations/init even with inputs?
- generalize add_input_driven_row to any set of initially fixed nodes (ie from param mutations or init)
	- just use all on self.pins?
'''

import util, net

CUPY, cp = util.import_cp_or_np(try_cupy=1) #should import numpy as cp if cupy not installed

####################################################################################################

class Drivers:
	def __init__(self, G, params, pinned_nodes, inputs = None):
		# note renamed drivers -> pins
		assert(isinstance(G,net.ParityNet))
		G.prep_augmented_ldoi()
		# builds G.num_exp_clauses, intdtype, clauses2nodes, nodes2clauses, composite_counts
		
		self.dtype=bool

		self.G = G 
		self.params = params
		self.inputs = inputs

		self.n2= G.n_neg          
		self.n = G.n 
		assert(self.n*2 == self.n2)
		self.N = G.num_exp_clauses

		self.dims = 2 # haven't implemented higher orders
		self.shape = tuple([self.n2 for _ in range(self.dims)])
		self.pins = cp.zeros(self.shape, dtype=self.dtype) 
		self.base_drivers = cp.zeros(self.shape, dtype=self.dtype)

		self.output_indices = cp.array(self.G.output_indices())
		self.diag = cp.arange(self.n2)
		for i in range(self.dims-1):
			ind = [slice(None) for _ in range(self.dims) ] # unforunately not in gpu
			ind[i]=ind[-1]=self.diag
			self.pins[tuple(ind)] = 1
			self.base_drivers[tuple(ind)] = 1
		self.base_pins = cp.zeros(self.n2,dtype=self.dtype) 

		# note that universally pinned nodes (base_pins) override the diagonally pinned nodes
		for nodeNum in pinned_nodes:
			self.pins[...,nodeNum] = 1 # i think elipses here is correct...
			self.pins[...,(nodeNum+self.n) % self.n2] = 0 # override in case fixed in diag
			self.base_pins[nodeNum] = 1
		self.pins_compls = cp.roll(self.pins,self.n,axis=-1)

		#self.block = self.pins_compls & ~(cp.roll(self.base_drivers,self.n,axis=1)) # inputs and mutated nodes, but not drivers
		self.block = cp.roll(self.base_pins,self.n)

		self.L = cp.zeros(self.shape,dtype=self.dtype) 
		self.L_next = cp.zeros(self.shape,dtype=self.dtype) 


	def update_to_next(self):
		self.L = self.L_next.copy()

	def check_continue(self):
		
		prev = [self.L]
		nexxt = [self.L_next]
		for i in range(len(prev)):
			assert(cp.all(nexxt[i] | ~prev[i])) # prev \subseteq nexxt
			if not cp.all(prev[i]==nexxt[i]):
				return True # ie continue search 
		return False # done with search 

	def remove_base_pinned_nodes(self):
		self.L[:,self.base_pins] = 0
