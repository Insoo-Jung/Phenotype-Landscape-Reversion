'''
contructs Net objects that include graph structure and the associated logic functions                        
last update 7/17/23 

lots of clean up todo, see below

'''

# 

import os, sys, yaml, util, math, itertools, pickle
import util
from copy import deepcopy

CUPY, cp = util.import_cp_or_np(try_cupy=1) #should import numpy as cp if cupy not installed
# in particular Fmapd should be cupy (if cupy is being used instead of numpy)


# POSS reorganization:
    # pass a F dict + node names? which can come from file OR code (2 branches)
    # see rdze.build_net_with_logic() for what a mess curr is
    # 2 step init? init is blank, then need to load via file, prev net, or F+nodes
        # unless building from prev G, need self.build_negative_nodes() 

# TODO: 
#   add explicit debug function for net construction
#       check that certain objs same, certain diff
#       ex G.F and node.F()
#   Gpar should have all parity nodes in G.nodes too..but worried will break other things
#       like with Gpar.n vs .n_neg problem is on the one had want to be able to Gpar in place of G without changing anything
#           on the other hand want Gpar to be a proper network where !x is a sep node altogether...(why?)

# TODO LATER:
#   node.F should be composed of other node object, rather than names jp
#   util: custom print that auto states the file/fn calling it
#   add time-reverse net (and inheritance with it)
#   Net._build_net() and build_Fmap() are still messy af
#   check for repeated clauses in build_Fmapd_and_A()


class Net:
    def __init__(self, params, model_file=None, G=None, init_for_sim=True, read_file_directly=False, reduce_net=False):
        # will construct from graph G if passed
        # otherwise will construct from the default model_file in params (unless otherwise specified)
        # generally use init_for_sim, unless prepare() will be called later

        if params is None:
            #print("WARNING net.py: building net with default params (ok for network reduction and parity construction)")
            params = {'debug':True,'model_file':model_file} 

        if model_file is None:
            model_file=params['model_file']

        self.params = params
        self.F= {} # logic function for each node
        self.F_orig = {} # back up to restore after a mutation 
        self.A=None # adjacency matrix, unsigned
        self.A_signed=None # signed adj matrix, not used by ParityNet
        self.A_mf = None # mean field expected pr that a node activates a neighbor (nodes + parity nodes, but regular net)
        self.Fmapd = None

        self.n = 0  # # regular nodes
        self.n_neg = 0 # num regular nodes + num negative nodes

        self.nodes = [] # object to hold name & number, only regular nodes included
        self.parityNodes = [] #includes complement nodes too, formerly "allNodes"

        self.nodeNames = [] # num -> name
        self.nodeNums = {} # name -> num

        self.encoding = 'bnet' # encoding format used by the input file, note bnet is default
        self.not_string = '!' # from encoding, called often enough that it is useful to store, again bnet by default
        
        self.num_clauses=0 # total number of clauses
        self.max_literals=0 # max literals per clause
        self.max_clauses=0 # max clauses per node

        self.debug = params['debug']
        self.was_reduced = False

        # CONSTRUCTION FROM PREVIOUS NET
        if G is not None: 
            self.copy_from_net(G)

        # messy...
        if (init_for_sim and G is None) or read_file_directly:
            # CONSTRUCTION FROM FILE
            self.read_from_file(model_file,reduce_net=reduce_net)


        # APPLY MUTATIONS AND BUILD Fmapd
        if init_for_sim:

            self.prepare(params)


            if util.istrue(params,['multicellular','active']):
                if util.istrue(params,['multicellular','external']):
                    self.external = cp.array(params['multicellular']['external']).astype(bool)
                else:
                    if 'external_pr' in params['multicellular']:
                        p_ext = params['multicellular']['external_pr']
                    else:
                        p_ext = .5
                    self.external = cp.random.choice([0,1],p=[1-p_ext,p_ext],size=(self.n,)).astype(bool)
                print("net.py: external=",self.external.astype(int))
                # for now, which nodes are external is just random

        # messy
        #params['input_indices'] = self.input_indices()
        #params['output_indices'] = self.output_indices()

    def __str__(self):
        # TODO: make this more complete
        return "Net:\nF =" + str(self.F)

    def prepare(self,params): 
        # applies setting mutations and builds Fmapd
        self.params=params # since sometimes alter the params applied to same net here
        self.restore_from_prev_mutations()
        self.add_self_loops()
        self.check_mutations()
        self.apply_mutations()
        self.build_Fmapd_and_A()

    def mutate(self,mutant,value):
        # mutates and returns copies to avoid altering the original
        # mutant is name, value is int

        Gcopy = deepcopy(self)
        params_copy = Gcopy.params  
        params_copy['mutations'][mutant] = value
        if isinstance(self,ParityNet):
            if self.not_string in mutant: 
                params_copy[mutant.replace(self.not_str,'')] = (value+1)%2
            else:
                params_copy[self.not_str + mutant] = (value+1)%2

        Gcopy.prepare(params_copy)

        #assert(0) # explicitly debug
        return Gcopy, params_copy

    def nonIO_nodes(self):
        io_indices =  self.input_indices() + self.output_indices()
        for node in self.nodes:
            if node.num not in io_indices:
                yield node

    def copy_from_net(self,G):
        self.__dict__ = deepcopy(G.__dict__) # blasphemy to some

    def add_self_loops(self):
        for node in self.nodes:
            # if function is 0 or 1 set to a self-loop and an init call
            if self.F[node.name] in [[['0']],[['1']]]:
                self.params['init'][node.name] = int(self.F[node.name][0][0])
                self.F[node.name] = [[node.name]]
            if node.name in self.params['inputs'] and util.istrue(self.params,'pin_inputs'):
                self.F[node.name] = [[node.name]]

    def _read_node_name(self,line_parts, node_fn_split, strip_from_node):
        node_name = line_parts[0].strip().replace(' ','').replace('\t','')
        for symbol in strip_from_node:
            node_name = node_name.replace(symbol,'')
        if self.debug:
            #print('net.py:',node_name)
            assert(node_name not in self.nodeNames)
        return node_name 

    def read_from_file(self, net_file, reduce_net=False):
        # inits network, except for F_mapd (which maps the logic to an actual matrix execution)
        # net file should be in DNF, see README for specifications

        if not os.path.isfile(net_file):
            sys.exit("Can't find network file: " + str(net_file)) 
        
        if reduce_net:
            # import statements here suck...but don't want someone to have to install pyeda if not used
            import espresso

        with open(net_file,'r') as file:
            format_name=self._get_encoding(file,net_file)
            node_fn_split, clause_split, literal_split, not_str, strip_from_clause, strip_from_node = get_file_format(format_name)

            self.encoding = format_name
            self.not_string = not_str 

            loop = 0
            while True:
                line = self._read_line(file, loop)
                if (line is None) or (line.isspace()): #i.e eof
                    break 
                line = line.strip().split(node_fn_split)
                node_name = self._read_node_name(line, node_fn_split, strip_from_node)
                self.add_node(node_name)
                if (len(line[1])==0) or (line[1].isspace()):
                    self.F[node_name] = [[node_name]]
                elif reduce_net:
                    #print("\nnode name",node_name,'\nfn=',line[1])
                    self.F[node_name] = espresso.reduce_nonDNF_string(line[1])
                    unreduced_fn = espresso.reformat_nonDNF_string(line[1])
                    self.was_reduced = self.check_if_reduced(unreduced_fn, self.F[node_name],node_name)
                else:
                    self.F[node_name] = self.read_in_function(line[1],format_name)

                loop += 1

        self.add_hidden_inputs() # any nodes without a function are assumed to be inputs (with a self-loop)
        self.build_negative_nodes() # this is just a pass function for Parity nets

    def read_in_function(self, line_part, format_name):
        fn = []
        node_fn_split, clause_split, literal_split, not_str, strip_from_clause, strip_from_node = get_file_format(format_name)

        clauses = line_part.split(clause_split)
        for clause in clauses:
            this_clause=[]
            for symbol in strip_from_clause:
                clause = clause.replace(symbol,'')

            literals = clause.split(literal_split)
            for j in range(len(literals)):
                literal_name = literals[j]
                for symbol in strip_from_node:
                    literal_name = literal_name.replace(symbol,'')
                literal_name = literal_name.replace(' ','').replace('\t','')
                this_clause += [literal_name]
            
            fn += [this_clause]
        return fn

    def check_if_reduced(self, unreduced_fn, reduced_fn, node_name):
        orig_num_ele = sum([len(s) for s in unreduced_fn])
        reduced_num_ele = sum([len(s) for s in reduced_fn])
        node_fn_was_reduced = (orig_num_ele != reduced_num_ele)
        was_reduced = node_fn_was_reduced | self.was_reduced
        if node_fn_was_reduced:
            print("\tfunction for",node_name,"was reduced:\n",unreduced_fn,'-->',reduced_fn)
        return was_reduced

    def add_hidden_inputs(self):
        to_add = []
        for node in self.nodeNames:
            for clause in self.F[node]:
                for input_node in clause:
                    input_node = input_node.replace('!','')
                    if (input_node not in self.nodeNames) and (input_node not in to_add):
                        to_add += [input_node]
        for node in to_add:
            #print('adding node',node)
            self.add_node(node)
            self.F[node] = [[node]] #ie a self loop

    def write_to_file(self, output_file):

        with open(output_file,'w') as ofile:
            if self.encoding != 'bnet':
                assert(0) # generally not supported
                ofile.write(self.encoding + '\n')
            else:
                pass # write it in case does not exist
            # note this overwrites
        node_fn_split, clause_split, literal_split, not_str, strip_from_clause, strip_from_node = get_file_format(self.encoding)
        if self.encoding == 'bnet':
            # these aren't necessary, but they do look nicer
            node_fn_split = ',\t'

        if isinstance(self,ParityNet):
            nodes = self.parityNodes
        else:
            nodes = self.nodes
        i=0
        for node in nodes:
            fn=node.F()
            with open(output_file,'a') as ofile:
                if i>0:
                    ofile.write('\n')
                ofile.write(node.name + node_fn_split + self._fn_to_str(fn))
                i+=1


    def _fn_to_str(self,fn):
        # fn = F['nodeName'] = [clause1,clause2,...] where clause = [lit1,lit2,...] 
        # input_nodes are names

        # maybe should store all of these in net obj itself? ex G.format.not_str
        node_fn_split, clause_split, literal_split, not_str, strip_from_clause, strip_from_node = get_file_format(self.encoding)
        if self.encoding == 'bnet':
            # these aren't necessary, but they do look nicer
            clause_split, literal_split = ' | ', ' & '

        fn_str = ''
        i=0
        for clause in fn:
            assert(len(clause)>0)
            if i!=0:
                fn_str += clause_split
            if len(strip_from_clause) > 0 and len(clause)>1:
                fn_str += strip_from_clause[0]
            j=0
            for lit in clause:
                if j!=0:
                    fn_str += literal_split
                fn_str += lit
                j+=1 
            if len(strip_from_clause) > 0 and len(clause)>1:
                fn_str += strip_from_clause[1]
            i+=1
        return fn_str

    
    def add_node(self,nodeName,isNegative=False):
        # note that F is not added for negatives, since simulation just negates their positive nodes directly

        newNode = Node(self,nodeName,self.n_neg,isNegative=isNegative)
        self.nodeNames += [nodeName]
        self.nodeNums[nodeName] = self.n_neg
        self.parityNodes += [newNode] 
        self.n_neg += 1
        
        if not isNegative:
            self.F[nodeName] = []
            self.nodes += [newNode]
            self.n+=1

    def rm_node(self,node):
        for node2 in self.parityNodes:
            k = node2.name
            if self.nodeNums[k] > node.num:
                self.nodeNums[k] -= 1 
                node2.num -= 1

        if not node.isNegative:
            if node not in self.nodes:
                print('name=',node.name,'not in',[n for n in node.name])
            self.nodes.remove(node)
            self.n -= 1
        del self.nodeNums[node.name]
        self.nodeNames.remove(node.name)
        self.parityNodes.remove(node)

        self.n_neg -= 1
        del self.F[node.name]
        del node

    def nodesByName(self,name):
        return self.nodes[self.nodeNums[name]]

    def build_negative_nodes(self):
        # put the negative nodes as 2nd half of node array
        if not isinstance(self,ParityNet):
            for i in range(self.n):
                self.add_node(self.not_string + self.nodeNames[i],isNegative=True)

    def check_mutations(self): 
        if 'mutations' in self.params.keys():
            for k in self.params['mutations'].keys():
                if k not in self.nodeNames:
                    sys.exit("\nSetting file specifies mutation on " + str(k) +", but this node is not in the network!\n") 

    def restore_from_prev_mutations(self):
        for node in self.nodeNames:
            if node in self.F_orig:
                self.F[node] = deepcopy(self.F_orig[node])

    def apply_mutations(self):
        if 'mutations' in self.params.keys() and len(self.params['mutations']) > 0:
            for node in self.nodeNames:
                if node in self.params['mutations']:
                    self.F_orig[node] = deepcopy(self.F[node])
                    lng = len(self.F[node]) 
                    # will be rendundant, but avoids possible issues with changing # of clauses
                    # otherwise need to rebuild Fmapd each time
                    self.F[node] = [[node] for _ in range(lng)]
                    self.params['init'][node] = int(self.params['mutations'][node])
                    if self.debug:
                        assert(self.params['init'][node] in [0,1]) 



    def build_Fmapd_and_A(self): 
        
        # Fmapd is a clause index with which to compress the clauses -> nodes
        #   nodes_to_clauses: computes which nodes are inputs (literals) of which clauses
        #   clauses_to_threads: maps which clauses will be processed by which threads (and #threads = #nodes)
        #   threads_to_nodes: maps which threads are functions of which nodes
        
        if util.istrue(self.params,'load_logic_pickle'):
            print("\nnet.build_Fmapd_and_A(): LOADING NET FROM PICKLE FILE\n")
            # load pickle instead of recalculating
            with open(self.params['load_logic_pickle'], 'rb') as f:
                data = pickle.load(f)
                self.Fmapd, self.A, self.A_signed, self.A_mf = data['Fmapd'], data['A'], data['A_signed'], data['A_mf']
                self.num_clauses, self.max_literals, self.max_clauses = data['num_clauses'], data['max_literals'], data['max_clauses']
            return

        n = self.n

        nodes_to_clauses = []        
        clauses_to_threads = []
        threads_to_nodes = [] 
        nodes_clause = {i:[] for i in range(n)}


        # ADJACENCY MATRIX, unsigned
        self.A = cp.zeros((n,n),dtype=bool)
        if not isinstance(self,ParityNet):
            self.A_signed = cp.zeros((n,n),dtype=cp.int8)
        self.num_clauses = 0

        # BUILDING NODES->CLAUSES
        for i in range(n):
            numbered_clauses = [self._clause_num_and_add_to_A(i,c) for c in self.nodes[i].F()]
            self.max_clauses = max(self.max_clauses, len(numbered_clauses))
            nodes_to_clauses += numbered_clauses
            nodes_clause[i] = [z for z in range(self.num_clauses, self.num_clauses+len(numbered_clauses))]
            self.num_clauses += len(numbered_clauses)

        self._make_square_clauses(nodes_to_clauses) 
        nodes_to_clauses = cp.array(nodes_to_clauses,dtype=self._get_index_dtype(self.n)) # the literals in each clause

        if self.num_clauses>0:
            assert(nodes_to_clauses.shape == (self.num_clauses,self.max_literals))

        # BUILDING CLAUSES->THREADS
        m=min(self.params['clause_bin_size'],self.max_clauses) 

        i=0
        while sum([len(nodes_clause[i2]) for i2 in range(n)]) > 0:
            # ie exit when have mapped all clauses

            # after each clauses_to_threads[i], need to add an index for where the new (compressed) clauses will go
            this_set=[[] for _ in range(n)]
            threads_to_nodes += [cp.zeros((n,n),dtype=bool)]
            sorted_keys = sorted(nodes_clause, key=lambda k: len(nodes_clause[k]), reverse=True)
            nodes_clause = {k:nodes_clause[k] for k in sorted_keys}
            node_indx = 0
            prev_take_from_node = take_from_node = sorted_keys[node_indx]

            for j in range(n):
                if sum([len(nodes_clause[i2]) for i2 in range(n)]) > 0:
                    take_from_node = sorted_keys[node_indx]
                    threads_to_nodes[i][j,take_from_node]=1 
                    if len(nodes_clause[take_from_node]) >= m:
                        this_set[j] = nodes_clause[take_from_node][:m]
                        if len(nodes_clause[take_from_node]) == m:
                            node_indx += 1
                        del nodes_clause[take_from_node][:m]
                    else:
                        top = len(nodes_clause[take_from_node])
                        this_set[j] = nodes_clause[take_from_node][:top] # why is [:top] nec?
                        rem = m-(top)
                        for k in range(rem):
                            this_set[j] += [this_set[j][-1]] #use a copy of prev clause to make matrix square (assuming DNF)
 
                        del nodes_clause[take_from_node][:top]
                        node_indx += 1

                else: #finished, just need to filll up the array
                    threads_to_nodes[i][j,prev_take_from_node]=1 
                    this_set[j] = this_set[j-1]

                prev_take_from_node = take_from_node

            clauses_to_threads += [this_set]    
            i+=1
            if i>1000000:
                sys.exit("ERROR: infinite loop likely in net.build_Fmapd_and_A()")
        
        if self.params['num_samples']<254 and self.num_clauses<254:
            thread_dtype = cp.uint8
        elif self.params['num_samples']<65533 and self.num_clauses<65533:
            thread_dtype = cp.uint16
        else:
            thread_dtype = cp.uint32
        clauses_to_threads = cp.array(clauses_to_threads,dtype=thread_dtype)
        threads_to_nodes = cp.array(threads_to_nodes,dtype=bool)
        # nodes_to_clauses already converted
        clause_mapping = {'nodes_to_clauses':nodes_to_clauses, 'clauses_to_threads':clauses_to_threads, 'threads_to_nodes':threads_to_nodes}

        #print('nodes->clauses\n',nodes_to_clauses,'\n\nclauses->threads\n',clauses_to_threads,'\n\nthreads->nodes\n',threads_to_nodes)
        self.Fmapd = clause_mapping

        if util.istrue(self.params,['save_logic_pickle']):
            pickle_file = self.params['model_file'].split('.')
            pickle_file = pickle_file[0] + 'Logic.pickle'
            data = {} 
            data['Fmapd'], data['A'], data['A_signed'], data['A_mf'] = self.Fmapd, self.A, self.A_signed, self.A_mf 
            data['num_clauses'], data['max_literals'], data['max_clauses'] = self.num_clauses, self.max_literals, self.max_clauses
            util.pickle_it(data, pickle_file)

    

    def to_networkx(self):
        import networkx as nx # move higher up if called more

        Gnx = nx.DiGraph()
        for node in self.nodes:
            Gnx.add_node(node.name)

        # get edges from functions
        for node in self.nodes:
            for clause in node.F():
                for source in clause:
                    source = source.replace('!','')
                    Gnx.add_edge(source, node.name) # redundant edges will be ignored anyways

        return Gnx  


    def _clause_num_and_add_to_A(self, source_node_num, clause):            
        clause_fn = [] 
        self.max_literals = max(self.max_literals, len(clause))
        
        for k in range(len(clause)):
            
            literal_node = self.nodeNums[clause[k]]
            clause_fn += [literal_node]

            if self.parityNodes[literal_node].isNegative:  
                literal_node-=self.n 
            self.A[source_node_num, literal_node] = 1

            if not isinstance(self,ParityNet):
                if self.not_string in clause[k]:
                    sign=-1
                else:
                    sign=1
                self.A_signed[source_node_num, literal_node] = sign

        return clause_fn

    def _make_square_clauses(self, nodes_to_clauses):
        for clause in nodes_to_clauses:
            for k in range( len(clause), self.max_literals ):
                clause += [clause[-1]]

    def _get_index_dtype(self, max_n):
        # redundant w a function in util
        if max_n<256: 
            return cp.uint8
        elif max_n<65536:
            return cp.uint16
        else:
            return cp.uint32

    def _read_line(self,file,loop):
        line = file.readline()
        if loop > 1000000:
                sys.exit("Hit an infinite loop, unless net is monstrously huge") 
        if line==None or len(line)==0:
            return None 
        return line


    def _get_encoding(self,file,net_file):
        extension = net_file.split('.')
        if extension[-1] == 'bnet':
            return 'bnet'
        else:
            return file.readline().replace('\n','')


    def input_indices(self):
        return [self.nodeNums[self.params['inputs'][i]] for i in range(len(self.params['inputs']))]
    
    def output_indices(self):
        return [self.nodeNums[self.params['outputs'][i]] for i in range(len(self.params['outputs']))]
    
    def mutant_indices(self):
        return [self.nodeNums[k] for k in self.params['mutations'].keys()]

    def get_input_sets(self):
        # assumes that 2^#inputs can fit in memory
        #input_indices = self.input_indices(params)
        return list(itertools.product([0,1],repeat=len(self.params['inputs'])))

    def print_matrix_names(self,X):
        # assumes X is a set of network states
        print('\nnetwork state:')
        for x in X: # one network state
            s=''
            for i in range(len(x)):
                s+=self.nodeNames[i] + ":" + str(x[i]) + " "
            print(s)
        print('\n')

    def node_vector_to_names(self,v):
        # where v is a list of node indices
        names = []
        for i in v:
            names+=[self.nodeNames[i]]
        return names

    def node_vector_bool_to_names(self,v):
        # where v[i]=1 implies node i should be included
        names = []
        for i in range(len(v)):
            if v[i]:
                names+=[self.nodeNames[i]]
        return names

    def node_matrix_bool_to_names(self, X):
        return [self.node_vector_bool_to_names(X[i]) for i in range(len(X))]

    def certain_nodes_to_names(self,s):
        # where s is a list of node states
        names = []
        for i in range(len(s)):
            if s[i] != 2:
                assert(s[i] in [0,1])
                names += [self.nodeNames[i]+'='+str(int(s[i]))]
        return names

    def build_parity(self, init_for_sim=True):
        # builds ParityNet from a regular net
        assert(not isinstance(self, ParityNet))
        G2 = ParityNet(self.params, G=self,init_for_sim=init_for_sim)
        if init_for_sim:
            G2.prepare(G2.params)
        return G2

    def build_time_reversal(self):
        # alter the logic functions then reprep for sim
        if (not isinstance(self, ParityNet)):
            G1 = self.build_parity()
        else:
            G1 = self

        G2 = deepcopy(G1)
        for node in G2.parityNodes: 
            if self.not_string in node.name: 
                compl = node.name.replace(self.not_string,'')
            else:
                compl = self.not_string + node.name 

            G2.F[node.name] = deepcopy(G1.F[compl])

            # swap any self loops
            for i in range(len(G2.F[node.name])):
                clause = G2.F[node.name][i]
                if node.name in G2.F[node.name][i]:
                    clause = ['compl_placeholder' if x==node.name else x for x in clause]
                if compl in clause: 
                    clause = [node.name if x==compl else x for x in clause]
                G2.F[node.name][i] = [compl if x=='compl_placeholder' else x for x in clause]

            assert(G2.F[node.name]==node.F())
        G2.build_Aexp() # need to rebuild
        G2.prepare(G2.params) # need to rebuild
        return G2 

##################################################################################################


class Node:
    def __init__(self,G,name,num,isNegative=False):
        self.name=name 
        self.num=num 
        self.isNegative = isNegative 
        self.G=G

    def F(self):
        #if self.name not in self.G.F:
        #   print("\nERROR: node ",self.name," not assigned a function in Net")
        #   sys.exit()
        return self.G.F[self.name]  # will raise error if node name not assigned in Net's F

    def setF(self,new_fn):
        self.G.F[self.name] = new_fn



##################################################################################################

class ParityNet(Net):
    def __init__(self,params, G=None, model_file=None, init_for_sim=True, read_file_directly=False):
    
        if model_file is None and G is None:
            model_file=params['parity_model_file'] # i.e. default

        self.parityNodes = [] # these include all regular nodes and their complements
        self.debug = params['debug']
        super().__init__(params, G=G, model_file=model_file,init_for_sim=init_for_sim, read_file_directly=read_file_directly)
        if G is not None:
            self.build_parity_functions(G)
        #    self.parityNodes = self.nodes # this seems very stupid, rmd and prob broke some other shit
        # assumes that negative nodes are already in parity form (i.e. have their logic)
        # note that composite nodes aren't formally added, since only A_exp is actually used (hence parity, not expanded net)

        if init_for_sim:
            self.build_Aexp()
            self.build_counts()

    def prep_augmented_ldoi(self):
        # clauses are DNF
        clauses = []
        maxlng = 0
        for node in self.parityNodes:
            for clause in node.F():
                if clause not in clauses:
                    clauses += [clause]
                    maxlng = max(maxlng,len(clause))
        assert(self.n_neg == len(self.parityNodes))

        self.num_exp_clauses = len(clauses)
        self.intdtype = util.get_int_dtype(maxlng)
        self.clauses2nodes = cp.zeros((self.num_exp_clauses, self.n_neg),dtype=bool)
        self.nodes2clauses = cp.zeros((self.n_neg, self.num_exp_clauses),dtype=self.intdtype)

        for i in range(self.num_exp_clauses):
            for j in range(self.n_neg):
                if clauses[i] in self.parityNodes[j].F():
                    self.clauses2nodes[i,j] = 1
                if self.parityNodes[j].name in clauses[i]:
                    self.nodes2clauses[j,i] = 1

        self.composite_counts = cp.sum(self.nodes2clauses,axis=0)
        '''
        print("nodes=",self.nodeNames)
        print("\nc2n\n",self.clauses2nodes.astype(int))
        print("\nn2c\n",self.nodes2clauses.astype(int))
        print("\ncounts=\n",self.composite_counts)
        '''


    def add_node(self,nodeName,isNegative=False):
        if self.not_string in nodeName:
            isNegative=True
        newNode = Node(self,nodeName,self.n_neg,isNegative=isNegative)
        self.parityNodes += [newNode] ##
        self.nodeNames += [nodeName]
        self.nodeNums[nodeName] = self.n_neg
        self.n_neg += 1
        
        self.F[nodeName] = [] ## vs only if nonNeg

        if not isNegative:
            self.nodes += [newNode]
            self.n += 1

        if self.debug:
            if self.not_string in nodeName:
                positive_name = nodeName.replace(self.not_string,'')
                #print("pos name=",positive_name,'vs',nodeName,self.nodesByName(nodeName).num,self.nodesByName(positive_name).num)
                assert(self.nodesByName(nodeName).num - self.n == self.nodesByName(positive_name).num)
                # for example, LDOI relies on this precise ordering



    def build_Aexp(self):
        self.n_exp = self.n_neg # build_Aexp will iterate this
        N = self.n_neg+self._num_and_clauses()
        self.A_exp = cp.zeros((N,N),dtype=bool) #adj for expanded net 
        composites = []

        for node in self.parityNodes:
            if self.debug:
                assert(node.F() == self.F[node.name]) 

            for clause in node.F():
                if len(clause)>1: # make a composite node
                    self.A_exp[self.n_exp,node.num]=1
                    for j in range(len(clause)):
                        self.A_exp[self.nodeNums[clause[j]],self.n_exp]=1
                    self.n_exp+=1
                elif clause not in ['0','1',['0'],['1']]: # ignore tautologies
                    self.A_exp[self.nodeNums[clause[0]],node.num]=1
                else:
                    print("net.py build_Aexp(): tautology ignored")
        
        if self.debug:
            assert(N==self.n_exp)

    def build_counts(self):
        index_dtype = util.get_int_dtype(max(cp.sum(self.A_exp,axis=0))) 

        self.counts = cp.sum(self.A_exp,axis=0,dtype=index_dtype)
        self.counts[:self.n_neg] = 1 # non composite nodes are OR gates
        #self.counts = cp.tile(self.counts,self.n_exp).reshape(self.n_exp,self.n_exp) # num input nodes needed to activate
        # broadcast might be more effic than tile/reshape, although this is not called often

    def build_parity_functions(self,originalG):
        from espresso import not_to_dnf # cant install pyeda on arch linux, so trying to avoid

        for node in self.parityNodes:
            if not node.isNegative: #should not already contain complement nodes jp
                complName = self.not_string + node.name 
                assert(complName not in self.F.keys())
                if complName not in self.nodeNames:
                    self.add_node(complName,isNegative=True)
                self.F[complName] = not_to_dnf(originalG.F[node.name],self)


    def nodesByName(self,name):
        return self.parityNodes[self.nodeNums[name]]


    def _num_and_clauses(self):
        count=0
        for node in self.parityNodes:
            for clause in node.F():
                if len(clause)>1: 
                    count+=1
        return count

    def count_clauses_and_lits(self):
        for node in self.parityNodes:
            self.num_clauses += len(node.F()) 
            self.max_clauses = max(self.max_clauses, len(node.F()))
            for clause in node.F():
                self.max_literals = max(self.max_literals, len(clause))

    def input_names_from_vector(self,v):
        inpt_ind = self.input_indices()
        #print(inpt_ind,v)
        names = []
        for i in range(len(v)):
            if v[i]==0:
                names+=[self.nodeNames[inpt_ind[i]+self.n]]
            else:
                assert(v[i]==1)
                names+=[self.nodeNames[inpt_ind[i]]]
        return names

    def extract_regular_net(self, init_for_sim=False):
        G2 = Net(self.params, G=self,init_for_sim=init_for_sim)
        for node in G2.nodes:
            del G2.nodeNums['!' + node.name]
            del G2.F['!' + node.name]
            G2.nodeNames.remove('!' + node.name)
            node.G = G2
        if init_for_sim:
            G2.prepare(G2.params)
        assert(not isinstance(G2, ParityNet))
        return G2


##################################################################################################


def get_file_format(format_name):
    recognized_formats = ['DNFwords','DNFsymbolic','bnet']
    if format_name not in recognized_formats:
        format_name = 'bnet' # assumed
        #sys.exit("ERROR net.py: first line of network file is the format, which must be one of" + str(recognized_formats))
    
    if format_name == 'DNFsymbolic':
        node_fn_split = '\t'
        clause_split = ' '
        literal_split = '&'
        not_str = '-'
        strip_from_clause = []
        strip_from_node = []
    elif format_name == 'DNFwords':
        node_fn_split = '*= '
        clause_split = ' or '
        literal_split = ' and '
        not_str = 'not '
        strip_from_clause = ['(',')']
        strip_from_node = ['(',')']
    elif format_name == 'bnet':
        node_fn_split = ','
        clause_split = '|'
        literal_split = '&'
        not_str = '!'
        strip_from_clause = ['(',')']
        strip_from_node = []

    return node_fn_split, clause_split, literal_split, not_str, strip_from_clause,strip_from_node

