
from networkx import DiGraph
import networkx as nx
import numpy as np
from easydict import EasyDict as edict
import operator
import uuid
import underworld as uw
from underworld import function as fn

class MatGraph(DiGraph):
    def __init__(self):
        DiGraph.__init__(self) #Call the parent class init function
        self.number_changed = 1
        self.condition_list = []


    def add_transition(self, nodes, function, FnOperator, value, combineby = 'and'):
        """
        Function that allows you to easily add material transitions in Underworld simulations.
        This function creates graph 'nodes', representing the two materials.
        It also provided a simple way of decribing the rules determining the transision process

        Parameters
        ----------
        nodes : Tuple
            Representing the material Indexes of materials which a transition is likely to occur, e.g crust to eclogite.
        function: underworld.function._function.Function
            (could also be a constant)
        nOperator: operator
            operators must be provided in function form through the operator package, eg. operator.gt(1., 2.)
        value: float
            the value will be compared to the providided function, given the provided operator
        combineby: string
            'and' or 'or', defaults to 'and'. If multiple rules are provided for a single edge in the graph (representing the material transition)
            then they be applied in the sense of any ('or'), or all ('and')
        """
        firstEdge = True
        try:
            self[nodes[0]][nodes[1]] #see if the node exists
            #get names of previous condition dict:
            prevdname = self[nodes[0]][nodes[1]].keys()[0]
            firstEdge = False
        except:
            self.add_node(nodes[0])
            self.add_node(nodes[1])
            self.add_edges_from([nodes])
        #create a random name for dictionary (we need to have a different key for each condition on the graph edge)
        dname = uuid.uuid4()
        self[nodes[0]][nodes[1]][dname] = {}
        self[nodes[0]][nodes[1]][dname]['function'] = function
        self[nodes[0]][nodes[1]][dname]['operator'] =  FnOperator
        self[nodes[0]][nodes[1]][dname]['value'] =  value
        self[nodes[0]][nodes[1]][dname]['combineby'] = 'and'
        if combineby == 'or':
            self[nodes[0]][nodes[1]][dname]['combineby'] =  'or'
        if not firstEdge:
            assert self[nodes[0]][nodes[1]][dname]['combineby'] == self[nodes[0]][nodes[1]][prevdname]['combineby'], "if the graph has multiple conditions on an edge, provided 'combineby' string must be identical to avoid ambiguity."
    def walk_update(self, swarm, materialVariable):
        """
        A function that allows you to update an underworld material swarm, given the directed graph (and rules) describing those transitions

        Parameters
        ----------
        swarm: underworld.swarm._swarm.Swarm
        materialVariable: underworld.swarm._swarmvariable.SwarmVariable

        Notes
        ----------
        A better way of doing this is the future would be to link the graph and conditions directly to the underworld fn.branching.conditional()
        This would amount to writing the graph conditions directly to a condition list.
        """
        self.number_changed = 0
        #Loop through particles
        for particleID in range(swarm.particleCoordinates.data.shape[0]):
            edgeFound = False
            partCoords = [swarm.particleCoordinates.data[particleID][0],swarm.particleCoordinates.data[particleID][1]]
            node = materialVariable.data[particleID][0] #each Mat is represented by a node in the graph
            #Loop through egdes
            for otherNode in self[node].keys(): #Now loop through connected nodes (graph edges)
                if edgeFound:
                    break
                #Loop through conditions on each edge
                for cond in self[node][otherNode].keys(): #Loop though conditions governing material transformation
                    if self[node][otherNode][cond]['operator'](self[node][otherNode][cond]['function'].evaluate(partCoords)[0][0], self[node][otherNode][cond]['value']):
                        edgeFound = otherNode
                        if self[node][otherNode][cond]['combineby'] == 'or':
                            break
                    else:
                        edgeFound = False
                        if self[node][otherNode][cond]['combineby'] == 'and':
                            break

            if edgeFound:
                materialVariable.data[particleID] = edgeFound
                self.number_changed += 1 #Utilising that the True + n = 1 + n

    def build_condition_list(self, materialVariable):
        self.condition_list = [] #empty the condition list
        dm = 1e-6
        for node in self.nodes(): #Loop through nodes of graph
            for otherNode in self[node].keys(): #loop through all egdes from a given node
                #if node < otherNode:
                #this returns true for all particles with materialIndex == node (direct comparison isn't supported)
                #checkFrom = ((materialVariable > (node-dm)) and (materialVariable < (node+dm)))  #this one wasn't working
                checkFrom = operator.and_((materialVariable > (node - dm) ),
                           (materialVariable < (node + dm) ))
                condIt = 0
                for cond in self[node][otherNode].keys(): #loop through all conditions attached to the graph edge
                    op = self[node][otherNode][cond]['operator']    #
                    fun = self[node][otherNode][cond]['function']   #{extract function, operator, value}
                    val = self[node][otherNode][cond]['value']      #
                    condExp = op(fun, val)  #Now provide the function & value to the operator, return result as a variable
                    if condIt == 0:
                        totCond = condExp #if this is the first condition, assign to totCond
                    else: #if this is the NOT first condition, combine conditin with previous totCond (using AND or OR)
                        if self[node][otherNode].values()[0]['combineby'] == 'or':
                            totCond = operator.or_(totCond, condExp)
                        else:
                            totCond = operator.and_(totCond, condExp)
                    condIt += 1

                #When we pass this on to fn.branching.conditional, we only want to apply it to paticles where
                # matIndex == node, which occurs where checkFrom == True, 1
                combCond = operator.and_(totCond, checkFrom)
                #combCond = totCond
                self.condition_list.append(((combCond), otherNode))
        self.condition_list.append((True ,          materialVariable)) #if no conditions are true, return current matId
