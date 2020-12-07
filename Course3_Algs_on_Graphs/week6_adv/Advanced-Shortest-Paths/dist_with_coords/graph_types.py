#!/usr/bin/python3
 
# graph_types.py

#=============================================
# python Type Aliases

# Edge is a connection between two vertices
class Edge:
	# Source      Vertex
	# Destination Vertex
	# Weight      float64
    pass # an empty class that can be used as a struct

from typing import NewType

UserId = NewType('UserId', int)
some_id = UserId(524313)

"""
# Path is a series of edges
Path = list[Edge]

# Previous tracks the previous vertex in the Path
Previous dict[Vertex]*Vertex

# Distance tracks our distance estimates
Distance = map[Vertex]float64

# Visited tracks visited vertices
Visited = map[Vertex]bool
"""
