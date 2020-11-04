#Uses python3

import sys


def negative_cycle(adj, cost):
    # Output 1 if the graph contains a cycle of negative weight
    #  and 0 otherwise.
    
    #write your code here
    return 0

nlr_ei = (−log(r_ei))

# Task. Given an directed graph with possibly negative edge weights
#  and with 𝑛 vertices and 𝑚 edges, check whether it contains a cycle
#  of negative weight.
# Input Format: A graph is given in the standard format.
# Constraints: 1 ≤ 𝑛 ≤ 10^3, 0 ≤ 𝑚 ≤ 10^4,
#  edge weights are integers of absolute value at most 103.
# Output Format: Output 1 if the graph contains a cycle of negative weight
#  and 0 otherwise.

if __name__ == '__main__':
    input = sys.stdin.read()
    data = list(map(int, input.split()))
    n, m = data[0:2]
    data = data[2:]
    edges = list(zip(zip(data[0:(3 * m):3], data[1:(3 * m):3]), data[2:(3 * m):3]))
    data = data[3 * m:]
    adj = [[] for _ in range(n)]
    cost = [[] for _ in range(n)]
    for ((a, b), w) in edges:
        adj[a - 1].append(b - 1)
        cost[a - 1].append(w)
    print(negative_cycle(adj, cost))
