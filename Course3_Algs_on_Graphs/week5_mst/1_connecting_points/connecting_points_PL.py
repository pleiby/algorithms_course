#Uses python3
import sys
import math

def minimum_distance(x, y):
    result = 0.
    #write your code here
    return result

# Task. 
# Given 𝑛 points on a plane, connect them with segments of minimum total length
#  such that there is a path between any two points. Recall that the length of
#  a segment with endpoints (𝑥1,𝑦1) and (𝑥2,𝑦2) is equal to √︀(𝑥1 − 𝑥2)^2 + (𝑦1 − 𝑦2)^2.

# Input Format. The first line contains the number 𝑛 of points. 
# Each of the following 𝑛 lines defines a point (𝑥𝑖, 𝑦𝑖).

def process_input_data():
    data = list(map(int, input.split()))
    n = data[0]
    x = data[1::2]
    y = data[2::2]
    return


if __name__ == '__main__':
    input = sys.stdin.read()

    print("{0:.9f}".format(minimum_distance(x, y)))
