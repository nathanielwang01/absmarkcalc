import argparse
import ast
from abs_markov import solution

'''
Simple driver script to demonstrate the function.

Takes command line argument as a 2d array in format where each row represents a state,
each item in a row / sum(row) represents the probability of moving to each corresponding
state from the given state, and absorbing states are represented by a row of 0's.

Returns a list of integers wherein the last value is a denominator and others are numerators
in fractions representing the probabilities of reaching each absorbing state. The probabilites
are ordered respective to each state, but they are not explicitly labeled.

e.g.
[[0, 1, 2], <-- 1/3 chance to go to state 1, 2/3 to state 2
 [1, 0, 1], <-- 1/2 chance to go to state 0, 1/2 to state 2
 [0, 0, 0]] <-- Absorbing state
as [[0,1,2],[1,0,1],[0,0,0]]
returns [1, 1]
'''

parser = argparse.ArgumentParser(description='Find the end state probabilities of an AMC.')
parser.add_argument('matrix', metavar='M', type=str, nargs=1, help='matrix as 2d list')

args = parser.parse_args()

mat = ast.literal_eval(args.matrix[0])

print(solution(mat))
