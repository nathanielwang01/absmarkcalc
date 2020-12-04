from fractions import Fraction
from fractions import gcd
from operator import sub

# Utility function for testing
def print_mat(matrix):
    for row in matrix:
        for item in row:
            print(item, end="\t")
        print()
    print()

# Turns list of Fractions into list of numerators append denominator
def format_fracs(frac_list):
    # Must LCM all values in list
    a = frac_list[0].denominator
    for i in range(1, len(frac_list)):
        b = frac_list[i].denominator

        a = (a * b) // gcd(a, b)

    # a = LCM
    formatted_list = [frac.numerator * (a // frac.denominator) for frac in frac_list]
    formatted_list.append(a)

    return formatted_list

# Adapted from Wikipedia pseudocode @https://en.wikipedia.org/wiki/Row_echelon_form
def rref(mat):
    lead = 0
    n = len(mat)

    for r in range(n):
        #print("looped")
        if n <= lead:
            #print("return 1")
            return
        
        i = r
        while mat[i][lead] == 0:
            i += 1
            if n == i:
                i = r
                lead += 1
                if n == lead:
                    #print("return 2")
                    return
        
        if i != r:
            # Swap rows i and r
            mat[r], mat[i] = mat[i], mat[r]

        if mat[r][lead] != 0:
            # Divide row r by mat[r][lead]
            mat[r] = [item / mat[r][lead] for item in mat[r]]

        for i in range(n):
            if i != r:
                # Subtract mat[i][lead] multiplied by row r from row i
                temp = [mat[i][lead] * item for item in mat[r]]
                mat[i] = list(map(sub, mat[i], temp))
            
        lead += 1
 
def solution(m):
    # Track indices of absorbing and transient rows
    a_rows = []
    t_rows = []

    prob_mat = []

    # Create P matrix
    n = len(m)
    for row_idx in range(n):
        row_denom = sum(m[row_idx])

        if row_denom:
            # Fractionalize row probs
            row_prob = [Fraction(item, row_denom) for item in m[row_idx]]

            t_rows.append(row_idx)
            prob_mat.append(row_prob)
        
        else:
            # 0 -> fractions
            row_prob = [Fraction(0) for item in m[row_idx]]
            row_prob[row_idx] = 1
            
            a_rows.append(row_idx)
            prob_mat.append(row_prob)

    header = [i for i in range(n)]

    # Put P in canonical form (minus '1' probabilities)
    mina_idx = min(a_rows)
    maxt_idx = max(t_rows)

    # If necessary, exchange rows & columns
    while maxt_idx > mina_idx:
        # Swap columns
        for row_idx in range(n):
            prob_mat[row_idx][maxt_idx], prob_mat[row_idx][mina_idx] = prob_mat[row_idx][mina_idx], prob_mat[row_idx][maxt_idx]

        # Swap rows
        prob_mat[maxt_idx], prob_mat[mina_idx] = prob_mat[mina_idx], prob_mat[maxt_idx]

        # Update tracking lists
        a_rows[a_rows.index(mina_idx)], t_rows[t_rows.index(maxt_idx)] = t_rows[t_rows.index(maxt_idx)], a_rows[a_rows.index(mina_idx)]

        # Update header
        header[mina_idx], header[maxt_idx] = header[maxt_idx], header[mina_idx]

        # Update max and min values
        mina_idx = min(a_rows)
        maxt_idx = max(t_rows)

    # Slice for Q
    q_mat = [prob_mat[i][:maxt_idx + 1] for i in range(maxt_idx + 1)]

    # N = (I - Q)^-1
    n_mat = []
    for i in range(len(q_mat)):
        n_mat.append([])
        for k in range(len(q_mat)):
            if i == k:
                n_mat[i].insert(k, 1 - q_mat[i][k])
                n_mat[i].append(1)
            else:
                n_mat[i].insert(k, 0 - q_mat[i][k])
                n_mat[i].append(0)

    # Perform row reduction for inversion
    rref(n_mat)
    n_mat = [n_mat[i][maxt_idx + 1:] for i in range(maxt_idx + 1)]

    # B = N * R
    r_mat = [prob_mat[i][maxt_idx + 1:] for i in range(maxt_idx + 1)]

    b_mat = []
    # Dot product
    for row in n_mat:
        # Multiply contents of row by each column in other matrix
        b_mat.append([])
        for col_idx in range(len(r_mat[0])):
            col = []
            for row_idx in range(len(r_mat)):
                col.append(r_mat[row_idx][col_idx] * row[row_idx])
            b_mat[-1].append(sum(col))

    # Row containing 0 of B contains desired probabilities
    probs = sorted(list(zip(header[len(n_mat):], b_mat[header.index(0)])))
    probs = [elem[1] for elem in probs]

    return format_fracs(probs)
