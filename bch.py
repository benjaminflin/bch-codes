# Ben Sherwin, Ben Flin, Max Helman
# Idea: interpret binary int as polynomial in F_2[x], e.g. 0b10110 = x^4 + x^2 + x
# The advantage of this is that adding two polynomials in F_2 can be calculated by XOR
#   and many other attributes/manipulations are easier/faster

# Returns poly1 * poly2 in F_2[x]
def mult_as_polys(poly1, poly2):
    result = 0
    for i in range(poly1.bit_length()):
        for j in range(poly2.bit_length()):
            # First, we get the first i and j bits of poly1 and poly2 by rightshifting (>>) by i and j respectively.
            # Then, we multiply the first i and j bits/coefficients in F_2 using bitwise and (&). 
            # Then, we (&) the result of the multiplication with 1 in order to set the first (i-1) bits to 0.
            # We leftshift (<<) the single bit back to (i + j) since now that single bit represents the coeffecient of x^(i+j).  
            # Finally, we add it to result using xor (^).
            # E.g.
            # poly1 = 0b101 = x^2 + 1,  poly2 = 0b011 = x + 1
            # result = 0b1111 = x^3 + x^2 + x + 1  
            result ^= ((poly1 >> i) & (poly2 >> j) & 1) << (i + j)
    return result

# Returns poly1 / poly2, poly1 mod poly2
def bitwise_long_divide(poly1, poly2):
    deg_dividend = poly1.bit_length()
    deg_divisor = poly2.bit_length()
    remainder = poly1
    quotient = 0
    # Implements polynomial long division through bit manipulations
    # E.g. for the first iteration of long division:  
    # 
    #        0b10 (quotient |= deg1 - deg2)
    #      ------- 
    # 0b11 | 0b101 
    #       -0b110 (poly2 << (deg1 - deg2))
    #      -------  
    #        0b011 (remainder ^= poly2 << (deg1 - deg2))

    # Continue dividing until deg(current dividend) < deg(divisor)   
    while deg_dividend >= deg_divisor:
        # We obtain the current remainder by rightshifting the divisior 
        # by the difference between the degrees of the dividend and remainder. This step allows
        # us to get the polynomial to subtract from the dividend. Then, to
        # subtract, we use bitwise xor (^).
        remainder ^= poly2 << (deg_dividend - deg_divisor)

        # The digit associated with the quotient can be thought of similarly
        # as above. We rightshift a single bit by the difference of the degrees
        # and add it to the quotient. This provides the necessary binary int x, such that 
        # mult_as_polys(x, poly2) = (poly2 << (deg1 - deg2))
        quotient |= 1 << (deg_dividend - deg_divisor)

        # The remainder becomes the new dividend for the next iteration
        deg_dividend = remainder.bit_length()

    return quotient, remainder

# calculates poly(alpha^i) in GF(2^(deg generator))
def eval_at_alpha_pow_bitwise(poly, i, generator):
    # for each power of poly(x^i), calculate and add up corresponding power of alpha:
    p = poly
    result = 0
    while p != 0:
        result ^= bitwise_long_divide(1 << ((p.bit_length() - 1) * i) % (2**(generator.bit_length() - 1) - 1), generator)[1]
        #result ^= bitwise_long_divide(1 << ((p.bit_length() - 1) * i), generator)
        p ^= 1 << (p.bit_length() - 1)

    return result


# Polynomials for a (255, 231) BCH code with t = 3
# Citation:
# https://link.springer.com/content/pdf/bbm%3A978-1-4899-2174-1%2F1.pdf 
primitive_poly = 0b100011101
# generator_poly = lcm(m_1, ..., m_6) = (100011101)*(101110111)*(111110011) for t = 3
generator_poly = 0b1101110111010000110110101
# max number of correctable errors
t = 3

# message is 255 - deg(generator_poly) + 1 bits of data
message = 0
for char in "Inbar is a great professor!":
    message = message << 8
    message += ord(char)

print("Original Message: " + str(message))

# encodes message
encoded_msg = mult_as_polys(message, generator_poly)

# introduces errors
added_errors = [51, 150, 218]
for err in added_errors:
    encoded_msg ^= 1 << err

print("Introducing errors at: " + str(added_errors))

# generates the syndromes of the 255 bit encoded polynomial
def get_syndromes(encoded, primitive):
    return [eval_at_alpha_pow_bitwise(encoded, i, primitive) for i in range(1, 2*t + 1)]

syndromes = get_syndromes(encoded_msg, primitive_poly)
#print([bin(x) for x in syndromes])

restored_msg = bitwise_long_divide(encoded_msg, generator_poly)[0]

corrected = False

# check if no errors
if syndromes == (2*t) * [0]:
    print("No errors!")
    corrected = True

# store the powers of alpha in a table
alpha_pows = [bitwise_long_divide(1 << i, primitive_poly)[1] for i in range(255)]

# check for 1 error
# find alpha^i = first syndrome and flip the ith bit and see if it's a valid message
if not corrected: 
    error_bit = alpha_pows.index(syndromes[0])
    potential_fix = encoded_msg ^ (1 << error_bit)

    # check for correctness
    if get_syndromes(potential_fix, primitive_poly) == (2*t) * [0]:
        restored_msg = bitwise_long_divide(potential_fix, generator_poly)[0]
        print("1 error at " + str(error_bit) + "th bit")
        corrected = True


# making a class for numpy matrix operations
# GF256 ~= F_2[x] / (p(x))
class GF256(object):
    def __init__(self, poly):
        self.poly = poly
    def __add__(self, x):
        return GF256(self.poly ^ x.poly)
    def __mul__(self, x):
        return GF256(bitwise_long_divide(mult_as_polys(self.poly, x.poly), primitive_poly)[1])
    def __sub__(self, x):
        return GF256(self.poly ^ x.poly)
    def __truediv__(self, x):
        return GF256(alpha_pows[(alpha_pows.index(self.poly) - alpha_pows.index(x.poly)) % (2**(primitive_poly.bit_length() - 1) - 1)])
    def __eq__(self, x):
        return self.poly == x.poly
    def __ne__(self, x):
        return self.poly != x.poly
    def __pow__(self, int_power):
        if int_power == 0:
            return GF256(1)
         
        val = GF256(self.poly)
        for _ in range(int_power - 1):
            val *= GF256(self.poly)
        return val
    def __repr__(self):
        return str(self.poly) #bin(self.poly)

import numpy as np
# Peterson–Gorenstein–Zierler algorithm to locate error correction polynomial
v = t
s_mat = np.array([[GF256(syndromes[i+j]) for i in range(v)] for j in range(v)], dtype=np.dtype(GF256))

# recursively calculate determinant of np array of GF256 objects using LaPlace extension
def rec_det(arr):
    dim = arr.shape[0]

    # do base cases manually
    if dim == 1:
        return arr[0][0]
    elif dim == 2:
        return (arr[0][0] * arr[1][1]) + (arr[1][0] * arr[0][1])
    
    # LaPlace extend and recursively calculate each smaller matrix
    d = GF256(0)
    for i in range(dim):
        smaller_arr = np.array([[arr[j][k] for j in range(1, dim)] for k in range(dim) if k != i], dtype=np.dtype(GF256))
        d = d + (arr[0][i] * rec_det(smaller_arr))
    return d

def get_minor(arr, i, j):
    dim = arr.shape[0]
    return np.array([[arr[k][l] for k in range(dim) if k != i] for l in range(dim) if l != j])

def get_inv(arr):
    dim = arr.shape[0]
    if dim == 1:
        return np.array([[GF256(1) / arr[0][0]]])

    minors = np.array([[rec_det(get_minor(s_mat, i, j)) for i in range(dim)] for j in range(dim)])
    det_inv = GF256(1) / rec_det(arr)
    return minors * det_inv

if not corrected:
    while v > 0:
        s_mat = np.array([[GF256(syndromes[i+j]) for i in range(v)] for j in range(v)])
        if rec_det(s_mat) != GF256(0):
            break
        else: 
            v -= 1
            if v == 0:
                print("empty error locator polynomial")

    c_mat = np.array([GF256(syndromes[v+i]) for i in range(v)], dtype=np.dtype(GF256)).T

    lambda_mat = []
    if v > 0:
        lambda_mat = np.matmul(get_inv(s_mat), c_mat) 
        lambda_mat = np.append(lambda_mat, [GF256(1)])
        #print(lambda_mat)


    errors = []
    # evaluate lambda polynomial at all powers of alpha looking for zeros
    for i in range(len(alpha_pows)):
        # lambda mat = [l1, l2, l3] --> l1x^3 + l2x^2 + l3x + 1
        # zeros of this poly are powers of alpha, and correspond to errors
        val = GF256(0)
        for j in range(len(lambda_mat)):
            val += lambda_mat[j] * ((GF256(alpha_pows[i]) ** (len(lambda_mat) - j - 1)))

        if val == GF256(0):
            errors.insert(0, alpha_pows.index((GF256(1) / GF256(alpha_pows[i])).poly))

    print("Found errors at: " + str(errors))
    potential_fix2 = encoded_msg
    for error in errors:
        potential_fix2 ^= 1 << error

    if get_syndromes(potential_fix2, primitive_poly) == (2*t) * [0]:
        restored_msg = bitwise_long_divide(potential_fix2, generator_poly)[0]
    else:
        print("Detected too many errors!")


print("Restored Message: " + str(restored_msg))

text = ""
recovered_text = restored_msg
while recovered_text != 0:
    text = chr(recovered_text & 0b11111111) + text
    recovered_text = recovered_text >> 8

print(text)