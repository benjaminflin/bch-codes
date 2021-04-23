
# Ben Sherwin, Ben Flin, Max Helman
# Implementation of BCH error correcting code
#
# Requires Python version >= 3.6
# 
# Interprets polynomials in F_2[x] as a binary int, e.g. 0b10110 = x^4 + x^2 + x
#   The advantage of this is that adding two polynomials in F_2 can be calculated by XOR
#   and many other attributes/manipulations are easier/faster.
# 
# Class GF2n represents GF(2^n), the Galois field with 2^n elements as a polynomial
#   f(x) + (p(x)) where p(x) is an irreducible polynomial of degree n+1 in F_2[x].
#   Represent each element of GF(2^n) as α^i for some generator α.
#   We take α = x + (p(x)), so α^i can be calculated as x^i mod p(x).
#
# A BCH generator polynomial g(x) is calculated as the lcm(m_1, ..., m_2t) where t is the
#   number of errors we want to correct and m_i is the minimal polynomial of alpha^i in F_2[x].
#   
# Then, we encode a binary message m by treating it as a polynomial in F_2[x] as above,
#   and multiplying e(x) = m(x) * g(x), again treating the encoded message e as both a 
#   polynomial and a binary number.  Note e must be exactly 2^n - 1 bits long, meaning higher
#   values of t mean larger generator polynomials g, and less useful data can be transmitted.
#
# We can check for errors as follows: since g(α^i) = 0 for 0 <= i <= 2t, syndrome[i] = 
#   R(α^i) = m(α^i)g(α^i) + err(α^i) = 0 + err(α^i) = err(α^i), where R is the received message.
#   Then, there is an error if at least one of the syndromes is nonzero.  Furthermore, since 
#   each syndrome is an element of GF(2^n), there are (2^n - 1)^t different nonzero syndrome
#   combinations corresponding exactly to the (2^n - 1)^t possible correctable errors.
#
# In order to find this correspondence between the syndromes and the errors, we calculate 
#   the error location polynomial which has factors (x * α^i - 1) using the PGZ algorithm.
#   Then, we factor the error location polynomial using brute force.  Since the roots are 
#   elements of GF(2^n), this is feasable for small n (we use n = 8 by default, but this can
#   be modified by replacing the primitive polynomial p(x)), however for much larger codewords,
#   it is better to factor the error location polynomial using the Forney algorithm, as
#   detailed in the writeup.

import random
import numpy as np

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

# Calculates poly(α^i) in GF(2^(deg(generator) - 1))
def eval_at_alpha_pow_bitwise(poly, i, generator):
    # For each power of poly(x^i), calculate and add up corresponding power of α:
    p = poly
    result = 0
    while p != 0:
        # Move coefficient x^k to x^(i+k) for highest k
        result ^= 1 << ((p.bit_length() - 1) * i)
        # Pop highest order coefficient of p (i.e. decrement degree of polynomial)
        p ^= 1 << (p.bit_length() - 1)
    # return poly(x^i) mod p(x)
    return bitwise_long_divide(result, generator)[1]

# Irreducible polynomial of degree 9 in F_2[x]
# Generates GF(2^8), so our default message length is q^m - 1 = 256 bits
primitive_poly = 0b100011101

# Store the powers of alpha in a table for future use
alpha_pows = [bitwise_long_divide(1 << i, primitive_poly)[1] for i in range(255)]

# Class for numpy matrix operations in GF(256)
# GF2n ~= F_2[x] / (primitive_poly)
class GF2n(object):
    def __init__(self, poly):
        self.poly = poly
    def __add__(self, x):
        return GF2n(self.poly ^ x.poly)
    def __mul__(self, x):
        return GF2n(bitwise_long_divide(mult_as_polys(self.poly, x.poly), primitive_poly)[1])
    def __sub__(self, x):
        return GF2n(self.poly ^ x.poly)
    def __truediv__(self, x):
        # For a / b, write a = α^i, b = α^j.  Then, calculate α^(i - j mod 2^n - 1)
        return GF2n(alpha_pows[(alpha_pows.index(self.poly) - alpha_pows.index(x.poly)) % (2**(primitive_poly.bit_length() - 1) - 1)])
    def __eq__(self, x):
        return self.poly == x.poly
    def __ne__(self, x):
        return self.poly != x.poly
    def __pow__(self, int_power):
        if int_power == 0:
            return GF2n(1)
         
        val = GF2n(self.poly)
        for _ in range(int_power - 1):
            val *= GF2n(self.poly)
        return val
    def __repr__(self):
        return str(self.poly)


# Gets the (i,j)th minor of a matrix
def get_minor(arr, i, j):
    dim = arr.shape[0]
    return np.array([[arr[k][l] for k in range(dim) if k != i] for l in range(dim) if l != j])

# Recursively calculate determinant of np array of GF2n objects using LaPlace extension
def rec_det(arr):
    dim = arr.shape[0]

    # do base cases manually
    if dim == 1:
        return arr[0][0]
    elif dim == 2:
        return (arr[0][0] * arr[1][1]) + (arr[1][0] * arr[0][1])
    
    # Laplace extend and recursively calculate each smaller matrix
    d = GF2n(0)
    for i in range(dim):
        smaller_arr = get_minor(arr, 0, i)
        d = d + (arr[0][i] * rec_det(smaller_arr))
    return d

# Gets the inverse of a matrix in GF2(n) by method of adjugate minors 
def get_inv(arr):
    dim = arr.shape[0]
    if dim == 1:
        return np.array([[GF2n(1) / arr[0][0]]])

    minors = np.array([[rec_det(get_minor(s_mat, i, j)) for i in range(dim)] for j in range(dim)])
    det_inv = GF2n(1) / rec_det(arr)
    return minors * det_inv


# Generator polynomials for a (255, 231) BCH code with t = 0 - 10
# Citation:
# https://link.springer.com/content/pdf/bbm%3A978-1-4899-2174-1%2F1.pdf 

# generator_poly = lcm(m_1, ..., m_6) = (100011101)*(101110111)*(111110011) for t = 3
generators = [1, \
              0b100011101, \
              0b10110111101100011, \
              0b1101110111010000110110101, \
              0b111101110010110110100001011111101, \
              0o23157564726421, \
              0o16176560567636227, \
              0o7633031270420722341, \
              0o2663470176115333714567, \
              0o52755313540001322236351, \
              0o22624710717340432416300455]

# t = max number of correctable errors
t = int(input(f"Enter maximum number of correctable errors (0-{len(generators) - 1}):")) % (len(generators) - 1)
generator_poly = generators[t]

# Ask user for encoded message
input_message = input(f"Enter a message to encode ({(256 - generator_poly.bit_length()) // 8} bytes/characters max): ") 

# Truncate message to 256 - deg(generator) bytes
input_message = input_message[:((256 - generator_poly.bit_length()) // 8)]

# The following encodes in ascii the message to be sent
# message is 255 - deg(generator_poly) + 1 bits of data
message = int.from_bytes(input_message.encode('utf-8'), byteorder='little')

print("Original Message: " + hex(message))

# Encodes message into 255 bit codeword by multiplying by the generator polynomial 
encoded_msg = mult_as_polys(message, generator_poly)

added_errors = [ random.randint(0, 255) for x in
    range(int(input(f"Enter number of errors (max {t}): "))) ]

# Introduces errors to the message at the 3 following positions in the message
for err in added_errors:
    encoded_msg ^= 1 << err

print("Introducing errors at: " + str(added_errors))

# Generates the syndromes of the 255 bit encoded polynomial
def get_syndromes(encoded, primitive):
    # Each syndrome s_j = R(alpha^i) where R is our encoded word for 2*t possible syndromes 
    # Since each power of alpha is a root of our generator polynomial, evaluating R(x) = G(x) + E(x) at 
    # alpha^i will simply get us our error polynomial E(alpha^i). 
    return [eval_at_alpha_pow_bitwise(encoded, i, primitive) for i in range(1, 2*t + 1)]
syndromes = get_syndromes(encoded_msg, primitive_poly)

# Get restored message by dividing by generator polynomial (may be mangled due to errors)
restored_msg = bitwise_long_divide(encoded_msg, generator_poly)[0]

corrected = False


# The following is the Peterson–Gorenstein–Zierler algorithm to locate 
# error correction polynomial 
# v represents is the hamming distance to the closest codeword
v = t
# Matrix S of size v x v
# This matrix is singular when no codeword is exactly hamming distance v away 
s_mat = np.array([[GF2n(syndromes[i+j]) for i in range(v)] for j in range(v)])

# Check if no errors
if syndromes == (2*t) * [0]:
    print("No errors!")
    corrected = True
else:
    # Procedure to check det(S) = 0 and decrement v if S is singular 
    # Note, this assumes that if det(S) = 0, v must be a lower value 
    while v > 0:
        s_mat = np.array([[GF2n(syndromes[i+j]) for i in range(v)] for j in range(v)])
        if rec_det(s_mat) != GF2n(0):
            break
        else: 
            v -= 1
            if v == 0:
                print("Empty error locator polynomial")

    # Build vector C of size (v x 1)
    c_mat = np.array([GF2n(syndromes[v+i]) for i in range(v)]).T

    # SΛ = -C so Λ is given by (S^-1)C
    # We ignore the negative sign because the coefficients are in Z_2
    lambda_mat = []
    if v > 0:
        lambda_mat = np.matmul(get_inv(s_mat), c_mat) 
        # Add one to the lambda mat to yield
        # our error locator polynomial: λ_1 x^n + ... + λ_n x + 1
        lambda_mat = np.append(lambda_mat, [GF2n(1)])

    errors = []
    # Evaluate lambda polynomial at all powers of alpha looking for zeros
    # Zeros of this poly are powers of alpha, and correspond to errors
    for i in range(len(alpha_pows)):
        # Add up powers of lambda(alpha^j) for deg(lambda_mat)
        val = GF2n(0)
        for j in range(len(lambda_mat)):
            val += lambda_mat[j] * ((GF2n(alpha_pows[i]) ** (len(lambda_mat) - j - 1)))

        # If there is a zero, then the power of alpha i gives us the position
        # of the error inside the codeword, i.e.
        # Λ(x) = (α^i_1 - 1)(α^i_2 - 1)...(α^i_n - 1)
        # yields zeros α^-i_1 ... α^-i_n.
        # i_1, ..., i_n are the positions of the errors.
        if val == GF2n(0):
            errors.insert(0, alpha_pows.index((GF2n(1) / GF2n(alpha_pows[i])).poly))

    print("Located errors at: " + str(errors))
    # Flip the bits at each of the errors 
    potential_fix = encoded_msg
    for error in errors:
        potential_fix ^= 1 << error

    # If there are no errors, then we divide by the generator polynomial to get
    # the restored message 
    if get_syndromes(potential_fix, primitive_poly) == (2*t) * [0]:
        restored_msg = bitwise_long_divide(potential_fix, generator_poly)[0]
    else:
        print("Detected too many errors!")

print("Restored Message: " + hex(restored_msg))

# Decode the restored message as text 
text = restored_msg.to_bytes(restored_msg.bit_length(), byteorder='little').decode('utf-8')

print(text)