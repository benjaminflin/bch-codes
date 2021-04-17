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

# print(bin(eval_at_alpha_pow_bitwise(0b1101110111010000110110101, 9, 0b100011101)))

# Polynomials for a (255, 231) BCH code with t = 3
# Citation:
# https://link.springer.com/content/pdf/bbm%3A978-1-4899-2174-1%2F1.pdf 
primitive_poly = 0b100011101
# generator_poly = lcm(m_1, ..., m_6) = (100011101)*(101110111)*(111110011)
generator_poly = 0b1101110111010000110110101

#message = 2**232 - 1      # all 1s
message = 0              # all 0s

encoded_msg = mult_as_polys(message, generator_poly)

encoded_msg ^= 1 << 50

# syndromes:
for i in range(1, 7):
    print(bin(eval_at_alpha_pow_bitwise(encoded_msg, i, primitive_poly)))

#print(bin(bitwise_long_divide(encoded_msg, generator_poly)[0]))

print(bin(bitwise_long_divide(1 << 50, primitive_poly)[1]))