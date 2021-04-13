import random
import functools
from sympy import *
from sympy.abc import x, alpha
from sympy.polys.galoistools import gf_rem, gf_compose

class Polynomial:
    def __init__(self, coeff = []):
        self.coeff = coeff
    def __mul__(self, other):
        poly1 = self
        poly2 = other
        result = Polynomial()
        result.coeff = (poly1.deg() + poly2.deg() + 1) * [0]
        for i in range(poly1.deg() + 1):
            for j in range(poly2.deg() + 1):
                result.coeff[i+j] ^= (poly1.coeff[i] & poly2.coeff[j])
        return result
        
    def deg(self):
        return len(self.coeff) - 1
    
    def __str__(self):
        string = ""
        n = self.deg()
        for a in self.coeff:
            string += str(a) + "x^" + str(n) + (" + " if n > 0 else "")
            n = n - 1
        return string 
    
    def pretty(self):
        string = ""
        for i in self.coeff:
            string += str(i)
        return string

    def evaluate(self, value):
        if type(value) == int:
            result = 0
            for i in range(self.deg() + 1):
                result += value**i * self.coeff[self.deg() - i]
            return result
        else:
            return self.evaluate_polynomial(value)
    
    def evaluate_polynomial(self, poly):
        # result = Polynomial([])
        # for i in range(self.deg() + 1):
        #     if(self.coeff[i] == 1):
        #         result += poly**(self.deg() - i)
        # return result

        return Polynomial(gf_compose(ZZ.map(self.coeff), ZZ.map(poly.coeff), 2, ZZ))
        

    def __add__(self, other):
        deg = max(self.deg(), other.deg()) 
        coeff = (deg + 1) * [0]

        poly1 = self.coeff[::-1]
        poly2 = other.coeff[::-1]

        for i in range(deg + 1):
            if i < len(poly1):
                coeff[i] += poly1[i] 
            if i < len(poly2):
                coeff[i] += poly2[i] 
            coeff[i] %= 2 

        result = Polynomial(coeff[::-1]) 
        return result

    def __pow__(self, exp):
        # power by repeated multiplication
        result = Polynomial(self.coeff)
        for i in range(exp - 1):
            result = result * self
        return result
    
    def remainder(self, poly):
        return Polynomial(gf_rem(ZZ.map(self.coeff), ZZ.map(poly.coeff), 2, ZZ))



# 1. Create a list (list1) of x^i mod generator_poly for 0 <= i < 256 (really alphas)
# 2. Another, parallel list (list2) of each of the above polynomials with coefficients 
#    listwise interpreted as a binary integer (x^2 + 1 = 101 = 5)
# 3. In order to calculate generator_poly(alpha^i), calculate generator_poly(x^i), 
#    and then for each term x^j, look up and xor list2[j % 255]. (adding in F_2[x])
# 4. Interpret final result as polynomial.


# (n,k) = (255, 239) generator polynomial 
# Citation: https://link.springer.com/content/pdf/bbm%3A978-1-4899-2174-1%2F1.pdf
# The octal notation was converted to binary
generator_poly = Poly.from_list([ int(i) for i in list("10110111101100011")], gens=[x], domain=GF(2)) 

print(generator_poly.all_coeffs())
basis_poly = Poly.from_list([ int(i) for i in list("100011101")], gens=[x], domain=GF(2)) # NOT 101110111

powers_of_alpha = [ Poly(x ** i, gens=[x], domain=GF(2)).rem(basis_poly) for i in range(256) ]
list2 = [ p.all_coeffs() for p in powers_of_alpha ]

def eval_at_alpha_pow(poly, i):
    result = [0, 0, 0, 0, 0, 0, 0, 0]
    poly = compose(poly, x**i)
    for j in poly.all_coeffs():
        for i in range(len(list2[j % 255])):
            result[i] ^= list2[j % 255][i]  
    return Poly.from_list(result, gens=[x], domain=GF(2))

print(eval_at_alpha_pow(generator_poly, 2))








# if __name__ == '__main__':  

#     err_rate = float(input("Enter desired error rate: "))

#     print("--- Message to be encoded ---")
#     print(msg_poly.pretty())

#     # generate encoded_msg 
#     encoded_msg = generator_poly * msg_poly

#     print("--- Encoded message ---")
#     print(encoded_msg.pretty())
#     # randomly permute encoded message 
#     for i in range(len(encoded_msg.coeff)):
#         if random.random() < err_rate:
#             encoded_msg.coeff[i] = random.randint(0, 1) 

#     print("--- Permuted message ---")
#     print(encoded_msg.pretty())

#     poly1 = Polynomial([1,0,1,1])
#     poly2 = Polynomial([1,1,0])
#     print(poly1)
#     print(poly2)
#     print(poly1*poly2)
#     print(generator_poly.evaluate(alpha).pretty())