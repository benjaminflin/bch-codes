import random

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
        result = Polynomial([])
        for i in range(self.deg() + 1):
            if(self.coeff[i] == 1):
                result += poly**(self.deg() - i)
        return result
        

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


# minimum distance
distance = 3 

# (n,k) = (255, 239) generator polynomial 
# Citation: https://link.springer.com/content/pdf/bbm%3A978-1-4899-2174-1%2F1.pdf
# The octal notation was converted to binary
generator_poly = Polynomial([ int(x) for x in list("10110111101100011")])
alpha = Polynomial([ int(x) for x in "101110111" ]) 

# 239 = bit message to send
msg = 239 * [0]
for i in range(239):
    msg[i] = random.randint(0, 1) 
msg_poly = Polynomial(msg) 



if __name__ == '__main__':  

    err_rate = float(input("Enter desired error rate: "))

    print("--- Message to be encoded ---")
    print(msg_poly.pretty())

    # generate encoded_msg 
    encoded_msg = generator_poly * msg_poly

    print("--- Encoded message ---")
    print(encoded_msg.pretty())
    # randomly permute encoded message 
    for i in range(len(encoded_msg.coeff)):
        if random.random() < err_rate:
            encoded_msg.coeff[i] = random.randint(0, 1) 

    print("--- Permuted message ---")
    print(encoded_msg.pretty())

    poly1 = Polynomial([1,0,1,1])
    poly2 = Polynomial([1,1,0])
    print(poly1)
    print(poly2)
    print(poly1*poly2)
    print(generator_poly.evaluate(alpha).pretty())





