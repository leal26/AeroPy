from functools import reduce

primes = [2, 3, 5, 7, 11, 13]

def product(numbers):
    print(numbers)
    p = reduce(lambda x, y: x * y, numbers)
    return p 

print(product(primes))
# 30030

print(product(primes))
# [2, 3, 5, 7, 11, 13]