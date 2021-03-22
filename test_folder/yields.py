
#
# Testing generators and yields
#

def gen1():
    yield "I am the bone of my sword"
    yield "Steel is my body and fire is my blood"
    yield "I have created over a thousand blades"


temp_gen1 = gen1()
print(next(temp_gen1))
print(next(temp_gen1))
print(next(temp_gen1))

#
# Testing runtime and execution time
#

import time

def seq(start,end):
    lst = []

    for i in range(start,end):
        lst.append(i)
    return lst


start_time = time.time()

data = seq(100,100_000_000)

end_time = time.time()


print(f"execution time: {end_time - start_time}")