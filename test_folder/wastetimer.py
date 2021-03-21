import time


def measuretime(func):
    def wrapper():
        starttime = time.perf_counter()
        func()
        endtime = time.perf_counter()
        print(f"Time needed: {endtime - starttime} seconds")
    return wrapper


@measuretime

def wastetime():
    sum([i**2 for i in range(1000000)])


wastetime()