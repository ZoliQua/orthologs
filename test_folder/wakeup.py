import time


def sleep(func):
    def wrapper():
        time.sleep(10)
        return func()
    return wrapper


@sleep

def wakeup():
    print("Get up! Your break is over.")


wakeup()