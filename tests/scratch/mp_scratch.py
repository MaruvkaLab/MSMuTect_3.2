import time
import multiprocessing
import random


def wait_t():
    time.sleep(random.randint(1, 5))


a=multiprocessing.Pool(processes=3)
p_count=0
active_processes = []
while p_count<10:
    time.sleep(0.5)
    num_active_proceses = sum([1 for p in active_processes if not p.ready()])
    if num_active_proceses < 3:
        print(int(time.time()))
        active_processes.append(a.apply_async(wait_t))
        p_count+=1
    else:
        print(num_active_proceses)
print("croc")
a.close()
a.join()
print("trap")