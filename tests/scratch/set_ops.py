





# class croc:
#     def __init__(self, x):
#         self.z = 5
#         self.x = x
#
#
# class trap:
#     def __init__(self, y):
#         self.z = 5
#         self.y = y
#
# zed_croc = croc(3)
# a={croc(3), croc(4), croc(3), zed_croc}
# b={croc(3), zed_croc}
# print(len(a.union(b)))

from queue import LifoQueue

a=LifoQueue()
a.put("croc")
a.put("trap")
print(list(a))
# print(a.get())
# print(a.get())
# print(a.get(False))
# print(a.qsize())