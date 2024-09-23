from queue import LifoQueue


class MSMuLIFOqueue(LifoQueue):
    def __init__(self):
        # not making a subclass since I don't want to deal with the constructor
        super().__init__()

    def access_last_element(self):
        last_element = self.get()
        self.put(last_element)
        return last_element

    def get(self):
        return super().get(False)

    def __len__(self):
        return super().qsize()

    def __iter__(self):
        # defines conversion to list
            while len(self) != 0:
                yield self.get()

# a=MSMuLIFOqueue()
# a.put(5)
# a.put(9)
# for i in range(3):
#     a.get()
# print(list(a))