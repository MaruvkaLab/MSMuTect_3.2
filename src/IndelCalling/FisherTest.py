from collections import defaultdict
import numpy as np


class Fisher:
    def __init__(self):
        self.already_computed = defaultdict(lambda: -1)
        self.already_computed[0] = 1
        self.already_computed[1] = 1

    def factorial(self, n: int) -> float:
        if self.already_computed[n] != -1:
            return self.already_computed[n]
        else:
            answer = n * self.factorial(n-1)
            self.already_computed[n] = answer
            return answer

    def choose(self, n: int, k: int) -> float:
        answer = self.factorial(n)/(self.factorial(k)*self.factorial(n-k))
        return answer

    def test(self, first_set: np.array, second_set: np.array):
        p_value = 1
        for i in range(first_set.size):
            p_value *= self.choose(first_set[i] + second_set[i], first_set[i])
        p_value /= self.choose(np.sum(first_set)+np.sum(second_set), np.sum(first_set))
        return p_value
