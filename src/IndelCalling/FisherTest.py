import math, sys
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
        elif n > sys.getrecursionlimit() - 1:  # to avoid recursion overload
            ans = math.factorial(n)
            self.already_computed[n] = ans
            return ans
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
            # casted to int, so if number is too large for numpy int 64 bits
            p_value *= self.choose(int(first_set[i] + second_set[i]), int(first_set[i]))
        p_value /= self.choose(int(np.sum(first_set)+np.sum(second_set)), int(np.sum(first_set)))
        return p_value
