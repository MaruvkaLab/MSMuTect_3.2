def equivalent_lists(l_a: list, l_b: list) -> bool:
    if len(l_a) != len(l_b):
        return False
    for i in range(len(l_a)):
        if l_a[i] != l_b[i]:
            return False
    return True