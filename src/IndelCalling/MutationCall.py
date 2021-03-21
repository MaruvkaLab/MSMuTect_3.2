

class MutationCall:
    NO_NORMAL_ALLELES = -4
    BORDERLINE_NONMUTATION = -3
    TOO_MANY_ALLELES = -2
    INSUFFICIENT = -1
    NOT_MUTATION = 0
    MUTATION = 1

    def __init__(self, call: int, p_value=-1):
        self.call = call
        self.p_value = p_value

    def format_pval(self):
        if self.p_value == -1:
            return 'NA'
        else:
            return str(self.p_value)

    @staticmethod
    def header():
        return "CALL\tP_VALUE"

    def __str__(self):
        return f"{self.call}\t{self.format_pval()}"