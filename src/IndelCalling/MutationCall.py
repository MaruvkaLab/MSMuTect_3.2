from src.IndelCalling.AlleleSet import AlleleSet
from src.IndelCalling.AICs import AICs


class MutationCall:
    # pseudo enum
    NO_NORMAL_ALLELES = -4
    BORDERLINE_NONMUTATION = -3
    TOO_MANY_ALLELES = -2
    INSUFFICIENT = -1
    NOT_MUTATION = 0
    MUTATION = 1

    def __init__(self, call: int, normal_alleles: AlleleSet, tumor_alleles: AlleleSet, aic_values: AICs, p_value=-1):
        self.call = call
        self.normal_alleles = normal_alleles
        self.tumor_alleles = tumor_alleles
        self.aic_values = aic_values
        self.p_value = p_value

