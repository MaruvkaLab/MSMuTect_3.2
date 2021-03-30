class AICs:
    def __init__(self, normal_normal: float = "NA", normal_tumor: float = "NA", tumor_tumor: float = "NA", tumor_normal: float = "NA"):
        self.normal_normal = normal_normal
        self.normal_tumor = normal_tumor
        self.tumor_tumor = tumor_tumor
        self.tumor_normal = tumor_normal

    @staticmethod
    def header():
        return "AIC_NORMAL_NORMAL\tAIC_NORMAL_TUMOR\tAIC_TUMOR_TUMOR\tAIC_TUMOR_NORMAL"

    def __str__(self):
        return f"{self.normal_normal}\t{self.normal_tumor}\t{self.tumor_tumor}\t{self.tumor_normal}"
