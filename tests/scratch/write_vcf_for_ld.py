from typing import List, Tuple


def formulate_vcf_line(chromosome: str, start: int) -> str:
    # ref_repeats = current_line.num_ref_repeats
    ref = "A"
    alt = "B"

    output_line = [
        chromosome,
        str(start),
        ".",  # id
        ref,
        alt,
        ".",  # quality,
        "PASS",  # FILTER
        "."  # info
    ]
    return "\t".join(output_line)



def create_vcf(output_vcf_fp: str, alternate_alleles: List[Tuple[str, int]]):
    vcf_lines = "\n".join([formulate_vcf_line(a[0], a[1]) for a in alternate_alleles])
    header = """
##fileformat=VCFv4.0
##<ID=ID,Number=,Type=,Description="NOTE: FOR ALL IMPURE LOCI, THE LOCATION OF THE MUTATION MAY BE INCORRECT BY UP TO THE LENGTH OF THE MICROSATELLITE LOCUS. THE EXACT REF AND ALT SEQUENCE MAY BE ERRANT AS WELL">
##fileformat=VCFv4.2
##source=MSMuTect
##phasing=partial
##contig=<ID=1>
##contig=<ID=1>
##contig=<ID=2>
##contig=<ID=3>
##contig=<ID=4>
##contig=<ID=5>
##contig=<ID=6>
##contig=<ID=7>
##contig=<ID=8>
##contig=<ID=9>
##contig=<ID=10>
##contig=<ID=11>
##contig=<ID=12>
##contig=<ID=13>
##contig=<ID=14>
##contig=<ID=15>
##contig=<ID=16>
##contig=<ID=17>
##contig=<ID=18>
##contig=<ID=19>
##contig=<ID=20>
##contig=<ID=21>
##contig=<ID=22>
##contig=<ID=X>
##contig=<ID=Y>
##FILTER=<ID=NM,Description="Not mutation">
##FILTER=<ID=AN,Description="Either tumor or normal lacks alleles">
##FILTER=<ID=RR,Description="Reversion to reference">
##FILTER=<ID=FFT,Description="Failed Fisher Test">
##FILTER=<ID=TMA,Description="Too many alleles">
##FILTER=<ID=INS,Description="Insufficient support">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
            """
    with open(output_vcf_fp, 'a') as ovf:
        ovf.write(header + "\n" + vcf_lines)


create_vcf("/home/avraham/MaruvkaLab/Magpie/linkage/plink/synthetic/synvcf.vcf", [])