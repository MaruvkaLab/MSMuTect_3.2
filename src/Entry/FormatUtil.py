from typing import List

from src.IndelCalling.AlleleSet import AlleleSet


def format_list(results: List[str], num_components: int) -> str:
    components = [str(n) for n in results]
    if len(components) < num_components:
        for i in range(num_components - len(components)):
            components.append("NA")
    else:
        components = components[:num_components]
    return "\t".join(components)


def format_alleles(allelic_data: List[AlleleSet]) -> List[str]:
    output_lines = [
        f"{allelic_datum.histogram.locus}\t{str(allelic_datum.histogram)}\t{str(allelic_datum)}"for allelic_datum in allelic_data]
    return output_lines