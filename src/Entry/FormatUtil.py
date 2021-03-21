from typing import List


def format_list(results: List[str], num_components: int) -> str:
    components = [str(n) for n in results]
    if len(components) < num_components:
        for i in range(num_components - len(components)):
            components.append("NA")
    else:
        components = components[:num_components]
    return "\t".join(components)


def format_alleles(allelic_data: list) -> List[str]: # List[AlleleSet] not declared to avoid circular import
    output_lines = [
        f"{allelic_datum.histogram.locus}\t{str(allelic_datum.histogram)}\t{str(allelic_datum)}"for allelic_datum in allelic_data]
    return output_lines