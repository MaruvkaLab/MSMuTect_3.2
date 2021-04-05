from typing import List


def format_list(results: List[str], num_components: int) -> str:
    components = [str(n) for n in results]
    if len(components) < num_components:
        for i in range(num_components - len(components)):
            components.append("NA")
    else:
        components = components[:num_components]
    return "\t".join(components)


