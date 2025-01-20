from itertools import combinations

class Fragment:
    """
    Represents a single fragment with a list of indices and a label.
    """
    def __init__(self, indices, label=None, index=None):
        """
        Initializes a Fragment with indices and a label.
        :param indices: List of indices in the fragment.
        :param label: Label describing the fragment's type.
        :param index: Monomer index (for tracking original fragments).
        """
        if len(indices) != len(set(indices)):
            raise ValueError(f"Duplicate indices in fragment: {indices}")
        self.indices = indices
        self.label = label or "monomer"
        self.index = index

        # If it's a monomer, include its index in the label
        if label == "monomer" and index is not None:
            self.label = f"monomer[{index}]"

    def union(self, other):
        """
        Creates a new Fragment that is the union of this fragment and another fragment.
        Removes duplicate indices in the process.
        """
        new_indices = sorted(set(self.indices + other.indices))
        all_indices = [self.index, other.index]
        all_indices = [idx for idx in all_indices if idx is not None]
        new_label = f"union[{','.join(map(str, all_indices))}]"
        return Fragment(new_indices, new_label)

    def intersection(self, other):
        """
        Creates a new Fragment that is the intersection of this fragment and another fragment.
        Keeps only indices that are common between the two fragments.
        """
        new_indices = sorted(set(self.indices) & set(other.indices))
        if not new_indices:
            return Fragment([], "empty intersection")
        all_indices = [self.index, other.index]
        all_indices = [idx for idx in all_indices if idx is not None]
        new_label = f"intersection[{','.join(map(str, all_indices))}]"
        return Fragment(new_indices, new_label)

    def __repr__(self):
        return f"Fragment(indices={self.indices}, {self.label})"




class FragmentList:
    """
    A collection of Fragment objects with operations for unions and intersections.
    """
    def __init__(self, fragments_data, allow_overlapping):
        """
        Initializes the FragmentList with data and a flag for overlapping behavior.
        :param fragments_data: List of lists where each sublist represents a fragment's indices.
        :param allow_overlapping: Boolean indicating if overlapping indices are allowed.
        """
        self.fragments = [Fragment(indices, label="monomer", index=i) for i, indices in enumerate(fragments_data)]
        self.allow_overlapping = allow_overlapping

        if not allow_overlapping:
            self._validate_no_overlaps(fragments_data)

    def _validate_no_overlaps(self, fragments_data):
        """
        Validates that no indices are repeated across fragments if overlapping is not allowed.
        :param fragments_data: List of lists of fragment indices.
        :raises ValueError: If overlapping indices are found.
        """
        all_indices = set()
        for fragment in fragments_data:
            for index in fragment:
                if index in all_indices:
                    raise ValueError(f"Overlapping index found: {index}. Overlapping indices are not allowed.")
                all_indices.add(index)
    def generate_combinations(self, order):
        """
        Generates all unique fragment combinations of the specified order.
        - For order 1, include all monomers and their intersections (if allow_overlapping).
        - For higher orders, generate unions and intersections as appropriate.
        :param order: The order of combinations to generate.
        :return: A list of Fragment objects (unions and intersections).
        """
        results = []

        if order == 1:
            # Add all fragments directly as monomers
            results.extend(self.fragments)

            # Add intersections if overlapping is allowed
            if self.allow_overlapping:
                for frag1, frag2 in combinations(self.fragments, 2):
                    intersection = frag1.intersection(frag2)
                    if intersection.indices:  # Only add non-empty intersections
                        results.append(intersection)
        else:
            # For higher orders, compute unions and intersections
            results.extend(self.fragments)
            for combo in combinations(self.fragments, order):
                try:
                    # Compute union
                    combined_union = combo[0]
                    for frag in combo[1:]:
                        combined_union = combined_union.union(frag)
                    results.append(combined_union)
                except ValueError:
                    # Skip invalid unions
                    pass

                if self.allow_overlapping:
                    # Compute intersection
                    combined_intersection = combo[0]
                    for frag in combo[1:]:
                        combined_intersection = combined_intersection.intersection(frag)
                    if combined_intersection.indices:  # Only add non-empty intersections
                        results.append(combined_intersection)

        return results

    def __repr__(self):
        return f"FragmentList(fragments={self.fragments}, allow_overlapping={self.allow_overlapping})"
