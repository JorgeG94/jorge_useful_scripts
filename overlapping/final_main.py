from fragment import FragmentList
from read_print import read_json

def main():
    file_path = "input.json"  # Replace with your JSON file path
    try:
        # Read JSON and initialize FragmentList
        fragments_data, allow_overlapping = read_json(file_path)
        fragment_list = FragmentList(fragments_data, allow_overlapping)

        # Specify the desired order N
        N = 2  # Change this value to test fragments of a specific order
        total_fragments = len(fragments_data)

        # Validate N <= total_fragments
        if N > total_fragments:
            print(f"Error: Requested order N ({N}) exceeds the number of available fragments ({total_fragments}).")
            return

        # Generate all fragments (monomers, unions, intersections)
        results = fragment_list.generate_combinations(N)

        # Print the entire list of fragments
        print(f"\nGenerated fragments for order {N}:")
        for frag in results:
            print(frag)

    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    main()

