from fragment import Fragment

# Create example fragments
frag1 = Fragment([0, 1, 2])
frag2 = Fragment([3, 4, 5])
frag3 = Fragment([6, 7, 8,1])

# Union operation
try:
    dimer = frag1.union(frag2)
    print(f"Union of Frag1 and Frag2: {dimer}")
except ValueError as e:
    print(f"Error during union: {e}")

# Intersection operation
intersection = frag1.intersection(frag3)
print(f"Intersection of Frag1 and Frag3: {intersection}")

# Invalid union example (overlapping indices)
frag4 = Fragment([2, 5])
try:
    invalid_union = frag1.union(frag4)
    print(f"Invalid Union Result: {invalid_union}")
except ValueError as e:
    print(f"Error during invalid union: {e}")

