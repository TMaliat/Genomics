import copy

width = 0
solutions = set()  # To store unique solutions
all_solutions = []  # To store all solutions, unique or not

def partial_digest(L, filename, step_filename):
    global width, solutions, all_solutions
    width = max(L)
    solutions = set()  # Reset unique solutions
    all_solutions = []  # Reset all solutions
    
    L.remove(width)
    X = {0, width}
    place(L, X, step_filename)
    
    # Write all solutions (unique or not) to the specified file
    with open(filename, "w") as file:
        # Only write unique solutions to file
        for sol in set(all_solutions):  # Using set to remove duplicates
            file.write(f"{sorted(sol)}\n")  # Write each solution to the file
    
    # Print unique solutions in the terminal
    print(f"Unique solutions for {filename[:-4]}:", [sorted(sol) for sol in solutions])

def is_subset_of(main, subset):
    mainList = copy.deepcopy(main)
    subList = copy.deepcopy(subset)
    
    # Handle single value main list
    if not isinstance(mainList, list):
        mainList = [mainList, 'JJ']
    if not isinstance(subList, list):
        subList = [subList, 'JJ']

    for x in subList:
        if x in mainList:
            mainList.remove(x)
        else:
            return False
    return True

def log_step(value, L, X, y_or_width_y, step_filename):
    """
    Log each step into step.txt for each test case
    """
    with open(step_filename, "a") as step_file:
        step_file.write(f"Inserted: {value}, L: {sorted(L)}, X: {sorted(X)}, Value used: {y_or_width_y}\n")

def place(L, X, step_filename):
    global solutions, all_solutions
    if len(L) == 0:
        solutions.add(tuple(sorted(X)))  # Add only unique solutions to solutions set
        all_solutions.append(tuple(sorted(X)))  # Append all solutions to all_solutions
        return
    
    L.sort()
    for y in L[::-1]:
        listX = list(X)
        
        # First option: Place y
        deltaYX = [abs(y - x_i) for x_i in listX]
        
        if is_subset_of(L, deltaYX):
            X.add(y)
            for l in deltaYX:
                L.remove(l)
            log_step(y, L, X, y, step_filename)  # Log the step
            place(L, X, step_filename)
            # Backtrack
            X.remove(y)
            for l in deltaYX:
                L.append(l)

        # Second option: Place width - y
        deltaWidthYX = [abs(x_i - abs(width - y)) for x_i in listX]
        
        if is_subset_of(L, deltaWidthYX):
            X.add(abs(width - y))
            for l in deltaWidthYX:
                L.remove(l)
            log_step(abs(width - y), L, X, width - y, step_filename)  # Log the step
            place(L, X, step_filename)
            # Backtrack
            X.remove(abs(width - y))
            for l in deltaWidthYX:
                L.append(l)

if __name__ == "__main__":

    # Given Case 1
    L1 = [2, 2, 3, 3, 4, 5, 6, 7, 8, 10]
    
    # Test Case for Backtracking
    L3 = [1, 1, 2, 2, 2, 3, 3, 4, 4, 5, 5, 5,6, 7, 7, 7, 8, 9, 10, 11, 12]

    # Call partial_digest for each case and output to specific files
    partial_digest(L1, "sol1.txt", "step1.txt")
    partial_digest(L3, "sol3.txt", "step3.txt")
