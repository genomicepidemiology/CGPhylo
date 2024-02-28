def load_phylip_matrix(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()
        size = int(lines[0].strip())
        matrix = [[0.0 for _ in range(size)] for _ in range(size)]
        for i, line in enumerate(lines[1:]):
            parts = line.strip().split()
            for j, distance in enumerate(parts[1:], start=1):
                matrix[i][j-1] = float(distance)
                matrix[j-1][i] = float(distance)  # Ensure symmetric assignment
    return matrix

def normalize_matrix(matrix):
    max_distance = max(max(row) for row in matrix)
    normalized_matrix = [[distance / max_distance for distance in row] for row in matrix]
    return normalized_matrix

def calculate_rmse(matrix1, matrix2):
    size = len(matrix1)
    error_sum = 0
    for i in range(size):
        for j in range(size):
            error_sum += (matrix1[i][j] - matrix2[i][j]) ** 2
    mean_error = error_sum / (size * size)
    rmse = mean_error ** 0.5
    return rmse

# Placeholder filenames with the provided paths
filename1 = '/home/people/malhal/papers/cgphylo/mashtree/aureus_ont/corrected_matrix'
filename2 = '/home/people/malhal/papers/cgphylo/cgphylo/results/aureus_ont/corrected_distance_matrix_1M.txt'

# Load the matrices
matrix1 = load_phylip_matrix(filename1)
matrix2 = load_phylip_matrix(filename2)

# Normalize the matrices
normalized_matrix1 = normalize_matrix(matrix1)
normalized_matrix2 = normalize_matrix(matrix2)

# Calculate RMSE
rmse = calculate_rmse(normalized_matrix1, normalized_matrix2)

print("RMSE between normalized matrices:", rmse)
