import scipy.io

# Replace 'path_to_file' with the path to your .mat file
file_path = 'Escherichia_coli_str_K_12_substr_MG1655.mat'
mat = scipy.io.loadmat(file_path)

# Remove '__header__', '__version__', and '__globals__' entries
mat = {key: val for key, val in mat.items() if not key.startswith('__')}

# Print each variable and its contents
for key in mat:
    print(f"{key}:")
    print(mat[key])
    print("\n")

