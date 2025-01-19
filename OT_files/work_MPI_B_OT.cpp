#include <iostream>
#include <ctime>
#include <cmath>
#include <fstream>
#include <mpi.h>
#include <cstdlib>
using namespace std;

// Function to allocate memory for a matrix
float** allocate_matrix(int N) {
    float** matrix = new float*[N];
    for (int i = 0; i < N; i++) {
        matrix[i] = new float[N];
    }
    return matrix;
}

// Function to free the memory of a matrix
void free_matrix(float** matrix, int N) {
    for (int i = 0; i < N; i++) {
        delete[] matrix[i];
    }
    delete[] matrix;
}

void init_matrix(float** M, int N) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            M[i][j] = static_cast<float>(rand()) / RAND_MAX; // Random value in [0.0, 1.0]
        }
    }
}

bool checkSymMPI(float** M, int N, int rank, int size) {
    bool isSymmetric = true;
    for (int i = rank; i < N; i += size) {
        for (int j = i + 1; j < N; ++j) {
            if (M[i][j] != M[j][i]) {
                isSymmetric = false;
                break;
            }
        }
    }

    // Combine results from all processes
    bool globalSymmetric;
    MPI_Allreduce(&isSymmetric, &globalSymmetric, 1, MPI_C_BOOL, MPI_LAND, MPI_COMM_WORLD);
    return globalSymmetric;
}

void matTransposeMPI(float** M, float** T, int N, int rank, int size) {
    // Check symmetry within matTransposeMPI
    bool isSymmetric = checkSymMPI(M, N, rank, size);
    if (rank == 0) {
        cout << "Symmetry check: " << (isSymmetric ? "True" : "False") << endl;
    }

    // Perform matrix transposition
    for (int i = rank; i < N; i += size) {
        for (int j = 0; j < N; ++j) {
            T[j][i] = M[i][j];
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
}

int compare(const void* a, const void* b) {
    return (*(double*)a > *(double*)b) - (*(double*)a < *(double*)b);
}

void printMatrix(float** matrix, int N, const char* label, FILE* file) {
    fprintf(file, "%s:\n", label);
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            fprintf(file, "%.4f ", matrix[i][j]);
        }
        fprintf(file, "\n");
    }
    fprintf(file, "\n");
}

void printValues(double* results, FILE* file, int n_times, int N, int size, float** M, float** T) {
    qsort(results, n_times, sizeof(double), [](const void* a, const void* b) {
        return (*(double*)a > *(double*)b) - (*(double*)a < (*(double*)b));
    });
    for (int i = (n_times / 10) * 3; i < (n_times / 10) * 7; i++) {
        fprintf(file, "MPI_B,%d,%d,%.9f\n", size, N, results[i]);
    }

    printMatrix(M, N, "Original Matrix (M) in Bcast:", file);
    printMatrix(T, N, "Transposed Matrix (T)in Bcast:", file);
}

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (argc != 3) {
        if (rank == 0) {
            cerr << "Usage: " << argv[0] << " <matrix_exponential> <num_times>" << endl;
        }
        MPI_Finalize();
        return 1;
    }

    int matrix_exponential = atoi(argv[1]);
    int N = pow(2, matrix_exponential);
    if (N < 2) {
        if (rank == 0) {
            cerr << "Matrix exponential must be a positive integer between 4 and 12" << endl;
        }
        MPI_Finalize();
        return 1;
    }

    int n_times = atoi(argv[2]);
    if (n_times < 1) {
        if (rank == 0) {
            cerr << "The number of times must be a positive integer" << endl;
        }
        MPI_Finalize();
        return 1;
    }

    if (size > N) {
        MPI_Finalize();
        return 0;
    }

    // Adjust size if 2^N is not divisible by size
    int adjusted_size = size;
    if (N % size != 0) {
        adjusted_size = N / size * size; // Largest multiple of size <= N
    }

    MPI_Comm new_comm;
    MPI_Comm_split(MPI_COMM_WORLD, rank < adjusted_size ? 0 : MPI_UNDEFINED, rank, &new_comm);
    if (rank >= adjusted_size) {
        MPI_Finalize();
        return 0;
    }
    MPI_Comm_rank(new_comm, &rank);
    MPI_Comm_size(new_comm, &size);

    srand(static_cast<unsigned int>(time(0) + rank));

    float** M = allocate_matrix(N);
    float** T = allocate_matrix(N);

    if (rank == 0) {
        init_matrix(M, N);
    }

    // Broadcast the matrix to all processes
    for (int i = 0; i < N; ++i) {
        MPI_Bcast(M[i], N, MPI_FLOAT, 0, new_comm);
    }

    double* uresults = new double[n_times];

    for (int i = 0; i < n_times; i++) {
        double start = MPI_Wtime();

        // Transpose the matrix (including symmetry check inside)
        matTransposeMPI(M, T, N, rank, size);

        double end = MPI_Wtime();
        uresults[i] = end - start;
    }

    if (rank == 0) {
        const char* result = "result.csv";
        FILE* file = fopen(result, "a");
        if (file == NULL) {
            cerr << "Error opening file" << endl;
            free_matrix(M, N);
            free_matrix(T, N);
            delete[] uresults;
            MPI_Finalize();
            return 1;
        }

        printValues(uresults, file, n_times, N, size, M, T);
        fclose(file);
    }

    free_matrix(M, N);
    free_matrix(T, N);
    delete[] uresults;

    MPI_Finalize();
    return 0;
}
