#include <iostream>
#include <ctime>
#include <cmath>
#include <fstream>
#include <mpi.h>
#include <cstdlib>
using namespace std;

// Function to allocate memory for a matrix
float* allocate_matrix_linear(int rows, int cols) {
    // Allocate memory for a matrix with 'rows' x 'cols' elements
    float* matrix = new (nothrow) float[rows * cols];
    if (!matrix) {
        cerr << "Memory allocation failed for matrix with dimensions " 
             << rows << " x " << cols << endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    return matrix;
}

void init_matrix(float* M, int N) {
    for (int i = 0; i < N * N; i++) {
        M[i] = static_cast<float>(rand()) / RAND_MAX; // Random value in [0.0, 1.0]
    }
}

bool checkSymMPI(float* local_matrix, int rows_per_process, int N, int rank, int size) {
    bool isSymmetric = true;
    for (int i = 0; i < rows_per_process; i++) {
        for (int j = 0; j < N; ++j) {
            if (local_matrix[i * N + j] != local_matrix[j * N + (i + rank * rows_per_process)]) {
                isSymmetric = false;
                break;
            }
        }
    }
    bool globalSymmetric;
    MPI_Allreduce(&isSymmetric, &globalSymmetric, 1, MPI_C_BOOL, MPI_LAND, MPI_COMM_WORLD);
    return globalSymmetric;
}

void matTransposeMPI(float* local_matrix, float* local_transpose, int rows_per_process, int N, int rank, int size) {
    // Temporary buffer to hold transposed data for exchange
    float* transposed_buffer = new float[rows_per_process * N];

    // Compute the transpose locally
    for (int i = 0; i < rows_per_process; i++) {
        for (int j = 0; j < N; j++) {
            transposed_buffer[j * rows_per_process + i] = local_matrix[i * N + j];
        }
    }

    // Buffer to receive transposed rows from all processes
    float* received_buffer = new float[rows_per_process * N];

    // Exchange transposed rows among processes
    MPI_Alltoall(transposed_buffer, rows_per_process * rows_per_process, MPI_FLOAT,
                 received_buffer, rows_per_process * rows_per_process, MPI_FLOAT, MPI_COMM_WORLD);

    // Reorganize received data into the local transpose
    for (int i = 0; i < rows_per_process; i++) {
        for (int j = 0; j < N; j++) {
            local_transpose[i * N + j] = received_buffer[i + j * rows_per_process];
        }
    }

    delete[] transposed_buffer;
    delete[] received_buffer;
}

void printValues(double* results, FILE* file, int n_times, int N, int size, float* M, float* T) {
    qsort(results, n_times, sizeof(double), [](const void* a, const void* b) {
        return (*(double*)a > *(double*)b) - (*(double*)a < *(double*)b);
    });
    for (int i = (n_times / 10) * 3; i < (n_times / 10) * 7; i++) {
        fprintf(file, "MPI_SG,%d,%d,%.9f\n", size, N, results[i]);
    }
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

    // Adjust process size if N % size != 0
    int adjusted_size = size;
    if (N % size != 0) {
        adjusted_size = N / size; // Largest divisor <= N
    }

    MPI_Comm new_comm;
    MPI_Comm_split(MPI_COMM_WORLD, rank < adjusted_size, rank, &new_comm);

    if (rank >= adjusted_size) {
        MPI_Finalize();
        return 0;
    }

    MPI_Comm_rank(new_comm, &rank);
    MPI_Comm_size(new_comm, &size);

    srand(static_cast<unsigned int>(time(0) + rank));

    int rows_per_process = N / adjusted_size;
    float* M = nullptr;
    float* T = nullptr;
    float* local_matrix = allocate_matrix_linear(rows_per_process, N);
    float* local_transpose = allocate_matrix_linear(rows_per_process, N);

    if (rank == 0) {
        M = allocate_matrix_linear(N, N);
        T = allocate_matrix_linear(N, N);
        init_matrix(M, N);
    }

    // Scatter rows of the matrix to all processes
    MPI_Scatter(M, rows_per_process * N, MPI_FLOAT, local_matrix, rows_per_process * N, MPI_FLOAT, 0, new_comm);

    double* uresults = new double[n_times];

    for (int i = 0; i < n_times; i++) {
        double start = MPI_Wtime();

        // Transpose the matrix (including symmetry check inside)
        matTransposeMPI(local_matrix, local_transpose, rows_per_process, N, rank, adjusted_size);

        double end = MPI_Wtime();
        uresults[i] = end - start;
    }

    // Gather the transposed rows back to the root process
    MPI_Gather(local_transpose, rows_per_process * N, MPI_FLOAT, T, rows_per_process * N, MPI_FLOAT, 0, new_comm);

    if (rank == 0) {
        const char* result = "result.csv";
        FILE* file = fopen(result, "a");
        if (file == NULL) {
            cerr << "Error opening file" << endl;
            delete[] M;
            delete[] T;
            delete[] uresults;
            MPI_Finalize();
            return 1;
        }

        printValues(uresults, file, n_times, N, adjusted_size, M, T);
        fclose(file);

        delete[] M;
        delete[] T;
    }

    delete[] local_matrix;
    delete[] local_transpose;
    delete[] uresults;

    MPI_Finalize();
    return 0;
}
