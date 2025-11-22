#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <chrono>
#include <iomanip>

using namespace std::chrono;

#define Vector std::vector
#define String std::string

// iterations = 4972
// time = 134 seconds

// iterations = 3213
// time = 116 seconds

const double Lx = 1;
const double Ly = 1;

const int N = 100;
const int M = 100;

const double Re = 50;
const double rho = 1;
const double dt = 0.001;
const double eps = 1e-5;
const double eps_P = 1e-5;
const int stop_iteration = 1e4;
const int stop_iteration_P = 1e4;


void save_1d_vector(const String path_to_file, const Vector<double> &arr,
                    int precision=6, String sep=" ", String end="\n") {
    
    std::ofstream output_file(path_to_file);

    if (!output_file.is_open()) {
        std::cout << "I can't save vector to the file " 
                  << path_to_file << " :( "<< std::endl;
        return;
    }

    int size = arr.size();
    for (int i = 0; i < size; ++i) {
        output_file << arr[i] << (i < size - 1 ? sep : end);
    }

    output_file.close();
    std::cout << "The vector succesfully saved to the file " 
              << path_to_file << " :) " << std::endl;
}

void save_2d_vector(const String path_to_file, const Vector<Vector<double>> &arr,
                    int precision=6, String sep=" ", String end="\n") {

    std::ofstream file(path_to_file);

    if (!file.is_open()) {
        std::cout << "I can't save matrix to the file " 
                  << path_to_file << " :( "<< std::endl;
        return;
    }

    int rows = arr.size();
    int cols = arr[0].size();
    for (int j = 0; j < rows; ++j) {
        for (int i = 0; i < cols; ++i) {
            file << std::fixed << std::setprecision(precision) 
                 << arr[j][i] << (i < cols - 1 ? sep : end);
        }
    }

    file.close();
    std::cout << "The matrix is succesfully saved to the file " 
              << path_to_file << " :) " << std::endl;
}

void print_1d_vector(const Vector<double> &U, String sep=" ", String end="\n") {
    int size = U.size();
    for (int i = 0; i < size; ++i) {
        std::cout << U[i] << (i < size - 1 ? sep : end);
    }
}

void print_2d_vector(const Vector<Vector<double>> &U, String sep=" ", String end="\n") {
    int rows = U.size();
    int cols = U[0].size();
    for (int j = 0; j < rows; ++j) {
        for (int i = 0; i < cols; ++i) {
            std::cout << U[j][i] << (i < cols - 1 ? sep : end);
        }
    }             
}

double max_abs_diff(const Vector<Vector<double>> &A, const Vector<Vector<double>> &B) {
    int rows = A.size();
    int cols = A[0].size();

    if (rows != B.size() or cols != B[0].size()) {
        return -1;
    }

    double maximum = 0;
    for (int j = 0; j < rows; ++j) {
        for (int i = 0; i < cols; ++i) {
            maximum = fmax(maximum, fabs(A[j][i] - B[j][i]));
        }
    }
    return maximum;
}


void Burger_boundary_condition(Vector<Vector<double>> &U, Vector<Vector<double>> &V) {
    // U(t, 0<x<0.4, y=0) = 0
    // V(t, 0<x<0.4, y=0) = 0
    for (int i = 0; i < (0.4 / Lx) * N; ++i) {
        U[0][i] = 0;
        V[0][i] = 0;
    }
    // Uy(t, 0.4<x<0.7, y=0) = 0
    // Vy(t, 0.4<x<0.7, y=0) = 0    
    for (int i = (0.4 / Lx) * N; i < (0.7 / Lx) * N; ++i) {
        U[0][i] = U[1][i];
        V[0][i] = V[1][i];
    }
    // U(t, 0.7<x<1, y=0) = 0
    // V(t, 0.7<x<1, y=0) = 0
    for (int i = (0.7 / Lx) * N; i < N; ++i) {
        U[0][i] = 0;
        V[0][i] = 0;
    }
    // Uy(t, 0<x<0.3, y=1) = 0
    // Vy(t, 0<x<0.3, y=1) = 0
    for (int i = 0; i < (0.3 / Lx) * N; ++i) {
        U[M-1][i] = U[M-2][i];
        V[M-1][i] = V[M-2][i];
    }
    // U(t, 0.3<x<1, y=1) = 0
    // V(t, 0.3<x<1, y=1) = 0
    for (int i = (0.3 / Lx) * N; i < N; ++i) {
        U[M-1][i] = 0;
        V[M-1][i] = 0;
    }
    // U(t, x=0, 0<y<1) = 0
    // V(t, x=0, 0<y<1) = 0
    for (int j = 0; j < M; ++j) {
        U[j][0] = 0;
        V[j][0] = 0;
    }
    // U(t, x=1, 0<y<0.7) = 0
    // V(t, x=1, 0<y<0.7) = 0
    for (int j = 0; j < (0.7 / Ly) * M; ++j) {
        U[j][N-1] = 0;
        V[j][N-1] = 0;
    }
    // U(t, x=1, 0.7<y<1) = -1
    // V(t, x=1, 0.7<y<1) = 0
    for (int j = (0.7 / Ly) * M; j < M; ++j) {
        U[j][N-1] = -1;
        V[j][N-1] = 0;
    }
}

void Poisson_boundary_condition(Vector<Vector<double>> &P) {
    // Py(t, 0<x<0.4, y=0) = 0
    for (int i = 0; i < (0.4 / Lx) * N; ++i) {
        P[0][i] = P[1][i];
    }
    // P(t, 0.4<x<0.7, y=0) = 0  
    for (int i = (0.4 / Lx) * N; i < (0.7 / Lx) * N; ++i) {
        P[0][i] = 0;
    }
    // Py(t, 0.7<x<1, y=0) = 0
    for (int i = (0.7 / Lx) * N; i < N; ++i) {
        P[0][i] = P[1][i];
    }
    // P(t, 0<x<0.3, y=1) = 0
    for (int i = 0; i < (0.3 / Lx) * N; ++i) {
        P[M-1][i] = 0;
    }
    // Py(t, 0.3<x<1, y=1) = 0
    for (int i = (0.3 / Lx) * N; i < N; ++i) {
        P[M-1][i] = P[M-2][i];
    }
    // Px(t, x=0, 0<y<1) = 0
    for (int j = 0; j < M; ++j) {
        P[j][0] = P[j][1];
    }
    // Px(t, x=1, 0<y<0.7) = 0
    for (int j = 0; j < (0.7 / Ly) * M; ++j) {
        P[j][N-1] = P[j][N-2];
    }
    // P(t, x=1, 0.7<y<1) = 1
    for (int j = (0.7 / Ly) * M; j < M; ++j) {
        P[j][N-1] = 1;
    }
}

void Explicit_Burger(Vector<Vector<double>> &S_new, Vector<Vector<double>> &S_old, 
                     Vector<Vector<double>> &U_old, Vector<Vector<double>> &V_old, 
                        int N, int M, double dx, double dy, double dt, double nu){
    
    for(int j = 1; j < M-1; ++j){
        for(int i = 1; i < N-1; ++i){
            // central scheme
            S_new[j][i] = S_old[j][i] + dt*(
                - U_old[j][i]*(S_old[j][i+1] - S_old[j][i-1]) 
                                / (2*dx)
                - V_old[j][i]*(S_old[j+1][i] - S_old[j-1][i]) 
                                / (2*dy)
                + nu*(
                    (S_old[j][i+1] - 2*S_old[j][i] + S_old[j][i-1]) 
                        / (dx*dx) 
                    + (S_old[j+1][i] - 2*S_old[j][i] + S_old[j-1][i])
                        / (dy*dy))
            );
        }
    }
}


void Jacobi_Poisson(Vector<Vector<double>> &P_new, Vector<Vector<double>> &P_old,
                    Vector<Vector<double>> &U_new, Vector<Vector<double>> &V_new, 
                    int N, int M, double dx, double dy, double dt, double rho) {
    
    for(int j = 1; j < M-1; ++j){
        for(int i = 1; i < N-1; ++i){
            // central scheme
            P_new[j][i] = (
                dy*dy*(P_old[j][i+1] + P_old[j][i-1])
                + dx*dx*(P_old[j+1][i] + P_old[j-1][i])) 
                    / (2*(dx*dx + dy*dy))
                - dx*dx*dy*dy*rho
                    / (2*dt*(dx*dx + dy*dy))
                * ((U_new[j][i+1] - U_new[j][i-1]) 
                        / (2*dx) 
                    + (V_new[j+1][i] - V_new[j-1][i]) 
                        / (2*dy)
            );
        }
    }
}



int main() {
    const double dx = Lx / (N - 1);
    const double dy = Ly / (M - 1);

    Vector<Vector<double>> U_old(M, Vector<double>(N, 0));
    Vector<Vector<double>> U_new(M, Vector<double>(N, 0));

    Vector<Vector<double>> V_old(M, Vector<double>(N, 0));
    Vector<Vector<double>> V_new(M, Vector<double>(N, 0));

    Vector<Vector<double>> P_old(M, Vector<double>(N, 0));
    Vector<Vector<double>> P_new(M, Vector<double>(N, 0));

    auto start = high_resolution_clock::now();

    Burger_boundary_condition(U_old, V_old);

    int iteration = 0;
    double maximum = 0;
    int iteration_P = 0;
    double maximum_P = 0;
    do {
        Explicit_Burger(U_new, U_old, U_old, V_old, N, M, dx, dy, dt, 1.0 / Re);
        Explicit_Burger(V_new, V_old, U_old, V_old, N, M, dx, dy, dt, 1.0 / Re);

        Burger_boundary_condition(U_new, V_new);


        Poisson_boundary_condition(P_old);

        iteration_P = 0;
        maximum_P = 0;
        do {
            Jacobi_Poisson(P_new, P_old, U_new, V_new, N, M, dx, dy, dt, rho);

            Poisson_boundary_condition(P_new);

            maximum_P = max_abs_diff(P_old, P_new);

            P_old = P_new;

            iteration_P++;
        } while (iteration_P < stop_iteration_P and maximum_P > eps_P);

        // Do not change
        // 3. Correction for U^(n+1), V^(n+1)
        for(int j = 1; j < M-1; ++j){
            for(int i = 1; i < N-1; ++i){
                U_new[j][i] = U_new[j][i] - dt / rho * (P_new[j][i+1] - P_new[j][i-1]) / (2*dx);
                V_new[j][i] = V_new[j][i] - dt / rho * (P_new[j+1][i] - P_new[j-1][i]) / (2*dy);
            }
        }

        Burger_boundary_condition(U_new, V_new);

        maximum = max_abs_diff(U_old, U_new);
        maximum = fmax(maximum, max_abs_diff(V_old, V_new));

        if (iteration % 10 == 0) {
            std::cout << iteration << " " << maximum << std::endl;
        }

        U_old = U_new;
        V_old = V_new;

        iteration++;
    } while (iteration < stop_iteration and maximum > eps);

    std::cout << iteration << " " << maximum << std::endl;

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    printf("calculation time: %.6f seconds\n", duration.count() / 1e6);

    Vector<Vector<double>> X(M, Vector<double>(N, 0));
    Vector<Vector<double>> Y(M, Vector<double>(N, 0));

    Vector<double> x(N, 0);
    Vector<double> y(M, 0);

    for (int i = 0; i < N; ++i) {
        x[i] = i*dx;
    }

    for (int j = 0; j < M; ++j) {
        y[j] = j*dy;
    }

    for (int j = 0; j < M; ++j) {
        for (int i = 0; i < N; ++i) {
            X[j][i] = i*dx;
            Y[j][i] = j*dy;
        }
    }

    save_2d_vector("U.txt", U_new, 6, "\t", "\n");
    save_2d_vector("V.txt", V_new, 6, "\t", "\n");
    save_2d_vector("P.txt", P_new, 6, "\t", "\n");
    // save_2d_vector("X.txt", X, 6, "\t", "\n");
    // save_2d_vector("Y.txt", Y, 6, "\t", "\n");

    save_1d_vector("x.txt", x, 6, "\n");
    save_1d_vector("y.txt", y, 6, "\n");

    return 0;
}