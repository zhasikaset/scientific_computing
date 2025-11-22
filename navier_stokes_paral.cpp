#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <chrono>
#include <iomanip>
#include <mpi.h>

using namespace std::chrono;

#define Vector std::vector
#define String std::string

// iterations = 4972
// time = 50 seconds

// iterations = 3213
// time = 46 seconds

const double Lx = 1;
const double Ly = 1;

const int N = 100;
const int M = 100;

const double Re = 50;
const double nu = 1.0 / Re;
const double rho = 1;
const double dt = 0.001;
const double eps = 1e-5;
const double eps_P = 1e-5;
const int stop_iteration = 1e4;
const int stop_iteration_P = 1e4;

// Количество блоков
const int Px = 2;
const int Py = 2;

double max_abs_diff(const Vector<double> &A, const Vector<double> &B,
                        int N_local, int M_local) {

    double maximum = 0;
    for (int j = 1; j < M_local + 1; ++j) {
        for (int i = 1; i < N_local + 1; ++i) {
            maximum = fmax(maximum, fabs(A[(j)*(N_local + 2) + (i)] - B[(j)*(N_local + 2) + (i)]));
        }
    }

    return maximum;
}

void Burger_boundary_condition(Vector<double> &U, Vector<double> &V, int N_local, int M_local, 
                                int x_start, int y_start, int coords[], 
                                int west, int east, int north, int south) {

    // Нижняя граница (по y=0)
    if (coords[0] == 0) {
        for (int i = 1; i < N_local + 1; ++i) {
            if (i - 1 + x_start < (0.4 / Lx) * N) {
                U[(1)*(N_local + 2) + (i)] = 0;
                V[(1)*(N_local + 2) + (i)] = 0;
            } else if (i - 1 + x_start < (0.7 / Lx) * N) {
                U[(1)*(N_local + 2) + (i)] = U[(2)*(N_local + 2) + (i)];
                V[(1)*(N_local + 2) + (i)] = V[(2)*(N_local + 2) + (i)];
            } else {
                U[(1)*(N_local + 2) + (i)] = 0;
                V[(1)*(N_local + 2) + (i)] = 0;
            }
        }
    }

    // Верхняя граница (по y=1)
    if (coords[0] == Py - 1) {
        for (int i = 1; i < N_local + 1; ++i) {
            if (i - 1 + x_start < (0.3 / Lx) * N) {
                U[(M_local)*(N_local + 2) + (i)] = U[(M_local - 1)*(N_local + 2) + (i)];
                V[(M_local)*(N_local + 2) + (i)] = V[(M_local - 1)*(N_local + 2) + (i)];
            } else {
                U[(M_local)*(N_local + 2) + (i)] = 0;
                V[(M_local)*(N_local + 2) + (i)] = 0;
            }
        }
    }

    // Левая граница (x=0)
    if (coords[1] == 0) {
        for (int j = 1; j < M_local + 1; ++j) {
            U[(j)*(N_local + 2) + (1)] = 0;
            V[(j)*(N_local + 2) + (1)] = 0;
        }
    }

    // Правая граница (x=1)
    if (coords[1] == Px - 1) {
        for (int j = 1; j < M_local + 1; ++j) {
            if (j - 1 + y_start < (0.7 / Ly) * M) {
                U[(j)*(N_local + 2) + (N_local)] = 0;
                V[(j)*(N_local + 2) + (N_local)] = 0;
            } else {
                U[(j)*(N_local + 2) + (N_local)] = -1;
                V[(j)*(N_local + 2) + (N_local)] = 0;
            }
        }
    }
}

void Poisson_boundary_condition(Vector<double> &P, int N_local, int M_local, 
                                int x_start, int y_start, int coords[], 
                                int west, int east, int north, int south) {

    if (coords[0] == 0) {
        for (int i = 1; i < N_local + 1; ++i) {
            if (i - 1 + x_start < (0.4 / Lx) * N) {
                P[(1)*(N_local + 2) + (i)] = P[(2)*(N_local + 2) + (i)];
            } 
            else if (i - 1 + x_start < (0.7 / Lx) * N) {
                P[(1)*(N_local + 2) + (i)] = 0;
            } else {
                P[(1)*(N_local + 2) + (i)] = P[(2)*(N_local + 2) + (i)];
            }
        }
    }

    // Верхняя граница (по y=1)
    if (coords[0] == Py - 1) {
        for (int i = 1; i < N_local + 1; ++i) {
            if (i - 1 + x_start < (0.3 / Lx) * N) {                
                P[(M_local)*(N_local + 2) + (i)] = 0;
            } else {
                P[(M_local)*(N_local + 2) + (i)] = P[(M_local - 1)*(N_local + 2) + (i)];
            }
        }
    }

    // Левая граница (x=0)
    if (coords[1] == 0) {
        for (int j = 1; j < M_local + 1; ++j) {
            P[(j)*(N_local + 2) + (1)] = P[(j)*(N_local + 2) + (2)];
        }
    }

    // Правая граница (x=1)
    if (coords[1] == Px - 1) {
        for (int j = 1; j < M_local + 1; ++j) {
            if (j - 1 + y_start < (0.7 / Ly) * M) {
                P[(j)*(N_local + 2) + (N_local)] = P[(j)*(N_local + 2) + (N_local - 1)];
            } else {
                P[(j)*(N_local + 2) + (N_local)] = 1;
            }
        }
    }            
}

void exchange_boundary_points(Vector<double> &U, int N_local, int M_local,
                                int west, int east, int north, int south, MPI_Comm cart_comm) {

    Vector<double> send_row(N_local, 0);
    Vector<double> recv_row(N_local, 0);

    Vector<double> send_col(M_local, 0);
    Vector<double> recv_col(M_local, 0);

    // отправка верхнему соседу
    if (north != -1) {
        for (int i = 1; i < N_local + 1; ++i) {
            send_row[i-1] = U[(1)*(N_local + 2) + (i)];
        }

        MPI_Sendrecv(send_row.data(), N_local, MPI_DOUBLE, north, 0, 
                        recv_row.data(), N_local, MPI_DOUBLE, north, 0, 
                            cart_comm, MPI_STATUS_IGNORE);

        for (int i = 1; i < N_local + 1; ++i) {
            U[(0)*(N_local + 2) + (i)] = recv_row[i-1];
        }
    }

    // отправка нижнему соседу
    if (south != -1) {
        for (int i = 1; i < N_local + 1; ++i) {
            send_row[i-1] = U[(M_local)*(N_local + 2) + (i)];
        }

        MPI_Sendrecv(send_row.data(), N_local, MPI_DOUBLE, south, 0, 
                        recv_row.data(), N_local, MPI_DOUBLE, south, 0, 
                            cart_comm, MPI_STATUS_IGNORE);
    
        for (int i = 1; i < N_local + 1; ++i) {
            U[(M_local + 1)*(N_local + 2) + (i)] = recv_row[i-1];
        }
    }

    // Отправка левому соседу
    if (west != -1) {
        for (int j = 1; j < M_local + 1; ++j) {
            send_col[j-1] = U[(j)*(N_local + 2) + (1)];
        }

        MPI_Sendrecv(send_col.data(), M_local, MPI_DOUBLE, west, 0, 
                        recv_col.data(), M_local, MPI_DOUBLE, west, 0, 
                            cart_comm, MPI_STATUS_IGNORE);

        for (int j = 1; j < M_local + 1; ++j) {
            U[(j)*(N_local + 2) + (0)] = recv_col[j-1];
        }
    }

    // Отправка правому соседу
    if (east != -1) {
        for (int j = 1; j < M_local + 1; ++j) {
            send_col[j-1] = U[(j)*(N_local + 2) + (N_local)];
        }

        MPI_Sendrecv(send_col.data(), M_local, MPI_DOUBLE, east, 0, 
                        recv_col.data(), M_local, MPI_DOUBLE, east, 0, 
                            cart_comm, MPI_STATUS_IGNORE);

        for (int j = 1; j < M_local + 1; ++j) {
            U[(j)*(N_local + 2) + (N_local + 1)] = recv_col[j-1];
        }
    }
}

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
        output_file << std::fixed << std::setprecision(precision) 
                    << arr[i] << (i < size - 1 ? sep : end);
    }

    output_file.close();
    std::cout << "The vector succesfully saved to the file " 
              << path_to_file << " :) " << std::endl;
}

void save_2d_vector(const String path_to_file, const Vector<double> &arr,
                    int rows, int cols, int precision=6, String sep=" ", String end="\n") {

    std::ofstream file(path_to_file);

    if (!file.is_open()) {
        std::cout << "I can't save matrix to the file " 
                  << path_to_file << " :( "<< std::endl;
        return;
    }

    for (int j = 0; j < rows; ++j) {
        for (int i = 0; i < cols; ++i) {
            file << std::fixed << std::setprecision(precision) 
                 << arr[j * cols + i] << (i < cols - 1 ? sep : end);
        }
    }

    file.close();
    std::cout << "The matrix is succesfully saved to the file " 
              << path_to_file << " :) " << std::endl;
}


int main() {
    const double dx = Lx / (N - 1);
    const double dy = Ly / (M - 1);

    MPI_Init(NULL, NULL);
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Status status;
    const int master = 0;

    int rank, numprocs;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &numprocs);

    if (numprocs != Px * Py) {
        if (rank == 0) {
            std::cout << "Error: size does not match Px * Py" << std::endl;
        }
        MPI_Abort(comm, -1);
    }

    // rows and colums
    int dims[2] = {Py, Px};
    int periods[2] = {0, 0};
    int reorder = true;

    MPI_Comm cart_comm;
    MPI_Cart_create(comm, 2, dims, periods, reorder, &cart_comm);

    int coords[2];
    MPI_Cart_coords(cart_comm, rank, 2, coords);

    int rank_local;
    MPI_Cart_rank(cart_comm, coords, &rank_local);

    int west, east, north, south;

    // нумерация как обычно слева направо
    MPI_Cart_shift(cart_comm, 1, 1, &west, &east);
    MPI_Cart_shift(cart_comm, 0, 1, &north, &south);

    int N_local = N / Px;
    int M_local = M / Py;

    int global_i_start = coords[1] * N_local;
    int global_j_start = coords[0] * M_local;

    Vector<double> U_old((M_local + 2) * (N_local + 2), 0);
    Vector<double> U_new((M_local + 2) * (N_local + 2), 0);

    Vector<double> V_old((M_local + 2) * (N_local + 2), 0);
    Vector<double> V_new((M_local + 2) * (N_local + 2), 0);

    Vector<double> P_old((M_local + 2) * (N_local + 2), 0);
    Vector<double> P_new((M_local + 2) * (N_local + 2), 0);

    double start_time, end_time;
    if (rank == master) {
        start_time = MPI_Wtime();
    }

    Burger_boundary_condition(U_old, V_old, N_local, M_local, global_i_start, global_j_start, coords, west, east, north, south);

    int iteration = 0;
    double maximum = 0;
    int iteration_P = 0;
    double maximum_P = 0;
    do {
        // 1. Решение Бюргерса
        // 1.1. Обмен сообщениями
        // Отправка верхнему соседу
        exchange_boundary_points(U_old, N_local, M_local, west, east, north, south, cart_comm);
        exchange_boundary_points(V_old, N_local, M_local, west, east, north, south, cart_comm);

        // 1.2. Явная схема Бюргерса
        for(int j = 1; j < M_local + 1; ++j){
            for(int i = 1; i < N_local + 1; ++i){
                // central scheme
                U_new[(j)*(N_local + 2) + (i)] = U_old[(j)*(N_local + 2) + (i)] + dt*(
                    - U_old[(j)*(N_local + 2) + (i)]*(U_old[(j)*(N_local + 2) + (i+1)] - U_old[(j)*(N_local + 2) + (i-1)]) 
                                    / (2*dx)
                    - V_old[(j)*(N_local + 2) + (i)]*(U_old[(j+1)*(N_local + 2) + (i)] - U_old[(j-1)*(N_local + 2) + (i)]) 
                                    / (2*dy)
                    + nu*(
                        (U_old[(j)*(N_local + 2) + (i+1)] - 2*U_old[(j)*(N_local + 2) + (i)] + U_old[(j)*(N_local + 2) + (i-1)]) 
                            / (dx*dx) 
                        + (U_old[(j+1)*(N_local + 2) + (i)] - 2*U_old[(j)*(N_local + 2) + (i)] + U_old[(j-1)*(N_local + 2) + (i)])
                            / (dy*dy))
                );

                V_new[(j)*(N_local + 2) + (i)] = V_old[(j)*(N_local + 2) + (i)] + dt*(
                    - U_old[(j)*(N_local + 2) + (i)]*(V_old[(j)*(N_local + 2) + (i+1)] - V_old[(j)*(N_local + 2) + (i-1)]) 
                                    / (2*dx)
                    - V_old[(j)*(N_local + 2) + (i)]*(V_old[(j+1)*(N_local + 2) + (i)] - V_old[(j-1)*(N_local + 2) + (i)]) 
                                    / (2*dy)
                    + nu*(
                        (V_old[(j)*(N_local + 2) + (i+1)] - 2*V_old[(j)*(N_local + 2) + (i)] + V_old[(j)*(N_local + 2) + (i-1)]) 
                            / (dx*dx) 
                        + (V_old[(j+1)*(N_local + 2) + (i)] - 2*V_old[(j)*(N_local + 2) + (i)] + V_old[(j-1)*(N_local + 2) + (i)])
                            / (dy*dy))
                );
            }
        }

        // 1.3. Ставим граничные условия
        Burger_boundary_condition(U_new, V_new, N_local, M_local, global_i_start, global_j_start, 
                                    coords, west, east, north, south);

        // 1.4. Обмен сообщениями
        // Отправка верхнему соседу
        exchange_boundary_points(U_new, N_local, M_local, west, east, north, south, cart_comm);
        exchange_boundary_points(V_new, N_local, M_local, west, east, north, south, cart_comm);

        // 2. Решение Пуассона
        // Граничные условия Пуассона
        Poisson_boundary_condition(P_old, N_local, M_local, global_i_start, global_j_start, 
                                    coords, west, east, north, south);

        iteration_P = 0;
        maximum_P = 0;
        do {
            // 2.1. Обмен сообщениями
            exchange_boundary_points(P_old, N_local, M_local, west, east, north, south, cart_comm);

            // 2.2. Метод Якоби
            for(int j = 1; j < M_local + 1; ++j){
                for(int i = 1; i < N_local + 1; ++i){
                    // central scheme
                    P_new[(j)*(N_local + 2) + (i)] = (
                        dy*dy*(P_old[(j)*(N_local + 2) + (i+1)] + P_old[(j)*(N_local + 2) + (i-1)])
                        + dx*dx*(P_old[(j+1)*(N_local + 2) + (i)] + P_old[(j-1)*(N_local + 2) + (i)])) 
                            / (2*(dx*dx + dy*dy))
                        - dx*dx*dy*dy*rho
                            / (2*dt*(dx*dx + dy*dy))
                        * ((U_new[(j)*(N_local + 2) + (i+1)] - U_new[(j)*(N_local + 2) + (i-1)]) 
                                / (2*dx) 
                            + (V_new[(j+1)*(N_local + 2) + (i)] - V_new[(j-1)*(N_local + 2) + (i)]) 
                                / (2*dy)
                    );
                }
            }

            // 2.3. Ставим граничные условия
            Poisson_boundary_condition(P_new, N_local, M_local, global_i_start, global_j_start, 
                                        coords, west, east, north, south);

            double maximum_local_P = max_abs_diff(P_new, P_old, N_local, M_local);

            MPI_Allreduce(&maximum_local_P, &maximum_P, 1, MPI_DOUBLE, MPI_MAX, cart_comm);

            P_old = P_new;

            iteration_P++;
        } while (iteration_P < stop_iteration_P and maximum_P > eps_P);

        // 3. Исправление скоростей
        // 3.1. Обмен сообщениями для Пуассона
        exchange_boundary_points(P_new, N_local, M_local, west, east, north, south, cart_comm);

        // 3.2. Явная схема
        for(int j = 1; j < M_local + 1; ++j){
            for(int i = 1; i < N_local + 1; ++i){
                U_new[(j)*(N_local + 2) + (i)] = U_new[(j)*(N_local + 2) + (i)] - dt / rho * (P_new[(j)*(N_local + 2) + (i+1)] - P_new[(j)*(N_local + 2) + (i-1)]) / (2*dx);
                V_new[(j)*(N_local + 2) + (i)] = V_new[(j)*(N_local + 2) + (i)] - dt / rho * (P_new[(j+1)*(N_local + 2) + (i)] - P_new[(j-1)*(N_local + 2) + (i)]) / (2*dy);
            }
        }

        // 3.3. Ставим граничные условия
        Burger_boundary_condition(U_new, V_new, N_local, M_local, global_i_start, global_j_start, coords, west, east, north, south);

        // 4. Нахождение максимума
        double maximum_local_U = max_abs_diff(U_new, U_old, N_local, M_local);
        double maximum_local_V = max_abs_diff(V_new, V_old, N_local, M_local);
        double maximum_local = fmax(maximum_local_U, maximum_local_V);

        MPI_Allreduce(&maximum_local, &maximum, 1, MPI_DOUBLE, MPI_MAX, cart_comm);

        // 5. Перезаписываем массивы
        U_old = U_new;
        V_old = V_new;

        if (rank == master and iteration % 10 == 0) {
            std::cout << iteration << " " << maximum << std::endl;
        }

        // 6. Увеличиваем счетчик итерации
        iteration++;
    } while (iteration < stop_iteration and maximum > eps);

    if (rank == master) {
        std::cout << iteration << " " << maximum << std::endl;
        end_time = MPI_Wtime();
        printf("Calculation time: %f seconds\n", end_time - start_time);
        fflush(stdout);
    }

    Vector<double> x(N, 0);
    Vector<double> y(M, 0);

    for (int i = 0; i < N; ++i) {
        x[i] = i*dx;
    }

    for (int j = 0; j < M; ++j) {
        y[j] = j*dy;
    }

    std::string s = "Source\\";

    Vector<double> Result(M_local * N_local, 0);
    for (int j = 1; j < M_local + 1; ++j) {
        for (int i = 1; i < N_local + 1; ++i) {
            Result[(j-1)*(N_local) + (i-1)] = U_new[(j)*(N_local + 2) + (i)];
        }
    }
    save_2d_vector(s + std::to_string(rank) + "U.txt", Result, M_local, N_local, 6, "\t");

    for (int j = 1; j < M_local + 1; ++j) {
        for (int i = 1; i < N_local + 1; ++i) {
            Result[(j-1)*(N_local) + (i)] = V_new[(j)*(N_local + 2) + (i)];
        }
    }
    save_2d_vector(s + std::to_string(rank) + "V.txt", Result, M_local, N_local, 6, "\t");

    save_1d_vector("x.txt", x, 6, "\n");
    save_1d_vector("y.txt", y, 6, "\n");

    MPI_Finalize();
    return 0;
}