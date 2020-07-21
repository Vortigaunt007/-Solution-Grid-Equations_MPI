
#include "Solver.h"


Solver::Solver(MatrixCSR &M_in, Vector &b_in): A(M_in), b(b_in), JacobiPre()
{

}

void Solver::InitJacobiPreconditioner()
{
    // Jacobi Preconditioner is diagonal matrix: JacobiPre[diag_elem] = A[diag_elem]
    int myranksize;
    int world_rank;

    MPI_Comm_size(MPI_COMM_WORLD, &myranksize);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    const int n = (int)A.AI.size() - 1;
    for(int i = 0; i < n; i++) { // string number
        for(int j_ = A.AI[i]; j_ < A.AI[i+1]; j_++) {
            int j = A.L2G[A.AJ[j_]]; // column number
            if (A.L2G[i] == j) {// diagonal element
                JacobiPre.AI.push_back(JacobiPre.AJ.size());
                JacobiPre.AJ.push_back(A.L2G[j]); // number of column of the matrix A diagonal
                JacobiPre.A.push_back(A.A[j_]); // diagonal element of A
            }
        }
    }

    JacobiPre.AI.push_back(JacobiPre.AJ.size());

    JacobiPre.L2G = A.L2G;
    JacobiPre.G2L = A.G2L;
    JacobiPre.Part = A.Part;
    JacobiPre.size_Halo = A.size_Halo;
    JacobiPre.size_Own = A.size_Own;
}


double Solver::ConjugateGradient(MPI_Exchange K)
{    
    int myranksize;
    int world_rank;

    MPI_Comm_size(MPI_COMM_WORLD, &myranksize);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    BasicOperation Basic_Operation;
    int size_with_Halo = A.size_Own + A.size_Halo;
    int N = size_with_Halo;
    Vector x_k(N), z_k(N), r_k(N), p_k(N), q_k(N), Ax_k(N);
    double alpha, beta, rho_k_1 = 1.0, rho_k = 1.0;
    int count_step; // k
    bool convergence;
    double time_startSolve;

    // initialization step


    A.Fill_Matrix(InitOffDiagElement); // Fill A in MatrixCSR A
    b.Fill_Vector(A); // Fill right part

    InitJacobiPreconditioner(); // Fill A in MatrixCSR Jacobi Preconditioner
    MatrixCSR M;
    M.L2G = JacobiPre.L2G;
    M.G2L = JacobiPre.G2L;
    M.Part = JacobiPre.Part;
    M.size_Halo = JacobiPre.size_Halo;
    M.size_Own = JacobiPre.size_Own;
    inverse_Matrix(JacobiPre, M); // M = JacobiPre^(-1)
    count_step = 1;
    convergence = false;

    long long count_mult = 0, count_summ = 0, count_dot = 0;
    double start_time, time_mult_Matrix_Vector, time_summ, time_dot;

    start_time = MPI_Wtime();
    for (int i = 0; i < 100; i++)
        Basic_Operation.mult_Matrix_Vector(A, x_k, Ax_k, K);
    time_mult_Matrix_Vector = MPI_Wtime() - start_time;
    time_mult_Matrix_Vector /= 100.0;
    count_mult = Basic_Operation.countMult;
    count_mult /= 100;

    start_time = MPI_Wtime();
    for (int i = 0; i < 100; i++)
        Basic_Operation.summAlpha(1.0, b, -1.0, Ax_k, r_k, size_with_Halo);
    time_summ = MPI_Wtime() - start_time;
    count_summ = Basic_Operation.countSummAlpha;

    start_time = MPI_Wtime();
    for (int i = 0; i < 100; i++)
        rho_k = Basic_Operation.dotProduct(r_k, z_k, A.size_Own);
    time_dot = MPI_Wtime() - start_time;
    count_dot = Basic_Operation.countDot;

    Basic_Operation.resetCounters();  // resetting counters: countDot, timeDot ...

    // first step of the conjugate gradient method

    Basic_Operation.summAlpha(1.0, b, -1.0, Ax_k, r_k, size_with_Halo); // r_0 = b - Ax_0 = b
    Basic_Operation.resetCounters();

    time_startSolve = MPI_Wtime();

    int k = 1;
    double alpha_prev = 0.0;

    while (!convergence) {
        Basic_Operation.update(A, r_k, K);
        Basic_Operation.mult_Matrix_Vector(M, r_k, z_k, K);

        alpha = Basic_Operation.dotProduct(r_k, z_k, A.size_Own);

        if(k == 1) {
            for (int i = 0; i < M.size_Own; i++)
                p_k[i] = z_k[i];
        } else {
            beta = alpha/alpha_prev;
            Basic_Operation.countSolver += 1;
            Basic_Operation.summAlpha(1.0, z_k, beta, p_k, p_k, size_with_Halo);
        }

        Basic_Operation.update(A, p_k, K);
        Basic_Operation.mult_Matrix_Vector(A, p_k, q_k, K);

        double dottmp, gamma;
        dottmp = Basic_Operation.dotProduct(p_k, q_k, A.size_Own);

        Basic_Operation.countSolver += 1;
        gamma = alpha / dottmp;

        Basic_Operation.summAlpha(1.0, x_k, gamma, p_k, x_k, size_with_Halo);
        Basic_Operation.summAlpha(1.0, r_k, -gamma, q_k, r_k, size_with_Halo);

        alpha_prev = alpha;

        //if (world_rank == 1)
            std::cout << "Iteration number = " << k << " " << "r_k = " << alpha << std::endl;
        if ((alpha < eps) || (k >= maxiter))
            convergence = true;
        else
            k++;
    }

/*
    while (!convergence) {
        Basic_Operation.mult_Matrix_Vector(A, p_k, q_k, K);
        Basic_Operation.update(A, q_k, K);
        alpha = rho_k / Basic_Operation.dotProduct(p_k, q_k, A.size_Own);
        Basic_Operation.countSolver += 1;

        Basic_Operation.summAlpha(1.0, x_k, alpha, p_k, x_k, size_with_Halo);
        Basic_Operation.summAlpha(1.0, r_k, -alpha, q_k, r_k, size_with_Halo);
       // if (world_rank == 0)
          std::cout << "Iteration number = " << count_step << " " << "r_k = " << sqrt(Basic_Operation.dotProduct(r_k, r_k, size_with_Halo)) << std::endl;
        if ((rho_k < eps) || (count_step >= maxiter))
            convergence = true;
        else
            count_step++;
        rho_k_1 = rho_k;
        Basic_Operation.mult_Matrix_Vector(M, r_k, z_k, K);
        Basic_Operation.update(A, z_k, K);
        rho_k = Basic_Operation.dotProduct(r_k, z_k, A.size_Own);
        beta = rho_k / rho_k_1;
        Basic_Operation.countSolver += 1;
        Basic_Operation.summAlpha(1.0, z_k, beta, p_k, p_k, size_with_Halo);
    }
*/
    timeSolveCG = MPI_Wtime() - time_startSolve;


    double maxT;
    MPI_Reduce(&timeSolveCG, &maxT, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    double gSolver = Basic_Operation.countSolver/(timeSolveCG * 1E9);
    double gAxpy = count_summ/(Basic_Operation.timeSummAlpha* 1E9);
    double gDot = count_dot/(Basic_Operation.timeDot * 1E9);
    double gSP_MV = count_mult/(Basic_Operation.timeMult * 1E9);

    MPI_Allreduce(&gSolver, &gSolver, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&gAxpy, &gAxpy, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&gDot, &gDot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&gSP_MV, &gSP_MV, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    Vector res(N);
    res = x_k;

    double controlValue;
    controlValue = Basic_Operation.dotProduct(x_k, x_k, A.size_Own);


    if(world_rank == 0) {
        std::cout << "Time working solver: " << maxT << std::endl;
        std::cout << "Iter: " << k << std::endl;
        std::cout << "Tol " << alpha << std::endl;
        std::cout << "Solver GFLOPS: " << gSolver << std::endl;
        std::cout << "SpMV_ELLPACK GFLOPS: " << gSP_MV << std::endl;
        std::cout << "AXPBY GFLOPS: " << gAxpy << std::endl;
        std::cout << "DOT GFLOPS: " << gDot << std::endl;
        std::cout << "L2 for result vector " << sqrt(controlValue) << std::endl;
    }

    return sqrt(controlValue);
}

/*
void Solver::saveToFileData(std::string filename) const
{
    std::ofstream fout(filename);

    fout << "A.A.size() = " << A.A.size() << std::endl;
    fout << "A.AI.size() = " << A.AI.size() << std::endl;
    fout << "A.AJ.size() = " << A.AJ.size() << std::endl;

    fout << "time_spent_filling = " << time_spent_filling << std::endl;
    fout << "time_spent_solving = " << time_spent_solving << std::endl;
    fout << "time_dotProduct = " << time_dotProduct << std::endl;
    fout << "time_summAlpha = " << time_summAlpha << std::endl;
    fout << "time_mult_Matrix_Vector = " << time_mult_Matrix_Vector << std::endl;

    fout.close();
}
*/
