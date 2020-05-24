#include "conjugate_gradient_parallel.h"
#include <stdlib.h>
#include "matrix.h"

process_data set_up_world(int Np, int N)
{
    process_data row;
    int period, size, large_count, col_cnt_dsp[3], buf[3 * Np];
    MPI_Aint lb, extent;
    MPI_Datatype row_t;

    // Store number of processes Np and dimension N
    row.N = N;
    row.Np = Np;

    // Create 1D communicator and save ranks and coordinates
    period = 1;
    MPI_Cart_create(MPI_COMM_WORLD, 1, &Np, &period, 0, &row.comm);
    MPI_Comm_rank(row.comm, &row.rank);
    MPI_Cart_coords(row.comm, row.rank, 1, &row.coord);

    // Calculate the number of rows handled by each process
    large_count = N % Np;
    row.count_min = N / Np;
    row.count_max = (large_count == 0) ? (row.count_min) : (row.count_min + 1);
    row.count = (row.coord < large_count) ? (row.count_max) : (row.count_min);
    row.displ = row.coord * (row.count_min) + ((row.coord <= large_count) ? (row.coord) : (large_count));

    // Create types for a block within a row, a transposed block, and a full row
    MPI_Type_vector(row.count, 1, row.N, MPI_DOUBLE, &row_t);
    MPI_Type_create_resized(row_t, 0, sizeof(double), &row.block_trans_t);

    MPI_Type_vector(1, row.N, 1, MPI_DOUBLE, &row.row_t);
    MPI_Type_commit(&row.block_trans_t);
    MPI_Type_commit(&row.row_t);

    // Gather rank, count, and displacement of each coordinate
    col_cnt_dsp[0] = row.coord;
    col_cnt_dsp[1] = row.count;
    col_cnt_dsp[2] = row.displ;

    MPI_Allgather(col_cnt_dsp, 3, MPI_INT, buf, 3, MPI_INT, row.comm);

    row.ranks = malloc(Np * sizeof(int));
    row.counts = malloc(Np * sizeof(int));
    row.displs = malloc(Np * sizeof(int));

    for (int i = 0; i < Np; ++i)
    {
        row.ranks[buf[3 * i]] = i;
        row.counts[i] = buf[3 * i + 1];
        row.displs[i] = buf[3 * i + 2];
    }

    return row;
}

// Exits program with error message if ptr is NULL
void malloc_test(void *ptr)
{
    if (ptr == NULL)
    {
        printf("malloc failed\n");
        MPI_Abort(MPI_COMM_WORLD, 0);
        exit(0);
    }
}

// Generate a random positive definite equation with solution by
// 1. Generating a random matrix A
// 2. Calculating A = A' + A
// 3. Adding N to each element on the diagonal of A
// 4. Generating a random solution x
// 5. Calculating the rhs in A * x = b
equation_data random_linear_system(process_data row)
{
    equation_data equation;
    double *B, *a, *b, *x_recv, *x_send, *x_tmp;
    int B_size, coord_work, coord_send, coord_recv, rank_work, rank_send,
            rank_recv, rank_up, rank_down, rank_block, displ_work, displ_send,
            count_work, count_send;
    MPI_Request send_req, send_reqs[row.Np];
    MPI_Status recv_stat;

    equation.N = row.N;
    equation.A = random_matrix(row.count, row.N);
    malloc_test(equation.A);

    if (row.Np > 1)
    {
        B_size = row.count_max * row.count_max;
        B = malloc(B_size * sizeof(double));
        malloc_test(B);
    }

    // Calculate A = 0.5*(A' + A) and add N to the diagonal
    for (int n = 0; n < row.Np; ++n)
    {
        coord_work =(row.coord + n) % row.Np;
        coord_recv =(row.Np + row.coord - 1 - n) % row.Np;
        coord_send =(coord_work + 1) % row.Np;
        rank_work = row.ranks[coord_work];
        rank_recv = row.ranks[coord_recv];
        rank_send = row.ranks[coord_send];
        displ_work = row.displs[rank_work];
        count_work = row.counts[rank_work];

        if ((n < row.Np - 1) &&(coord_work != coord_recv))
            MPI_Isend(&equation.A[row.displs[rank_recv]], row.counts[rank_recv],
                    row.block_trans_t, rank_recv, n, row.comm, &send_req);

        if (n == 0)
        {
            // Don't use the buffer, just calculate the diagonal block addition in-place
            for (int i = 0; i < row.count; ++i)
            {
                a = equation.A +(i * row.N + displ_work);
                for (int j = 0; j < count_work; ++j)
                {
                    b = equation.A +(j * row.N + displ_work);
                    if (j < i)
                        a[j] = b[i];
                    else
                        a[j] = 0.5 *(a[j] + b[i]);
                    if (j == i)
                        a[j] += row.N;
                }
            }
        }
        else if (n > row.Np / 2)
        {
            // Just copy B
            for (int i = 0; i < row.count; ++i)
            {
                a = equation.A +(i * row.N + displ_work);
                b = B +(i * count_work);
                for (int j = 0; j < count_work; ++j)
                    a[j] = b[j];
            }
        }
        else
        {
            // Add B to A
            for (int i = 0; i < row.count; ++i)
            {
                a = equation.A +(i * row.N + displ_work);
                b = B +(i * count_work);
                for (int j = 0; j < count_work; ++j)
                    a[j] = 0.5 *(a[j] + b[j]);
            }
        }

        if (n < row.Np - 1)
        {
            if (coord_work != coord_recv)
            {
                MPI_Recv(&(B[0]), B_size, MPI_DOUBLE, rank_send, n, row.comm,
                        MPI_STATUS_IGNORE);
                MPI_Wait(&send_req, MPI_STATUS_IGNORE);
            }
            else
            {
                MPI_Sendrecv(&equation.A[row.displs[rank_recv]],
                        row.counts[rank_recv], row.block_trans_t,
                        rank_recv, n, &(B[0]), B_size, MPI_DOUBLE,
                        rank_send, n, row.comm, MPI_STATUS_IGNORE);
            }
        }
    }

    if (row.Np > 1)
        free(B);			// Useless memory, free before allocating more

    // Generate random solution x, zero matrix b, and memory for x*
    equation.x = random_matrix(row.count, 1);
    malloc_test(equation.x);
    equation.x_star = calloc(row.count, sizeof(double));
    equation.b = calloc(row.count, sizeof(double));
    malloc_test(equation.x_star);
    malloc_test(equation.b);

    // Create tempory x vectors for sending and receiving
    x_recv = malloc(row.count_max * sizeof(double));
    x_send = malloc(row.count_max * sizeof(double));
    malloc_test(x_recv);
    malloc_test(x_send);

    rank_up = row.ranks[(row.Np + row.coord - 1) % row.Np];
    rank_down = row.ranks[(row.coord + 1) % row.Np];

    // Initially store local x in x_recv
    memcpy(x_recv, equation.x, row.count * sizeof(double));

    // Perform matrix-vector multiplication to calculate rhs b
    for (int n = 0; n < row.Np; ++n)
    {
        x_tmp = x_recv;
        x_recv = x_send;
        x_send = x_tmp;

        if (n < row.Np - 1)
            MPI_Isend(x_send, row.count_max, MPI_DOUBLE, rank_up, 111, row.comm,
                      &send_req);

        rank_block = row.ranks[(row.coord + n) % row.Np];
        displ_send = row.displs[rank_block];
        count_send = row.counts[rank_block];
        for (int i = 0; i < row.count; ++i)
        {
            a = equation.A +(i * row.N + displ_send);
            for (int j = 0; j < count_send; ++j)
                equation.b[i] += x_send[j] * a[j];
        }

        if (n < row.Np - 1)
        {
            MPI_Recv(x_recv, row.count_max, MPI_DOUBLE, rank_down, 111,
                     row.comm, MPI_STATUS_IGNORE);
            MPI_Wait(&send_req, MPI_STATUS_IGNORE);
        }
    }

    free(x_recv);
    free(x_send);

    return equation;
}

void solve_conjugate_gradient_par(process_data row, equation_data equation, int max_steps, double tol)
{
    double *a, *p, *p_recv, *p_send, *p_tmp, *r, *z, *x;
    double alpha, alpha_tmp, beta, gamma, gamma_new, gamma_tmp;
    int rank_up, rank_down, rank_block, displ_send, count_send;
    MPI_Request send_req;

    x = &(equation.x_star[0]);

    p = malloc(row.count * sizeof(double));
    r = malloc(row.count * sizeof(double));
    z = malloc(row.count * sizeof(double));

    malloc_test(p);
    malloc_test(r);
    malloc_test(z);

    p_recv = malloc(row.count_max * sizeof(double));
    p_send = malloc(row.count_max * sizeof(double));
    malloc_test(p_recv);
    malloc_test(p_send);

    rank_up = row.ranks[(row.Np + row.coord - 1) % row.Np];
    rank_down = row.ranks[(row.coord + 1) % row.Np];

    // x = [0 ... 0] - initial guess
    // r = b
    // p = r
    // gamma = r' * r
    gamma_tmp = 0.0;
    for (int i = 0; i < row.count; ++i)
    {
        x[i] = 0.0;
        r[i] = equation.b[i];
        p[i] = r[i];
        p_recv[i] = p[i];
        z[i] = 0.0;
        gamma_tmp += r[i] * r[i];
    }
    MPI_Allreduce(&gamma_tmp, &gamma, 1, MPI_DOUBLE, MPI_SUM, row.comm);

    for (int n = 0; n < max_steps; ++n)
    {
        // z = A * p
        for (int m = 0; m < row.Np; ++m)
        {
            p_tmp = p_recv;
            p_recv = p_send;
            p_send = p_tmp;

            if (m < row.Np - 1)
                MPI_Isend(p_send, row.count_max, MPI_DOUBLE, rank_up, 222,
                          row.comm, &send_req);

            rank_block = row.ranks[(row.coord + m) % row.Np];
            displ_send = row.displs[rank_block];
            count_send = row.counts[rank_block];

            for (int i = 0; i < row.count; ++i)
            {
                a = equation.A +(i * row.N + displ_send);
                for (int j = 0; j < count_send; ++j)
                    z[i] += p_send[j] * a[j];
            }

            if (m < row.Np - 1)
            {
                MPI_Recv(p_recv, row.count_max, MPI_DOUBLE, rank_down, 222,
                         row.comm, MPI_STATUS_IGNORE);
                MPI_Wait(&send_req, MPI_STATUS_IGNORE);
            }
        }

        // alpha = gamma /(p' * z)
        alpha_tmp = 0.0;
        for (int i = 0; i < row.count; ++i)
            alpha_tmp += p[i] * z[i];

        MPI_Allreduce(&alpha_tmp, &alpha, 1, MPI_DOUBLE, MPI_SUM, row.comm);

        alpha = gamma / alpha;

        // x = x + alpha * p
        // r = r - alpha * z
        // gamma_new = r' * r
        gamma_tmp = 0.0;
        for (int i = 0; i < row.count; ++i)
        {
            x[i] += alpha * p[i];
            r[i] -= alpha * z[i];
            gamma_tmp += r[i] * r[i];
        }

        MPI_Allreduce(&gamma_tmp, &gamma_new, 1, MPI_DOUBLE, MPI_SUM,
                      row.comm);

        if (sqrt(gamma_new) < tol)
            break;

        beta = gamma_new / gamma;

        // p = r + beta * p;
        for (int i = 0; i < row.count; ++i)
        {
            p[i] = r[i] + beta * p[i];
            p_recv[i] = p[i];
            z[i] = 0.0;
        }

        // gamma = gamma_new
        gamma = gamma_new;
    }

    free(p);
    free(r);
    free(z);
    free(p_send);
    free(p_recv);

    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    return;
}
