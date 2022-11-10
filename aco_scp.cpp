#include "aco_scp.hpp"
#include "set_covering_problem.hpp"
#include "pheromone.hpp"
#include <cstdio>
#include <fstream>
#include <cmath>
#include <map>
#include <algorithm>
#include <set>
#include <chrono>
#include <deque>
#include <cfloat>

const bool DEBUG = false;

// コンストラクタ
ACO_SCP::ACO_SCP(int irow, int icolumn, int max_col, int max_row,
                 const vector<int>& vcost, const vector< vector<int> >& vrow,
                 const vector< vector<int> >& vcolumn,
                 int alpha, int beta, double rho,
                 RandomGenerator& rg,
                 int nants)
  : row(irow), column(icolumn), \
    max_num_of_rows_covered_one_column(max_row), max_num_of_column_covers_one_row(max_col), \
    cost(vcost), row_covered(vrow), column_covers(vcolumn), \
    alpha(alpha), beta(beta), rho(rho), \
    random_generator01(rg), \
    nants(nants)
{
  known_optimal_value = 0;
  time_limit = DEFAULT_TIME_LIM;
}

ACO_SCP::ACO_SCP(Set_covering_problem& scp,
                 int alpha, int beta, double rho,
                 RandomGenerator& rg,
                 int nants,
                 int time_limit, int maxiter)
  : row(scp.get_row()), column(scp.get_column()),           \
    max_num_of_rows_covered_one_column(scp.get_max_rows()), \
    max_num_of_column_covers_one_row(scp.get_max_columns()), \
    cost(scp.get_costs()), \
    row_covered(scp.get_row_covered()), column_covers(scp.get_column_covers()), \
    alpha(alpha), beta(beta), rho(rho), \
    random_generator01(rg), \
    nants(nants),
    time_limit(time_limit), maxiter(maxiter)
{
  known_optimal_value = 0;
  time_limit = DEFAULT_TIME_LIM;
}

// main
void ACO_SCP::main()
{
  // ACOの実行
  try { do_acoscp();}
  catch (const char* error_message) {}
}

// 初期解を与える
void ACO_SCP::set_initial_solution(const vector<int> &sol)
{
  Family_of_subsets init(row, column);
  for (const int j : sol)
    init.add_column(j);
  initial_solution = init;
}

// ACO
void ACO_SCP::do_acoscp()
{
  if (DEBUG)
  {
    // Problem
    for (int j=0; j<column; j++) printf("%d,", cost[j]);
    printf("\n");

    for (int j=0; j<column; j++) printf("%lu,", column_covers[j].size());
    printf("\n");

    for (int i=0; i<row; i++) printf("%lu,", row_covered[i].size());
    printf("\n");
  }
  // Greedy
  // if (initial_solution.get_size() == 0) {
  //   fprintf(stderr, "start greedy method\n");
  //   initial_solution = greedy_method();

  //   initial_solution.display();
  //   throw "initial_solution is obtained";
  // }

  initial_solution = greedy_method();
  if (DEBUG) initial_solution.display();

  incumbent_solution = initial_solution;
  Q = initial_solution.get_cost();

  Pheromone pheromone(column, Q, alpha, beta, rho);
  pheromone.init_pheromone_trails( 1.0 );
  if (DEBUG)
  {
    pheromone.print_pheromone();
    pheromone.print_total();
  }

  clock_t iteration_start = clock();       // スタート時間
  clock_t iteration_end = iteration_start; // 終了時間


  /*** 終了条件 ***/
  int iteration = 0;
  // 時間が LIMIT_TIME秒 を超えたら終了する
  while ((double)(iteration_end - iteration_start) / CLOCKS_PER_SEC < time_limit)
  {
    // Move ants
    if (DEBUG) printf("Move Ants\n");
    for (int ant = 0; ant < nants; ant++)
    {
      //   const int INF = 1e6;
      Family_of_subsets solution(row, column);

      // Repeat the following until a feasible solution is obtained
      while (!solution.is_feasible())
      {
        // Select an uncovered row choosed_row;
        int ci = (int)(solution.uncovered_rows.size() * random_generator01.get());
        int choosed_row = solution.uncovered_rows[ci];
        if (DEBUG) printf("Choose row %d\n", choosed_row);

        // Pickup columns that cover choosed_row
        vector<int> cov;
        for (int j : row_covered[choosed_row])
        {
          if (!solution.contains(j)) cov.push_back(j);
        }

        // Choose randomly a column cj
        // cov: choosed_rowを含む列の一覧
        if (DEBUG)
        {
          for (int j : cov)
          {
            if (solution.contains(j)) { printf("Error\n"); exit(1); }
            printf("%d ", j);
          }
          printf("\n");
        }
        //uniform_int_distribution<int> urndcov(0, (int)cov.size()-1);
        //int idj = urndcov(mt);
        int idx = (int)(cov.size() * random_generator01.get());
        int cj = cov[idx];

        pheromone.compute_total_information(solution);
        if (DEBUG) pheromone.print_pheromone();
        if (DEBUG)
        {
          for (int j : cov)
          {
            printf("%f ", pheromone.get_total(j));
          }
          printf("\n");
        }


        double total_sum = 0.0;
        for (int j : cov)
        {
          total_sum += pheromone.get_total(j);
        }
        // 列を選択
        double rnd = total_sum * random_generator01.get();
        double partial_sum = 0.0;
        for (int cj : cov)
        {
          partial_sum += pheromone.get_total(cj);
          if (partial_sum > rnd) break;
        }
        // column cj is chosen

        // Add the column cj to the solution
        if (DEBUG) printf("Add column %d\n", cj);
        solution.add_column(cj);

      } // Ene while !solution.feasible()

      if (DEBUG) printf("End Ant %d\n", ant);

      // Update the iteration best solution
      if (ant == 0)
      {
        iteration_best_solution = solution;
      }
      else if (iteration_best_solution.get_cost() > solution.get_cost())
      {
        if (DEBUG) printf("Iteration_Best_Solution is updated.\n");
        iteration_best_solution = solution;
      }
    } // End Move all ants

    /// Update Best-So-Far solution
    if (incumbent_solution.get_cost() > iteration_best_solution.get_cost())
    {
      if (DEBUG) printf("Incumbent_Solution is updated.\n");
      incumbent_solution = iteration_best_solution;
    }

    // Pheromone update
    pheromone.evaporation();
    if (DEBUG) pheromone.print_pheromone();

    // incumbent_solution のcolumn一覧
    pheromone.update_pheromone_trail(incumbent_solution);
    if (DEBUG) { printf("Add pheromone "); pheromone.print_pheromone(); }

    if (DEBUG) printf("End %d-iteration\n", iteration);
    iteration++;

    // 終了時間
    iteration_end = clock();

  } // End while Iteration


  printf("Ants Finished\n");
  incumbent_solution.display();

  printf("End\n");
}

// 貪欲法
Family_of_subsets ACO_SCP::greedy_method()
{
  Family_of_subsets solution(row, column);

  // それぞれの列の[コスト / |新規にカバーできる行数|]の値
  double *cost_per_rows = new double[column];
  for (int j = 0; j < column; j++)
    cost_per_rows[j] = cost[j] / static_cast<double>(column_covers[j].size());

  while (!solution.is_feasible()) {
    // 現在最も一行辺りのコストが低い列を選び解に追加する
    int min_j = -1;
    double min_gamma_j = INF;
    for (int j = 0; j < column; j++) {
      if (cost_per_rows[j] < min_gamma_j) {
        min_j = j;
        min_gamma_j = cost_per_rows[j];
      }
    }
    solution.add_column(min_j);

    // 変動する列だけc/rの値を更新
    for (const int i : column_covers[min_j]) {
      for (const int j : row_covered[i]) {
        if (solution.contains(j)) {
          cost_per_rows[j] = INF;
        }
        else {
          int row_not_covered = 0;
          for (const int c : column_covers[j])
            if (solution.get_theta(c,0))
              row_not_covered++;
          cost_per_rows[j] = cost[j] / static_cast<double>(row_not_covered);
        }
      }
    }
  }

  delete[] cost_per_rows;

  return solution;
}
