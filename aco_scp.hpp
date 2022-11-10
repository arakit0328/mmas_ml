#pragma once

#include "SCP.hpp"
#include <string>
#include <vector>
#include <random>

// 0-1乱数生成クラス
class RandomGenerator
{
public:
  // コンストラクタ
  RandomGenerator(int seed = 1)
    : generator_(seed), distribution_(0.0, 1.0) {}
  ~RandomGenerator(){}

  // 0-1乱数を一つ返す
  double get() { return distribution_(generator_); }
private:
  std::mt19937_64 generator_;
  std::uniform_real_distribution<double> distribution_;
};


class ACO_SCP
{
public:
  ACO_SCP(int irow, int icolumn, int max_col, int max_row,
          const vector<int>& vcost, const vector< vector<int> >& vrow,
          const vector< vector<int> >& vcolumn,
          int alpha, int beta, double rho,
          RandomGenerator& rg,
          int nants);

  ACO_SCP(Set_covering_problem& scp,
          int alpha, int beta, double rho,
          RandomGenerator& rg,
          int nants,
          int time_limit, int maxiter);

  void main();

  // 初期解を与える
  void set_initial_solution(const vector<int> &sol);

  void do_acoscp();

  // Greedy Method
  Family_of_subsets greedy_method();

private:
  const int row, column;   // row行column列の行列
  const vector<int>& cost; // 各列のコスト |column|個
  const vector< vector<int> >& row_covered; // 各行が何列目によってカバーされているか |row|個
  const vector< vector<int> >& column_covers; // 各列が何行目をカバーするのか |column|個

  /// Parameters for ACO
  int nants;                    // num of ants
  double Q = 10;

  /// End: parameters for ACO
  vector<double> pheromone;     // pheromone
  vector<double> heuristic;     // pheromone
  int alpha;
  int beta;
  double rho;

  time_t time_limit;  // 探索の制限時間 単位は秒
  int maxiter;
  int known_optimal_value; // 既知の最適解のコスト 小規模なインスタンスの実験用

  int upper_bound; // 探索中に発見した実行可能解の目的関数値の内最小のもの
  int max_num_of_rows_covered_one_column; // ある一列によってカバーされている最大の行数 t
  int max_num_of_column_covers_one_row; // ある一行をカバーしている最大の列数 l

  Family_of_subsets initial_solution;    // 初期解
  Family_of_subsets incumbent_solution;  // 現段階での最適解
  Family_of_subsets iteration_best_solution; // 1回の探索での最適解
  Family_of_subsets *searching_solution; // 現在探索中の解

  RandomGenerator random_generator01;       // 0-1乱数生成

  const bool debug = false;
  const int INF = 1e6;

  const int DEFAULT_TIME_LIM = 30;  // default is different for each instance, min180
  const int MINUMUM_TRIAL_LOCAL_SEARCH = 100;  // default 100

};                              // end ACO_SCP


#endif
