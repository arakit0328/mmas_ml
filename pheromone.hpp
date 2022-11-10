#pragma once

#include "SCP.hpp"

class Pheromone
{
public:
  Pheromone(int column, double alpha, double beta, double rho, double tau_max, double tau_min);
  ~Pheromone();
  void init_pheromone_trails(double initial_trail);
  void update_pheromone_trail(SCPsolution& F); // フェロモンを解Fを使って更新
  void compute_total_information(SCPsolution& F);
  double HEURISTIC(int j, SCPsolution& F);
  double get_total(int j);      // 列jのtotal information
  void print_pheromone();
  void print_total();

private:
  int col;                      // num of columns
  double alpha, beta;
  double rho;                   // evaporation ratio
  double tau_max, tau_min;      // for max-min aco
  double *pheromone;            // Pheromone array
  double *total;                // Total infomation
};


