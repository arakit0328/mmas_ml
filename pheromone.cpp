#include <cstdio>
#include <cmath>

#include "pheromone.hpp"

Pheromone::Pheromone(int column,
                     double alpha,
                     double beta,
                     double rho,
                     double tau_max,
                     double tau_min)
  : alpha(alpha), beta(beta), rho(rho)
{
  pheromone = new double [column];
  total = new double [column];
  col = column;
  tau_max = tau_max;
  tau_min = tau_min;
}

Pheromone::~Pheromone()
{
  delete [] pheromone;
  delete [] total;
}


double Pheromone::HEURISTIC(int j, SCPsolution &cs)
{
  return cs.SCORE[j];
}


void Pheromone::init_pheromone_trails(double initial_trail = 1.0)
{
  /* Initialize pheromone trails */
  for (int j = 0; j < col; j++)
  {
    pheromone[j] = initial_trail;
    total[j] = initial_trail;
  }
}


// 解Fに含まれる列のフェロモンを増加
void Pheromone::update_pheromone_trail(SCPsolution& global_best_cs)
{
  for (int j = 1; j <= col; ++j)
  {
    pheromone[j] = rho * pheromone[j] ;
  }

  double t = 0;

  for (int j = 1; j <= global_best_cs.K; ++j)
  {
    t = pheromone[j] + (1 - rho);

    if (t < tau_min)
      pheromone[j] = tau_min;
    else if (t > tau_max)
      pheromone[j] = tau_max;
    else
      pheromone[j] = t;
  }
}


// Computes total infomation and stores in total[1:col]
void Pheromone::compute_total_information(SCPsolution& cs)
{
  for (int j = 0; j < col; j++)
  {
    total[j] = pow(pheromone[j], alpha) * pow(HEURISTIC(j, cs), beta);
  }
}


// 列jのtotal information
double Pheromone::get_total(int j)
{
  return total[j];
}


void Pheromone::print_pheromone()
{
  for (int j=0; j < col; j++)
  {
    printf("%.2f,", pheromone[j]);
  }
  printf("\n");
}


void Pheromone::print_total()
{
  for (int j=0; j < col; j++)
  {
    printf("%.2f,", total[j]);
  }
  printf("\n");
}
