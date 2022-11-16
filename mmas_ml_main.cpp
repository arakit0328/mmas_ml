#include "SCPv.hpp"
#include "Random.hpp"
#include <cstdlib>
#include <iostream>
#include <cstdio>
#include <cassert>
#include <vector>
#include <algorithm>
#include <random>
#include <limits>
using namespace std;

// 乱数発生クラス Random.hpp を参照
// rnd() で整数乱数
// rnd(a, b) とするとa以上b以下の一様乱数を生成できる
// Rand rnd;
// mt19937_64 rnd(seed);

// greedy_neighborhood_search で使う繰り返しの回数
int numIteration = 500;

// graspで使う alpha の値
double alpha = 0.95;

// row_weight の最大値
const int Weight_Threshold = 1000;
const double Oblivious_Ratio = 0.3;

// BMS
const double BMS_ratio = 0.66;



// 貪欲法：スコア最大の列をK列選ぶ
// 引数の cs に結果が入る
void greedy_construction(SCPinstance &pData, SCPsolution &cs, Rand &rnd)
{
  int maxc;

  cs.initialize(pData); 
  for (int k = 0; k < cs.K; k++)
  {
    maxc = cs.get_column_maxscore(pData, rnd);
    cs.add_column(pData, maxc, false);
  } // End for k
}


// GRASP法：スコアが alpha * (最大値 - 最小値) 以上である列からランダムに一つ選ぶ
// 引数の cs に結果が入る
void grasp_construction(SCPinstance &pData, SCPsolution &cs, double alpha, Rand &rnd)
{
  cs.initialize(pData);

  int c;
  for (int k = 0; k < cs.K; k++)
  {
    c = cs.get_column_grasp(pData, alpha, rnd);
    cs.add_column(pData, c, false);
  } // End for k
}


// 配列の順序をランダムに入れ替える
void random_permutation(vector<int>& A, Rand& rnd)
{
  int j;
  int n = A.size();
  for (int i = 0; i < n-1; ++i)
  {
    j = rnd(i, n-1);
    swap(A[i], A[j]);
  }
}


// 単純な改善法
// 引数の cs に結果が入る
void simple_neighborhood_search(SCPinstance &pData,
                                SCPsolution &cs,
                                Rand& rnd)
{
  int K = cs.K;
  vector<int> idx;
  int c1, cov1;
  int c2, cov2;

  idx = cs.CS;
  random_permutation(idx, rnd);
  cov1 = cs.num_Cover;

  for (int i=0; i < K; ++i)
  {
    c1 = idx[i];
    cs.remove_column(pData, c1);

    // 最大スコアの列
    c2 = cs.get_column_maxscore(pData, rnd);
    cs.add_column(pData, c2, false);
    cov2 = cs.num_Cover;

    if (cov1 > cov2)
    {
      cs.remove_column(pData, c2);
      cs.add_column(pData, c1, false);
    }
    else
    {
      cov1 = cs.num_Cover;
    }
  } // End for i
}


// 貪欲初期解＋単純局所探索を niter 回繰り返し
// 引数の cs と best_cs は変更される
// 引数の best_cs に結果が入る
void greedy_neighborhood_search(SCPinstance &pData,
                                int K,
                                int niter,
                                SCPsolution &cs,
                                SCPsolution &best_cs,
                                Rand& rnd)
{
  // 貪欲法
  for (int iter = 0; iter < niter; ++iter)
  {
    // 貪欲法
    greedy_construction(pData, cs, rnd);

    // 局所探索
    simple_neighborhood_search(pData, cs, rnd);

    if (best_cs.num_Cover < cs.num_Cover)
    {
      best_cs = cs;
    }
  } // End for iter
}


// GRASP初期解＋単純局所探索を niter 回繰り返し
// 引数の cs と best_cs は変更される
// 引数の best_cs に結果が入る
void grasp_neighborhood_search(SCPinstance &pData,
                               int K,
                               int alpha,
                               int niter,
                               SCPsolution &cs,
                               SCPsolution &best_cs,
                               Rand& rnd)
{
  for (int iter = 0; iter < niter; ++iter)
  {
    // 初期解を生成
    grasp_construction(pData, cs, alpha, rnd);

    // 局所探索
    simple_neighborhood_search(pData, cs, rnd);

    if (best_cs.num_Cover < cs.num_Cover)
    {
      best_cs = cs;
    }
  } // End for iter
}



int get_column_BMS(SCPinstance& instance,
                   SCPsolution& cs,
                   vector<int>& score,
                   vector<int>& subscore,
                   vector<int>& conf,
                   double bms_ratio,
                   Rand& rnd)
{
  std::vector<int> maxCols;
  int maxScore = numeric_limits<int>::min(), maxc = 0;

  vector<int> bms = cs.CS;

  random_permutation(bms, rnd);

  for (int i = 0; i < bms_ratio * cs.K; i++)
  {
    int c = bms[i];

    if (conf[c] == 0) { continue; }

    // 最大スコアの列をチェック
    if (maxScore < score[c])
    {
      maxScore = score[c];
      maxCols.clear();
      maxCols.push_back(c);
    }
    else if (maxScore == score[c])
      maxCols.push_back(c);
  } // End for c

  if (maxCols.size() == 1) maxc = maxCols[0];
  else
  {
    vector<int> maxSub;
    int max_subscore = numeric_limits<int>::min();

    for (int c : maxCols)
    {
      if (subscore[c] > max_subscore)
      {
        max_subscore = subscore[c];
        maxSub.clear();
        maxSub.push_back(c);
      }
      else if (subscore[c] == max_subscore)
        maxSub.push_back(c);
    } // end for c

    if (maxSub.size() == 1) maxc = maxSub[0];
    else
    {
      int j = rnd(0, maxSub.size() - 1);
      maxc = maxSub[j];
    }
  }

  return maxc;
}


// r_select行をカバーする列からスコア最大のものを選択する
int get_column_maxscore(SCPinstance& instance,
                        SCPsolution& cs,
                        int r_select,
                        vector<int>& score,
                        vector<int>& subscore,
                        vector<int>& conf,
                        Rand& rnd)
{
  std::vector<int> maxCols;
  int maxScore = -instance.numColumns, maxc = 0;

  for (int c : instance.RowCovers[r_select])
  {
    if (cs.SOLUTION[c] == 1) { printf("Column %d should not in CS.\n", c); exit(1); }
    if (conf[c] == 0) { continue; }

    // 最大スコアの列をチェック
    if (maxScore < cs.SCORE[c])
    {
      maxScore = cs.SCORE[c];
      maxCols.clear();
      maxCols.push_back(c);
    }
    else if (maxScore == cs.SCORE[c])
      maxCols.push_back(c);
  } // End for c

  if (maxCols.size() == 1) maxc = maxCols[0];
  else
  {
    vector<int> maxSub;
    int max_subscore = numeric_limits<int>::min();

    for (int c : maxCols)
    {
      if (subscore[c] > max_subscore)
      {
        max_subscore = subscore[c];
        maxSub.clear();
        maxSub.push_back(c);
      }
      else if (subscore[c] == max_subscore)
        maxSub.push_back(c);
    } // end for c

    if (maxSub.size() == 1) maxc = maxSub[0];
    else
    {
      int j = rnd(0, maxSub.size() - 1);
      maxc = maxSub[j];
    }
  }

  return maxc;
}


// Oblivious local search
SCPsolution ObliviousLS(SCPinstance& instance,
                        const SCPsolution& CS,
                        //SCPsolution& CSimproved,
                        int weight_threshold,
                        double oblivious_ratio,
                        double bms_ratio,
                        int niter,
                        Rand& rnd)
{
  // conf, score, subscore, weight はここで定義すればいい
  vector<int> score(instance.numColumns, 0);
  vector<int> subscore(instance.numColumns, 0);
  vector<int> row_weight(instance.numRows, 1);
  vector<int> conf(instance.numColumns, 1);
  vector<int> time(instance.numColumns, 0);

  // cs からスコアの初期値を求める
  for (int j = 0; j < instance.numColumns; j++)
  {
    if (CS.SOLUTION[j])
    {
      for (int r : instance.ColEntries[j])
      {
        if (CS.COVERED[r] == 1) score[j]--;
        if (CS.COVERED[r] == 2) subscore[j]--; // twoToone
      }
    }
    else
    {
      // 列jが解に入っていない
      for (int r : instance.ColEntries[j])
      {
        if (CS.COVERED[r] == 0) score[j]++;
        if (CS.COVERED[r] == 1) subscore[j]++; // oneTotwo
      }
    }
  }

  SCPsolution CSbest(instance, CS.K, weight_threshold, oblivious_ratio);
  SCPsolution CSw(instance, CS.K, weight_threshold, oblivious_ratio);

  CSw = CS;
  CSbest = CSw;

  for (int iter = 0; iter < niter; iter++)
  {
    // 列を削除
    int c = get_column_BMS(instance, CSw, score, subscore, conf, bms_ratio, rnd);
    CSw.remove_column(instance, c);

    time[c] = iter;
    // スコア更新
    score[c] = -score[c];
    subscore[c] = -subscore[c];

    for (int r : instance.ColEntries[c]) // 列cがカバーする行
    {
      // r行がカバーされなくなったら，rを含む行のスコアを増加
      if (CSw.COVERED[r] == 0)
      {
        for (int rc : instance.RowCovers[r]) // r行をカバーする列
        {
          if (rc != c) score[rc] += row_weight[r];
        } // End: for ri
      } // End if covered[r] == 0

      // r行が1回カバーされるようになったら，rを含むCSの要素のスコアを減少
      if (CSw.COVERED[r] == 1)
      {
        for (int rc : instance.RowCovers[r]) // r行をカバーする列
        {
          if (CSw.SOLUTION[rc])
          {
            score[rc] -= row_weight[r];
          }
          else
          {
            subscore[rc]++;     // oneTotwo
          }
        }
      } // End if covered[r] == 1

      // r行が2回カバーされるようになったら，rを含むCSの要素のsubscoreを減少
      if (CSw.COVERED[r] == 2)
      {
        for (int rc : instance.RowCovers[r]) // r行をカバーする列
        {
          if (CSw.SOLUTION[rc])
          {
            subscore[rc]--;        // twoToone
          }
        }
      } // End if covered[r] == 1
    }
    // スコア更新終了
    // confの更新
    conf[c] = 0;
    for (int cn : instance.Neighborhood[c]) conf[cn] = 1;
    // confの更新終了

    // 列を追加
    int x = rnd(1, instance.numRows - CSw.num_Cover); // カバーされていない行数までの乱数
    int r_select = 0;           // 選択する行番号 カバーされていない行のx番目になる
    while (x > 0)
    {
      if (CSw.COVERED[r_select] == 0) x--;
      r_select++;
    }
    r_select--;

    assert(CSw.COVERED[r_select] == 0);

    c = get_column_maxscore(instance, CSw, r_select, score, subscore, conf, rnd);
    CSw.add_column(instance, c, true);

    time[c] = iter;
    // スコアを更新
    score[c] = -score[c];
    subscore[c] = -subscore[c];
    // スコア更新
    for (int r : instance.ColEntries[c]) // 列cがカバーする行
    {
      // r行が初めてカバーされたら，rを含む列のスコアを減少
      if (CSw.COVERED[r] == 1)
      {
        for (int rc : instance.RowCovers[r])
        {
          if (rc != c)
          {
            score[rc] -= row_weight[r];
            subscore[rc]++;     // oneTotwo
          }
        } // End: for ri
      } // End if covered[r] == 1

      // r行が2回カバーされたら，rを含むCSwの要素のスコアを増加
      if (CSw.COVERED[r] == 2)
      {
        for (int rc : instance.RowCovers[r]) // r行をカバーする列
        {
          if (CSw.SOLUTION[rc])
          {
            score[rc] += row_weight[r];
            subscore[rc]--;     // twoToone
          }
        }
      } // End if covered[r] == 2
    }

    // confの更新
    // 候補解に含まれる列のconfはc以外は1になる
    for (int c : CSw.CS) conf[c] = 1;
    conf[c] = 0;
    for (int cn : instance.Neighborhood[c]) conf[cn] = 1;
    // confの更新終了

    // row_weightとscoreの更新
    bool flag_reduce = false;
    for (int i = 0; i < instance.numRows; i++)
    {
      if (CSw.COVERED[i] == 0)
      {
        row_weight[i]++;
        if (row_weight[i] > weight_threshold) flag_reduce = true;
        for (int c : instance.RowCovers[i])
        {
          if (CSw.SOLUTION[c]) score[c]--;
          else score[c]++;
        }
      }
    }

    if (flag_reduce)
    {
      for (int i = 0; i < instance.numRows; i++)
      {
        row_weight[i] = (int)(oblivious_ratio * (double)row_weight[i]);
        if (row_weight[i] ==0) row_weight[i] = 1;
      }
      for (int j = 0; j < instance.numColumns; j++)
      {
        score[j] = (int)(oblivious_ratio * (double)score[j]);
      }
    }

    if (CSbest.num_Cover < CSw.num_Cover)
    {
      CSbest = CSw;
    }

    cout << "Iter: " << iter << " " << CSbest.num_Cover << " " << CSw.num_Cover << endl;
  } // end for iter

  return CSbest;

} // End ObliviousLS





// メイン関数
int main(int argc, char** argv)
{
  //コマンドライン引数の数が少なければ強制終了
  if (argc < 3){
    cout << "Usage: ./command filename K(int)" << endl;
    return 0;
  }
  char *FileName = argv[1];
  FILE *SourceFile = fopen(FileName,"r");
  int K = atoi(argv[2]);

  // SCPのインスタンスを読み込む
  SCPinstance  pData(SourceFile);

  // initialize
  SCPsolution CS(pData, K, 10, 0.6);
  SCPsolution CSimp(pData, K, 10, 0.6);
  // End Initialize;

  int seed = 0;
  Rand rnd;
  rnd.seed(seed);


  greedy_construction(pData, CS, rnd);
  cout << "Greedy: " << CS.num_Cover << ": ";
  CS.print_solution();

  CSimp = ObliviousLS(pData, CS, Weight_Threshold, Oblivious_Ratio, BMS_ratio, numIteration, rnd);
  cout << "LS: " << CSimp.num_Cover << ": ";
  CSimp.print_solution();


  // // //CS.double_layer_selection(pData, rnd);

  // // 結果
  printf("Covers %d rows.\n", CS.num_Cover);
  printf("Covers %d rows.\n", CSimp.num_Cover);



  return 0;
}
