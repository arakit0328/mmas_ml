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


// ObliviousLS で使う繰り返しの回数
int numIteration = 5000;

// row_weight の最大値
const int Weight_Threshold = 1000;
const double Oblivious_Ratio = 0.3;

// BMS
const double BMS_ratio = 0.66;




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



// BMS selection, the selected column is removed from CS
int get_column_BMS(SCPinstance& instance,
                   SCPsolution& cs,
                   vector<int>& score,
                   vector<int>& subscore,
                   vector<int>& conf,
                   vector<int>& time,
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
      vector<int> minTime;
      int mt = numeric_limits<int>::max();

      for (int c : maxSub)
      {
        if (time[c] < mt)
        {
          mt = time[c];
          minTime.clear();
          minTime.push_back(c);
        }
        else if (time[c] == mt)
          minTime.push_back(c);
      }

      int j = rnd(0, minTime.size() - 1);
      maxc = minTime[j];
    }
  }

  return maxc;
}


// スコア最大のものを選択する
int get_column_maxscore(SCPinstance& inst,
                        SCPsolution& cs,
                        vector<int>& score,
                        vector<int>& subscore,
                        vector<int>& time,
                        Rand& rnd)
{
  std::vector<int> maxCols;
  int maxScore = -inst.numColumns, maxc = 0;

  for (int c = 0; c < inst.numColumns; c++)
  {
    if (cs.SOLUTION[c] == 1) { printf("Column %d should not in CS.\n", c); exit(1); }

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
    int j = rnd() % maxCols.size();
    maxc = maxCols[j];
  }

  return maxc;
}



// r_select行をカバーする列からスコア最大のものを選択する
int get_column_maxscore(SCPinstance& inst,
                        SCPsolution& cs,
                        int r_select,  // the selected row
                        vector<int>& score,
                        vector<int>& subscore,
                        vector<int>& time,
                        vector<int>& conf,
                        Rand& rnd)
{
  std::vector<int> maxCols;
  int maxScore = -inst.numColumns, maxc = 0;

  for (int c : inst.RowCovers[r_select])
  {
    if (cs.SOLUTION[c] == 1) { printf("Column %d should not in CS.\n", c); exit(1); }
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
      vector<int> minTime;
      int mt = numeric_limits<int>::max();

      for (int c : maxSub)
      {
        if (time[c] < mt)
        {
          mt = time[c];
          minTime.clear();
          minTime.push_back(c);
        }
        else if (time[c] == mt)
          minTime.push_back(c);
      }

      int j = rnd(0, minTime.size() - 1);
      maxc = minTime[j];
    }
  }

  return maxc;
}


// Update score when column c is added
void add_update_score(const SCPinstance& inst,
                      const int c,
                      const SCPsolution& cs,
                      vector<int>& score,
                      vector<int>& subscore,
                      vector<int>& row_weight)
{
  // スコアを更新
  score[c] = -score[c];
  subscore[c] = -subscore[c];

  // スコア更新
  for (int r : inst.ColEntries[c]) // 列cがカバーする行
  {
    // r行が初めてカバーされたら，rを含む列のスコアを減少
    if (cs.COVERED[r] == 1)
    {
      for (int rc : inst.RowCovers[r])
      {
        if (rc != c)
        {
          score[rc] -= row_weight[r];
          subscore[rc]++;     // oneTotwo
        }
      } // End: for ri
    } // End if covered[r] == 1

      // r行が2回カバーされたら，rを含むCSwの要素のスコアを増加
    if (cs.COVERED[r] == 2)
    {
      for (int rc : inst.RowCovers[r]) // r行をカバーする列
      {
        if (cs.SOLUTION[rc])
        {
          score[rc] += row_weight[r];
          subscore[rc]--;     // twoToone
        }
      }
    } // End if covered[r] == 2
  }
}


// Update score when column c is removed
void remove_update_score(const SCPinstance& inst,
                         const int c,
                         const SCPsolution& cs,
                         vector<int>& score,
                         vector<int>& subscore,
                         vector<int>& row_weight)
{
  score[c] = -score[c];
  subscore[c] = -subscore[c];

  // スコア更新
  for (int r : inst.ColEntries[c]) // 列cがカバーする行
  {
    // r行がカバーされなくなったら，rを含む行のスコアを増加
    if (cs.COVERED[r] == 0)
    {
      for (int rc : inst.RowCovers[r]) // r行をカバーする列
      {
          if (rc != c) score[rc] += row_weight[r];
      } // End: for ri
    } // End if covered[r] == 0

    // r行が1回カバーされたら，rを含むCSの要素のスコアを減少
    if (cs.COVERED[r] == 1)
    {
      for (int rc : inst.RowCovers[r]) // r行をカバーする列
      {
        if (cs.SOLUTION[rc])
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
      if (cs.COVERED[r] == 2)
      {
        for (int rc : inst.RowCovers[r]) // r行をカバーする列
        {
          if (cs.SOLUTION[rc])
          {
            subscore[rc]--;        // twoToone
          }
        }
      } // End if covered[r] == 1
  }
}


// CSに含まれない列から最大スコアのものを選んで返す
int simple_column_maxscore(SCPinstance &inst,
                           SCPsolution& cs,
                           vector<int>& score,
                           Rand& rnd)
{
  std::vector<int> maxCols;
  int maxScore = 0, maxc = 0;

  for (int c = 0; c < inst.numColumns; c++)
  {
    if (cs.SOLUTION[c]) { continue; }

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
    int j = rnd(0, maxCols.size() - 1);
    maxc = maxCols[j];
  }

  return maxc;
}


// Select K columns at random.
// 引数の cs に結果が入る
SCPsolution random_construction(SCPinstance& inst,
                                int k,
                                Rand &rnd)
{
  SCPsolution cs(inst, k);
  vector<int> cols(inst.numColumns);

  for (int j = 0; j < inst.numColumns; ++j) cols[j] = j;
  random_permutation(cols, rnd);
  for (int k = 0; k < cs.K; k++) cs.add_column(inst, cols[k]);

  return cs;
}


// conf = 1 for all neighbors of c
void update_conf(SCPinstance& inst, int c, vector<int>& conf)
{
  for (int cn : inst.Neighborhood[c]) conf[cn] = 1;
}


SCPsolution DoubleLayerSelection(SCPinstance& inst,
                                 int k,
                                 Rand& rnd)
{
  SCPsolution cs(inst, k);
  vector<int> score(inst.numColumns, 0);
  int c;


  for (int c = 0; c < inst.numColumns; c++)
  {
    score.push_back(inst.ColEntries[c].size());
  }

  for (int i = 0; i < k; i++)
  {
    std::vector<int> maxCols;
    int maxScore = 0;

    for (int c = 0; c < inst.numColumns; c++)
    {
      if (cs.SOLUTION[c] == 1) { continue; }

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

    if (maxCols.size() == 1) c = maxCols[0];
    else
    {
      int j = rnd() % maxCols.size();
      c = maxCols[j];
    }
    // cの選択終了

    cs.add_column(inst, c);

    // スコアを更新
    for (int r : inst.ColEntries[c]) // 列cがカバーする行
    {
      if (cs.COVERED[r] == 1)
        for (int rc : inst.RowCovers[r])
          if (rc != c) score[rc]--;
    } // End if covered[r] == 1
  }

  return cs;
}




// Oblivious local search
SCPsolution ObliviousLS(SCPinstance& inst,
                        SCPsolution& CS,
                        int weight_threshold,
                        double oblivious_ratio,
                        double bms_ratio,
                        int niter,
                        Rand& rnd)
{
  // conf, score, subscore, weight はここで定義すればいい
  vector<int> score(inst.numColumns, 0);
  vector<int> subscore(inst.numColumns, 0);
  vector<int> row_weight(inst.numRows, 1);
  vector<int> conf(inst.numColumns, 1);
  vector<int> time(inst.numColumns, 0);

  // cs からスコアの初期値を求める
  for (int j = 0; j < inst.numColumns; j++)
  {
    if (CS.SOLUTION[j])
    {
      for (int r : inst.ColEntries[j])
      {
        if (CS.COVERED[r] == 1) score[j]--;
        if (CS.COVERED[r] == 2) subscore[j]--; // twoToone
      }
    }
    else
    {
      // 列jが解に入っていない
      for (int r : inst.ColEntries[j])
      {
        if (CS.COVERED[r] == 0) score[j]++;
        if (CS.COVERED[r] == 1) subscore[j]++; // oneTotwo
      }
    }
  }

  SCPsolution CSbest(inst, CS.K);
  SCPsolution CSw(inst, CS.K);

  CSw = CS;
  CSbest = CSw;

  for (int iter = 0; iter < niter; iter++)
  {
    //cout << "Iter: " << iter << " ";
    // 列を削除
    int c = get_column_BMS(inst, CSw, score, subscore, conf, time, bms_ratio, rnd);
    CSw.remove_column(inst, c);
    remove_update_score(inst, c, CSw, score, subscore, row_weight);
    // スコア更新終了

    time[c] = iter;

    // confの更新
    conf[c] = 0;
    update_conf(inst, c, conf);
    // End: remove a column

    //cout << "remove " << c << " ";


    // 列を追加
    int x = rnd(1, inst.numRows - CSw.num_Cover); // カバーされていない行数までの乱数
    int r_select = 0;           // 選択する行番号 カバーされていない行のx番目になる
    while (x > 0)
    {
      if (CSw.COVERED[r_select] == 0) x--;
      r_select++;
    }
    r_select--;

    assert(CSw.COVERED[r_select] == 0);

    c = get_column_maxscore(inst, CSw, r_select, score, subscore, time, conf, rnd);
    CSw.add_column(inst, c);
    add_update_score(inst, c, CSw, score, subscore, row_weight);

    time[c] = iter;

    // confの更新
    for (int ci : CSw.CS) conf[ci] = 1; // 候補解に含まれる列のconfはc以外は1になる
    conf[c] = 0;
    update_conf(inst, c, conf);
    // confの更新終了

    //cout << " add " << c << " ";

    // row_weightとscoreの更新
    bool flag_reduce = false;
    for (int i = 0; i < inst.numRows; i++)
    {
      if (CSw.COVERED[i] == 0)
      {
        row_weight[i]++;
        if (row_weight[i] > weight_threshold) flag_reduce = true;
        for (int c : inst.RowCovers[i])
        {
          if (CSw.SOLUTION[c]) score[c]--;
          else score[c]++;
        }
      }
    }

    if (flag_reduce)
    {
      for (int i = 0; i < inst.numRows; i++)
      {
        row_weight[i] = (int)(oblivious_ratio * (double)row_weight[i]);
        if (row_weight[i] ==0) row_weight[i] = 1;
      }
      for (int j = 0; j < inst.numColumns; j++)
      {
        score[j] = (int)(oblivious_ratio * (double)score[j]);
      }
    }

    if (CSbest.num_Cover < CSw.num_Cover)
    {
      CSbest = CSw;
    }

    //cout << CSbest.num_Cover << " " << CSw.num_Cover << endl;
  } // end for iter

  return CSbest;

} // End ObliviousLS


// Oblivious local search
SCPsolution SimpleLS(SCPinstance& inst,
                     SCPsolution& CS,
                     double bms_ratio,
                     int niter,
                     Rand& rnd)
{
  // conf, score, subscore, weight はここで定義すればいい
  vector<int> score(inst.numColumns, 0);
  vector<int> subscore(inst.numColumns, 0);
  vector<int> row_weight(inst.numRows, 1);
  vector<int> conf(inst.numColumns, 1);
  vector<int> time(inst.numColumns, 0);

  // cs からスコアの初期値を求める
  for (int j = 0; j < inst.numColumns; j++)
  {
    if (CS.SOLUTION[j])
    {
      for (int r : inst.ColEntries[j])
      {
        if (CS.COVERED[r] == 1) score[j]--;
        if (CS.COVERED[r] == 2) subscore[j]--; // twoToone
      }
    }
    else
    {
      // 列jが解に入っていない
      for (int r : inst.ColEntries[j])
      {
        if (CS.COVERED[r] == 0) score[j]++;
        if (CS.COVERED[r] == 1) subscore[j]++; // oneTotwo
      }
    }
  }

  SCPsolution CSbest(inst, CS.K);
  SCPsolution CSw(inst, CS.K);

  CSw = CS;
  CSbest = CSw;

  for (int iter = 0; iter < niter; iter++)
  {
    //cout << "Iter: " << iter << " ";
    // 列を削除
    int c = get_column_BMS(inst, CSw, score, subscore, conf, time, bms_ratio, rnd);
    CSw.remove_column(inst, c);
    remove_update_score(inst, c, CSw, score, subscore, row_weight);
    // スコア更新終了

    time[c] = iter;

    // confの更新
    conf[c] = 0;
    update_conf(inst, c, conf);
    // End: remove a column

    //cout << "remove " << c << " ";


    // 列を追加
    int x = rnd(1, inst.numRows - CSw.num_Cover); // カバーされていない行数までの乱数
    int r_select = 0;           // 選択する行番号 カバーされていない行のx番目になる
    while (x > 0)
    {
      if (CSw.COVERED[r_select] == 0) x--;
      r_select++;
    }
    r_select--;

    assert(CSw.COVERED[r_select] == 0);

    c = get_column_maxscore(inst, CSw, r_select, score, subscore, time, conf, rnd);
    CSw.add_column(inst, c);
    add_update_score(inst, c, CSw, score, subscore, row_weight);

    time[c] = iter;

    // confの更新
    for (int ci : CSw.CS) conf[ci] = 1; // 候補解に含まれる列のconfはc以外は1になる
    conf[c] = 0;
    update_conf(inst, c, conf);
    // confの更新終了

    //cout << " add " << c << " ";


    if (CSbest.num_Cover < CSw.num_Cover)
    {
      CSbest = CSw;
    }

    //cout << CSbest.num_Cover << " " << CSw.num_Cover << endl;
  } // end for iter

  return CSbest;

} // End SimpleLS



SCPsolution MMAS_ML(SCPinstance& inst,
                    int k,
                    Rand& rnd)
{
  SCPsolution cs_best(inst, k);
  SCPsolution cs = DoubleLayerSelection(inst, k, rnd);
  cout << "End DoubleLayerSelection" << endl;

  SCPsolution cs_0 = ObliviousLS(inst, cs, Weight_Threshold, Oblivious_Ratio, BMS_ratio, numIteration, rnd);
  SCPsolution cs_1 = SimpleLS(inst, cs, BMS_ratio, numIteration, rnd);

  cout << "CS0 covers " << cs_0.num_Cover << endl;
  cout << "CS1 covers " << cs_1.num_Cover << endl;

  if (cs_0.num_Cover > cs_1.num_Cover) cs_best = cs_0;
  else cs_best = cs_1;

  return cs_best;
}




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
  SCPinstance  instance(SourceFile);

  // initialize
  SCPsolution CS(instance, K);
  SCPsolution CSimp(instance, K);
  // End Initialize;

  int seed = 0;
  Rand rnd;
  rnd.seed(seed);

  CSimp = MMAS_ML(instance, K, rnd);

  // 結果
  printf("Covers %d rows.\n", CSimp.num_Cover);



  return 0;
}
