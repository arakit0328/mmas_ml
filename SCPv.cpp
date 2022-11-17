#include "SCPv.hpp"
#include <vector>
#include <iostream>
#include <algorithm>
#include "Random.hpp"

//extern Rand rnd;

//
//
//  Class SCPinstance
//
//

// コンストラクタ
SCPinstance::SCPinstance(FILE *SourceFile)
{
  if (SourceFile == NULL) throw DataException();
  else
  {
    int R, C, cost;
    fscanf(SourceFile, "%d", &R);
    fscanf(SourceFile, "%d", &C);
    numRows = R;
    numColumns = C;

    int* nCov = new int [numColumns];
    int* idx = new int [numColumns];
    for (int j = 0; j < numColumns; j++) {
      nCov[j] = 0;
      idx[j] = 0;
    }

    for (int j = 0; j < numColumns; j++)
    {
      std::vector<int> nn;
      Neighborhood.push_back(nn);
    }

    // read costs
    for(int j = 0; j < C; j++)
    {
      fscanf(SourceFile, "%d", &cost);
      Weight.push_back(cost);
    }

    // ファイルから各行の情報を読む
    int CoverNo, CoverID;
    std::vector<int> cov;


    for(int i = 0;  i < numRows; i++)
    {
      if (fscanf(SourceFile, "%d", &CoverNo) == EOF)
        throw (DataException());

      cov.clear();
      for(int j = 0; j < CoverNo; j++)
      {
        if (fscanf(SourceFile,"%d", &CoverID) == EOF)
          throw (DataException());

        if (CoverID >= 1 && CoverID <= numColumns)
        {
          cov.push_back(CoverID - 1);
          nCov[CoverID - 1]++;
        }
        else
          throw (DataException());
      }

      RowCovers.push_back(cov);
    }
    // ファイルの読み込み終了

    // 列の情報を作成
    for (int j = 0; j < numColumns; ++j) {
      std::vector<int> c(nCov[j]);
      ColEntries.push_back(c);
    }

    for (int i = 0; i < numRows; i++)
    {
      for (int c : RowCovers[i]) {
        ColEntries[c][idx[c]] = i;
        idx[c]++;
      }
    }
    // 列の情報の作成終了
    delete [] nCov;
    delete [] idx;
  }


  // 近傍を作る
  for (int j = 0; j < numColumns; ++j)
  {
    for (int r : ColEntries[j])
    {
      for (int c : RowCovers[r])
      {
        if (c != j) Neighborhood[j].push_back(c);
      }
    }
  }

  for (int j = 0; j < numColumns; ++j)
  {
    std::sort(Neighborhood[j].begin(), Neighborhood[j].end());
    Neighborhood[j].erase(std::unique(Neighborhood[j].begin(), Neighborhood[j].end()), Neighborhood[j].end());
  }
  // 近傍の生成終了

  

  
  // 簡単な方法でデータの正しさを確認
  long Sum1 = 0, Sum2 = 0;
  for (int i = 0; i < numRows; i++) Sum1 += RowCovers[i].size();
  for (int j = 0; j < numColumns; j++) Sum2 += ColEntries[j].size();

  //Incorrect Source File!
  if (Sum1 != Sum2)  throw (DataException());

  // 密度の計算
  Density = (float)Sum1/(numColumns * numRows);
}

// End: コンストラクタ


// デストラクタ
//SCPinstance::~SCPinstance()
// End: デストラクタ

// End SCPinstance

//
//
// Class SCPsolution
//
//

// コンストラクタ
SCPsolution::SCPsolution(SCPinstance &instance, int k)
{
  //instance = instance;

  nRow = instance.numRows;
  nCol = instance.numColumns;
  K = k;
  num_Cover = 0;
  
  for (int j = 0; j < nCol; ++j)
  {
    SOLUTION.push_back(0);
  }

  for (int i = 0; i < nRow; i++)
  {
    COVERED.push_back(0);
  }
}


// デストラクタ
SCPsolution::~SCPsolution()
{
}


// 候補解を初期化
void SCPsolution::initialize(SCPinstance &instance)
{
  num_Cover = 0;

  for (int j = 0; j < nCol; ++j)
  {
    SOLUTION[j] = 0;
  }

  for (int i = 0; i < nRow; ++i)
  {
    COVERED[i] = 0;
  }

  // cs には最初は大きい値を詰めておく
  CS.clear();
}




// CSに含まれない列から最大スコアのものを選んで返す
// int SCPsolution::get_column_BMS(SCPinstance &instance, Rand& rnd)
// {
//   std::vector<int> maxCols;
//   int maxScore = -instance.numColumns, maxc = 0;

//   for (int c = 0; c < instance.numColumns; c++)
//   {
//     if (SOLUTION[c] == 0) { continue; }
//     if (conf[c] == 0) { continue; }

//     // 最大スコアの列をチェック
//     if (maxScore < SCORE[c])
//     {
//       maxScore = SCORE[c];
//       maxCols.clear();
//       maxCols.push_back(c);
//     }
//     else if (maxScore == SCORE[c])
//       maxCols.push_back(c);
//   } // End for c

//   if (maxCols.size() == 1) maxc = maxCols[0];
//   else
//   {
//     int j = rnd(0, maxCols.size() - 1);
//     maxc = maxCols[j];
//   }

//   return maxc;
// }


// CSに含まれない列から alpha * 最大スコア 以上のものを選んで返す
// int SCPsolution::get_column_grasp(SCPinstance &instance, double alpha, Rand& rnd)
// {
//   std::vector<int> Cols;
//   int maxScore = 0, minScore = instance.numRows;

//   for (int c = 0; c < instance.numColumns; c++)
//   {
//     if (SOLUTION[c]) { continue; }

//     // 最大スコアと最小スコアを確認
//     if (maxScore < SCORE[c]) maxScore = SCORE[c];
//     else if (minScore > SCORE[c]) minScore = SCORE[c];
//   } // End for c

//   for (int c = 0; c < instance.numColumns; c++)
//   {
//     if (SOLUTION[c]) { continue; }

//     // スコアが大きいものをColsへ格納
//     if (SCORE[c] >= minScore + alpha * (maxScore - minScore))
//       Cols.push_back(c);
//   } // End for c

//   int j = rnd(0, Cols.size() - 1);
//   int col = Cols[j];

//   return col;
// }


// CSに列cを追加する
void SCPsolution::add_column(SCPinstance &instance, int c)
{
  if (SOLUTION[c])
  {
    printf("Column %d has already contained in CS\n", c);
    exit(1);
  }

  SOLUTION[c] = 1;

  // CSに列cを追加
  CS.push_back(c);

  for (int r : instance.ColEntries[c]) {
    COVERED[r]++;
    if (COVERED[r] == 1) {
      num_Cover++;		// カバーされる行の数が増える
    }
  }
} // End add_column


// CSから列cを削除する
void SCPsolution::remove_column(SCPinstance &instance, int c)
{
  if (SOLUTION[c] == 0)
  {
    printf("Column %d is not contained in CS\n", c);
    exit(1);
  }

  SOLUTION[c] = 0;

  // 列cを削除
  CS.erase(remove(CS.begin(), CS.end(), c), CS.end());

  for (int r : instance.ColEntries[c]) {
    COVERED[r]--;
    if (COVERED[r] == 0) {
      num_Cover--;        // カバーされる行の数が増える
    }
  }
} // End remove_column


// CSの中身を表示
void SCPsolution::print_solution()
{
  for (int c : CS)
  {
    printf("%d ", c + 1);
  }
  printf("\n");
} // End print_solution
