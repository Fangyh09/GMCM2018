#include <iostream>
// matlab


//#include <math.h>
//#include <vector>
#include <algorithm>
//#include <stdio.h>
//#include <cstring>
#include "Genetic.h"

//using namespace std;
//
extern Individual Population[GROUP_SCALE + 1];

extern vector<Edge> selected_edge;
//vector<Edge> selected_edge;
extern vector<Edge> all_edges;

extern int edgemapd2c[61 - 11];
extern int edgemapc2d[61 - 11];

extern string CITY[12];
//int PC[12] = { 3735, 233, 957, 1093, 962, 1604, 90, 3048, 1091, 2420, 3268, 673 };
//int PC[12] = { 100, 1, 60, 100, 962, 1604, 90, 3048, 1091, 2420, 3268, 673 };

extern int PC[12];
extern int DIS[12][12];

extern double  EDGE_PROB[61 - 11];
// for ok function
extern double dp[MAXN][MAXN];
extern double cap[MAXN][MAXN];

extern int xy2idx[MAXN][MAXN];
// for dfs
extern int vis[MAXN];
extern vector<int> G[MAXN];
extern vector<vector<int>> dfs_ans;

extern double B[MAXM];

extern int tot_avail_edge;


//
//const int MAXN = 12 + 1;
//const int MAXM = 12 * 11 + 5;
//int F[MAXM];
//string CITY[12] = { "BeiJing", "WuLuMuQi", "ZhengZhou ", "HaErBin", "XiAn", "ChengDu", "LaSa", "ChongQing", "wuHan", "ShangHai", "GuangZhou", "QunMing" };
//int PC[12] = { 3735, 233, 957, 1093, 962, 1604, 90, 3048, 1091, 2420, 3268, 673 };
//int PC[12] = { 100, 1, 60, 100, 962, 1604, 90, 3048, 1091, 2420, 3268, 673 };
//
//int PC[12] = { 3735, 2398, 9532, 3799, 3813, 8262, 331, 3048, 5885, 2420, 11000, 4771 };
//int DIS[12][12] = { { 0, 2400, 620, 1060, 900, 1500, 2500, 1450, 1050, 1100, 1870, 2100 },
//{ 2400, 0, 2440, 3050, 2110, 2060, 1600, 2340, 2760, 3270, 3280, 2520 },
//{ 620, 2440, 0, 1650, 500, 1000, 2180, 860, 460, 830, 1300, 1500 },
//{ 1060, 3050, 1650, 0, 1960, 2580, 3560, 2500, 2000, 1680, 2800, 3150 },
//{ 900, 2110, 500, 1960, 0, 620, 1750, 580, 650, 1220, 1320, 640 },
//{ 1500, 2060, 1000, 2580, 620, 0, 1250, 300, 980, 1660, 1230, 640 },
//{ 2500, 1600, 2180, 3560, 1750, 1250, 0, 1520, 2230, 2900, 2310, 1270 },
//{ 1450, 2340, 860, 2500, 580, 300, 1520, 0, 720, 1410, 940, 650 },
//{ 1050, 2760, 460, 2000, 650, 980, 2230, 720, 0, 690, 830, 1300 },
//{ 1100, 3270, 830, 1680, 1220, 1660, 2900, 1410, 690, 0, 1220, 1960 },
//{ 1870, 3280, 1300, 2800, 1320, 1230, 2310, 940, 830, 1220, 0, 1080 },
//{ 2100, 2520, 1500, 3150, 1200, 640, 1270, 650, 1300, 1960, 1080, 0 } };
//
//struct Edge
//{
//	int x, y;
//	int idx;
//	double c;
//	double v;
//	Edge() {}
//	Edge(int xin, int yin) {
//       x = xin;
//       y = yin;
//	}
//	void print() {
//		printf("x=%d,y=%d,c=%f\n", x, y, v);
//	}
//};
//
//double getValue(double x, double c1, double c2, double w) {
//	return w * sqrt(c1 * c2) * 0.01 * x;
//}
//
//double mapping(double x) {
//	if (x <= 600) {
//		return 32;
//	}
//	else if (x <= 1200) {
//		return 16;
//	}
//	else if (x <= 3000) {
//		return 8;
//	}
//	else {
//		return 0;
//	}
//}
//
//
int cmp(const Edge& a, const Edge& b) {
	return a.v > b.v;
}
//
//int find(int x) {
//	return F[x] == -1 ? x : F[x] = find(F[x]);
//}
//
//void merge(int x, int y) {
//	int fx = find(x);
//	int fy = find(y);
//	if (fx != fy) {
//		F[fx] = fy;
//	}
//}
//
//vector<int> MaxST(vector<Edge> edges, int n) {
//	vector<int> ans;
//	int len = edges.size();
//	int nums = 0;
//	for (int i = 0; i < len; i++) {
//		cout << "in it" << endl;
//		int x = edges[i].x;
//		int y = edges[i].y;
//		double v = edges[i].v;
//		if (find(x) != find(y)) {
//			merge(x, y);
//			ans.push_back(edges[i].idx);
//			cout << "has it" << endl;
//			nums++;
//		}
//		if (nums >= n - 1) {
//			break;
//		}
//	}
//	return ans;
//}
//
//vector<int> G[MAXN];
//int vis[MAXN];
//double dp[MAXN][MAXN];
//double cap[MAXN][MAXN];
//
//
//vector<vector<int> > dfs_ans;
// beg is not changed
// s is current begin
//
//
extern int A[66][1800];
extern double V[1800];


int F[MAXM];
//
//
//bool ok(vector<int> vec) {
//   for (int i = 0;i < MAXN; i ++) {
//        for (int j = 0;j < MAXN; j ++) {
//            dp[i][j] = -0x3f3f3f3f;
//            cap[i][j] = 0;
//        }
//   }
//   int n = vec.size();
//   int p_s = vec[0];
//   int p_e = vec[n - 1];
//   for (int i = 0;i < n;i ++) {
//       for (int j = 0;j < n;j ++) {
//           if (i == j) continue;
//           int p1 = vec[i];
//           int p2 = vec[j];
//           dp[i][j] = sqrt(PC[p1] * PC[p2]) * mapping(DIS[p1][p2]);
//           cap[i][j] = mapping(DIS[p1][p2]);
//       }
//   }
//   for (int len = 2; len < n; len ++) {
//       for (int i = 0; i + len < n; i ++) {
//           int j = i + len;
//           for (int k = i + 1; k < j; k ++) {
//               dp[i][j] = max(dp[i][j], dp[i][k] + dp[k][j]);
//               if (dp[i][j] < dp[i][k] + dp[k][j]) {
//                  dp[i][j] = dp[i][k] + dp[k][j];
//               }
//            cap[i][j] = min(cap[i][k], cap[k][j]);
//
//           }
//       }
//   }
//   return sqrt(PC[p_s] * PC[p_e]) * cap[0][n-1] >= dp[0][n - 1];
//}
//
//
//void dfs(int beg, int s, int e, vector<int> pres, double sum) {
//	if (sum > sqrt(PC[beg] * PC[e]) * 8) {
//	    return;
//	}
//	if (pres.size() > 4) {
//		return ;
//	}
//	if (pres.size() >= 2 && pres[0] == 0 && pres[1] == 6) {
//		int a = 1;
//	}
//	printf("s = %d, e = %d\n", s, e);
//	if (s == e) {
//		cout << "sum" << sum << endl;
//		cout << "s-e" << sqrt(PC[beg] * PC[e]) << endl;
//		double se_sum = getValue(1, PC[beg], PC[e], 1);
//		if (ok(pres)) {
//			dfs_ans.push_back(pres);
//		}
//		return;
//	}
//	vector<int> nxt = G[s];
//	int n = nxt.size();
//
//	for (int i = 0; i < n; i++) {
//		int p = nxt[i];
//		if (vis[p] == 0) {
//			vis[p] = 1;
//			double cur_val = getValue(mapping(DIS[s][p]), PC[s], PC[p], 1);
//			double nsum = sum + cur_val;
//			vector<int> npres = pres;
//			npres.push_back(p);
//			dfs(beg, p, e, npres, nsum);
//			vis[p] = 0;
//		}
//	}
//}
//
//
int main()
{
   // for matlab
   if(!libMyLinprogInitialize()) {
    std::cout << "Could not initialize libMySin!" << std::endl;
    return -1;
   }

  //freopen("out.txt", "w", stdout);
  // edges
   int N = MAXN;
   int M = NUM_EDGES;
   int cur_m = 0;
   int idx = 0;
   int mark=0;
   for (int i = 0;i < 66; i ++) {B[i] = 0;}
   // for (int i = 0;i < 1800;i ++) {V[i] = 0;}
   // for (int i = 0;i < 66; i ++) {LOW[i] = 0;}
   memset(xy2idx, -1, sizeof(xy2idx));
     tot_avail_edge = 0;
	for (int i = 0; i < N; i++) {
		for (int j = i + 1; j < N; j++) {
			Edge e;
			e.x = i;
			e.y = j;
			e.c = mapping(DIS[i][j]);
			e.v = getValue(e.c, PC[i], PC[j], 1);
      // if (i == 1 || i == 6 || i == 11 || j == 1 || j == 6 || j == 11) {
      //   e.v *= 10;
      // } 

			if (e.c > 0) {
        if (i == 1 && j == 6) {
          mark = e.idx;
        }
				e.idx = idx;
				all_edges.push_back(e);
               xy2idx[e.x][e.y] = e.idx;
               xy2idx[e.y][e.x] = e.idx;
               B[e.idx] = e.c;
               idx += 1;
               tot_avail_edge += 1;
			}
		}
	}

   memset(F, -1, sizeof(F));
   vector<Edge> sorted_edges = all_edges;
   sort(sorted_edges.begin(), sorted_edges.end(), cmp);
   UnionFind uf;
   vector<int> max_st = uf.max_st(sorted_edges, N);
   int max_st_len = max_st.size();
   if (max_st.size() != N - 1) {
     cout << "error in max_st_len" << endl;
     return 1;
   }
   //int max_st_len = max_st.size();
    cur_m += max_st.size();

    int flag[MAXM];
    memset(flag, 0 , sizeof(flag));

   for (int i = 0;i < max_st_len;i ++){
       selected_edge.push_back(all_edges[max_st[i]]);
       flag[max_st[i]] = 1;
   }
   // // add edge
   // for (int i = 0,sz = sorted_edges.size();i < sz && cur_m < M; i ++) {
   //      idx = sorted_edges[i].idx;
   //      if (flag[idx] == 0) {
   //          selected_edge.push_back(all_edges[idx]);
   //          flag[idx] = 1;
   //          cur_m += 1;
   //      }
   //  }


   int cur_pos = 0;
   int cur_nums = 0;
   double sum_prob = 0;

   while (cur_pos < tot_avail_edge) {
      if (flag[cur_pos] == 0) {
        double tmp_val = all_edges[cur_pos].v/ sqrt(PC[all_edges[cur_pos].x] * PC[all_edges[cur_pos].y]);
        tmp_val = tmp_val * max(PC[all_edges[cur_pos].x], PC[all_edges[cur_pos].y]);
        EDGE_PROB[cur_nums] = tmp_val;
        //EDGE_PROB[cur_nums] = 1;
        edgemapd2c[cur_pos] = cur_nums;
        edgemapc2d[cur_nums] = cur_pos;
        sum_prob += all_edges[cur_pos].v;
        cur_nums += 1;
      }
      cur_pos += 1;
   }
   EDGE_PROB[mark] = 11000;
   // for prob
   for (int i = 0;i < cur_nums; i ++) {
     EDGE_PROB[i] /= sum_prob;
   }



//   vector<Edge> can_edge = selected_edge;
//   for (int i = 0,sz = sorted_edges.size();i < sz && cur_m < M; i ++) {
//        idx = sorted_edges[i].idx;
//        if (flag[idx] == 0) {
//            can_edge.push_back(edges[idx]);
//            flag[idx] = 1;
//            cur_m += 1;
//        }
//    }
//    cout << cur_m << endl;
//    for (int i = 0, sz = can_edge.size(); i < sz; i ++) {
//            cout << i << endl;
//        can_edge[i].print();
//    }
    //return  0;

//    selected_edge.push_back(Edge(0, 1));
//    selected_edge.push_back(Edge(0, 6));
//    selected_edge.push_back(Edge(0, 10));
//    selected_edge.push_back(Edge(0, 9));
//    selected_edge.push_back(Edge(0, 3));
//
//    selected_edge.push_back(Edge(0, 7));
//    selected_edge.push_back(Edge(0, 8));
//    selected_edge.push_back(Edge(7, 10));
//    selected_edge.push_back(Edge(8, 10));
//
//    selected_edge.push_back(Edge(0, 11));
//    selected_edge.push_back(Edge(0, 7));
//    selected_edge.push_back(Edge(9, 10));
//    selected_edge.push_back(Edge(7, 9));
//    selected_edge.push_back(Edge(0, 5));

//	for (int i = 0, esz = can_edge.size(); i < esz; i++) {
//		Edge e = can_edge[i];
//		G[e.x].push_back(e.y);
//		G[e.y].push_back(e.x);
//		printf("e.x=%d, e.y=%d\n", e.x, e.y);
//	}
   // return 0;


//   vector<vector<int> > all_dfs_ans;
//   for (int s = 0; s < N; s ++) {
//       for (int e = s + 1; e < N; e ++) {
//       	memset(vis, 0, sizeof(vis));
//       	vis[s] = 1;
//       	vector<int> pres;
//       	pres.push_back(s);
//           dfs_ans = vector<vector<int> > ();
//       	dfs(s, s, e, pres, 0);
//           printf("s=%d, e=%d\n", s, e);
//       	cout << dfs_ans.size() << endl;
//       	for (int i = 0, isz = dfs_ans.size(); i < isz; i++) {
//       		printf("%d\n", i);
//       		vector<int> pres = dfs_ans[i];
//               all_dfs_ans.push_back(pres);
//       		 for (int j = 0, jsz = pres.size(); j < jsz; j++) {
//       		 	printf("%d ", pres[j]);
//       		 }
//       		 printf("\n");
//       	}
//       	 cout << "finish" << endl;
//       }
//   }
  // return 0;
//   cout << "now print" << endl;
//
 //  int tot = 0;
 //  vector<pair<double,int>> vpp;
 //  vector<int> best_sol;
 //  int sol1[16][2] = {{0,3},{0,1},{0,9},{0,2},{0,8}, {0,4},{0,6}, {0,10},
 //                     {2,4}, {2,8}, {4,7}, {5,7}, {7,8}, {7,10}, {8,10},{10, 11}};

 //  for (int i = 0,isz = all_dfs_ans.size(); i < isz; i ++) {
 //       if (all_dfs_ans[i].size() == 2) {
 //           int x1 = all_dfs_ans[i][0];
 //           int x2 = all_dfs_ans[i][1];
 //           for (int j = 0;j < 16;j ++) {
 //               int y1 = sol1[j][0];
 //               int y2 = sol1[j][1];
 //               if ((x1 == y1 && x2 == y2) || (x2 == y1 && x1 == y2)) {
 //                   best_sol.push_back(i);
 //               }
 //           }
 //       }
 //      for (int j = 0,jsz = all_dfs_ans[i].size(); j < jsz;j ++) {
 //          printf("%d ", all_dfs_ans[i][j]);
 //      }
 //      tot += 1;
 //      pair<double,int> pr;
 //      int cur_len = all_dfs_ans[i].size();
 //      pr.first = sqrt(PC[all_dfs_ans[i][0]] * PC[all_dfs_ans[i][cur_len - 1]]);
 //      pr.second = i;
 //      cout << pr.first <<  " " << pr.second << endl;
 //      vpp.push_back(pr);
 //      cout << endl;
 //  }
 // // return 0;
 //  sort(vpp.begin(), vpp.end());
 //  for (int i = 0,sz = vpp.size(); i < sz; i ++) {
 //      int idx = vpp[i].second;
 //      cout << "path ";
 //      for (int j = 0,jsz = all_dfs_ans[idx].size(); j < jsz;j ++) {
 //          printf("%d ", all_dfs_ans[idx][j]);
 //      }
 //      cout << "  cost = " << vpp[i].first <<  endl;
 //  }
 //  //return 0;
 //  cout << "tot" << tot << endl;

 //  memset(A, 0, sizeof(A));
 //  for (int i = 0,isz = all_dfs_ans.size(); i < isz; i ++) {
 //      for (int j = 0,jsz = all_dfs_ans[i].size(); j + 1 < jsz;j ++) {
 //          int p1 = all_dfs_ans[i][j];
 //          int p2 = all_dfs_ans[i][j + 1];
 //          int cur_idx = xy2idx[p1][p2];
 //          A[cur_idx][i] += 1;
 //      }
 //      int link_len = all_dfs_ans[i].size();
 //      V[i] = sqrt(PC[all_dfs_ans[i][0]] * PC[all_dfs_ans[i][link_len - 1]]);
 //      cout << endl;
 //  }

 //  for (int i = 0, sz = selected_edge.size(); i < sz; i ++) {
 //      Edge e = selected_edge[i];
 //      int p1 = e.x;
 //      int p2 = e.y;
 //      LOW[xy2idx[p1][p2]] = 0.01;
 //  }

 // // return  0;

 //  int num_lines = all_dfs_ans.size();
 //  cout << "num_lines" << num_lines << "\n";
 //  FILE *fp1 = freopen("A.txt", "w", stdout);


 //  int nrows = 66 - 5;
 //  int ncols = num_lines;
 //  for (int i = 0; i < nrows; i ++ ){
 //      for (int j = 0; j < ncols; j ++) {
 //          printf("%d%s", A[i][j], j==num_lines-1?"\n":",");
 //      }
 //  }
 //  fclose(fp1);

 //  fp1 = freopen("UPP.txt", "w", stdout);
 //  for (int i = 0; i < nrows; i ++) {
 //      printf("%f%s", B[i], i == nrows - 1? "\n" : ",");
 //  }
 //  fclose(fp1);

 //  fp1 = freopen("V.txt", "w", stdout);
 //  for (int i = 0; i < ncols; i ++) {
 //      printf("%f%s", V[i], i == ncols - 1? "\n" : ",");
 //  }
 //  fclose(fp1);

 //  fp1 = freopen("LOW.txt", "w", stdout);
 //  for (int i = 0; i < nrows; i ++) {
 //      printf("%f%s", LOW[i], i == nrows - 1? "\n" : ",");
 //  }
 //    cout << "result " << endl;
 //   for (int i = 0;i < best_sol.size(); i ++) {
 //       printf("%d%s", best_sol[i], i == best_sol.size() - 1 ? "\n" : " ");
 //   }

 //   fclose(fp1);

    cout << "I'm coming" << endl;
    int Xnration;
    int i;
    int seed = 123456789;

    showTime();
    initGroup(seed);
    evaluate();

    selectBest();

    for (Xnration = 0; Xnration < MAX_GENS; Xnration++)
    {
        selector(seed);
        crossover(seed);
        mutate(seed);
        report(Xnration);
        evaluate();
        elitist();
    }

    cout << "\n";
    cout << "  Best member after " << MAX_GENS << " Xnrations:\n";
    cout << "\n";

    for (i = 0; i < N_VARS; i++)
    {
        cout << "  X(" << i + 1 << ") = " << Population[GROUP_SCALE].Xn[i] << "\n";
    }
    cout << "\n";
    cout << "  Best Fitness = " << Population[GROUP_SCALE].Fitness << "\n";

    showTime();
    libMyLinprogTerminate();
    mclTerminateApplication();
    while (1);
    return 0;
   return 0;
}




//int main() {
//    std::cout << "hello" << std::endl;
//    Flow flow;
//    int x[10] = {1,3};
//    flow.compute_fitness(x, 10);
//}




