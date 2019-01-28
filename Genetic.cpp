/********************************
@author: fangyh09
********************************/

# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>
# include "Genetic.h"
#include <unistd.h>
#include <stdio.h>


string CITY[12] = { "BeiJing", "WuLuMuQi", "ZhengZhou ", "HaErBin", "XiAn", "ChengDu", "LaSa", "ChongQing", "wuHan", "ShangHai", "GuangZhou", "QunMing" };
//int PC[12] = { 3735, 233, 957, 1093, 962, 1604, 90, 3048, 1091, 2420, 3268, 673 };
//int PC[12] = { 100, 1, 60, 100, 962, 1604, 90, 3048, 1091, 2420, 3268, 673 };

int PC[12] = { 3735, 2398, 9532, 3799, 3813, 8262, 331, 3048, 5885, 2420, 11000, 4771 };
int DIS[12][12] = { { 0, 2400, 620, 1060, 900, 1500, 2500, 1450, 1050, 1100, 1870, 2100 },
{ 2400, 0, 2440, 3050, 2110, 2060, 1600, 2340, 2760, 3270, 3280, 2520 },
{ 620, 2440, 0, 1650, 500, 1000, 2180, 860, 460, 830, 1300, 1500 },
{ 1060, 3050, 1650, 0, 1960, 2580, 3560, 2500, 2000, 1680, 2800, 3150 },
{ 900, 2110, 500, 1960, 0, 620, 1750, 580, 650, 1220, 1320, 640 },
{ 1500, 2060, 1000, 2580, 620, 0, 1250, 300, 980, 1660, 1230, 640 },
{ 2500, 1600, 2180, 3560, 1750, 1250, 0, 1520, 2230, 2900, 2310, 1270 },
{ 1450, 2340, 860, 2500, 580, 300, 1520, 0, 720, 1410, 940, 650 },
{ 1050, 2760, 460, 2000, 650, 980, 2230, 720, 0, 690, 830, 1300 },
{ 1100, 3270, 830, 1680, 1220, 1660, 2900, 1410, 690, 0, 1220, 1960 },
{ 1870, 3280, 1300, 2800, 1320, 1230, 2310, 940, 830, 1220, 0, 1080 },
{ 2100, 2520, 1500, 3150, 1200, 640, 1270, 650, 1300, 1960, 1080, 0 } };

double  EDGE_PROB[61 - 11] = {0};
int edgemapd2c[61 - 11];
int edgemapc2d[61 - 11];
// for ok function
double dp[MAXN][MAXN];
double cap[MAXN][MAXN];

int xy2idx[MAXN][MAXN];

vector<Edge> selected_edge;
vector<Edge> all_edges;

double B[MAXM];

// for dfs
int vis[MAXN];
vector<int> G[MAXN];
vector<vector<int>> dfs_ans;
vector<double> dfs_ans_val;

int tot_avail_edge = 0;

// max path
int A[66][1800];
double V[1800];
double LOW[66];
double ZERO[1800];
double DA[66 * 2 * 1800];
double UPPER_LOW[66 * 2];
double matlab_ret[1800];



//申请种群内存，其中多加1个是放置上一代中最优秀个体
struct Individual Population[GROUP_SCALE + 1];

X_Range  XnRange[N_VARS] = { { -3.0,12.1}, {4.1,5.8} };

//有交配权的所有父代进行交叉
void crossover(int &seed)
{
    const double a = 0.0;
    const double b = 1.0;
    int mem;
    int one;
    int first = 0;
    double x;

    for (mem = 0; mem < GROUP_SCALE; ++mem)
    {
        x = randT(0.0,1.0);
        //x = r8_uniform_ab(a, b, seed);//产生交配概率

        if (x < P_MATING)
        {
            ++first;

            if (first % 2 == 0)//交配
            {
                Xover(one, mem, seed);
            }
            else
            {
                one = mem;
            }

        }
    }
    return;
}

//对最差的一代和最优的一代的处理，起到优化的目的
void elitist()
{
    int i;
    double best;
    int best_mem;
    double worst;
    int worst_mem;

    best = Population[0].Fitness;
    worst = Population[0].Fitness;

    for (i = 0; i < GROUP_SCALE - 1; ++i)
    {
        if (Population[i + 1].Fitness < Population[i].Fitness)
        {

            if (best <= Population[i].Fitness)
            {
                best = Population[i].Fitness;
                best_mem = i;
            }

            if (Population[i + 1].Fitness <= worst)
            {
                worst = Population[i + 1].Fitness;
                worst_mem = i + 1;
            }

        }
        else
        {

            if (Population[i].Fitness <= worst)
            {
                worst = Population[i].Fitness;
                worst_mem = i;
            }

            if (best <= Population[i + 1].Fitness)
            {
                best = Population[i + 1].Fitness;
                best_mem = i + 1;
            }

        }

    }

//对于当前代的最优值的处理，如果当前的最优值小于上一代则将上一代的值最优个体取代当前的最弱个体
//基因保留
    if (Population[GROUP_SCALE].Fitness <= best)
    {
        for (i = 0; i < N_VARS; i++)
        {
            Population[GROUP_SCALE].Xn[i] = Population[best_mem].Xn[i];
        }
        Population[GROUP_SCALE].Fitness = Population[best_mem].Fitness;
    }
    else
    {
        for (i = 0; i < N_VARS; i++)
        {
            Population[worst_mem].Xn[i] = Population[GROUP_SCALE].Xn[i];
        }
        Population[worst_mem].Fitness = Population[GROUP_SCALE].Fitness;
    }
    return;
}


//计算适应度值
void evaluate()
{
    int member;
    int i;
    int x[N_VARS + 1];

    for (member = 0; member < GROUP_SCALE; member++)
    {
        for (i = 0; i < N_VARS; i++)
        {
            //x[i + 1] = Population[member].Xn[i];
            x[i] = Population[member].Xn[i];
        }

        //todo
        //Population[member].Fitness = 21.5 + x[1] * sin(4 * PI*x[1]) + x[2] * sin(20 * PI*x[2]);
        // todo
        Flow flow;
        Population[member].Fitness = -flow.compute_fitness(x, N_VARS);
    }
    return;
}


//产生整形的随机数
int i4_uniform_ab(int a, int b, int &seed)
{
    int c;
    const int i4_huge = 2147483647;
    int k;
    float r;
    int value;

    if (seed == 0)
    {
        cerr << "\n";
        cerr << "I4_UNIFORM_AB - Fatal error!\n";
        cerr << "  Input value of SEED = 0.\n";
        exit(1);
    }
    //保证a小于b
    if (b < a)
    {
        c = a;
        a = b;
        b = c;
    }

    k = seed / 127773;
    seed = 16807 * (seed - k * 127773) - k * 2836;

    if (seed < 0)
    {
        seed = seed + i4_huge;
    }

    r = (float)(seed)* 4.656612875E-10;
    //
    //  Scale R to lie between A-0.5 and B+0.5.
    //
    r = (1.0 - r) * ((float)a - 0.5)
        + r   * ((float)b + 0.5);
    //
    //  Use rounding to convert R to an integer between A and B.
    //
    value = round(r);//四舍五入
    //保证取值不越界
    if (value < a)
    {
        value = a;
    }
    if (b < value)
    {
        value = b;
    }

    return value;
}

//初始化种群个体
void initGroup(int &seed)

{
    int i;
    int j;
    double lbound;
    double ubound;
    //
    //  initGroup variables within the bounds
    //
    for (i = 0; i < N_VARS; i++)
    {
        //input >> lbound >> ubound;

        for (j = 0; j < GROUP_SCALE; j++)
        {
            Population[j].Fitness = 0;
            Population[j].ReFitness = 0;
            Population[j].SumFitness = 0;
            // init
            Population[j].Xn[i] = 0;
            //Population[j].Xn[i] = randT(XnRange[i].Lower, XnRange[i].Upper);
            //Population[j].Xn[i] = r8_uniform_ab(XnRange[i].Lower, XnRange[i].Upper, seed);
        }
    }

    for (j = 0; j < GROUP_SCALE;j ++) {
        int cnt1 = 0;
        for (i = 0; i < N_VARS; i ++) {
            double p = randT(0.0, 1.0);
            if (p < EDGE_PROB[i]) {
                Population[j].Xn[i] = 1;
                cnt1 += 1;
            }
            if (cnt1 >= NUM_ONES) {
                break;
            }
            else if (N_VARS - (i + 1) == NUM_ONES - cnt1) {
                for (int k = i + 1; k < N_VARS; k ++) {
                    Population[j].Xn[k] = 1;
                }
                break;
            }
        }
        // check
        int sum = 0;
        for (int i = 0; i < N_VARS; i ++) {
            sum += Population[j].Xn[i];
        }
        if (sum != NUM_ONES) {
            cout << "error !!!" << endl;
            return ;
        }
    }


    return;
}


//挑选出最大值，保存在种群数组的最后一个位置
void selectBest()
{
    int cur_best;
    int mem;
    int i;

    cur_best = 0;

    for (mem = 0; mem < GROUP_SCALE; mem++)
    {
        if (Population[GROUP_SCALE].Fitness < Population[mem].Fitness)
        {
            cur_best = mem;
            Population[GROUP_SCALE].Fitness = Population[mem].Fitness;
        }
    }

    for (i = 0; i < N_VARS; i++)
    {
        Population[GROUP_SCALE].Xn[i] = Population[cur_best].Xn[i];
    }

    return;
}

//个体变异
void mutate(int &seed)
{
    const double a = 0.0;
    const double b = 1.0;
    int i;
    int j;
    double lbound;
    double ubound;
    double x;

    for (i = 0; i < GROUP_SCALE; i++)
    {
        for (j = 0; j < N_VARS; j++)
        {
            //x = r8_uniform_ab(a, b, seed);
            x = randT(a, b);//突变概率
            if (x < P_MUTATION * EDGE_PROB[j])
            {
                if (Population[i].Xn[j] == 1) continue; 
                // int v1 = int(randT(0.0, 0.9999999) * N_VARS);
                int v1 = j;
                int v2 = v1;
                while (v1 == v2 || Population[i].Xn[v1] == Population[i].Xn[v2]) {
                    v2 = int(randT(0.0, 0.9999999) * N_VARS);
                }
                int tmp = Population[i].Xn[v1];
                Population[i].Xn[v1] = Population[i].Xn[v2];
                Population[i].Xn[v2] = tmp;
                //lbound = XnRange[j].Lower;
                //ubound = XnRange[j].Upper;
                //Population[i].Xn[j] = randT(lbound, ubound);
                //Population[i].Xn[j] = r8_uniform_ab(lbound, ubound, seed);

            }
        }

         // check
        int sum = 0;
        for (j = 0; j < N_VARS; j++) {
            sum += Population[i].Xn[j];
        }
        if (sum != NUM_ONES) {
            cout << "error !!!" << endl;
            int a = 1/0;
            return ;
        }
    }

    return;
}

//模板函数，用于生成各种区间上的数据类型
template<typename T>
T randT(T Lower, T Upper)
{
    return rand() / (double)RAND_MAX *(Upper - Lower) + Lower;
}

//产生小数随机数
double r8_uniform_ab(double a, double b, int &seed)

{
    int i4_huge = 2147483647;
    int k;
    double value;

    if (seed == 0)
    {
        cerr << "\n";
        cerr << "R8_UNIFORM_AB - Fatal error!\n";
        cerr << "  Input value of SEED = 0.\n";
        exit(1);
    }

    k = seed / 127773;
    seed = 16807 * (seed - k * 127773) - k * 2836;

    if (seed < 0)
    {
        seed = seed + i4_huge;
    }

    value = (double)(seed)* 4.656612875E-10;

    value = a + (b - a) * value;

    return value;
}

//输出每一代进化的结果
void report(int Xnration)
{
    double avg;
    double best_val;
    int i;
    double square_sum;
    double stddev;
    double sum;
    double sum_square;

    if (Xnration == 0)
    {
        cout << "\n";
        cout << "  Xnration       Best            Average       Standard \n";
        cout << "  number           value           Fitness       deviation \n";
        cout << "\n";
    }
    sum = 0.0;
    sum_square = 0.0;

    for (i = 0; i < GROUP_SCALE; i++)
    {
        sum = sum + Population[i].Fitness;
        sum_square = sum_square + Population[i].Fitness * Population[i].Fitness;
    }

    avg = sum / (double)GROUP_SCALE;
    square_sum = avg * avg * GROUP_SCALE;
    stddev = sqrt((sum_square - square_sum) / (GROUP_SCALE - 1));
    best_val = Population[GROUP_SCALE].Fitness;

    cout << "  " << setw(8) << Xnration
        << "  " << setw(14) << best_val
        << "  " << setw(14) << avg
        << "  " << setw(14) << stddev << "\n";

    return;
}

//选择有交配权的父代
void selector(int &seed)
{
    struct Individual NewPopulation[GROUP_SCALE + 1];//临时存放挑选的后代个体
    const double a = 0.0;
    const double b = 1.0;
    int i;
    int j;
    int mem;
    double p;
    double sum;

    sum = 0.0;
    for (mem = 0; mem < GROUP_SCALE; mem++)
    {
        sum = sum + Population[mem].Fitness;
    }
    //计算概率密度
    for (mem = 0; mem < GROUP_SCALE; mem++)
    {
        Population[mem].ReFitness = Population[mem].Fitness / sum;
    }
    // 计算累加分布，思想是轮盘法
    Population[0].SumFitness = Population[0].ReFitness;
    for (mem = 1; mem < GROUP_SCALE; mem++)
    {
        Population[mem].SumFitness = Population[mem - 1].SumFitness +
            Population[mem].ReFitness;
    }
    // 选择个体为下一代繁殖，选择优秀的可能性大，这是轮盘法的奥秘之处
    for (i = 0; i < GROUP_SCALE; i++)
    {
        p = r8_uniform_ab(a, b, seed);
        if (p < Population[0].SumFitness)
        {
            NewPopulation[i] = Population[0];
        }
        else
        {
            for (j = 0; j < GROUP_SCALE; j++)
            {
                if (Population[j].SumFitness <= p && p < Population[j + 1].SumFitness)
                {
                    NewPopulation[i] = Population[j + 1];
                }
            }
        }
    }
    //更新后代个体
    for (i = 0; i < GROUP_SCALE; i++)
    {
        Population[i] = NewPopulation[i];
    }
    return;
}

//显示系统时间
void showTime()
{
# define TIME_SIZE 40

    static char time_buffer[TIME_SIZE];
    const struct tm *tm;
    size_t len;
    time_t now;

    now = time(NULL);
    tm = localtime(&now);

    len = strftime(time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm);

    cout << time_buffer << "\n";

    return;
# undef TIME_SIZE
}

//交叉产生子代
void Xover(int one, int two, int &seed)
{
    int i;
    int point;
    double t;
    //随机选择交叉点，这里的点是以变量的整个长度为单位
    //point = randT<int>(0, N_VARS - 1);
    int dpos = 0;
    vector<int> v1;
    vector<int> v2;
    for (i = 0;i < N_VARS; i ++) {
        if (Population[one].Xn[i] + Population[two].Xn[i] >= 1 && Population[one].Xn[i] != Population[two].Xn[i]) {
            dpos += 1;
            if (Population[one].Xn[i] == 1) {
                v1.push_back(i);
            }
            else {
                v2.push_back(i);
            }
        }
    }
    if (dpos == 0) return;
    dpos /= 2;
    point = int(randT(0.0, 0.999999) * dpos) + 1;
    //point = randT<int> (0, dpos - 1);
    //point = i4_uniform_ab(0, N_VARS - 1, seed);
    //交叉
//    for (i = 0; i < point; i++)
//    {
//        t = Population[one].Xn[i];
//        Population[one].Xn[i] = Population[two].Xn[i];
//        Population[two].Xn[i] = t;
//    }
    //for one
    for (int j = 0; j < 1; j ++) {
        int i1 = int(randT(0.0, 0.999999) * point);
        int i2 = int(randT(0.0, 0.999999) * point);
        int pos1 = v1[i1];
        int pos2 = v2[i2];
        Population[one].Xn[pos1] = 1 - Population[one].Xn[pos1];
        Population[one].Xn[pos2] = 1 - Population[one].Xn[pos2];
        Population[two].Xn[pos1] = 1 - Population[two].Xn[pos1];
        Population[two].Xn[pos2] = 1 - Population[two].Xn[pos2];
    }

    int sum = 0;
    for (int j = 0; j < N_VARS; j++) {
        sum += Population[one].Xn[j];
    }
    if (sum != NUM_ONES) {
        cout << "error ~~~~~" << endl;
        int a = 1/0;
        return ;
    }

    sum = 0;
    for (int j = 0; j < N_VARS; j++) {
        sum += Population[two].Xn[j];
    }
    if (sum != NUM_ONES) {
        cout << "error ~~~~" << endl;
        int a = 1/0;
        return ;
    }


    return;
}

// **************************************
// ***************flow*******************
// **************************************

vector<Edge> build_edges(int x[], int num) {
    vector<Edge> ans = selected_edge;
    cout << "selected_edge size:" << ans.size() << endl;
    int cnt = 0;
    for (int i = 0;i < num;i ++) {
        if (x[i] == 1) {
            int mi = edgemapc2d[i];
            ans.push_back(all_edges[mi]);
            cnt ++;
        }
    }
    cout << "totoal 1 cnt" << cnt << endl;
    if (ans.size() != NUM_EDGES) {
        cout << "error !!! in NUM_EDGES" << endl;
        int a = 1/0;
    }
    return ans;
}

Flow::Flow() {
    N = 12;
}

Flow::Flow(int xin[], int num_in) {
    x = xin;
    num = num_in;
}

double mapping(double x) {
	if (x <= 600) {
		return 32;
	}
	else if (x <= 1200) {
		return 16;
	}
	else if (x <= 3000) {
		return 8;
	}
	else {
		return 0;
	}
}

double getValue(double x, double c1, double c2, double w) {
	return w * sqrt(c1 * c2) * 0.01 * x;
}

// for dfs

double ok(vector<int> vec) {
   for (int i = 0;i < MAXN; i ++) {
        for (int j = 0;j < MAXN; j ++) {
            dp[i][j] = -0x3f3f3f3f;
            cap[i][j] = 0;
        }
   }
   int n = vec.size();
   int p_s = vec[0];
   int p_e = vec[n - 1];
   for (int i = 0;i < n;i ++) {
       for (int j = 0;j < n;j ++) {
           if (i == j) continue;
           int p1 = vec[i];
           int p2 = vec[j];
           dp[i][j] = sqrt(PC[p1] * PC[p2]) * mapping(DIS[p1][p2]);
           cap[i][j] = mapping(DIS[p1][p2]);
       }
   }
   for (int len = 2; len < n; len ++) {
       for (int i = 0; i + len < n; i ++) {
           int j = i + len;
           for (int k = i + 1; k < j; k ++) {
               dp[i][j] = max(dp[i][j], dp[i][k] + dp[k][j]);
               if (dp[i][j] < dp[i][k] + dp[k][j]) {
                  dp[i][j] = dp[i][k] + dp[k][j];
               }
            cap[i][j] = min(cap[i][k], cap[k][j]);

           }
       }
   }
   if (sqrt(PC[p_s] * PC[p_e]) * cap[0][n-1] >= dp[0][n - 1]) {
    // return sqrt(PC[p_s] * PC[p_e]) * cap[0][n-1];
    return sqrt(PC[p_s] * PC[p_e]);
   }
   return -1;
}

//for dfs

void dfs(int beg, int s, int e, vector<int> pres, double sum) {
	//if (sum > sqrt(PC[beg] * PC[e]) * 8) {
	  //  return;
	//}
    if (beg == 9 and e == 10) {
        int c = 1;
    }
	if (pres.size() > 6) {
		return ;
	}
	if (pres.size() >= 2 && pres[0] == 0 && pres[1] == 6) {
		int a = 1;
	}
	//printf("s = %d, e = %d\n", s, e);
	if (s == e) {
		//cout << "sum" << sum << endl;
		//cout << "s-e" << sqrt(PC[beg] * PC[e]) << endl;
		//double se_sum = getValue(1, PC[beg], PC[e], 1);
        double res = ok(pres);
		if (res > 0) {
			dfs_ans.push_back(pres);
            dfs_ans_val.push_back(res);
		}
		return;
	}
	vector<int> nxt = G[s];
	int n = nxt.size();

	for (int i = 0; i < n; i++) {
		int p = nxt[i];
		if (vis[p] == 0) {
			vis[p] = 1;
			double cur_val = getValue(mapping(DIS[s][p]), PC[s], PC[p], 1);
            // if (s == 1 || s == 6 || s == 11 || p == 1 || p == 6 || p == 11) {
            //     cur_val *= 10;
            // } 
			double nsum = sum + cur_val;
			vector<int> npres = pres;
			npres.push_back(p);
			dfs(beg, p, e, npres, nsum);
			vis[p] = 0;
		}
	}
}

pair<vector<vector<int> >, vector<double> >build_path(int N) {
   vector<vector<int> > all_dfs_ans;
   vector<double> all_dfs_ans_val;
   for (int s = 0; s < N; s ++) {
       for (int e = s + 1; e < N; e ++) {
       	memset(vis, 0, sizeof(vis));
       	vis[s] = 1;
       	vector<int> pres;
       	pres.push_back(s);
        dfs_ans = vector<vector<int> > ();
        dfs_ans_val = vector<double> ();
       	dfs(s, s, e, pres, 0);
       // printf("s=%d, e=%d\n", s, e);
       	//cout << dfs_ans.size() << endl;
       	for (int i = 0, isz = dfs_ans.size(); i < isz; i++) {
       		//printf("%d\n", i);
       		vector<int> pres = dfs_ans[i];
            double pres_val = dfs_ans_val[i];
               all_dfs_ans.push_back(pres);
               all_dfs_ans_val.push_back(pres_val);
       		 // for (int j = 0, jsz = pres.size(); j < jsz; j++) {
       		 // 	printf("%d ", pres[j]);
       		 // }
       		 // printf("\n");
       	}
       	// cout << "finish" << endl;
       }
   }
   cout << "all_dfs_ans.size: " << all_dfs_ans.size() << endl;
   pair<vector<vector<int> >, vector<double> > ret = make_pair(all_dfs_ans, all_dfs_ans_val);
   // return all_dfs_ans;
   return ret;
}


double Flow::compute_fitness(int x[], int num) {
   // cout << "coming" << endl;
    vector<Edge> can_edge = build_edges(x, N_VARS);
   // cout << "ok!!!!" << endl;
    for (int i = 0; i < MAXN; i ++) {
        G[i].clear();
    }
    for (int i = 0, esz = can_edge.size(); i < esz; i++) {
		Edge e = can_edge[i];
		G[e.x].push_back(e.y);
		G[e.y].push_back(e.x);
		//printf("e.x=%d, e.y=%d\n", e.x, e.y);
	}

	pair<vector<vector<int> >, vector<double> > bp_ret = build_path(this->N); // use global G

    vector<vector<int> > all_dfs_ans = bp_ret.first;
    vector<double> all_dfs_ans_val = bp_ret.second;

     cout << "print path" << endl;
     int dim_edge = tot_avail_edge;
     int dim_path = all_dfs_ans.size();
     for (int i = 0, isz = all_dfs_ans.size(); i < isz; i ++) {
        printf("[");
        for (int j = 0, jsz = all_dfs_ans[i].size();j < jsz; j ++) {
            printf("%d%s", all_dfs_ans[i][j], j==jsz-1? "]\n" : ",");
        }
     }
     FILE *fp1 = fopen("path.txt", "a+");  // Open new stdout
     if (fp1 == NULL)
     {
         printf("Error opening file!\n");
         exit(1);
     }
     for (int i = 0, isz = all_dfs_ans.size(); i < isz; i ++) {
        for (int j = 0, jsz = all_dfs_ans[i].size();j < jsz; j ++) {
            fprintf(fp1, "%d%s", all_dfs_ans[i][j], j==jsz-1? ";" : ",");
        }
     }
     fprintf(fp1, "\n");
     fclose(fp1);

    // 61 x
    // build A and V
      memset(A, 0, sizeof(A));

   for (int i = 0;i < 1800;i ++) {V[i] = 0;}
   for (int i = 0;i < 66; i ++) {LOW[i] = 0;}

      for (int i = 0,isz = all_dfs_ans.size(); i < isz; i ++) {
          double cur_val = all_dfs_ans_val[i];
          for (int j = 0,jsz = all_dfs_ans[i].size(); j + 1 < jsz;j ++) {
              int p1 = all_dfs_ans[i][j];
              int p2 = all_dfs_ans[i][j + 1];
              int cur_idx = xy2idx[p1][p2];
              A[cur_idx][i] += 1;
          }
          

          V[i] = cur_val;
          int jsz_sz = all_dfs_ans[i].size();

          // if (all_dfs_ans[i][0] == 1 || all_dfs_ans[i][0] == 6 || all_dfs_ans[i][0] == 11) {
          //   V[i] *= 10;
          // }
          // else if (all_dfs_ans[i][jsz_sz - 1] == 1 || all_dfs_ans[i][jsz_sz - 1] == 6 || all_dfs_ans[i][jsz_sz - 1] == 11) {
          //   V[i] *= 10;
          // }
      }

    // build low
    // 1x61
    for (int i = 0;i < dim_edge; i ++) { LOW[i] = 0; }
    for (int i = 0, sz = selected_edge.size(); i < sz; i ++) {
        Edge e = selected_edge[i];
        int p1 = e.x;
        int p2 = e.y;
        LOW[xy2idx[p1][p2]] = 0.01;
    }

    // B
     int cur_pos = 0;
     for (int j = 0; j < dim_path; j ++) {
        for (int i = 0; i < dim_edge; i ++) {
            DA[cur_pos] = A[i][j];
            cur_pos ++;
        }
        for (int i = 0; i < dim_edge; i ++) {
            DA[cur_pos] = -A[i][j];
            cur_pos ++;
        }
     }

     double *a = DA;
     mwArray ma(dim_edge*2, dim_path, mxDOUBLE_CLASS);
     ma.SetData(a, dim_edge *2* dim_path);


     for (int i = 0;i < dim_edge;i ++) {
        UPPER_LOW[i] = B[i];
     }
     for (int i = 0; i < dim_edge;i ++) {
        UPPER_LOW[dim_edge + i] = -LOW[i];
     }
     cout << "UPPER_LOW" << endl;
     for (int i = 0;i < 2 * dim_edge; i++) {
        cout << UPPER_LOW[i] << " ";
     }
     cout << endl;
     double *b = UPPER_LOW;

     mwArray mb(dim_edge*2, 1, mxDOUBLE_CLASS);
     mb.SetData(b, dim_edge * 2);

     double* f = V;
     //cout << "print f" << endl;
     for (int i = 0;i < dim_path; i ++) {
        f[i] = -f[i] * 0.01;
        // cout << f[i] << " ";
     }
     //cout << endl;
     mwArray mf(1, dim_path, mxDOUBLE_CLASS);
     mf.SetData(f, 1 * dim_path);

     for (int i = 0;i < dim_path;i ++) {ZERO[i] = 0;}
     double* zero = ZERO;
     mwArray mx(1, dim_path, mxDOUBLE_CLASS);
     mx.SetData(zero, dim_path);

     // double y[] = {500, 500, 500};
     // mwArray my(1, 3, mxDOUBLE_CLASS);
     // my.SetData(y, 3);
     mwArray tmp;
     mwArray xx;
     mwArray fval;
     mylinprog(2, xx,fval,mf,ma,mb,mx,tmp);
     cout << xx.NumberOfElements() << endl;
     cout << "fval:" << fval << endl;
     cout << "xx: " << xx << endl;
     double matlab_score[1];
     fval.GetData(matlab_score, 1);
     int num_ret = xx.NumberOfElements();
     xx.GetData(matlab_ret, num_ret);

     FILE* fp2 = fopen("weight.txt", "a+");  // Open new stdout
     for (int i = 0; i < num_ret; i ++) {
        fprintf(fp2, "%f%s", matlab_ret[i], i == num_ret-1?"\n":";");
     }
     fclose(fp2);

     fp2 = fopen("score.txt", "a+");  // Open new stdout
     fprintf(fp2, "%f\n", matlab_score[0]);
     fclose(fp2);

    if (xx.NumberOfElements() > 0) {
        return fval;
    }
    else {
        return 0;
    }

  // return 1.0;
}
















