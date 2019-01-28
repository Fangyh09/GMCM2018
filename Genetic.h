/********************************
@author: fangyh09
********************************/
	
#ifndef _GENETIC_H_
#define _GENETIC_H_
#include <vector>
#include <string>
#include <cstring>
#include "libMyLinprog.h"
#include "mclmcrrt.h"
#include "mclmcr.h"
#include "mclcppclass.h"
#include "matrix.h"
using namespace std;

#define  PI  3.14159265358979323846

//遗传算法参数，种群规模（0~100）、繁殖代数、函数变量个数、交叉概率、编译概率
// # define GROUP_SCALE    50
// # define MAX_GENS       300
# define GROUP_SCALE    5
# define MAX_GENS       500
# define N_VARS         50
# define P_MATING       0.8
# define P_MUTATION     0.5

// # define P_MUTATION     0.15
//# define NUM_ONES  5
const int NUM_EDGES = 33;
const int NUM_ONES = NUM_EDGES - 11;
// flow
const int MAXN = 12 ;
const int MAXM = 12 * 11 / 2;


struct Individual
{
  //  double Xn[N_VARS];      //存放变量值
    int Xn[N_VARS];      //存放变量值
    double Fitness;         //适应值
    double ReFitness;       //适应值概率密度
    double SumFitness;      //累加分布，为轮盘转
};
struct X_Range
{
    double Upper;           //变量的上界取值
    double Lower;           //变量的下界取值
};

template<typename T>
T randT(T Lower, T Upper); //产生任意类型随机数函数

void crossover(int &seed);
void elitist();        //基因保留
void evaluate();

void initGroup(int &seed);

void selectBest();
void mutate(int &seed);

double r8_uniform_ab(double a, double b, int &seed);
int i4_uniform_ab(int a, int b, int &seed);

void report(int Xnration);
void selector(int &seed);
void showTime();
void Xover(int one, int two, int &seed);




struct Edge
{
	int x, y;
	int idx;
	double c;
	double v;
	Edge() {}
	Edge(int xin, int yin) {
       x = xin;
       y = yin;
	}
	void print() {
		printf("x=%d,y=%d,c=%f\n", x, y, v);
	}
};


class UnionFind {
public:
    UnionFind() {}
    // return edge idx
    vector<int> max_st(vector<Edge> edges, int N) {
        memset(F, -1, sizeof(F));
        // compute
        vector<int> ans;
        int len = edges.size();
        int nums = 0;
        for (int i = 0; i < len; i++) {
            int x = edges[i].x;
            int y = edges[i].y;
            if (find(x) != find(y)) {
                merge(x, y);
                ans.push_back(edges[i].idx);
                nums++;
            }
            if (nums >= N - 1) {
                break;
            }
        }
        return ans;
    }

private:
    int F[12 * 11];
    int find(int x) {
        return F[x] == -1 ? x : F[x] = find(F[x]);
    }

    void merge(int x, int y) {
        int fx = find(x);
        int fy = find(y);
        if (fx != fy) {
            F[fx] = fy;
        }
    }
};


class Flow
{
public:
    Flow();
    Flow(int x[], int num);
    double compute_fitness(int x[], int num);
private:
    UnionFind uf;
    int* x;
    int num;
    int N, M;
};



vector<Edge> build_edges(int x[], int num);

double mapping(double x);
double getValue(double x, double c1, double c2, double w);

// for dfs

double ok(vector<int> vec);
//for dfs

void dfs(int beg, int s, int e, vector<int> pres, double sum);
pair<vector<vector<int> >, vector<double> > build_path(int N);
//double Flow::compute_fitness(int x[], int num);


#endif // !_GENETIC_H_
