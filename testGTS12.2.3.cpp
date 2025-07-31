#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include <cmath>
#include <set>
#include <list>
#include <deque>  // 用于支持双端队列
//#include <Eigen/Dense> // 需要Eigen库，用于矩阵秩计算
#include <queue>  // 引入队列用于循环禁忌表
//#include "Eigen/Dense"
#include <stdexcept>
#include <chrono>
#include <numeric> // 添加头文件
//#include <fstream>
#include <fstream>
#include <sstream>

using namespace std;

#define GENE_LENGTH 352        // 旋转对称布尔函数的等价类数量
#define FULL_LENGTH 4096       // 13变量布尔函数的全集大小

ofstream process("/Users/mac/Desktop/study/科研/3.18汇报/testGTS/testGTS12.2/改成旋转对称、单点均匀、单点洗牌/process.txt");
ofstream result("/Users/mac/Desktop/study/科研/3.18汇报/testGTS/testGTS12.2/改成旋转对称、单点均匀、单点洗牌/result.txt");

// 循环左移函数：将整数 x 的二进制表示向左循环移位 1 位
// 输入: x - 要移位的整数；n - 二进制向量的位数
// 输出: 左移后的整数，最高位移到最低位（完成循环移位）
int cyclicShiftLeft(int x, int n) {
    return ((x << 1) & ((1 << n) - 1)) | (x >> (n - 1));
    // (x << 1): 将 x 的二进制表示整体左移 1 位，右侧补 0。
    // ((1 << n) - 1): 生成一个低 n 位全为 1 的掩码，限制结果在 n 位范围内。
    // (x >> (n - 1)): 将原始 x 的最高位移到最低位。
    // |: 将左移后的结果与右移的结果结合，完成循环移位。
}

// 生成旋转对称布尔函数的等价类代表元素
// 输入: representatives - 存储代表元素的数组；n - 布尔函数的变量个数
// 输出: 代表元素数组填充了等价类的代表
void computeCycleRepresentatives(int *representatives, int n) {
    int count = 0;  // 当前找到的等价类代表数量
    char visited[FULL_LENGTH] = {0};  // 标记数组，记录每个输入向量是否已被访问

    for (int i = 0; i < FULL_LENGTH; i++) { // 遍历所有可能的输入向量 (0 到 2^n - 1)
        if (visited[i]) continue;  // 如果当前输入向量已经被处理过，跳过

        int current = i;  // 当前输入向量
        int minVal = i;   // 当前等价类的最小值，初始化为自身

        do {
            visited[current] = 1;  // 标记当前向量为已访问
            current = cyclicShiftLeft(current, n);  // 对当前向量执行循环左移
            if (current < minVal) minVal = current; // 更新当前等价类的最小值
        } while (current != i);  // 如果循环移位返回到起始值，结束循环

        representatives[count++] = minVal; // 保存当前等价类的最小值作为代表元素

        if (count == GENE_LENGTH) break; // 如果已找到所有等价类的代表元素，提前退出
    }
}

void expandToFull8192(int *reduced, int *full, int *representatives,int n) {
    for (int i = 0; i < FULL_LENGTH; i++) { // 遍历 8192 个输入向量
        int current = i;
        int minVal = i;

        // 找到 i 的最小循环移位
        for (int j = 0; j < GENE_LENGTH; j++) {
            current = cyclicShiftLeft(current, n);
            if (current < minVal) minVal = current;
        }

        // 检查最小循环移位是否是代表元素
        for (int k = 0; k < GENE_LENGTH; k++) {
            if (minVal == representatives[k]) {
                full[i] = reduced[k]; // 将代表元素对应的值赋给 full[i]
                break;
            }
        }
    }
}


class BooleanFunction {
public:
    int n;
    int genes[GENE_LENGTH];
    int fullGenes[FULL_LENGTH];
    int representatives[GENE_LENGTH];
    vector<int> table;


    BooleanFunction(int n) : n(n) {

    	computeCycleRepresentatives(representatives, n);
        // 初始化基因型（保证平衡）
        do {
            for(int i=0; i<GENE_LENGTH; i++){
                genes[i] = rand()%2;
                //if(genes[i]) ones++;
            }
            // 将 632 位扩展到完整的 8192 位布尔函数
    		expandToFull8192(genes, fullGenes, representatives,n);
            // 将 int 数组转换为 vector<int>
    		vector<int> foff(begin(fullGenes), end(fullGenes));
    		table=foff;
        } while(!isBalanced(table));

	}

	// 检查布尔函数的真值表是否平衡
bool isBalanced(const vector<int>& f) {
    int ones_count = 0;
    int zeroes_count = 0;
    for (int i = 0; i < f.size(); i++) {
        if (f[i] == 1) ones_count++;
        else zeroes_count++;
    }
    return ones_count == zeroes_count;  // 检查 1 和 0 的数量是否相等
}

/*
    bool isBalanced(int []fullGenes,FULL_LENGTH) {
    	int ones_count = 0;
    	int zeroes_count = 0;
    	for (int i = 0; i < FULL_LENGTH; i++) {
        	if (f[i] == 1)
				ones_count++;
        	else
				zeroes_count++;
    	}
    	return ones_count == zeroes_count;  // 检查 1 和 0 的数量是否相等

    }
    */
	// 从文件加载布尔函数真值表
static BooleanFunction LoadFromFile(const string& filename, int n) {
    BooleanFunction bf(n);
    ifstream file(filename);
    if (!file) {
        cerr << "无法打开文件: " << filename << endl;
        exit(1);
    }

    int bit;
    vector<int> truth_table;
    while (file >> bit) {
        truth_table.push_back(bit);
    }

    if (truth_table.size() != (1 << n)) {
        cerr << "文件内容长度错误，应该是 " << (1 << n) << " 但读取到 " << truth_table.size() << endl;
        exit(1);
    }

    bf.table = truth_table;
    return bf;
}
	static BooleanFunction BentFunction(int n) {
        if (n % 2 != 0) {
            throw invalid_argument("n must be even for Bent function.");
        }

        BooleanFunction bf(n);
        int size = 1 << n; // 2^n

        // 递归生成Sylvester Hadamard矩阵的第一行
        vector<int> hadamard_row = {1};
        for (int k = 1; k <= n; ++k) {
            vector<int> new_row;
            for (int val : hadamard_row) {
                new_row.push_back(val);
                new_row.push_back(val);
            }
            hadamard_row = new_row;
            for (size_t i = hadamard_row.size() / 2; i < hadamard_row.size(); ++i) {
                hadamard_row[i] *= -1;
            }
        }

        // 选择第一行作为Bent函数（映射±1到0/1）
        for (int i = 0; i < size; ++i) {
            bf.table[i] = (hadamard_row[i] == 1) ? 0 : 1; // Bent函数真值表
        }

        return bf;
    }
	/*
    double fitness() {
        return non_linearity() - auto_correlation();
    }

    double non_linearity() {
        return accumulate(table.begin(), table.end(), 0); // 伪实现
    }

    double auto_correlation() {
        double sum = 0;
        for (size_t i = 0; i < table.size(); ++i) {
            sum += abs(table[i]);
        }
        return sum; // 伪实现
    }
    */
    vector<int> compute_walsh_spectrum() const {
        vector<int> W(table.begin(), table.end());
        for (int &v : W) v = (v == 1)? -1 : 1;
        int N = W.size();
        for (int len = 1; len < N; len *= 2) {
            for (int i = 0; i < N; i += 2 * len) {
                for (int j = 0; j < len; ++j) {
                    int a = W[i + j];
                    int b = W[i + j + len];
                    W[i + j] = a + b;
                    W[i + j + len] = a - b;
                }
            }
        }
        return W;
    }

    double cost1() const {
        vector<int> W = compute_walsh_spectrum();
        double sum = 0;
        for (int w : W) {
            sum += pow(abs(w), 4);
        }
        return sum / pow(2, 4 * n);
    }
    /*
    double cost2() {
        vector<int> W = compute_walsh_spectrum();
        int count = 0;
        for (int w : W) {
            if (w == 0) count++;
        }
        return count;
    }
	*/

	// GF(2)高斯消元计算矩阵的秩
static int gf2_rank(vector<vector<int>> matrix) {
    int rows = matrix.size();
    if (rows == 0) return 0;
    int cols = matrix[0].size();
    int rank = 0;

    for (int col = 0; col < cols; ++col) {
        // 寻找主元（首个该列为1的行）
        int pivot = -1;
        for (int i = rank; i < rows; ++i) {
            if (matrix[i][col] == 1) {
                pivot = i;
                break;
            }
        }
        if (pivot == -1) continue; // 该列全0

        // 交换行
        swap(matrix[rank], matrix[pivot]);

        // 消去其他行
        for (int i = 0; i < rows; ++i) {
            if (i != rank && matrix[i][col] == 1) {
                for (int j = col; j < cols; ++j) {
                    matrix[i][j] ^= matrix[rank][j]; // 异或操作（GF(2)减法）
                }
            }
        }
        rank++;
    }
    return rank;
}

int cost2() const {
    vector<int> W = compute_walsh_spectrum();

    // 收集满足 W_f(ω)=0 的向量（二进制表示）
    vector<vector<int>> vectors;
    for (int i = 0; i < W.size(); ++i) {
        if (W[i] == 0) {
            vector<int> vec(n);
            for (int j = 0; j < n; ++j) {
                vec[j] = (i >> j) & 1; // 转换为二进制向量
            }
            vectors.push_back(vec);
        }
    }

    // 计算GF(2)下的秩
    return BooleanFunction::gf2_rank(vectors);
}

void ensure_balance() {
    int ones_count = count(table.begin(), table.end(), 1);
    int half_size = (1 << n) / 2;

    if (ones_count != half_size) {
        int diff = abs(ones_count - half_size);
        random_device rd;
        mt19937 gen(static_cast<unsigned int>(chrono::steady_clock::now().time_since_epoch().count()));
        uniform_int_distribution<int> dist(0, table.size() - 1);

        if (ones_count > half_size) {
            // 1 过多，随机挑选部分 1 变为 0
            while (diff > 0) {
                int idx = dist(gen);
                if (table[idx] == 1) {
                    table[idx] = 0;
                    diff-=2;
                }
            }
        } else {
            // 0 过多，随机挑选部分 0 变为 1
            while (diff > 0) {
                int idx = dist(gen);
                if (table[idx] == 0) {
                    table[idx] = 1;
                    diff-=2;
                }
            }
        }
    }
}
	/*
	int cost2() const {
    	vector<int> W = compute_walsh_spectrum();

    	// 收集所有满足W_f(ω)=0的向量ω（二进制表示）
    	vector<Eigen::VectorXd> vectors;
    	for (int i = 0; i < W.size(); ++i) {
        	if (W[i] == 0) {
            	Eigen::VectorXd vec(n);  // n为变量数
            	for (int j = 0; j < n; ++j) {
                	vec(j) = (i >> j) & 1;  // 将i转换为二进制向量
            	}
            	vectors.push_back(vec);
        	}
    	}

    	// 构造矩阵并计算秩
    	Eigen::MatrixXd matrix(vectors.size(), n);
    	for (int i = 0; i < vectors.size(); ++i) {
        	matrix.row(i) = vectors[i];
    	}
    	return matrix.fullPivLu().rank();
	}
	*/
	//突变
    void mutate(double mutation_rate) {
    	if (rand() % 3 == 0) {
        bool ph=false;
    	while(ph==false){
    		// 单点突变：随机选择一个位置，将该位置的值取反
        	int idx = rand() % GENE_LENGTH;  // 随机选择突变位置
        	genes[idx] = 1 - genes[idx];  // 取反操作
        	// 将 632 位扩展到完整的 8192 位布尔函数
    		expandToFull8192(genes, fullGenes, representatives,n);
			//cout << "已将 632 位扩展到完整的 8192 位布尔函数\n";
			// 将 int 数组转换为 vector<int>
    		vector<int> foff(begin(fullGenes), end(fullGenes));
    		table=foff;
    		ph=isBalanced(table);
		}


    } else {
        bool ph=false;
    	while(ph==false){
    		// 洗牌突变：随机选择子串并随机打乱子串顺序
        	int start = rand() % GENE_LENGTH;  // 子串起始位置
        	int end = rand() % GENE_LENGTH;    // 子串结束位置
        	if (start > end) {  // 确保 start <= end
            	int temp = start;
            	start = end;
            	end = temp;
        	}

        	// 提取子串到临时数组中
        	int length = end - start + 1;
        	int subArray[length];
        	for (int i = 0; i < length; i++) {
            	subArray[i] = genes[start + i];
        	}

        	// 随机打乱子串
        	for (int i = length - 1; i > 0; i--) {
            	int j = rand() % (i + 1);
            	int temp = subArray[i];
            	subArray[i] = subArray[j];
            	subArray[j] = temp;
        	}

        	// 将打乱后的子串写回到原基因数组
        	for (int i = 0; i < length; i++) {
            	genes[start + i] = subArray[i];
        	}
        	// 将 632 位扩展到完整的 8192 位布尔函数
    		expandToFull8192(genes, fullGenes, representatives,n);
			//cout << "已将 632 位扩展到完整的 8192 位布尔函数\n";
			// 将 int 数组转换为 vector<int>
    		//vector<int> foff(begin(fullGenes), end(fullGenes));
    		vector<int> foff(fullGenes, fullGenes+FULL_LENGTH);
    		table=foff;
    		ph=isBalanced(table);
		}

    }
    	/*
        random_device rd;
        mt19937 gen(static_cast<unsigned int>(chrono::steady_clock::now().time_since_epoch().count()));
        uniform_real_distribution<double> dist(0.0, 1.0);
        for (int &bit : table) {
            if (mutation_rate >= dist(gen)  ) {
                bit = 1 - bit;
            }
        }
        // **保持平衡**
    	ensure_balance();
    	*/
    }
    //均匀交叉
    //交叉操作，使用已生成的Bent函数
    BooleanFunction crossover_with_Bent(const BooleanFunction &parent1, const BooleanFunction &parent2) {
       BooleanFunction child(n);  // 创建子代布尔函数
	   if (rand() % 3 == 0) {
    	bool ph=false;
    	while(ph==false){
    		// 单点交叉：随机选择一个位置，将父母的前半段和后半段交叉
        	int point = rand() % GENE_LENGTH;  // 随机选择交叉点
        	for (int i = 0; i < GENE_LENGTH; i++) {
            	child.genes[i] = (i < point) ? parent1.genes[i] : parent2.genes[i];
        	}
        	// 将 632 位扩展到完整的 8192 位布尔函数
    		expandToFull8192(child.genes, child.fullGenes, child.representatives,n);
			//cout << "已将 632 位扩展到完整的 8192 位布尔函数\n";
			// 将 int 数组转换为 vector<int>
    		vector<int> foff(begin(child.fullGenes), end(child.fullGenes));
    		child.table=foff;
    		ph=isBalanced(child.table);
		}

    } else {
        bool ph=false;
    	while(ph==false){
    		// 均匀交叉：对子代的每个基因位随机选择父母的基因
        	for (int i = 0; i < GENE_LENGTH; i++) {
            	child.genes[i] = (rand() % 2 == 0) ? parent1.genes[i] : parent2.genes[i];
        	}
        	// 将 632 位扩展到完整的 8192 位布尔函数
    		expandToFull8192(child.genes, child.fullGenes, child.representatives,n);
			//cout << "已将 632 位扩展到完整的 8192 位布尔函数\n";
			// 将 int 数组转换为 vector<int>
    		vector<int> foff(begin(child.fullGenes), end(child.fullGenes));
    		child.table=foff;
    		ph=isBalanced(child.table);
		}

    } /*
		random_device rd;
        mt19937 gen(static_cast<unsigned int>(chrono::steady_clock::now().time_since_epoch().count()));
        uniform_int_distribution<int> dist(0, 1);  // 随机选择基因位

        BooleanFunction child(n);  // 创建子代布尔函数
        for (int i = 0; i < (1 << n); ++i) {
            // 对每个基因位，随机从父代或Bent函数中选择一个基因
            child.table[i] = (dist(gen) == 0) ? bent.table[i] : other.table[i];
        }

		// **保持平衡**
    	child.ensure_balance();
    */
        return child;
    }
    /*
    BooleanFunction crossover(const BooleanFunction &other) {
        random_device rd;
        mt19937 gen(rd());
        uniform_int_distribution<int> dist(0, table.size() - 1);
        int point = dist(gen);
        BooleanFunction child(n);
        copy(table.begin(), table.begin() + point, child.table.begin());
        copy(other.table.begin() + point, other.table.end(), child.table.begin() + point);
        return child;
    }
    */
};

// 锦标赛选择方法
BooleanFunction tournament_selection(const vector<BooleanFunction>& population, int tournament_size) {
    random_device rd;
    mt19937 gen(static_cast<unsigned int>(chrono::steady_clock::now().time_since_epoch().count()));
    uniform_int_distribution<int> dist(0, population.size() - 1);

    vector<BooleanFunction> tournament;
    for (int i = 0; i < tournament_size; ++i) {
        tournament.push_back(population[dist(gen)]);
    }

    return *min_element(tournament.begin(), tournament.end(), [](const BooleanFunction &a, const BooleanFunction &b) {
        return a.cost1() < b.cost1();  // 选择 cost1 最小的个体
    });
}

// 遗传算法（GA）
BooleanFunction genetic_algorithm(int n, int pop_size, int generations, int tournament_size, double crossover_rate, double mutation_rate ) {
    vector<BooleanFunction> population(pop_size, BooleanFunction(n));  // 初始化种群大小为 500
    //BooleanFunction bent = BooleanFunction::BentFunction(n);  // 固定的 Bent 函数
	//BooleanFunction bent = BooleanFunction::LoadFromFile("GTS91.txt", n);
	// 初始化全局最优解，设置为初始化种群中的最优个体
    BooleanFunction global_best = *min_element(population.begin(), population.end(), [](const BooleanFunction &a, const BooleanFunction &b) {
        return a.cost1() < b.cost1();  // 选择 cost1 最小的个体
    });

    int gen=0;
    // 输出当前种群的最优个体适应度
    cout << "Generation " << gen  << " - Best individual cost1: " << global_best.cost1() << endl;
    result << "Generation " << gen  << " - Best individual cost1: " << global_best.cost1() << endl;
    for (gen = 0; gen < generations; ++gen) {
        vector<BooleanFunction> new_population;

        while (new_population.size() < pop_size) {
            // 锦标赛选择，选择 5 个个体中的最优者
            BooleanFunction parent = tournament_selection(population, tournament_size);

			// 随机决定是否进行交叉
			random_device rd;
			mt19937 gen(static_cast<unsigned int>(chrono::steady_clock::now().time_since_epoch().count()));
			uniform_real_distribution<double> dist(0.0, 1.0);  // 生成 0 到 1 之间的随机数
			BooleanFunction child(n);  // 创建子代布尔函数
			if (dist(gen) <= crossover_rate) {  // 如果随机数小于等于交叉率，进行交叉
    			child = parent.crossover_with_Bent(parent, global_best);
			} else {
   				child = parent;  // 如果不进行交叉，直接保留父代
			}
            /*
			// 与 Bent 函数进行均匀交叉
            BooleanFunction child = parent.crossover_with_Bent(parent,bent);
			*/
            // 突变
            //if (random_device{}() % 1000 < mutation_rate * 1000) {
            child.mutate(mutation_rate);
            //}

            new_population.push_back(child);
        }

        population = new_population;  // 更新种群

        // 当前新生成种群中的最优个体
        BooleanFunction best_current = *min_element(population.begin(), population.end(), [](const BooleanFunction &a, const BooleanFunction &b) {
            return a.cost1() < b.cost1();  // 选择 cost1 最小的个体
        });

        // 输出当前种群的最优个体适应度
        cout << "Generation " << gen + 1 << " - Best individual cost1: " << best_current.cost1() << endl;
		result << "Generation " << gen + 1 << " - Best individual cost1: " << best_current.cost1() << endl;
        // 比较并更新全局最优解
        if (best_current.cost1() < global_best.cost1()) {
            global_best = best_current;  // 更新全局最优解
        }
    }

	// 输出当前种群的最优个体适应度
    cout << "Generation " << gen  << " - Best individual cost1: " << global_best.cost1() << endl;
    result << "Generation " << gen  << " - Best individual cost1: " << global_best.cost1() << endl;
    // 返回迭代后的最优解（全局最优个体）
    return global_best;
}

//用2-opt方法在邻域内生成候选解列表
vector<BooleanFunction> generate_2opt_neighbors(const BooleanFunction &current_solution, int num_neighbors) {
    vector<BooleanFunction> neighbors;
    int size = current_solution.table.size();
    random_device rd;
    mt19937 gen(static_cast<unsigned int>(chrono::steady_clock::now().time_since_epoch().count()));
    uniform_int_distribution<int> dist(0, size - 1);

    while (neighbors.size() < num_neighbors) {
        BooleanFunction neighbor = current_solution;
        /*
        int i = dist(gen);
        int j = dist(gen);
        if (i != j) {  // 确保交换两个不同位置的值
            swap(neighbor.table[i], neighbor.table[j]);
            neighbors.push_back(neighbor);
        }
        */
        int i, j;

        // 确保选择的两个位置 i 和 j 的值是不同的（即一个是 0，一个是 1）
        do {
            i = dist(gen);
            j = dist(gen);
        } while (i == j || neighbor.table[i] == neighbor.table[j]);

        // 交换 0 和 1
        swap(neighbor.table[i], neighbor.table[j]);

        neighbors.push_back(neighbor);
    }

    return neighbors;
}

// 修改后的禁忌搜索
BooleanFunction tabu_search(BooleanFunction initial_solution, int iterations, int candidate_list_size, int tabu_size) {
    BooleanFunction best_solution = initial_solution;
    BooleanFunction current_solution = initial_solution;
    //vector<BooleanFunction> tabu_list;
    list<BooleanFunction> tabu_list;  // 先进先出的循环禁忌表
    set<vector<int>> tabu_set;  // 额外的哈希表用于快速查找

    int best_rank = current_solution.cost2(); // 记录最大秩
    int i=0;
	// 输出当前最优个体适应度
    cout << "Generation " << i  << " - Best individual cost2: " << best_solution.cost2() << endl;
    result << "Generation " << i  << " - Best individual cost2: " << best_solution.cost2() << endl;
    // 终止条件：如果 `cost2()` 充分大，表示可以构造非奇异矩阵
        if (best_rank >= current_solution.n) {
            cout << "找到满足条件的解，提前终止禁忌搜索" << endl;
            result << "找到满足条件的解，提前终止禁忌搜索" << endl;
            // 输出当前最优个体适应度
    		cout << " Generation " << i  << " - Best individual cost2: " << best_rank << endl;
    		result << " Generation " << i  << " - Best individual cost2: " << best_rank << endl;
            return best_solution;
        }else{

    for (i = 0; i < iterations; ++i) {
       	vector<BooleanFunction> neighbors = generate_2opt_neighbors(current_solution, candidate_list_size);

		//// **Step 1: 选择 cost2() 最大的邻居 S'**
        BooleanFunction best_neighbor = *max_element(neighbors.begin(), neighbors.end(), [](const BooleanFunction &a, const BooleanFunction &b) {
            return a.cost2() < b.cost2();
        });
        /*
        // **Step 2: 先和当前解 S 进行对比**
		if (best_neighbor.cost2() >= current_solution.cost2()) {
    		current_solution = best_neighbor;
			}

		// **Step 3: 再检查 S' 是否在禁忌表**
		if (tabu_set.find(best_neighbor.table) != tabu_set.end()) {
    		if (best_neighbor.cost2() > best_rank) {  // 期望准则，允许突破禁忌
    			// 覆盖禁忌状态
				tabu_set.erase(best_neighbor.table);
                auto it = find(tabu_list.begin(), tabu_list.end(), best_neighbor);
                if (it != tabu_list.end()) {
                    tabu_list.erase(it);  // 从队列中移除旧记录
                }
                // 更新解
        		best_rank = best_neighbor.cost2();
        		best_solution = best_neighbor;
    		} else {
        		continue;  // 如果 S' 被禁忌，且不满足期望准则，则跳过
    		}
		}
		*/

		// **Step 2: 检查禁忌状态**
        bool is_tabu = (tabu_set.find(best_neighbor.table) != tabu_set.end());

		// **Step 3: 应用期望准则**
        if (is_tabu) {
            if (best_neighbor.cost2() > best_rank) {
				// 覆盖禁忌状态
				/*
				tabu_set.erase(best_neighbor.table);
                auto it = find(tabu_list.begin(), tabu_list.end(), best_neighbor);
                if (it != tabu_list.end()) {
                    tabu_list.erase(it);  // 从队列中移除旧记录
                }
                */
                tabu_list.remove_if([&](const BooleanFunction& x) {
            		return x.table == best_neighbor.table;
        		});
                tabu_set.erase(best_neighbor.table);
                // 更新解
                best_rank = best_neighbor.cost2();
                best_solution = best_neighbor;
                current_solution = best_neighbor;  // 更新当前解为 S'
                // 输出当前最优个体适应度
    			cout << "S' is tabu, Generation " << i+1  << " - Best individual cost2: " << best_solution.cost2() << endl;
    			result << "S' is tabu, Generation " << i+1  << " - Best individual cost2: " << best_solution.cost2() << endl;
            } else {
                continue;  // 跳过禁忌解
            }
        } else {
            // **Step 4: 接受候选解**
            current_solution = best_neighbor;
            if (best_neighbor.cost2() > best_rank) {
                best_rank = best_neighbor.cost2();
                best_solution = best_neighbor;
                // 输出当前最优个体适应度
    			cout << "S' is not tabu, Generation " << i+1  << " - Best individual cost2: " << best_solution.cost2() << endl;
    			result << "S' is not tabu, Generation " << i+1  << " - Best individual cost2: " << best_solution.cost2() << endl;
            }
        }

        // **Step 5: 更新禁忌表（添加候选解 S'）**
        tabu_list.push_back(best_neighbor);
        tabu_set.insert(best_neighbor.table);

        // **Step 6: 维护禁忌表大小**
        if (tabu_list.size() > tabu_size) {
            BooleanFunction oldest = tabu_list.front();
            tabu_list.pop_front();
            tabu_set.erase(oldest.table);
        }

		/*
		// **Step 1: 与当前解 S 进行比较，决定是否更新**
        if (best_neighbor.cost2() >= current_solution.cost2()) {
            current_solution = best_neighbor;  // 更新当前解 S
        }

        // **Step 2: 更新禁忌表**
        tabu_list.push_back(current_solution);
        if (tabu_list.size() > tabu_size) {
            tabu_list.erase(tabu_list.begin());  // 保持禁忌表大小不超过 tabu_size
        }

         // **Step 3: 期望准则 (Aspiration Criterion)**
        // 如果 S' 比历史最优 S* 还要好，则覆盖禁忌状态，更新 S*
        if (best_neighbor.cost2() > best_rank) {
            best_rank = best_neighbor.cost2();
            best_solution = best_neighbor;
            // **从禁忌表中移除 S'，确保它不会被误禁**
            tabu_list.erase(remove(tabu_list.begin(), tabu_list.end(), best_neighbor), tabu_list.end());
        }
       */

        // 终止条件：如果 `cost2()` 充分大，表示可以构造非奇异矩阵
        if (best_rank >= current_solution.n) {
            cout << "找到满足条件的解，提前终止禁忌搜索" << endl;
            result << "找到满足条件的解，提前终止禁忌搜索" << endl;
            // 输出当前最优个体适应度
    		cout << " Generation " << i+1  << " - Best individual cost2: " << best_rank << endl;
    		result << " Generation " << i+1  << " - Best individual cost2: " << best_rank << endl;
            break;
        }
    }

    return best_solution;}
}
/*
BooleanFunction tabu_search(BooleanFunction initial_solution, int iterations, int tabu_size = 10) {
    BooleanFunction best_solution = initial_solution;
    BooleanFunction current_solution = initial_solution;
    vector<BooleanFunction> tabu_list;

    for (int i = 0; i < iterations; ++i) {
        vector<BooleanFunction> neighbors;
        for (int j = 0; j < 5; ++j) {
            BooleanFunction neighbor = current_solution;
            neighbor.mutate(0.05);
            neighbors.push_back(neighbor);
        }

        auto best_neighbor = *max_element(neighbors.begin(), neighbors.end(), [](const BooleanFunction &a, const BooleanFunction &b) {
            return a.fitness() < b.fitness();
        });

        if (best_neighbor.fitness() > best_solution.fitness()) {
            best_solution = best_neighbor;
        }

        tabu_list.push_back(current_solution);
        if (tabu_list.size() > tabu_size) {
            tabu_list.erase(tabu_list.begin());
        }

        current_solution = best_neighbor;
    }

    return best_solution;
}
*/
// GF(2)矩阵求逆（高斯-若尔当消元法）
vector<vector<int>> gf2_invert(const vector<vector<int>>& mat) {
    int n = mat.size();
    vector<vector<int>> aug(n, vector<int>(2*n, 0));

    // 构造增广矩阵 [mat | I]
    for (int i = 0; i < n; ++i) {
        copy(mat[i].begin(), mat[i].end(), aug[i].begin());
        aug[i][n + i] = 1;
    }

    // 高斯消元
    for (int col = 0; col < n; ++col) {
        int pivot = -1;
        for (int i = col; i < n; ++i) {
            if (aug[i][col] == 1) {
                pivot = i;
                break;
            }
        }
        if (pivot == -1) return {}; // 不可逆

        swap(aug[col], aug[pivot]);

        for (int i = 0; i < n; ++i) {
            if (i != col && aug[i][col] == 1) {
                for (int j = col; j < 2*n; ++j) {
                    aug[i][j] ^= aug[col][j];
                }
            }
        }
    }

    // 提取逆矩阵
    vector<vector<int>> inv(n, vector<int>(n));
    for (int i = 0; i < n; ++i) {
        copy(aug[i].begin() + n, aug[i].end(), inv[i].begin());
    }
    return inv;
}

BooleanFunction apply_linear_transformation(const BooleanFunction &f) {
    int n = f.n;
    vector<int> W = f.compute_walsh_spectrum();

    // 收集满足 W_f(ω)=0 的向量
    vector<vector<int>> vectors;
    for (int i = 0; i < (1 << n); ++i) {
        if (W[i] == 0) {
            vector<int> vec(n);
            for (int j = 0; j < n; ++j) {
                vec[j] = (i >> j) & 1;
            }
            vectors.push_back(vec);
        }
    }

    // 构造GF(2)矩阵并选择n个线性无关行
    vector<vector<int>> B;
    vector<vector<int>> current_mat;//存储当前已选向量
    for (auto& vec : vectors) {
        vector<vector<int>> tmp_mat = current_mat;//创建临时副本
		tmp_mat.push_back(vec);
        int r = BooleanFunction::gf2_rank(tmp_mat);//对副本计算秩
        if (r > B.size()) {
            B.push_back(vec);
            current_mat.push_back(vec);//更新当前矩阵
            if (B.size() == n) break;
        }
    }
    if (B.size() < n) {
		cout<< "不足n个线性无关行构造矩阵Bf，返回原函数" <<endl;
		result<< "不足n个线性无关行构造矩阵Bf，返回原函数" <<endl;
		return f; // 无法构造可逆矩阵
	}
    // 计算逆矩阵
    vector<vector<int>> B_inv = gf2_invert(B);
    if (B_inv.empty()){
    	cout<< "矩阵不可逆，返回原函数" <<endl;
    	result<< "矩阵不可逆，返回原函数" <<endl;
    	return f;
	}

    // 应用变换：f'(X) = f(B^{-1}X)
    BooleanFunction transformed(n);
    for (int x = 0; x < (1 << n); ++x) {
        vector<int> X(n);
        for (int j = 0; j < n; ++j) {
            X[j] = (x >> j) & 1;
        }

        // 计算 B^{-1}X (GF(2)矩阵乘法)
        vector<int> X_transformed(n, 0);
        // 修改矩阵乘法部分
		for (int i = 0; i < n; ++i) {
    		X_transformed[i] = inner_product(
        		B_inv[i].begin(), B_inv[i].end(),
        		X.begin(), 0,
        		[](int a, int b) { return a ^ b; }, // 异或累加
        		[](int a, int b) { return a & b; }  // 与操作
    		);
		}
        /*
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                X_transformed[i] ^= (B_inv[i][j] & X[j]); // 异或累加
            }
        }
        */

        // 转换为真值表索引
        int index = 0;
        for (int j = 0; j < n; ++j) {
            index |= (X_transformed[j] << j);
        }
        transformed.table[x] = f.table[index];
    }
    return transformed;
}


int main() {
	// 生成一个Bent函数（只生成一次）
    //BooleanFunction bent = BooleanFunction::BentFunction(8);//n要替换为具体的
    //ofstream result("result14.0.txt");
    //if (!result) { // 检查文件是否打开成功
        //cerr << "无法打开文件!" << endl;
        //return 1;
    //}

    auto start = chrono::high_resolution_clock::now();  // 开始计时

	cout << "运行遗传算法..." << endl;
	result << "运行遗传算法..." << endl;
    int n=12;
	int pop_size=500;//种群大小
	int generations=50;//迭代次数
	int tournament_size=5;//锦标赛个体数
	double crossover_rate=0.5;//交叉概率
	double mutation_rate=0.001;//突变概率
    BooleanFunction best_ga_solution = genetic_algorithm(n, pop_size, generations, tournament_size, crossover_rate, mutation_rate);
    cout << "Best GA Boolean Function found: ";
    result << "Best GA Boolean Function found: ";
    for (int bit : best_ga_solution.table) {
        cout << bit << " ";
        result << bit << " ";
    }
    cout << endl;
    result << endl;

    cout << "运行禁忌搜索..." << endl;
    result << "运行禁忌搜索..." << endl;
    int iterations=50;//迭代次数
	int candidate_list_size=4000;//候选列表长度
	int tabu_size=30;//禁忌列表长度
    BooleanFunction best_tabu_solution = tabu_search(best_ga_solution, iterations, candidate_list_size, tabu_size);
    cout << "Best tabu Boolean Function found: ";
    result << "Best tabu Boolean Function found: ";
    for (int bit : best_tabu_solution.table) {
        cout << bit << " ";
        result << bit << " ";
    }
    cout << endl;
    result << endl;

    cout << "应用线性变换..." << endl;
    result << "应用线性变换..." << endl;
    BooleanFunction best_solution = apply_linear_transformation(best_tabu_solution);

    cout << "Best Boolean Function found: ";
    result << "Best Boolean Function found: ";
    for (int bit : best_solution.table) {
        cout << bit << " ";
        result << bit << " ";
    }
    cout << endl;
    result << endl;

    auto end = chrono::high_resolution_clock::now();    // 结束计时
    chrono::duration<double> elapsed = end - start;

    cout << "程序计算耗时: " << elapsed.count() << " 秒" << endl;
    result << "程序计算耗时: " << elapsed.count() << " 秒" << endl;
    
    return 0;
}
