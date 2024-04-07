#include "tomita.h"
#include "tool.h"
#include <algorithm>
#include <cassert>
#include <climits>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <queue>

extern int p;
unsigned long long opt_mat_calls = 0;
unsigned long long pivot_mat_calls = 0;
unsigned long long rev_mat_calls = 0;
unsigned long long tplex = 0;
unsigned long long tplex_terminate = 0;

void Tomita_t::print_maximal() {
	int i;
	for (i = 0; i < R.top; i++) printf("%d ", g.new2old[R.v_list[i]]);
	printf("\n");
}

Tomita_t::Tomita_t(const char *r_file) {
	g.read_graph(r_file);
	maximal_clique = 0;

	int i, j, u, v;

	d = (int*) calloc(g.v_size, sizeof(int));
	stat_v = (int*) calloc(g.v_size, sizeof(int));
	rank = (int*) malloc(g.v_size * sizeof(int));
	lab = (int*) calloc(g.v_size, sizeof(int));
	adj = (int**) malloc(g.v_size * sizeof(int*));

	for (i = 0; i < g.e_size; i++) {
		u = g.edges[i].u;
		v = g.edges[i].v;
		d[u]++;
		d[v]++;
	}

	for (i = 0; i < g.v_size; i++) {
		adj[i] = (int*) malloc(d[i] * sizeof(int));
		d_max = d_max > d[i] ? d_max : d[i];
		d[i] = 0;
	}

	for (i = 0; i < g.e_size; i++) {
		u = g.edges[i].u;
		v = g.edges[i].v;
		adj[u][d[u]++] = v;
		adj[v][d[v]++] = u;
	}

	g.cuhash.resize(g.v_size);
	for (i = 0; i < g.v_size; i++) {
		g.cuhash[i].reserve(d[i] + 1);
		for (int32_t j = 0; j < d[i]; ++j)
			g.cuhash[i].insert(adj[i][j]);
	}

	P = (int**) malloc((d_max + 1) * sizeof(int*));
	for (i = 0; i <= d_max; i++) P[i] = (int*) malloc(g.v_size * sizeof(int));
	P_size = (int*) malloc((d_max + 1) * sizeof(int));

	X = (int**) malloc((d_max + 1) * sizeof(int*));
	for (i = 0; i <= d_max; i++) X[i] = (int*) malloc(g.v_size * sizeof(int));
	X_size = (int*) malloc((d_max + 1) * sizeof(int));

	//    d = (int**) malloc((d_max + 1) * sizeof(int*));
	//    d[0] = d0;
	//    for (i = 1; i <= d_max; i++) d[i] = (int*) malloc(g.v_size * sizeof(int));
}

Tomita_t::~Tomita_t() {
	int i;

	if (d) {
		//        for (i = 0; i <= d_max; i++) free(d[i]);
		free(d);
		d = nullptr;
	}

	if (rank) {
		free(rank);
		rank = nullptr;
	}

	if (lab) {
		free(lab);
		lab = nullptr;
	}

	if (stat_v) {
		free(stat_v);
		stat_v = nullptr;
	}

	if (adj) {
		for (i = 0; i < g.v_size; i++) free(adj[i]);
		free(adj);
		adj = nullptr;
	}

	if (P_size) {
		free(P_size);
		P_size = nullptr;

		for (i = 0; i <= d_max; i++) free(P[i]);
		free(P);
		P = nullptr;
	}

	if (X_size) {
		free(X_size);
		X_size = nullptr;

		for (i = 0; i <= d_max; i++) free(X[i]);
		free(X);
		X = nullptr;
	}

}
//
//void Tomita_t::Tomita_pivot_mat() {
//    int i, j, u, v, w, end;
//    Heap_t heap;
//    KeyVal_t kv;
//
//    heap.make_heap(d, g.v_size);
//    for (i = 0; i < g.v_size; i++) {
//        kv = heap.pop();
//        core_num = kv.val > core_num ? kv.val : core_num;
//        rank[kv.key] = i;
//        for (j = 0; j < d[kv.key]; j++)
//            heap.update(adj[kv.key][j]);
//    }
//
//    printf("Core number = %d\n", core_num);
//
//    mat = (int**) malloc((core_num + 1) * sizeof(int*));
//    for (i = 0; i <= core_num; i++) mat[i] = (int*) calloc(g.v_size, sizeof(int));
//    v2d = (int*) malloc(g.v_size * sizeof(int));
//
//    R.make_stack(core_num + 2);
//
//    for (u = 0; u < g.v_size; u++) {
//        P_size[1] = 0;
//        X_size[1] = 0;
//
//        for (i = 0; i < d[u]; i++) {
//            v = adj[u][i];
//            if (rank[v] < rank[u]) {
//                X[1][X_size[1]++] = v;
//                lab[v] = -1;
//            }
//            else {
//                P[1][P_size[1]++] = v;
//                lab[v] = 1;
//            }
//        }
//
////        for (i = 0; i < P_size[1]; i++) {
////            v = P[1][i];
////            v2d[v] = i;
////
////            for (j = 0; j < d[v]; j++) {
////                w = adj[v][j];
////
////                if (lab[w]) mat[i][w] = 1;
////            }
////        }
//
//        if (u % 10000 == 0)
//            printf("vertex %d: |P| = %d, |X| = %d\n", u, P_size[1], X_size[1]);
//
//        R.push(u);
//
//        Tomita_pivot_mat_rec(1);
//
//        R.pop();
//
////        for (i = 0; i < P_size[1]; i++) {
////            v = P[1][i];
////
////            for (j = 0; j < d[v]; j++) {
////                w = adj[v][j];
////                mat[v2d[v]][w] = 0;
////            }
////        }
//
//        for (i = 0; i < d[u]; i++) {
//            v = adj[u][i];
//            lab[v] = 0;
//        }
//
//    }
//
//    printf("#mc = %llu\n", maximal_clique);
//
//    for (i = 0; i <= core_num; i++) free(mat[i]);
//    free(mat);
//
//    free(v2d);
//}
//
//void Tomita_t::Tomita_pivot_mat_rec(int l) {
//    pivot_mat_calls++;
//
//    if (P_size[l] == 0) {
//        if (X_size[l] == 0) {
//            maximal_clique++;
////            print_maximal();
//        }
//        return;
//    }
//
//    int i, j, u, v;
//
//    // choose pivot ---
//
////    int pivot, cur_max, nb_cnt, pivot_P_size = 0;
////    int *pivot_P = (int*) malloc(P_size[l] * sizeof(int));
////
////    for (i = 0; i < P_size[l] + X_size[l]; i++) {
////        u = i < P_size[l] ? P[l][i] : X[l][i - P_size[l]];
////
////        nb_cnt = 0;
////        for (j = 0; j < P_size[l]; j++) {
////            v = P[l][j];
////            if (mat[v2d[v]][u]) nb_cnt++;
////        }
////
////        if (i == 0) {
////            pivot = u;
////            cur_max = nb_cnt;
////        }
////
////        else if (nb_cnt >= cur_max) {
////            pivot = u;
////            cur_max = nb_cnt;
////        }
////    }
////
////    for (i = 0; i < P_size[l]; i++) {
////        u = P[l][i];
////        if (!mat[v2d[u]][pivot]) pivot_P[pivot_P_size++] = u;
////    }
//
//    // ----------------
//
//    for (i = 0; i < P_size[l]; i++) {
//        u = P[l][i];
//
//        P_size[l + 1] = 0;
//        X_size[l + 1] = 0;
//
//        for (j = 0; j < P_size[l] + X_size[l]; j++) {
//            v = j < P_size[l] ? P[l][j] : X[l][j - P_size[l]];
//
//            if (mat[v2d[u]][v]) {
//
//                if (lab[v] == l) {
//                    P[l + 1][P_size[l + 1]++] = v;
//                    lab[v] = l + 1;
//                }
//                else if (lab[v] == -l) {
//                    X[l + 1][X_size[l + 1]++] = v;
//                    lab[v] = -(l + 1);
//                }
//            }
//        }
//
//        R.push(u);
//
//        Tomita_pivot_mat_rec(l + 1);
//
//        R.pop();
//
//        for (j = 0; j < P_size[l + 1]; j++) {
//            v = P[l + 1][j];
//            lab[v] = l;
//        }
//
//        for (j = 0; j < X_size[l + 1]; j++) {
//            v = X[l + 1][j];
//            lab[v] = -l;
//        }
//
//        lab[u] = -l;
//    }
//
//    for (i = 0; i < P_size[l]; i++) {
//        u = P[l][i];
//        lab[u] = l;
//    }
//
////    free(pivot_P);
//}

void Tomita_t::Tomita_opt_mat() {
	int i, j, u, v, w, end;
	Heap_t heap;
	KeyVal_t kv;

	heap.make_heap(d, g.v_size);
	for (i = 0; i < g.v_size; i++) {
		kv = heap.pop();
		core_num = kv.val > core_num ? kv.val : core_num;
		rank[kv.key] = i;
		for (j = 0; j < d[kv.key]; j++)
			heap.update(adj[kv.key][j]);
	}

	printf("Core number = %d\n", core_num);

	//    mat = (int**) malloc((core_num + 1) * sizeof(int*));
	//    for (i = 0; i <= core_num; i++) mat[i] = (int*) calloc(g.v_size, sizeof(int));
	//    v2d = (int*) malloc(g.v_size * sizeof(int));

	R.make_stack(core_num + 2);

	for (u = 0; u < g.v_size; u++) {
		P_size[1] = 0;
		X_size[1] = 0;

		for (i = 0; i < d[u]; i++) {
			v = adj[u][i];
			if (rank[v] < rank[u]) {
				X[1][X_size[1]++] = v;
				lab[v] = -1;
			}
			else {
				P[1][P_size[1]++] = v;
				lab[v] = 1;
			}
		}

		//        for (i = 0; i < P_size[1]; i++) {
		//            v = P[1][i];
		//            v2d[v] = i;
		//
		//            end = d[v];
		//            for (j = 0; j < end; j++) {
		//                w = adj[v][j];
		//
		//                if (lab[w]) mat[i][w] = 1;
		//            }
		//        }

		if (u % 10000 == 0)
			printf("vertex %d: |P| = %d, |X| = %d\n", u, P_size[1], X_size[1]);


		R.push(u);

		Tomita_opt_mat_rec(1);

		R.pop();


		//        for (i = 0; i < P_size[1]; i++) {
		//            v = P[1][i];
		//
		//            for (j = 0; j < d[v]; j++) {
		//                w = adj[v][j];
		//                mat[v2d[v]][w] = 0;
		//            }
		//        }

		for (i = 0; i < d[u]; i++) {
			v = adj[u][i];
			lab[v] = 0;
		}
	}

	printf("#mc = %llu\n", maximal_clique);

	//    for (i = 0; i <= core_num; i++) free(mat[i]);
	//    free(mat);
	//
	//    free(v2d);
}

void Tomita_t::Tomita_opt_mat_rec(int l) {
	opt_mat_calls++;
	if (P_size[l] == 0) {
		if (X_size[l] == 0) {
			maximal_clique++;
			//            print_maximal();
		}
		return;
	}

	int i, j, k, u, v, w, end;

	// choose pivots ---

	int pivot, cur_max, nb_cnt, opt_P_size = 0, cur_min = INT32_MAX, LR = 0;

	for (i = 0; i < X_size[l] + P_size[l]; i++) {

		u = i < X_size[l] ? X[l][i] : P[l][i - X_size[l]];

		nb_cnt = 0;

		if (d[u] > P_size[l]) {
			for (j = 0; j < P_size[l]; j++) {
				v = P[l][j];
				if (g.connect(u, v)) nb_cnt++;
			}
		}
		else {
			for (j = 0; j < d[u]; j++) {
				v = adj[u][j];
				if (lab[v] == l) nb_cnt++;
			}
		}

		if (i == 0) {
			pivot = u;
			cur_max = nb_cnt;
		}

		else if (nb_cnt > cur_max) {
			pivot = u;
			cur_max = nb_cnt;
		}

		if (i >= X_size[l]) { 

			if (nb_cnt < cur_min) cur_min = nb_cnt;

			if (nb_cnt == P_size[l] - 2) LR++;
		}
	}

	if (p > 0 && cur_min >= P_size[l] - p && X_size[l] == 0) {

		if (p == 1) maximal_clique++;

		else if (p == 2) maximal_clique += (1 << (LR >> 1));

		//        else if (p == 3) list_in_3plex(l);
		//
		//        else list_in_plex(l);

		return;
	}

	int* opt_P = (int*) malloc(P_size[l] * sizeof(int));

	for (i = 0; i < P_size[l]; i++) {
		u = P[l][i];
		if (!g.connect(pivot, u)) opt_P[opt_P_size++] = u;
	}

	// ----------------

	for (i = 0; i < opt_P_size; i++) {
		u = opt_P[i];

		P_size[l + 1] = 0;
		X_size[l + 1] = 0;


		if (d[u] > P_size[l] + X_size[l]) {
			for (j = 0; j < P_size[l] + X_size[l]; j++) {
				v = j < P_size[l] ? P[l][j] : X[l][j - P_size[l]];

				if (g.connect(u, v)) {

					if (lab[v] == l) {
						P[l + 1][P_size[l + 1]++] = v;
						lab[v] = l + 1;
					}
					else if (lab[v] == -l) {
						X[l + 1][X_size[l + 1]++] = v;
						lab[v] = -(l + 1);
					}
				}
			}
		}

		else {
			for (j = 0; j < d[u]; j++) {
				v = adj[u][j];
				if (lab[v] == l) {
					P[l + 1][P_size[l + 1]++] = v;
					lab[v] = l + 1;
				}
				else if (lab[v] == -l) {
					X[l + 1][X_size[l + 1]++] = v;
					lab[v] = -(l + 1);
				}

			}
		}


		//        if (d[u] > P_size[l]) {
		//            for (j = 0; j < P_size[l]; j++) {
		//                v = P[l][j];
		//                if (lab[v] == l && g.connect(u, v)) {
		//                    P[l + 1][P_size[l + 1]++] = v;
		//                    lab[v] = l + 1;
		//                }
		//            }
		//        }
		//        else {
		//            for (j = 0; j < d[u]; j++) {
		//                v = adj[u][j];
		//                if (lab[v] == l) {
		//                    P[l + 1][P_size[l + 1]++] = v;
		//                    lab[v] = l + 1;
		//                }
		//            }
		//        }

		//        if (d[u] > X_size[l]) {
		//            for (j = 0; j < X_size[l]; j++) {
		//                v = X[l][j];
		//                if (g.connect(u, v)) {
		//                    X[l + 1][X_size[l + 1]++] = v;
		//                    lab[v] = -(l + 1);
		//                }
		//            }
		//        }
		//        else {
		//            for (j = 0; j < d[u]; j++) {
		//                v = adj[u][j];
		//                if (lab[v] == -l) {
		//                    X[l + 1][X_size[l + 1]++] = v;
		//                    lab[v] = -(l + 1);
		//                }
		//            }
		//        }


		R.push(u);

		Tomita_opt_mat_rec(l + 1);

		R.pop();

		for (j = 0; j < P_size[l + 1]; j++) {
			v = P[l + 1][j];
			lab[v] = l;
		}

		for (j = 0; j < X_size[l + 1]; j++) {
			v = X[l + 1][j];
			lab[v] = -l;
		}

		lab[u] = -l;
	}

	for (i = 0; i < opt_P_size; i++) {
		u = opt_P[i];
		lab[u] = l;
	}

	free(opt_P);
}

//void Tomita_t::list_in_1plex(int l) {
//    int i, j, u, v, flg;
//
//    for (i = 0; i < X_size[l]; i++) {
//        u = X[l][i];
//
//        flg = 1;
//        for (j = 0; j < P_size[l]; j++) {
//            v = P[l][j];
//
//            flg &= mat[v2d[v]][u];
//        }
//        if (flg) return;
//    }
//
//    maximal_clique++;
//}
//
//void Tomita_t::list_in_2plex(int l) {
//    int i, j, k, u, v, F_size = 0, I_size = 0, flg, nums, is_maximal;
//    vector<int> I_l(P_size[l]), I_r(P_size[l]), cliques(P_size[l]);
//    vector<bool> vis(P_size[l], false);
//
//    for (i = 0; i < P_size[l]; i++) {
//
//        if (!vis[i]) {
//            u = P[l][i];
//            flg = 1;
//
//            for (j = i + 1; j < P_size[l]; j++) {
//                v = P[l][j];
//
//                if (!mat[v2d[u]][v]) {
//                    I_l[I_size] = u;
//                    I_r[I_size++] = v;
//                    flg = 0;
//                    vis[j] = true;
//                    break;
//                }
//            }
//            if (flg) cliques[F_size++] = u;
//        }
//    }
//
//    for (i = 0; i < (1 << I_size); i++) {
//
//        for (j = 0; j < I_size; j++) {
//
//            if ((i >> j) & 1) cliques[F_size + j] = I_l[j];
//            else cliques[F_size + j] = I_r[j];
//        }
//
//        is_maximal = 1;
//
//        for (j = 0; j < X_size[l]; j++) {
//            u = X[l][j];
//
//            flg = 1;
//            for (k = 0; k < F_size + I_size; k++) {
//                v = cliques[k];
//
//                flg &= mat[v2d[v]][u];
//            }
//            if (flg) {
//                is_maximal = 0;
//                break;
//            }
//        }
//
//        if (is_maximal)
//            maximal_clique++;
//    }
//
//}
//
//void Tomita_t::list_in_3plex(int l) {
//    int i, j, k, idx, u, v, w, F_size = 0, I_size = 0, flg, nums, is_maximal;
//    vector<int> cliques(P_size[l]), r_d(P_size[l], 0);
//    vector<int> I(P_size[l]);
//    vector<vector<int>> r_adj(P_size[l], vector<int>(2));
//    vector<bool> vis(P_size[l], false);
//    vector<vector<int>> chains;
//    queue<int> q;
//
//
//
//
//}
//
//void Tomita_t::list_in_plex(int l) {
//
//}

//void Tomita_t::Tomit_rev_mat() {
//    int i, j, u, v, w, end;
//    Heap_t heap;
//    KeyVal_t kv;
//
//    heap.make_heap(d, g.v_size);
//    for (i = 0; i < g.v_size; i++) {
//        kv = heap.pop();
//        core_num = kv.val > core_num ? kv.val : core_num;
//        rank[kv.key] = i;
//        for (j = 0; j < d[kv.key]; j++)
//            heap.update(adj[kv.key][j]);
//    }
//
//    printf("Core number = %d\n", core_num);
//
//    R.make_stack(core_num + 2);
//
//    for (u = 0; u < g.v_size; u++) {
//        P_size[1] = 0;
//        X_size[1] = 0;
//
//        for (i = 0; i < d[u]; i++) {
//            v = adj[u][i];
//            if (rank[v] < rank[u]) {
//                X[1][X_size[1]++] = v;
//                lab[v] = -1;
//            } else {
//                P[1][P_size[1]++] = v;
//                lab[v] = 1;
//            }
//        }
//
//        if (u % 10000 == 0)
//            printf("vertex %d: |P| = %d, |X| = %d\n", u, P_size[1], X_size[1]);
//
//
//        R.push(u);
//
//        Tomita_rev_mat_rec(1);
//
//        R.pop();
//
//        for (i = 0; i < d[u]; i++) {
//            v = adj[u][i];
//            lab[v] = 0;
//        }
//    }
//
//    printf("#mc = %llu\n", maximal_clique);
//
//}
//
//void Tomita_t::Tomita_rev_mat_rec(int l) {
//    rev_mat_calls++;
//    if (P_size[l] == 0) {
//        if (X_size[l] == 0) {
//            maximal_clique++;
////            print_maximal();
//        }
//        return;
//    }
//
//    int i, j, k, u, v, w, end, pivot, nb_cnt, cur_min, sizeP = P_size[l], sizeX = X_size[l];
//
//    while (true) {
//
//        for (i = 0; i < P_size[l]; i++) {
//            u = P[l][i];
//
//            nb_cnt = 0;
//            for (j = 0; j < P_size[l]; j++) {
//                v = P[l][j];
//                if (g.connect(u, v)) nb_cnt++;
//            }
//
//            if (i == 0) {
//                pivot = u;
//                k = i;
//                cur_min = nb_cnt;
//            }
//
//            else if (nb_cnt < cur_min) {
//                pivot = u;
//                k = i;
//                cur_min = nb_cnt;
//            }
//        }
//
//        if (cur_min >= P_size[l] - p) break;
//
//        P_size[l + 1] = 0;
//        X_size[l + 1] = 0;
//
//        for (j = 0; j < P_size[l] + X_size[l]; j++) {
//            v = j < P_size[l] ? P[l][j] : X[l][j - P_size[l]];
//
//            if (g.connect(u, v)) {
//
//                if (lab[v] == l) {
//                    P[l + 1][P_size[l + 1]++] = v;
//                    lab[v] = l + 1;
//                }
//                else if (lab[v] == -l) {
//                    X[l + 1][X_size[l + 1]++] = v;
//                    lab[v] = -(l + 1);
//                }
//            }
//        }
//
//        R.push(pivot);
//
//        Tomita_rev_mat_rec(l + 1);
//
//        R.pop();
//
//        for (j = 0; j < P_size[l + 1]; j++) {
//            v = P[l + 1][j];
//            lab[v] = l;
//        }
//
//        for (j = 0; j < X_size[l + 1]; j++) {
//            v = X[l + 1][j];
//            lab[v] = -l;
//        }
//
//        P[l][k] = P[l][--P_size[l]];
//        P[l][P_size[l]] = pivot;
//        X[l][X_size[l]++] = pivot;
//        lab[pivot] = -l;
//    }
//
//
//    if (P_size[l] > 0) {
//
//        if (p == 1) list_in_1plex(l);
//
//        else if (p == 2) list_in_2plex(l);
//
//        else if (p == 3) list_in_3plex(l);
//
//        else list_in_plex(l);
//    }
//
//    P_size[l] = sizeP;
//    X_size[l] = sizeX;
//
//    for (i - 0; i < P_size[l]; i++) {
//        u = P[l][i];
//        lab[u] = l;
//    }
//
//}
