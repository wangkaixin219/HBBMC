#include <iostream>

#include "tomita.h"

extern unsigned long long opt_mat_calls, pivot_mat_calls, rev_mat_calls, tplex, tplex_terminate;
int p;
//extern int pivot_pivot_num, pivot_non_pivot_num, opt_pivot_num, opt_non_pivot_num;

int main(int argc, char* argv[]) {

    Tomita_t alg(argv[1]);
    Stat_t res;

    p = atoi(argv[2]);

    int type = atoi(argv[3]);

//    res.setup();
//    pivot.Tomita_pivot();
//    res.collapse();
//    printf("(pivot) runtime = %.3lf ms (%d calls)\n\n", res.runtime, pivot_calls);
//    res.setup();
//    opt.Tomita_opt();
//    res.collapse();
//    printf("(opt) runtime = %.3lf ms (%d calls)\n\n", res.runtime, opt_calls);

    if (type == 0) {
        res.setup();
        alg.Tomita_opt_mat();
        res.collapse();
        printf("(opt_mat) runtime = %.3lf ms (%lld calls)\n", res.runtime, opt_mat_calls);
        printf("(terminate / all) = (%lld / %lld) = %.3lf%%\n\n", tplex_terminate, tplex, (double) tplex_terminate / tplex * 100);
    }

//    else if (type == 1) {
//        res.setup();
//        alg.Tomita_pivot_mat();
//        res.collapse();
//        printf("(pivot_mat) runtime = %.3lf ms (%lld calls)\n", res.runtime, pivot_mat_calls);
//    }
//
//    else if (type == 2) {
//        res.setup();
//        alg.Tomit_rev_mat();
//        res.collapse();
//        printf("(rev_mat) runtime = %.3lf ms (%lld calls)\n", res.runtime, rev_mat_calls);
//    }
//
//    else {
//        Graph_t g;
//        int n = atoi(argv[4]), rho = atoi(argv[5]);
//        g.er_model(n, rho);
//    }

    return 0;
}
