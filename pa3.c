#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/*
 * pa3.c
 * Karp's algorithm to compute minimum cycle mean and return one minimum-mean cycle.
 * Usage: ./pa3 in_file out_file1 out_file2
 * out_file1: binary float with minimum cycle mean
 * out_file2: text file with cycle nodes as: v0 v_{-1} v_{-2} ... (no trailing space)
 */

typedef struct {
    int to;
    int from;
    double w;
    int next; // for linked list of incoming edges
} Edge;

int main(int argc, char **argv) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s in_file out_file1 out_file2\n", argv[0]);
        return EXIT_FAILURE;
    }
    const char *infile = argv[1];
    const char *out1 = argv[2];
    const char *out2 = argv[3];

    FILE *f = fopen(infile, "r");
    if (!f) { perror(infile); return EXIT_FAILURE; }

    int n = 0, m = 0;
    if (fscanf(f, "V %d\n", &n) != 1) { fclose(f); fprintf(stderr, "bad V line\n"); return EXIT_FAILURE; }
    if (fscanf(f, "E %d\n", &m) != 1) { fclose(f); fprintf(stderr, "bad E line\n"); return EXIT_FAILURE; }

    // incoming adjacency via edge list
    int *head = calloc(n+1, sizeof(int));
    if (!head) { fclose(f); fprintf(stderr, "memory failure\n"); return EXIT_FAILURE; }
    for (int i=0;i<=n;i++) head[i] = -1;
    Edge *edges = calloc(m, sizeof(Edge));
    if (!edges) { free(head); fclose(f); fprintf(stderr, "memory failure\n"); return EXIT_FAILURE; }

    for (int i = 0; i < m; ++i) {
        int dest, src; double w;
        if (fscanf(f, "%d %d %lf\n", &dest, &src, &w) != 3) {
            // try without newline
            fseek(f, 0, SEEK_SET);
            break;
        }
        edges[i].to = dest;
        edges[i].from = src;
        edges[i].w = w;
        edges[i].next = head[dest];
        head[dest] = i;
    }
    fclose(f);

    // Karp DP: D[k][v], k = 0..n
    // allocate as flat arrays of size (n+1)*(n+1)
    long long rows = (long long)(n + 1);
    long long total = rows * rows;
    // guard against crazy sizes
    if (n < 1) { free(head); free(edges); fprintf(stderr, "empty graph\n"); return EXIT_FAILURE; }

    double *D = malloc(sizeof(double) * total);
    int *pred = malloc(sizeof(int) * total);
    if (!D || !pred) { free(head); free(edges); free(D); free(pred); fprintf(stderr, "memory allocation failure\n"); return EXIT_FAILURE; }

    const double INF = INFINITY;
    // initialize
    for (long long k = 0; k <= n; ++k) {
        for (int v = 0; v <= n; ++v) {
            D[k*rows + v] = INF;
            pred[k*rows + v] = -1;
        }
    }
    int s = 1; // arbitrary source 1
    D[0*rows + s] = 0.0;

    // DP iterations
    for (int k = 1; k <= n; ++k) {
        for (int v = 1; v <= n; ++v) {
            double best = INF;
            int best_u = -1;
            // for each incoming edge u->v
            for (int ei = head[v]; ei != -1; ei = edges[ei].next) {
                int u = edges[ei].from;
                double prev = D[(k-1)*rows + u];
                if (!isfinite(prev)) continue;
                double val = prev + edges[ei].w;
                if (val < best) { best = val; best_u = u; }
            }
            D[k*rows + v] = best;
            pred[k*rows + v] = best_u;
        }
    }

    // compute lambda* via Karp formula
    double lambda = INFINITY;
    int vstar = -1;
    for (int v = 1; v <= n; ++v) {
        double Dn = D[n*rows + v];
        if (!isfinite(Dn)) continue;
        double maxval = -INFINITY;
        int any = 0;
        for (int k = 0; k <= n-1; ++k) {
            double Dk = D[k*rows + v];
            if (!isfinite(Dk)) continue;
            any = 1;
            double val = (Dn - Dk) / (double)(n - k);
            if (val > maxval) maxval = val;
        }
        if (!any) continue;
        if (maxval < lambda) { lambda = maxval; vstar = v; }
    }

    if (vstar == -1 || !isfinite(lambda)) {
        // failed to compute
        free(head); free(edges); free(D); free(pred);
        fprintf(stderr, "failed to compute minimum cycle mean\n");
        return EXIT_FAILURE;
    }

    // Find a minimum-mean cycle in the modified graph (weights w - lambda).
    // We'll search for a zero-weight cycle using iterative DFS on the modified weights.
    int *out_head = calloc(n+1, sizeof(int));
    int *out_next = malloc(sizeof(int) * m);
    if (!out_head || !out_next) { free(head); free(edges); free(D); free(pred); free(out_head); free(out_next); fprintf(stderr, "memory failure\n"); return EXIT_FAILURE; }
    for (int i = 0; i <= n; ++i) out_head[i] = -1;
    for (int i = 0; i < m; ++i) { out_next[i] = out_head[edges[i].from]; out_head[edges[i].from] = i; }

    int *color = calloc(n+1, sizeof(int)); // 0 unvisited, 1 in stack, 2 done
    double *cum = malloc(sizeof(double) * (n+1));
    int *parent = malloc(sizeof(int) * (n+1));
    int *cycle = malloc(sizeof(int) * (n+1));
    if (!color || !cum || !parent || !cycle) { free(head); free(edges); free(D); free(pred); free(out_head); free(out_next); free(color); free(cum); free(parent); free(cycle); fprintf(stderr, "memory failure\n"); return EXIT_FAILURE; }

    for (int i = 1; i <= n; ++i) { color[i] = 0; parent[i] = -1; }
    int found = 0; int clen = 0;
    // tolerance for zero-weight cycle detection: scale with |lambda|
    // use an absolute tolerance that grows with lambda to account for numerical error
    const double TOL = fmax(1e-12, 1e-8 * (1.0 + fabs(lambda)));

    int *stack = malloc(sizeof(int) * (n+5));
    int *iter = malloc(sizeof(int) * (n+5));
    if (!stack || !iter) { free(head); free(edges); free(D); free(pred); free(out_head); free(out_next); free(color); free(cum); free(parent); free(cycle); free(stack); free(iter); fprintf(stderr, "memory failure\n"); return EXIT_FAILURE; }

    for (int s_n = 1; s_n <= n && !found; ++s_n) {
        if (color[s_n] != 0) continue;
        int sp = 0;
        stack[sp] = s_n; iter[sp] = out_head[s_n]; color[s_n] = 1; cum[s_n] = 0.0; parent[s_n] = -1; sp++;
        while (sp > 0 && !found) {
            int u = stack[sp-1]; int ei = iter[sp-1];
            if (ei != -1) {
                iter[sp-1] = out_next[ei];
                int v = edges[ei].to;
                double wprime = edges[ei].w - lambda;
                if (color[v] == 0) {
                    parent[v] = u;
                    cum[v] = cum[u] + wprime;
                    color[v] = 1;
                    stack[sp] = v; iter[sp] = out_head[v]; sp++;
                } else if (color[v] == 1) {
                    double cyc_w = cum[u] + wprime - cum[v];
                    if (fabs(cyc_w) <= TOL) {
                        clen = 0;
                        cycle[clen++] = v;
                        int cur = u;
                        while (cur != v && cur != -1 && clen <= n) {
                            cycle[clen++] = cur;
                            cur = parent[cur];
                        }
                        found = 1;
                        break;
                    }
                }
            } else {
                color[u] = 2;
                sp--;
            }
        }
    }

    free(out_head); free(out_next); free(color); free(cum); free(parent); free(stack); free(iter);

    // If DFS did not find a zero-weight cycle (due to numerical issues),
    // fall back to predecessor-based reconstruction from Karp DP tables.
    if (!found) {
        // walk n steps from vstar using pred[n][.] to reach a node on a cycle
        int y2 = vstar;
        for (int t = 0; t < n; ++t) {
            int p = pred[n*rows + y2];
            if (p == -1) break;
            y2 = p;
        }
        if (y2 != -1) {
            // collect cycle by following predecessors until repetition
            int *mark = calloc(n+1, sizeof(int));
            if (mark) {
                int cur = y2;
                while (cur != -1 && !mark[cur]) { mark[cur] = 1; cur = pred[n*rows + cur]; }
                if (cur != -1) {
                    // cur is in cycle; collect cycle nodes
                    int idx = 0; int startc = cur;
                    do {
                        cycle[idx++] = cur;
                        cur = pred[n*rows + cur];
                        if (idx > n) break;
                    } while (cur != startc && cur != -1);
                    clen = idx;
                    found = (clen > 0);
                }
                free(mark);
            }
        }
    }

    // write binary float lambda to out1
    FILE *fout1 = fopen(out1, "wb");
    if (!fout1) { perror(out1); free(head); free(edges); free(D); free(pred); free(cycle); return EXIT_FAILURE; }
    float lambda_f = (float)lambda;
    fwrite(&lambda_f, sizeof(float), 1, fout1);
    fclose(fout1);

    // write cycle to out2 in required format: dest then source order by following predecessors
    FILE *fout2 = fopen(out2, "w");
    if (!fout2) { perror(out2); free(head); free(edges); free(D); free(pred); free(cycle); return EXIT_FAILURE; }
    for (int i = 0; i < clen; ++i) {
        if (i) fprintf(fout2, " ");
        fprintf(fout2, "%d", cycle[i]);
    }
    fprintf(fout2, "\n");
    fclose(fout2);

    free(cycle);
    free(head);
    free(edges);
    free(D);
    free(pred);

    return EXIT_SUCCESS;
}
