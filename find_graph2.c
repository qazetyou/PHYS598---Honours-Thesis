#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#define MAXN 20
#include "../nauty.h"
#define MAXGRAPHS 100000
#define WORDSIZE 64

// A struct to store a discovered local complement plus the used path.
typedef struct {
    graph *gPerm;
    int pathPerm[MAXN];
    int depthPerm;
    int pathGraph[MAXN];
    int depthGraph;
} LCResult;

/*--------------------------------------------------------------------------------
 * are_isomorphic:
 *   Check if two graphs are isomorphic by comparing their canonical forms.
 *
 * Parameters:
 *   g1, g2 - The two graphs to compare.
 *   m      - The number of words per row in the graph's adjacency matrix.
 *   n      - The number of vertices in the graph.
 *
 * Returns:
 *   1 if the graphs are isomorphic, 0 otherwise.
 *-------------------------------------------------------------------------------*/
int are_isomorphic(graph *g1, graph *g2, int m, int n)
{
    int *lab1 = malloc(n * sizeof(int));
    int *ptn1 = malloc(n * sizeof(int));
    int *orbits1 = malloc(n * sizeof(int));
    int *lab2 = malloc(n * sizeof(int));
    int *ptn2 = malloc(n * sizeof(int));
    int *orbits2 = malloc(n * sizeof(int));
    static DEFAULTOPTIONS_GRAPH(options);
    statsblk stats;

    options.getcanon = TRUE;

    graph *canon_g1 = malloc(sizeof(graph) * m * n);
    graph *canon_g2 = malloc(sizeof(graph) * m * n);

    densenauty(g1, lab1, ptn1, orbits1, &options, &stats, m, n, canon_g1);
    densenauty(g2, lab2, ptn2, orbits2, &options, &stats, m, n, canon_g2);

    int result = memcmp(canon_g1, canon_g2, sizeof(graph) * m * n) == 0;

    free(lab1);
    free(ptn1);
    free(orbits1);
    free(lab2);
    free(ptn2);
    free(orbits2);
    free(canon_g1);
    free(canon_g2);
    return result;
}

/*--------------------------------------------------------------------------------
 * graph_seen:
 *   Check if the graph 'g' is isomorphic to any graph already stored in 'graphs'.
 *
 * Parameters:
 *   graphs - Array of pointers to previously discovered graphs.
 *   g      - The graph to check.
 *   m      - The number of words per row in the graph's adjacency matrix.
 *   n      - The number of vertices in the graph.
 *   count  - The number of graphs currently stored in 'graphs'.
 *
 * Returns:
 *   1 if an isomorphic graph is found in 'graphs', 0 otherwise.
 *-------------------------------------------------------------------------------*/
int graph_seen(graph *graphs[], graph *g, int m, int n, int count)
{
    for (int k = 0; k < count; ++k) {
        if (are_isomorphic(g, graphs[k], m, n)) {
            return 1; // Found an isomorphic match
        }
    }
    return 0; // No match found
}

/*--------------------------------------------------------------------------------
 * apply_permutation:
 *   Builds a permuted copy of the graph 'base' in 'dest' using the permutation
 *   array 'p'. For each edge (i, j) in 'base', an edge (p[i], p[j]) is added
 *   to 'dest'.
 *
 * Parameters:
 *   base - The original graph to permute.
 *   m    - The number of words per row in the graph's adjacency matrix.
 *   n    - The number of vertices in the graph.
 *   p    - The permutation array specifying the new vertex order.
 *
 * Returns:
 *   A pointer to the newly allocated permuted graph.
 *-------------------------------------------------------------------------------*/
graph *apply_permutation(graph *base, int m, int n, const int *p)
{
    size_t gSize = (size_t)m * n * sizeof(graph);

    graph *dest = malloc(gSize);
    if (!dest) {
        fprintf(stderr, "Memory allocation failed in apply_permutation()\n");
        exit(EXIT_FAILURE);
    }

    EMPTYGRAPH(dest, m, n);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (ISELEMENT(GRAPHROW(base, i, m), j)) {
                ADDELEMENT(GRAPHROW(dest, p[i], m), p[j]);
                ADDELEMENT(GRAPHROW(dest, p[j], m), p[i]);
            }
        }
    }
    return dest;
}

/*--------------------------------------------------------------------------------
 * are_equal:
 *   Check if two graphs have identical adjacency matrices.
 *
 * Parameters:
 *   g1, g2 - The two graphs to compare.
 *   m      - The number of words per row in the graph's adjacency matrix.
 *   n      - The number of vertices in the graph.
 *
 * Returns:
 *   1 if the graphs are identical, 0 otherwise.
 *-------------------------------------------------------------------------------*/
int are_equal(graph *g1, graph *g2, int m, int n)
{
    size_t gSize = (size_t)m * n * sizeof(graph);
    return (memcmp(g1, g2, gSize) == 0);
}

/*--------------------------------------------------------------------------------
 * local_complement:
 *   Toggles the edges among all pairs of neighbors of vertex 'v' in the graph 'g'.
 *   This operation is known as the local complement of 'v'.
 *
 * Parameters:
 *   g - The original graph.
 *   v - The vertex on which to perform the local complement.
 *   m - The number of words per row in the graph's adjacency matrix.
 *   n - The number of vertices in the graph.
 *
 * Returns:
 *   A pointer to the newly allocated graph with the local complement applied.
 *-------------------------------------------------------------------------------*/
graph *local_complement(graph *g, int v, int m, int n)
{
    size_t gSize = (size_t)m * n * sizeof(graph);

    graph *g_new = malloc(gSize);
    if (!g_new) {
        fprintf(stderr, "Memory allocation failed in local_complement()\n");
        exit(EXIT_FAILURE);
    }
    memcpy(g_new, g, gSize);

    for (int i = 0; i < n; ++i) {
        if (i == v || !ISELEMENT(GRAPHROW(g, v, m), i)) {
            continue;
        }
        for (int j = i + 1; j < n; ++j) {
            if (j == v || !ISELEMENT(GRAPHROW(g, v, m), j)) {
                continue;
            }
            if (ISELEMENT(GRAPHROW(g, i, m), j)) {
                DELELEMENT(GRAPHROW(g_new, i, m), j);
                DELELEMENT(GRAPHROW(g_new, j, m), i);
            } else {
                ADDELEMENT(GRAPHROW(g_new, i, m), j);
                ADDELEMENT(GRAPHROW(g_new, j, m), i);
            }
        }
    }
    return g_new;
}

/*--------------------------------------------------------------------------------
 * satisfiesAllPerms:
 *   Check if the graph 'Gprime' is invariant under all user-provided permutations.
 *   This means that applying each permutation to 'Gprime' results in a graph
 *   identical to 'Gprime'.
 *
 * Parameters:
 *   Gprime - The graph to check.
 *   perms  - Array of user-provided permutations.
 *   pcount - The number of permutations in 'perms'.
 *   m      - The number of words per row in the graph's adjacency matrix.
 *   n      - The number of vertices in the graph.
 *
 * Returns:
 *   1 if 'Gprime' is invariant under all permutations, 0 otherwise.
 *-------------------------------------------------------------------------------*/
int satisfiesAllPerms(graph *Gprime, int **perms, int pcount, int m, int n)
{
    for (int i = 0; i < pcount; i++) {
        graph *Gprime_p = apply_permutation(Gprime, m, n, perms[i]);
        if (!are_equal(Gprime, Gprime_p, m, n)) {
            free(Gprime_p);
            return 0;
        }
        free(Gprime_p);
    }
    return 1;
}

/*--------------------------------------------------------------------------------
 * print_graph:
 *   Print all edges of the graph 'g' in the format "u v", where 'u' and 'v'
 *   are the vertices connected by an edge.
 *
 * Parameters:
 *   g - The graph to print.
 *   m - The number of words per row in the graph's adjacency matrix.
 *   n - The number of vertices in the graph.
 *-------------------------------------------------------------------------------*/
void print_graph(graph *g, int m, int n)
{
    printf("Graph edges:\n");
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            if (ISELEMENT(GRAPHROW(g, i, m), j)) {
                printf("%d %d\n", i, j);
            }
        }
    }
}

/*--------------------------------------------------------------------------------
 * explore:
 *   Perform a depth-first search to find a local complement that satisfies
 *   all permutations, and optionally matches a comparison graph.
 *-------------------------------------------------------------------------------*/
int explore(graph *g, int m, int n,
    graph *graphs[], int *count,
    int *path, int depth, int maxDepth,
    LCResult *result,
    int **perms, int pcount, graph *gCompare,
    int permMatched, int graphMatched)
{
    // Determine how many bytes each graph needs
    size_t gSize = (size_t)m * n * sizeof(graph);

    /*
     * 1. Check if 'g' satisfies all permutations.
     *    If it does, store it in 'graphs' (if not the base graph)
     *    and record it in the result structure.
     */
    if (satisfiesAllPerms(g, perms, pcount, m, n)) {
        
        //*************************************************************************
        //print_graph(g, m, n);
        //*************************************************************************

        // If we haven't seen the graph yet,
        // store a copy in graphs[].
        if (!graph_seen(graphs, g, m, n, *count)) {
            graph *g_copy = malloc(gSize);
            if (!g_copy) {
                fprintf(stderr, "Memory allocation failed in explore()\n");
                exit(EXIT_FAILURE);
            }
            memcpy(g_copy, g, gSize);
            graphs[(*count)++] = g_copy;
        }

        // Record this as a graph that satisfies permutations
        result->gPerm = malloc(gSize);
        if (!result->gPerm) {
            fprintf(stderr, "Memory allocation failed for result->gPerm\n");
            exit(EXIT_FAILURE);
        }
        memcpy(result->gPerm, g, gSize);

        // Also record the vertex path (the sequence of toggles) used
        memcpy(result->pathPerm, path, depth * sizeof(int));
        result->depthPerm = depth;

        // Set the bit indicating "permutation match" is found
        permMatched = 1;
    }

    /*
     * 2. If gCompare is provided and g == gCompare (identical adjacency),
     *    record that we found a match for the comparison graph.
     */
    if (gCompare && are_equal(g, gCompare, m, n)) {
        memcpy(result->pathGraph, path, depth * sizeof(int));
        result->depthGraph = depth;
        graphMatched = 1;
    }

    /*
     * 3. If we haven't reached the maximum depth, perform local-complement
     *    toggles on each valid vertex and recurse.
     */
    if (depth < maxDepth) {
        for (int v = 0; v < n; ++v) {
            // (A) Skip toggling the same vertex consecutively
            if (depth > 0 && path[depth - 1] == v) {
                continue;
            }

            // (B) Check if vertex v has at least two neighbors
            int neighbor_count = 0;
            for (int w = 0; w < m; w++) {
                setword row_chunk = GRAPHROW(g, v, m)[w];
                neighbor_count += __builtin_popcountll(row_chunk);
                if (neighbor_count >= 2) break;
            }
            if (neighbor_count < 2) {
                continue;
            }

            // (C) Perform local complement on vertex v
            graph *g_next = local_complement(g, v, m, n);

            // (D) Record this vertex in the path
            path[depth] = v;

            // (E) Save the old flags
            int oldPermMatched = permMatched;
            int oldGraphMatched = graphMatched;

            // (F) Recurse with the new graph
            int res = explore(g_next, m, n,
                              graphs, count,
                              path, depth + 1, maxDepth,
                              result,
                              perms, pcount, gCompare,
                              permMatched, graphMatched);

            // (G) Update flags based on the recursive result
            permMatched = oldPermMatched || (res >> 1);
            graphMatched = oldGraphMatched || (res & 1);

            // (H) Free the toggled graph
            free(g_next);
        }
    }

    /*
     * 4. Combine the two boolean flags (permMatched, graphMatched)
     *    into a single integer: the high bit for permMatched,
     *    and the low bit for graphMatched.
     */
    return (permMatched << 1) + graphMatched;
}


/*--------------------------------------------------------------------------------
 * construct_graph:
 *   Construct a graph by reading edges from the user. The input format is
 *   pairs of integers "u v", where 'u' and 'v' are vertices connected by an edge.
 *   Input terminates when the pair "-1 -1" is entered.
 *
 * Parameters:
 *   g - The graph to construct.
 *   m - The number of words per row in the graph's adjacency matrix.
 *   n - The number of vertices in the graph.
 *-------------------------------------------------------------------------------*/
void construct_graph(graph *g, int m, int n)
{
    EMPTYGRAPH(g, m, n);

    printf("Enter edges as pairs (terminate with -1 -1):\n");
    int u, v;
    while (scanf("%d %d", &u, &v) == 2 && u != -1 && v != -1) {
        ADDONEEDGE(g, u, v, m);
    }
}

/*--------------------------------------------------------------------------------
 * main:
 *   The entry point of the program. Reads input from the user, constructs the
 *   base graph and (optionally) a comparison graph, and performs a depth-first
 *   search to explore local complements.
 *
 * Steps:
 *   1) Build the base graph 'g' and optionally the comparison graph 'gCompare'.
 *   2) Read user-provided permutations.
 *   3) Prompt the user for the maximum depth of exploration.
 *   4) Perform DFS to find local complements that satisfy permutations or match
 *      the comparison graph.
 *   5) Print the results and the total number of unique graphs discovered.
 *
 * Returns:
 *   The total number of unique graphs (non-isomorphic) discovered during DFS.
 *-------------------------------------------------------------------------------*/
int main(void)
{
    int n, m;
    printf("Number of vertices: ");
    if (scanf("%d", &n) != 1 || n <= 0 || n > MAXN) {
        fprintf(stderr, "Invalid number of vertices. Must be between 1 and %d.\n", MAXN);
        return EXIT_FAILURE;
    }
    m = SETWORDSNEEDED(n);

    graph g[MAXN * MAXM], gCompare[MAXN * MAXM];
    printf("Construct the base graph:\n");
    construct_graph(g, m, n);

    // Ask the user if they want to compare graphs
    int compareGraphs = 0;
    printf("\nDo you want to compare the local complements with another graph? (1 for Yes, 0 for No): ");
    if (scanf("%d", &compareGraphs) != 1 || (compareGraphs != 0 && compareGraphs != 1)) {
        fprintf(stderr, "Invalid input. Must be 1 (Yes) or 0 (No).\n");
        return EXIT_FAILURE;
    }

    if (compareGraphs) {
        printf("Construct the graph for comparison:\n");
        construct_graph(gCompare, m, n);
    }

    int pcount;
    printf("\nHow many permutations do you want to check? ");
    if (scanf("%d", &pcount) != 1 || pcount < 0) {
        fprintf(stderr, "Invalid number of permutations.\n");
        return EXIT_FAILURE;
    }

    // Read permutations
    int **perms = malloc(pcount * sizeof(int *));
    if (!perms) {
        fprintf(stderr, "Memory allocation failed for perms[]\n");
        return EXIT_FAILURE;
    }
    for (int i = 0; i < pcount; i++) {
        perms[i] = malloc(n * sizeof(int));
        if (!perms[i]) {
            fprintf(stderr, "Memory allocation failed for perms[%d]\n", i);
            return EXIT_FAILURE;
        }
        printf("Permutation %d (enter %d numbers):\n", i + 1, n);
        for (int j = 0; j < n; j++) {
            if (scanf("%d", &perms[i][j]) != 1) {
                fprintf(stderr, "Invalid input for permutation %d.\n", i + 1);
                return EXIT_FAILURE;
            }
        }
    }

    // Prompt for maxDepth
    int maxDepth;
    printf("Enter the maximum depth for exploration: ");
    if (scanf("%d", &maxDepth) != 1 || maxDepth < 0) {
        fprintf(stderr, "Invalid depth value.\n");
        return EXIT_FAILURE;
    }

    // Optionally track discovered graphs
    graph *graphs[MAXGRAPHS];
    int count = 0;

    // Path array
    int path[MAXN] = {0};

    // Prepare LCResult
    LCResult result;
    result.gPerm = NULL;
    memset(result.pathPerm, 0, sizeof(result.pathPerm));
    result.depthPerm = 0;
    memset(result.pathGraph, 0, sizeof(result.pathGraph));
    result.depthGraph = 0;

    // DFS over local complements
    int res = explore(g, m, n,
                      graphs, &count,
                      path, 0, maxDepth,
                      &result,
                      perms, pcount, compareGraphs ? gCompare : NULL,
                      0, 0);

    // Print results
    if (compareGraphs && (res % 2)) {
        printf("\nLocal complement found that is equal to the designated graph!\n");
        printf("Vertex path used (local complements on these vertices): ");
        for (int i = 0; i < result.depthGraph; i++) {
            printf("%d ", result.pathGraph[i]);
        }
        printf("\n");
    } else if (compareGraphs) {
        printf("\nNo local complement found that is equal to the designated graph.\n");
    }

    if (res >> 1) {
        printf("\nLocal complement found that is equal to all permutations!\n");
        printf("Vertex path used (local complements on these vertices): ");
        for (int i = 0; i < result.depthPerm; i++) {
            printf("%d ", result.pathPerm[i]);
        }
        printf("\nEdges of final G':\n");
        print_graph(result.gPerm, m, n);
        free(result.gPerm);
    } else {
        printf("\nNo local complement found that is equal to all permutations.\n");
    }

    // Print the total number of unique graphs
    printf("\nTotal unique graphs (non-isomorphic): %d\n", count);

    // Now print all discovered graphs
    if (count > 1) {
        printf("\nAll discovered graphs:\n");
        for (int idx = 0; idx < count; idx++) {
            printf("Graph #%d:\n", idx);
            print_graph(graphs[idx], m, n);
            printf("\n");
        }
    }

    // Cleanup
    for (int idx = 0; idx < count; idx++) {
        free(graphs[idx]);
    }
    for (int i = 0; i < pcount; i++) {
        free(perms[i]);
    }
    free(perms);

    return count; // Return the size of 'graphs'
}
