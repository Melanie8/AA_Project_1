/*
 * Written by :
 * Melanie Sedda <sedda@kth.se>
 * Catarina Vaz  <acvaz@kth.se>
 * October 2015
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <vector>
#include <utility>
#include <algorithm>
using namespace std;

int print = 0;
int print_time = 0;

// Initialize a matrix AxB
int** init_matrix_int(int A, int B) {
    int **matrix = (int **)calloc(A, sizeof(int *));
    int i;
    for (i=0; i<A; i++) {
        matrix[i] = (int *)calloc(B, sizeof(int));
    }
    return matrix;
}

double** init_matrix_double(int A, int B) {
    double **matrix = (double **)calloc(A, sizeof(double *));
    int i;
    for (i=0; i<A; i++) {
        matrix[i] = (double *)calloc(B, sizeof(double));
    }
    return matrix;
}

long** init_matrix_long(int A, int B) {
    long **matrix = (long **)calloc(A, sizeof(long *));
    int i;
    for (i=0; i<A; i++) {
        matrix[i] = (long *)calloc(B, sizeof(long));
    }
    return matrix;
}

long int** init_matrix_long_int(int A, int B) {
    long int **matrix = (long int **)calloc(A, sizeof(long int *));
    int i;
    for (i=0; i<A; i++) {
        matrix[i] = (long int *)calloc(B, sizeof(long int));
    }
    return matrix;
}


long timediff(clock_t t1, clock_t t2) {
    long elapsed;
    elapsed = ((double)t2 - t1) / CLOCKS_PER_SEC * 1000000;
    return elapsed;
}

// Compute the distances between all the pairs of points
long int** compute_distances(int N, double **points) {
    if (print) {
        printf("compute_distances\n");
    }
    clock_t start = clock();
    long long square = 0;
    double d = 0;
    long int **distances = init_matrix_long(N, N);
    if (!distances)
        return NULL;
    int i,j;
    for (i=0; i<N; i++) {
        distances[i][i] = 0;
        for (j=i+1; j<N; j++) {
            square = pow(points[i][0]-points[j][0],2) + pow(points[i][1]-points[j][1],2);
            d = round(sqrt(square));
            distances[i][j] = d;
            distances[j][i] = d;
            if (print)
                printf("Distance from %d to %d: %f\n", i, j, d);
        }
    }
    clock_t end = clock();
    long elapsed = timediff(start, end);
    if (print_time)
        printf("Compute distances took %ld microseconds\n", elapsed);
    return distances;
}

// Nearest neighbor tour
int* nearest_neighbor(int N, double** points, long int** distances) {
    if (print)
        printf("greedy_tour\n");
    clock_t start = clock();
    int *tour = (int *)malloc(N*sizeof(int));
    int *used = (int *)calloc(N, sizeof(int));
    long int length_tour = 0;
    tour[0] = 0;
    printf("0\n");
    if (print)
        printf("Point 0 of the tour: 0\n");
    used[0] = 1;
    int best,i,j;
    for (i=1; i<N; i++) {
       best = -1;
       for (j=1; j<N; j++) {
          if (used[j]==0 && (best == -1 || distances[tour[i-1]][j] < distances[tour[i-1]][best])){
             best = j;
          }
       }
       tour[i] = best;
       length_tour += distances[tour[i-1]][tour[i]];
       printf("%d\n", best);
       used[best] = 1;
    }
    length_tour += distances[tour[N-1]][tour[0]];
    if (print)
        printf("Length tour = %ld\n", length_tour);
    if (print_time) {
        clock_t end = clock();
        long elapsed = timediff(start, end);
        printf("Nearest neighbor took %ld microseconds\n", elapsed);
    }
    return tour;
}

// Greedy tour
int* greedy(int N, double** points, long int** distances) {
    if (print)
        printf("greedy_tour\n");
    clock_t start_greedy = clock();
    int *tour = (int *)malloc(N*sizeof(int));
    int *degrees = (int *)calloc(N, sizeof(int));
    long int length_tour = 0;
    int i,j,k,r;
    /*int k, i, j, r, besti, bestj, prev, current, proceed, connected;
    long int bestdist;
    int **used = init_matrix_int(N, N);*/

    // Create vector with all the distances
    vector<pair<int, pair<int, int> > > v;
    int L = (N*(N-1))/2;
    v.resize(L);
    k = 0;
    for (i=0; i<N; i++) {
        for (j=i+1; j<N; j++) {
            v[k].first = distances[i][j];
            v[k].second.first = i;
            v[k].second.second = j;
            k++;
        }
    }

    // Sort
    sort(v.begin(), v.end());

    // Union find
    int *sets = (int *)calloc(N, sizeof(int));
    vector <int> neighbor[1000];
    for (i=0; i<N; i++) {
        sets[i] = i;
    }
    int nedges  = 0;
    int s;
    for (k=0; k<L && nedges<N; k++){
        i = v[k].second.first;
        j = v[k].second.second;
        if (degrees[i]<=1 && degrees[j]<=1 && (sets[i]!=sets[j] || nedges==N-1)) {
            printf("i = %d, j = %d\n", i, j);
            neighbor[i].push_back(j);
            neighbor[j].push_back(i);
            degrees[i]++;
            degrees[j]++;
            s = sets[i];
            if (nedges<N-1){
                for (r=0; r<L; r++){
                    if (sets[r] == s)
                        sets[r] = sets[j];
                }
            }
            nedges++;
            length_tour += distances[i][j];
        }
    }
    printf("nedges = %d, length_tour = %ld\n", nedges, length_tour);
    for (i=0; i<N; i++){
        printf("i = %d, neighbor[i][0] = %d, neighbor[i][1] = %d", i, neighbor[i][0], neighbor[i][1]);
    }

    // Deduce the tour
    int current = 0;
    int prev = 0;
    int next;
    tour[0] = 0;
    printf("0\n");
    for (k = 1; k<N; k++){
        printf("0 : %d, 1 : %d\n", neighbor[current][0], neighbor[current][1]);
        if (neighbor[current][0] == prev)
            next = neighbor[current][1];
        else
            next = neighbor[current][0];
        printf("%d\n", next);
        prev = current;
        current = next;
        tour[k] = next;
    }

    /*for (k=0; k<N; k++) {
        // Find the best edge
        besti = 0;
        bestj = 0;
        bestdist = 0;
        for (i=0; i<N; i++) {
            for (j=i+1; j<N; j++) {
                // Update the best edge
                if (used[i][j]==0 && degrees[i]<=1 && degrees[j]<=1 && (bestdist == 0 || distances[i][j] < bestdist)){
                    // Check if it would create a cycle of length < N
                    connected = 0;
                    if (k < N-1){
                        proceed = 1;
                        prev = i;
                        current = i;
                        while(proceed){
                            proceed = 0;
                            for (r=0; r<N && proceed == 0; r++) {
                                if (r!=prev && r!=current && used[current][r]) {
                                    if (r==j) {
                                        connected = 1;
                                    } else {
                                        proceed = 1;
                                        prev = current;
                                        current = r;
                                    }
                                }
                            }
                        }
                    }
                    if (connected == 0) {
                        besti = i;
                        bestj = j;
                        bestdist = distances[besti][bestj];
                    }
                }
            }
        }
        // Add the best edge
        used[besti][bestj] = 1;
        used[bestj][besti] = 1;
        degrees[besti] += 1;
        degrees[bestj] += 1;
        length_tour += distances[besti][bestj];
        if (print)
            printf("\n \n Best edge : (%d,%d)\n \n", besti, bestj);
    }
    if (print)
        printf("length_tour = %ld\n", length_tour);
    // Deduce tour
    tour[0] = 0;
    printf("0\n");
    current = tour[0];
    prev = tour[0];
    int found_next;
    for (k=1; k<N; k++){
        found_next = 0;
        for (r=0; r<N && found_next == 0; r++) {
            if (r!=prev && r!=current && used[current][r]) {
                prev = tour[k-1];
                tour[k] = r;
                current = r;
                printf("%d\n", r);
                found_next = 1;
            }
        }
    }*/
    if (print_time) {
        clock_t end_greedy = clock();
        long elapsed = timediff(start_greedy, end_greedy);
        if (print_time)
            printf("Greedy took %ld microseconds\n", elapsed);
    }
    return tour;
}

// Clarke Wright tour
int* clarke_wright(int N, double** points, long int** distances){
    int* tour = (int *)malloc(N*sizeof(int));
    int hub = rand() % N;
    if (print)
        printf("Hub : %d\n", hub);
    long int **saving = init_matrix_long_int(N, N);
    int i, j, k;
    for (i=0; i<N; i++) {
        for (j=i+1; j<N; j++) {
            saving[i][j] = distances[i][hub]+distances[j][hub]-distances[i][j];
            if (print) {
                printf("Saving %d %d: %ld\n", i, j, saving[i][j]);
            }
        }
    }
    //sort(saving, &saving[(N-1)*N*0.5-1], );
    for (k=0; k<N; k++) {
        // Find the best edge
        for (i=0; i<N; i++) {
            for (j=i+1; j<N; j++) {

            }
        }
    }
    return tour;
}

int main(int argc, char *argv[]) {
    /* Get arguments */
    int N;
    scanf("%d", &N);
    double **points = init_matrix_double(N, 2);
    long i;
    for (i=0; i<N; i++) {
        scanf("%lf", points[i]);
        scanf("%lf", points[i]+1);
        if (print)
            printf("Node %ld: (%f, %f)\n", i, points[i][0], points[i][1]);
    }

    /* Compute distances between all the points */
    long **distances = compute_distances(N, points);

    /* Find a tour */
    int* tour_greedy = greedy(N, points, distances);
    //int* tour_nn = nearest_neighbor(N, points, distances);
    //int* tour = clarke_wright(points, distances);

    return EXIT_SUCCESS;
}
