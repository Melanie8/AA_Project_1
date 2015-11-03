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
int print_length = 0;

// Initialize a matrix AxB
double** init_matrix_double(int A, int B) {
    double **matrix = (double **)calloc(A, sizeof(double *));
    int i;
    for (i=0; i<A; i++) {
        matrix[i] = (double *)calloc(B, sizeof(double));
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
    double square = 0;
    double d = 0;
    long int **distances = init_matrix_long_int(N, N);
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

// Enhance existent tour
pair<long int, int *> enhance(int N, long int** distances, long int length_tour, int* tour){
    int cont = 1;
    int r, k, s, enhancement, temp;
    while(cont){
        cont = 0;
        for (r=1; r < (N+1)/2; r++) {
            for (k = 0; k < N; k++) {
                enhancement = distances[tour[k]][tour[(k+1)%N]] + distances[tour[(k+1+r)%N]][tour[(k+2+r)%N]] - distances[tour[k]][tour[(k+1+r)%N]] - distances[tour[(k+1)%N]][tour[(k+2+r)%N]];
                if (enhancement > 0) {
                    cont = 1;
                    // printf("Enhance for %d (r = %d)\n", enhancement, r);
                    // We switch k+1 and k+1+r, k+2 and k+r, etc
                    for (s = 1; k+s < k+2+r-s; s++) {
                        temp = tour[(k+s)%N];
                        tour[(k+s)%N] = tour[(k+2+r-s)%N];
                        tour[(k+2+r-s)%N] = temp;
                    }
                    length_tour -= enhancement;
                }
            }
        }
    }
    if (print) {
        for (k = 0; k < N; k++){
            printf("%d\n", tour[k]);
        }
    }
    if (print_length)
        printf("Length tour = %ld\n", length_tour);
    return make_pair(length_tour, tour);
}

// Nearest neighbor tour
pair<long int, int *> nearest_neighbor(int N, double** points, long int** distances) {
    if (print)
        printf("greedy_tour\n");
    clock_t start = clock();
    int *tour = (int *)malloc(N*sizeof(int));
    int *used = (int *)calloc(N, sizeof(int));
    long int length_tour = 0;
    tour[0] = 0;
    if (print)
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
       if (print)
          printf("%d\n", best);
       used[best] = 1;
    }
    length_tour += distances[tour[N-1]][tour[0]];
    if (print_length)
        printf("Length tour = %ld\n", length_tour);
    if (print_time) {
        clock_t end = clock();
        long elapsed = timediff(start, end);
        printf("Nearest neighbor took %ld microseconds\n", elapsed);
    }
    return enhance(N, distances, length_tour, tour);
}

// Greedy tour
pair<long int, int *> greedy(int N, double** points, long int** distances) {
    if (print)
        printf("greedy_tour\n");
    clock_t start_greedy = clock();
    int *tour = (int *)malloc(N*sizeof(int));
    int *degrees = (int *)calloc(N, sizeof(int));
    long int length_tour = 0;
    int i,j,k,r;

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
            neighbor[i].push_back(j);
            neighbor[j].push_back(i);
            degrees[i]++;
            degrees[j]++;
            s = sets[i];
            if (nedges<N-1){
                for (r=0; r<N; r++){
                    if (sets[r] == s)
                        sets[r] = sets[j];
                }
            }
            nedges++;
            length_tour += distances[i][j];
        }
    }
    // Deduce the tour
    int current = 0;
    int prev = 0;
    int next;
    tour[0] = 0;
    if (print)
        printf("0\n");
    for (k = 1; k<N; k++){
        if (neighbor[current][0] == prev)
            next = neighbor[current][1];
        else
            next = neighbor[current][0];
        if (print)
            printf("%d\n", next);
        prev = current;
        current = next;
        tour[k] = next;
    }
    if (print_length)
        printf("Length tour = %ld\n", length_tour);
    if (print_time) {
        clock_t end_greedy = clock();
        long elapsed = timediff(start_greedy, end_greedy);
        if (print_time)
            printf("Greedy took %ld microseconds\n", elapsed);
    }
    return enhance(N, distances, length_tour, tour);
}

// Clarke Wright tour
pair<long int, int *> clarke_wright(int N, double** points, long int** distances){
    if (print)
        printf("clarke-write\n");
    clock_t start_greedy = clock();

    int* tour = (int *)malloc(N*sizeof(int));
    int *degrees = (int *)calloc(N, sizeof(int));
    long int length_tour = 0;

    int hub = rand() % N;
    if (print)
        printf("Hub : %d\n", hub);
    int i, j, k, r, s, nedges;

    // Create vector with all the shortcuts
    vector<pair<int, pair<int, int> > > v;
    int L = ((N-1)*(N-2))/2;

    v.resize(L);
    k = 0;
    for (i=0; i<N; i++) {
        for (j=i+1; j<N; j++) {
            if (i!=hub && j!= hub) {
                v[k].first = distances[i][j]-distances[i][hub]-distances[j][hub];
                v[k].second.first = i;
                v[k].second.second = j;
                if (print)
                    printf("k = %d, shortcut = %d, i = %d, j = %d\n", k, v[k].first, v[k].second.first, v[k].second.second);
                k++;
            }
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
    nedges  = 0;
    for (k=0; k<L && nedges<N-2; k++){
        i = v[k].second.first;
        j = v[k].second.second;
        if (degrees[i]<=1 && degrees[j]<=1 && sets[i]!=sets[j]) {
            if (print)
                printf("Edge taken : (%d, %d), shortcut = %d\n", i, j, v[k].first);
            neighbor[i].push_back(j);
            neighbor[j].push_back(i);
            degrees[i]++;
            degrees[j]++;
            s = sets[i];
            if (nedges<N-1){
                for (r=0; r<N; r++){
                    if (sets[r] == s)
                        sets[r] = sets[j];
                }
            }
            nedges++;
            length_tour += distances[i][j];
        }
    }

    if(N <= 2) return make_pair(1000000000, tour);

    // Deduce the tour
    tour[0] = hub;
    if (print)
        printf("%d\n", hub);
    int first = -1;
    for(i=0; i<N && first == -1; i++) {
        if (degrees[i]==1) {
            first = i;
        }
    }
    tour[1] = first;
    if (print)
        printf("%d\n", first);
    length_tour += distances[first][hub];
    int current = first;
    int prev = hub;
    int next;
    for (k = 2; k<N; k++){
        if (neighbor[current][0] == prev)
            next = neighbor[current][1];
        else
            next = neighbor[current][0];
        if (print)
            printf("%d\n", next);
        prev = current;
        current = next;
        tour[k] = next;
    }
    length_tour += distances[current][hub];
    if (print_length)
        printf("Length tour = %ld\n", length_tour);
    if (print_time) {
        clock_t end_greedy = clock();
        long elapsed = timediff(start_greedy, end_greedy);
        if (print_time)
            printf("Clarke-Wright took %ld microseconds\n", elapsed);
    }
    return enhance(N, distances, length_tour, tour);
}

//Christofides
pair<long int,int *> christofides(int N, double** points, long int** distances){
	if(print)
		printf("christofides\n");
	clock_t start_greedy = clock();
	int* tour = (int *)malloc(N*sizeof(int));
	int *degrees = (int *)calloc(N, sizeof(int));
	long int length_tour = 0;
    int i,j,k,r;
    
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
    //sort
    sort(v.begin(),v.end());
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
        if (sets[i]!=sets[j] || nedges==N-1) {
            neighbor[i].push_back(j);
            neighbor[j].push_back(i);
            degrees[i]++;
            degrees[j]++;
            s = sets[i];
            if (nedges<N-1){
                for (r=0; r<N; r++){
                    if (sets[r] == s)
                        sets[r] = sets[j];
                }
            }
            nedges++;
            length_tour += distances[i][j];
        }
    }


// Optimal tour : brute force
int* brute_force(int N, double** points, long int** distances){
    if (print)
        printf("clarke-write\n");
    clock_t start_greedy = clock();

    int* tour = (int *)malloc(N*sizeof(int));

    //??

    if (print_time) {
        clock_t end_greedy = clock();
        long elapsed = timediff(start_greedy, end_greedy);
        if (print_time)
            printf("Clarke-Wright took %ld microseconds\n", elapsed);
    }
    return tour;

}

int main(int argc, char *argv[]) {
    /* Get arguments */
    int N, i;
    scanf("%d", &N);
    //N = 1000;
    double **points = init_matrix_double(N, 2);

    /*
    int x = time(NULL) % 1000;
    x = 38;
    srand(x);
    printf("Seed = %d\n", x);

    for(i = 0; i < N; i++) {
        points[i][0] = rand() % 1000;
        points[i][1] = rand() % 1000;
    }

    */

    srand(time(NULL));

    for (i=0; i<N; i++) {
        scanf("%lf", points[i]);
        scanf("%lf", points[i]+1);
        if (print)
            printf("Node %d: (%f, %f)\n", i, points[i][0], points[i][1]);
    }


    /* Compute distances between all the points */
    long **distances = compute_distances(N, points);

    /* Find tours */
    pair <long int, int*> tour_nn = nearest_neighbor(N, points, distances);
    pair <long int, int*> tour_greedy = greedy(N, points, distances);

    long int test_length = 1000000000;
    pair <long int, int*> tour_cw;
    for(i = 0; i < 5; i++){
        pair <long int, int*> tour_cw_test = clarke_wright(N, points, distances);
        if(tour_cw_test.first < test_length) {
            tour_cw = tour_cw_test;
            test_length = tour_cw_test.first;
        }
    }

    //int* tour_brute = brute_force(N, points, distances);

    /* Find the best one */
    //long int best_length = min(tour_nn.first, tour_greedy.first);
    long int best_length = min(tour_nn.first, min(tour_greedy.first, tour_cw.first));
    int *best_tour;
    if (best_length == tour_nn.first)
        best_tour = tour_nn.second;
    else if (best_length == tour_greedy.first)
        best_tour = tour_greedy.second;
    else
        best_tour = tour_cw.second;

    /* Print the best tour */
    for (i=0; i<N; i++){
        printf("%d\n", best_tour[i]);
    }
    if (print_length)
        printf("Best length = %ld\n", best_length);

    return EXIT_SUCCESS;
}
