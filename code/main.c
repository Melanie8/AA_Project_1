/*
 * Written by :
 * Melanie Sedda <sedda@kth.se>
 * Catarina Vaz  <acvaz@kth.se>
 * October 2015
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <dirent.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <fcntl.h>
#include <errno.h>
#include <stdbool.h>
#include <math.h>

int print = 0;

// Initialize a matrix AxB
double** init_matrix(long A, long B) {
    double **matrix = (double **)malloc(A*sizeof(double *));
    long i;
    for (i=0; i<A; i++) {
        matrix[i] = (double *)malloc(B*sizeof(double));
    }
    return matrix;
}

// Compute the distances between all the pairs of points
double** compute_distances(int N, double **points) {
    if (print)
        printf("compute_distances\n");
    long long square = 0;
    double d = 0;
    double **distances = init_matrix(N, N);
    if (!distances)
        return NULL;
    long i,j;
    for (i=0; i<N; i++) {
        distances[i][i] = 0;
        for (j=i+1; j<N; j++) {
            square = pow(points[i][0]-points[j][0],2) + pow(points[i][1]-points[j][1],2);
            d = round(sqrt(square));
            distances[i][j] = d;
            distances[j][i] = d;
            if (print)
                printf("Distance from %ld to %ld: %f\n", i, j, d);
        }
    }
    return distances;
}

// Nearest neighbor
long* nearest_neighbor(int N, double** points, double** distances) {
    if (print)
        printf("greedy_tour\n");
    long *tour = (long *)malloc(N*sizeof(long));
    long *used = (long *)calloc(N, sizeof(long));
    long length_tour = 0;
    tour[0] = 0;
    printf("0\n");
    if (print)
        printf("Point 0 of the tour: 0\n");
    used[0] = 1;
    long best,i,j;
    for (i=1; i<N; i++) {
       best = -1;
       for (j=1; j<N; j++) {
          if (print)
            printf("j = %ld, used[j] = %ld, distances[tour[i-1]][j] = %f, distances[tour[i-1]][best] = %f\n", j, used[j], distances[tour[i-1]][j], distances[tour[i-1]][best]);
          if (used[j]==0 && (best == -1 || distances[tour[i-1]][j] < distances[tour[i-1]][best])){
             if (print)
                printf("Update best = %ld\n", j);
             best = j;
          }
       }
       tour[i] = best;
       length_tour += distances[tour[i-1]][tour[i]];
       printf("%ld\n", best);
       if (print)
        printf("Point %ld of the tour: %ld\n", i, best);
       used[best] = 1;
    }
    length_tour += distances[tour[N-1]][tour[0]];
    if (print)
        printf("Length tour = %ld", length_tour);
    return tour;
}

// Clarke Wright tour
long* clarke_wright(int N, double** points, double** distances){
    long *tour = (long *)malloc(N*sizeof(long));
    long hub = rand() % N;
    if (print)
        printf("Hub : %ld\n", hub);
    long length = round(N*(N-1)*0.5);
    double *saving = (double *)malloc(length*sizeof(double));
    double **pairs = init_matrix(length, 2);
    long i, j;
    long count = 0;
    for (i=0; i<N; i++){
        if (print)
            printf("i = %ld\n", i);
        for (j=i+1; j<N; j++){
            if (print){
                printf("count / total = %ld/%ld\n", count,length);
                printf("distances[i][hub] = %lf\n", distances[i][hub]);
                printf("distances[j][hub] = %lf\n", distances[j][hub]);
                printf("distances[i][j] = %lf\n", distances[i][j]);
            }
            saving[count] = distances[i][hub]+distances[j][hub]-distances[i][j];
            if (print)
                printf("Saving %ld %ld: %lf\n", i, j, saving[count]);
            pairs[count][0] = i;
            pairs[count][1] = j;
            count +=1;
        }
    }
    //sort(saving, &saving[(N-1)*N*0.5-1], );
    return tour;
}

int main(int argc, char *argv[]) {
    /* Get arguments */
    int N;
    scanf("%d", &N);
    double **points = init_matrix(N, 2);
    long i;
    for (i=0; i<N; i++) {
        scanf("%lf", points[i]);
        scanf("%lf", points[i]+1);
        if (print)
            printf("Node %ld: (%f, %f)\n", i, points[i][0], points[i][1]);
    }

    /* Compute distances between all the points */
    double **distances = compute_distances(N, points);

    /* Find a tour */
    long *tour = nearest_neighbor(N, points, distances);
    //long *tour = clarke_wright(points, distances);

    return EXIT_SUCCESS;
}
