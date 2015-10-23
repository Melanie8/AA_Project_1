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

long N;


void free_matrix(long M, double **matrix) {
    if (!matrix) {
        return;
    }
    long i;
    for (i=0; i<M; i++) {
        if (matrix[i]) {
            free(matrix[i]);
            matrix[i] = NULL;
        }
    }
    free(matrix);
}

double** init_matrix(long M, long N) {
    double **matrix = (double **)malloc(M*sizeof(double *));
    if (!matrix) {
        fprintf(stderr, "malloc failed\n");
        return NULL;
    }
    long i;
    for (i=0; i<N; i++) {
        matrix[i] = (double *)malloc(N*sizeof(double));
        if (!matrix[i]) {
            fprintf(stderr, "malloc failed\n");
            free_matrix(M, matrix);
            return NULL;
        }
    }
    return matrix;
}


double** read_points(char *file){
    long i, j;
    FILE *fp;
    char *line = NULL;
    size_t len = 0;
    ssize_t read;

    fp = fopen(file, "r");
    if (fp == NULL)
        return NULL;

    read = getline(&line, &len, fp);
    N = strtol(line, (char **)NULL, 10);
    printf("N = %ld\n", N);
    double **points = init_matrix(N, 2);
    if (!points)
        return NULL;
    for (i=0; i<N; i++) {
        read = getline(&line, &len, fp);
        if (read <= 0) {
            perror("convert_line(getline failed)\n");
            if (line) {
                free(line); line = NULL;
            }
            return NULL;
        }
        for (j=0; j<2; j++) {
            sscanf(line, "%lf", &points[i][j]);
            printf("Line %ld column %ld: %lf\n", i, j, points[i][j]);

        }
    }

    fclose(fp);
    if (line)
        free(line);

    /*N = (long) strtol(line[i], (char **)NULL, 10);
    points = init_matrix(N, 2);
    if (!points)
        return NULL;

    for (i=0; i<N; i++) {
        read = getline(&line, &len, f);
        if (read <= 0) {
            perror("convert_line(getline failed)\n");
            if (line) {
                free(line); line = NULL;
            }
            return NULL;
        }
        printf("Line = %s", line)
        for (j=0; j<2; j++) {
            points[i][j] = strtod(line[i], NULL);
            printf("Line %ld column %ld: %ld", i, j, points[i][j]);

        }
    }*/
    return points;

}

double** compute_distances(double **points) {
    double **distances = init_matrix(N, N);
    if (!distances)
        return NULL;
    long i,j;
    for (i=0; i<N; i++) {
        distances[i][i] = 0;
        for (j=i+1; j<N; j++) {
            double d = sqrt(pow(points[i][0]-points[j][0],2) + pow(points[i][1]-points[j][1],2));
            distances[i][j] = d;
            distances[j][i] = d;
            printf("Distance from %ld to %ld: %f\n", i, j, d);
        }
    }
    return distances;
}

long* greedy_tour(double** points, double** distances) {
    long *tour = (long *)malloc(N*sizeof(long));
    long *used = (long *)calloc(N, sizeof(long));
    if (!tour || !used) {
        fprintf(stderr, "malloc failed\n");
        return NULL;
    }
    tour[0] = 0;
    printf("Point 0 of the tour: 0\n");
    used[0] = 1;
    long best,i,j;
    for (i=1; i<N; i++) {
       best = -1;
       for (j=0; j<N; j++) {
          if (used[j]==0 && (best = -1 || distances[tour[i-1]][j] < distances[tour[i-1]][best]))
             best = j;
       }
       tour[i] = best;
       printf("Point %ld of the tour: %ld\n", i, best);
       used[best] = 1;
    }
   return tour;
}

long* clarke_wright(double** points, double** distances){
    long *tour = (long *)malloc(N*sizeof(long));
    long hub = rand() % N;
    printf("Hub : %ld\n", hub);
    long length = round(N*(N-1)*0.5);
    double *saving = (double *)malloc(length*sizeof(double));
    double **pairs = init_matrix(length, 2);

    if (!saving || !pairs) {
        fprintf(stderr, "malloc failed\n");
        return NULL;
    }
    long i, j;
    long count = 0;
    for (i=0; i<N; i++){
        printf("i = %ld\n", i);
        for (j=i+1; j<N; j++){
            printf("count / total = %ld/%ld\n", count,length);
            printf("distances[i][hub] = %lf\n", distances[i][hub]);
            printf("distances[j][hub] = %lf\n", distances[j][hub]);
            printf("distances[i][j] = %lf\n", distances[i][j]);
            saving[count] = distances[i][hub]+distances[j][hub]-distances[i][j];
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
    /* int N = (int) strtol(argv[1], (char **)NULL, 10);
    printf("Number of nodes: %d\n", N);
    double **points = init_matrix(N, 2);
    if (!points)
        return EXIT_FAILURE;
    long i;
    for (i=0; i<N; i++) {
        points[i][0] = strtod(argv[2+2*i], (char **)NULL);
        points[i][1] = strtod(argv[3+2*i], (char **)NULL);
        printf("Node %ld: (%f, %f)\n", i, points[i][0], points[i][1]);
    }*/
    double **points = read_points(argv[1]);

    /* Compute distances between all the points */
    double **distances = compute_distances(points);

    /* Find a tour */
    //long *tour = greedy_tour(points, distances);
    long *tour = clarke_wright(points, distances);

    return EXIT_SUCCESS;
}
