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
#include <stdbool.h>
using namespace std;

int print = 0;
int print_time = 0;
int print_length = 0;
vector<pair<int, int> > closeto[1000];

int version = 0;
int nbEnhancement = 25;

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
    double dist_prec = 0;
    int d = 0;
    long int **distances = init_matrix_long_int(N, N);
    if (!distances)
        return NULL;
    int i,j;
    for (i=0; i<N; i++) {
        distances[i][i] = 0;
        for (j=i+1; j<N; j++) {
            dist_prec = sqrt(pow(points[i][0]-points[j][0],2) + pow(points[i][1]-points[j][1],2));
            d = round(dist_prec);
            distances[i][j] = d;
            distances[j][i] = d;
            /*if (print)
                printf("Distance from %d to %d: %f\n", i, j, d);*/
        }
    }
    clock_t end = clock();
    long elapsed = timediff(start, end);
    if (print_time)
        printf("Compute distances took %ld microseconds\n", elapsed);
    return distances;
}

int *neighbors_to_tour(int N, vector <int> neighbor[1000]){
    int *tour = (int *)malloc(N*sizeof(int));
    int current = 0;
    int prev = 0;
    int next;
    tour[0] = 0;
    if (print)
        printf("0\n");
    for (int k = 1; k<N; k++){
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
    return tour;
}

// Enhance existent tour (with 2-opt)
pair<long int, int *> enhance(int N, long int** distances, long int length_tour, int* tour){
	 clock_t start = clock();
    int cont = 1;
    int r, k, s, enhancement, temp;
    while(cont){
        cont = 0;
        for (r = 1; r < (N+1)/2; r++) {
            for (k = 0; k < N; k++) {
                enhancement = distances[tour[k]][tour[(k+1)%N]] + distances[tour[(k+1+r)%N]][tour[(k+2+r)%N]] - distances[tour[k]][tour[(k+1+r)%N]] - distances[tour[(k+1)%N]][tour[(k+2+r)%N]];
                if (enhancement > 0) {
                    cont = 1;
                    //printf("Enhance for %d : (%d -> %d) and (%d -> %d)\n", enhancement, tour[k], tour[(k+1)%N], tour[(k+1+r)%N], tour[(k+2+r)%N]);
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
        printf("Length tour after enhance = %ld\n", length_tour);

    if (print_time) {
        clock_t end = clock();
        long elapsed = timediff(start, end);
        printf("Enhance took %ld microseconds\n", elapsed);
    }

    return make_pair(length_tour, tour);
}


// Enhance existent tour (with 2-opt, quicker version)
pair<long int, int *> enhance2(int N, long int** distances, long int length_tour, int* tour){
    clock_t start = clock();
    int cont = 1;
    int pos[1000];
    int r, k, l, m, enhancement, temp;

    // pos[x] is the position of x in the actual tour
    for(k = 0; k < N; k++) pos[tour[k]] = k;

    while(cont){
        cont = 0;
        // We will try to change edge (a,b) with edge (c,d)
        for(int a = 0; a < N; a++) {
        	int pos_a = pos[a];
        	int pos_b = (pos_a + 1) % N;
        	int b = tour[pos_b];
            // We only consider c such that dist(b,c) < dist(b,a)
            for(l = 0; l < closeto[b].size() && closeto[b][l].first < distances[a][b]; l++) {
            	int c = closeto[b][l].second;
            	int pos_c = pos[c];
            	int pos_d = (pos_c + N-1) % N;
            	int d = tour[pos_d];
            	if(a != c && a != d && b != c && b != d){
            		enhancement = distances[a][b] + distances[c][d] - distances[a][d] - distances[b][c];
            		if(enhancement > 0){
            			// Instead of a -> b and d -> c, we put a -> d and b -> c, i.e. we replace the path (b ... d) by the path (d ... b)
            			int L = pos_d - pos_b + 1;
            			if(L < 0) L += N;
            			for(k = 0; k < L/2; k++){
            				int toreplace1 = (pos_b + k)%N;
            				int toreplace2 = (pos_d + N - k)%N;
            				int remember = tour[toreplace1];
            				tour[toreplace1] = tour[toreplace2];
            				tour[toreplace2] = remember;
            				pos[tour[toreplace1]] = toreplace1;
            				pos[tour[toreplace2]] = toreplace2;
            			}
            			length_tour -= enhancement;
            			cont = 1;
            			break; // We break the loop, because if we continue it then we re-use a -> b, which is not an edge anymore!
            		}
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
        printf("Length tour after enhance 2 = %ld\n", length_tour);

    if (print_time) {
        clock_t end = clock();
        long elapsed = timediff(start, end);
        printf("Enhance 2 took %ld microseconds\n", elapsed);
    }
    return make_pair(length_tour, tour);
}

bool isbetween(int x, int a, int b){ // Is a between x and y?
	if(b > a) return (x >= a && x <= b);
	else return (x >= a || x <= b);
}

void recompute(int N, int* tour, int* pos, vector<int> *neighbor){
	int current = 0;
   int prev = 0;
   int next;
   tour[0] = 0;
   pos[0] = 0;
   for (int k = 1; k<N; k++){
       if (neighbor[current][0] == prev)
           next = neighbor[current][1];
       else
           next = neighbor[current][0];
       prev = current;
       current = next;
       tour[k] = next;
       pos[tour[k]] = k;
   }
}

// Enhance existent tour (with 3-opt)
pair<long int, int *> enhance3(int N, long int** distances, long int length_tour, int* tour){
    clock_t start = clock();
    int cont = 1;
    int r, k, l, m, enhancement, temp;

    // pos[x] is the position of x in the actual tour

	 int *remembertour = (int *)malloc(N*sizeof(int));
	 int *besttour = (int *)malloc(N*sizeof(int));
	 long int rememberlength_tour = length_tour;
	 long int bestlength_tour = 1000000000;

	 for(int i = 0; i < N; i++)
	    remembertour[i] = tour[i];

	 vector<int> p;
	 p.resize(N);

	 for(int i = 0; i < N; i++)
	    p[i] = i;

	 for(int z = 0; z < nbEnhancement; z++){
    	random_shuffle(p.begin(), p.end());

    	if(z > 0){
    		for(int i = 0; i < N; i++)
    		    tour[i] = remembertour[i];
    		length_tour = rememberlength_tour;
    	}

		 int pos[1000];
		 vector<int> neighbor[1000];

		 for(k = 0; k < N; k++){
		 	pos[tour[k]] = k;
		 	neighbor[tour[k]].push_back(tour[(k+1)%N]);
		 	neighbor[tour[(k+1)%N]].push_back(tour[k]);
		 }

		 cont = 1;
		 while(cont){
		     cont = 0;
		     // We will try to change edge (a,b) with edge (c,d)
		     for(int q = 0; q < N; q++) {
		     	 int a = p[q];
		         int pos_a = pos[a];
		         int pos_b = (pos_a + 1) % N;
		     	 int b = tour[pos_b];
		     	 // We only consider c such that dist(b,c) < dist(b,a)
		         for(l = 0; l < closeto[b].size() && closeto[b][l].first < distances[a][b]; l++) {
		         	int c = closeto[b][l].second;
		         	int pos_c = pos[c];
		         	int pos_d = (pos_c + N-1) % N;
		         	int d = tour[pos_d];
		         	if(a != c && a != d && b != c && b != d){
		         		enhancement = distances[a][b] + distances[c][d] - distances[a][d] - distances[b][c];

		         		if(enhancement > 0){
		         			// Instead of a -> b and d -> c, we put a -> d and b -> c, i.e. we replace the path (b ... d) by the path (d ... b)
		         			//printf("Switch (%d %d) et (%d %d)\n", a, b, d, c);

		         			int L = pos_d-pos_b+1;
								if(L < 0) L += N;
								for(k = 0; k < L/2; k++)
								{
									int toreplace1 = (pos_b + k)%N;
									int toreplace2 = (pos_d + N - k)%N;
									int remember = tour[toreplace1];
									tour[toreplace1] = tour[toreplace2];
									tour[toreplace2] = remember;
									pos[tour[toreplace1]] = toreplace1;
									pos[tour[toreplace2]] = toreplace2;
								}

		         			if(neighbor[a][0] == b) neighbor[a][0] = d;
		         			else neighbor[a][1] = d;
		         			if(neighbor[b][0] == a) neighbor[b][0] = c;
		         			else neighbor[b][1] = c;
		         			if(neighbor[c][0] == d) neighbor[c][0] = b;
		         			else neighbor[c][1] = b;
		         			if(neighbor[d][0] == c) neighbor[d][0] = a;
		         			else neighbor[d][1] = a;
		         			length_tour -= enhancement;
		         			cont = 1;
		         			break; // We break the loop, because if we continue it then we re-use a -> b, which is not an edge anymore!
		         		}
		         		else {// 3-opt
		         			bool mustbreak = false;
		         			for(m = 0; m < closeto[d].size() && closeto[d][m].first + distances[b][c] < distances[a][b] + distances[c][d]; m++){
		         				int e = closeto[d][m].second;
		         				int pos_e = pos[e];
		         				int pos_f, whichcase;
		         				if(isbetween(pos_e, pos_b, pos_d))
		         				    pos_f = (pos_e + 1)%N;
		         				else
		         				    pos_f = (pos_e + N-1)%N;
		         				int f = tour[pos_f];
		         				if(f != a && f != b && f != c && f != d && e != a && e != b && e != c && e != d){
		         					enhancement = distances[a][b]+distances[c][d]+distances[e][f]-distances[b][c]-distances[d][e]-distances[f][a];
		         					if(enhancement > 0){
		         						//printf("3-Opt (%d %d) (%d %d) (%d %d)\n", a, b, c, d, e, f);
		         						if(neighbor[a][0] == b) neighbor[a][0] = f;
		         						else neighbor[a][1] = f;
		         						if(neighbor[b][0] == a) neighbor[b][0] = c;
		         						else neighbor[b][1] = c;
		         						if(neighbor[c][0] == d) neighbor[c][0] = b;
		         						else neighbor[c][1] = b;
		         						if(neighbor[d][0] == c) neighbor[d][0] = e;
		         						else neighbor[d][1] = e;
		         						if(neighbor[e][0] == f) neighbor[e][0] = d;
		         						else neighbor[e][1] = d;
		         						if(neighbor[f][0] == e) neighbor[f][0] = a;
		         						else neighbor[f][1] = a;
		         						recompute(N, tour, pos, neighbor);
		         						length_tour -= enhancement;
		         						cont = 1;
		         						mustbreak = true;
		         						break;
		         					}
		         				}
		         			}
		         			if(mustbreak) break;
		         		}

		         	}
		         }
		     }

		 }

		 if(length_tour < bestlength_tour){
		 	bestlength_tour = length_tour;
		 	for(int i = 0; i < N; i++) besttour[i] = tour[i];
		 }

    }


    if (print) {
        for (k = 0; k < N; k++){
            printf("%d\n", besttour[k]);
        }
    }
    if (print_length)
        printf("Length tour after enhance 3 = %ld\n", bestlength_tour);

    if (print_time) {
        clock_t end = clock();
        long elapsed = timediff(start, end);
        printf("Enhance 3 took %ld microseconds\n", elapsed);
    }

    return make_pair(bestlength_tour, besttour);
}

// Lin-Kernighan
/*pair <long int, vector <pair <int, int> > > lin_kernighan_sub(int N, long int** distances, long int length_tour, vector <pair <int, int> > neighbors){
    clock_t start = clock();

    int alternatives_tried[1000];
    int x_taken[1000][1000];
    int y_taken[1000][1000];
    long int new_length_tour = 0;

    // choose t1
    int first = 0;
    while(first < N){
        int i = 1;
        vector <int> t;
        t.push_back(-1); // to make it really start at index 1
        t.push_back(tour[first]);
        // choose x1 = (t1, t2) in T
        int second;
        if (alternatives_tried[first] == 0) {
            second = (first+1) % N;
            t.push_back(tour[second]);
        } else {
            second = (first-1) % N;
            t.push_back(tour[second]);
        }
        x_taken[t[1]][t[2]] = 1;
        x_taken[t[2]][t[1]] = 1;
        new_length_tour -= distances[t[1]][t[2]];

        // choose y1 = (t2, t3) not in T
        int G1 = 0;
        for (int t3=0; t3<N && G1 <=0; t3++){
            if (t3 != tour[(second+1) % N] && t3 != tour[(second-1) % N] && x_taken[t[2]][t3]==0 && y_taken[t[2]][t3]==0){
                G1 = distances[t[1]][t[2]]-distances[t[2]][t3];
            }
        }
        if (G1 > 0){
            t.push_back(t3);
            y_taken[t[2]][t[3]] = 1;
            y_taken[t[3]][t[2]] = 1;
            new_length_tour += distances[t[2][t[3]];
            int G = 1;
            while(G > 0){
                i++;
                // choose x_i = (t_{2i-1}, t_{2i}) in T
                t3 = t[2*i-1];
                bool cont = true;
                for (int n = 0; n < 2 && cont; n++){
                    t4 = neighbors[t3][n];
                    // if x_i â‰  y_s for all s < i and if t_{2i} is joined to t_1, the resulting configuration is a tour ??????????????
                    // !!!! allow unfeasible choice for i = 2
                    if (y_taken[t[2*i-1]][t[2*i]]==0 && ){
                        if (different){
                            cont = false;
                            t.push_back(t4);
                            taken[t[3]][t[4]] = 1;
                            taken[t[4]][t[3]] = 1;
                            new_length_tour += distances[t[3][t[4]];
                            // If T' is a better tour than T, let T = T' and redo the heuristic
                            if (new_length_tour < length_tour){
                                // remove the edges x
                                for (int j = 1; j < i; j++){
                                    if (neighbors[t[2*j-1]][0] == t[2*j])
                                        neighbors[t[2*j-1]][0] = -1;
                                    else
                                        neighbors[t[2*j-1]][1] = -1;
                                    if (neighbors[t[2*j]][0] == t[2*j-1])
                                        neighbors[t[2*j]][0] = -1;
                                    else
                                        neighbors[t[2*j]][1] = -1;
                                }
                                // add the edges y
                                for (int j = 1; j < i; j++){
                                    if (neighbors[t[2*j]][0] == -1
                                        neighbors[t[2*j]][0] = t[2*j+1];
                                    else
                                        neighbors[t[2*j]][1] = t[2*j+1];
                                    if (neighbors[t[2*j+1]][0] == -1
                                        neighbors[t[2*j+1]][0] = t[2*j];
                                    else
                                        neighbors[t[2*j+1]][1] = t[2*j];
                                }
                                return lin_kernighan_sub(N, distances, new_length_tour, neighbors);
                            }
                            // choose y_i = (t_{2i), t_{2i+1}) not in T
                            int G;
                            for (int t5 = 0; t5 < N; t5++){
                                if (x_taken[t4][t5] == 0 && y_taken[t4][t5] == 0) {
                                    G = 0;
                                    for (int j = 0; j < i; j++){
                                        G += distances[t[2*i-1]][t[2*i]] - distances[t[2*i]][t[2*i+1]];
                                    }
                                }
                            }
                        }
                    }

                }
            }
        }
        // choose untried alternative for t1
        if (alternatives_tried[first] == 0)
            alternatives_tried[first] = 1;
        else
            first++;
    }

    if (print_length)
        printf("Length tour after Lin-Kernighan = %ld\n", length_tour);

    if (print_time) {
        clock_t end = clock();
        long elapsed = timediff(start, end);
        printf("Lin-Kernighan took %ld microseconds\n", elapsed);
    }
    return make_pair(length_tour, new_neighbors);
}

pair<long int, int *> lin_kernighan(int N, long int** distances, long int length_tour, int* tour){
    vector <pair <int, int> > neighbors;
    neighbors.resize(N);
    for (i=0; i<N; i++){
        neighbors[i][0] = tour[(i-1) % N];
        neighbors[i][1] = tour[(i+1) % N];
    }
    pair <long int, vector <pair <int, int> > new_neighbors = lin_kernighan_sub(N, distances, length_tour, neighbors);

    return make_pair(new_neighbors[0], neighbors_to_tour(N, new_neighbors[1]));
}*/

// Nearest neighbor tour
pair<long int, int *> nearest_neighbor(int N, double** points, long int** distances, int nattempts) {
    vector <pair<long int, int *> > tours;
    tours.resize(nattempts);
    long int best_length = 1000000000;
    int best_tour = -1;
    for (int a = 0; a < nattempts; a++){
        if (print)
            printf("greedy_tour\n");
        clock_t start = clock();
        int *tour = (int *)malloc(N*sizeof(int));
        int *used = (int *)calloc(N, sizeof(int));
        long int length_tour = 0;
        int first = rand() % N;
        tour[0] = first;
        if (print)
            printf("Point 0 of the tour: %d\n", first);
        used[first] = 1;
        int best,i,j;
        for (i=1; i<N; i++){
           best = -1;
           for (j=0; j<N; j++) {
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
            printf("Length tour for nearest_neighbour = %ld\n", length_tour);
        if (print_time) {
            clock_t end = clock();
            long elapsed = timediff(start, end);
            printf("Nearest neighbor took %ld microseconds\n", elapsed);
        }
        if(version == 3) tours[a] = enhance3(N, distances, length_tour, tour);
        else if(version == 2) tours[a] = enhance2(N, distances, length_tour, tour);
        else if(version==1) tours[a] = enhance(N, distances, length_tour, tour);
        else return make_pair(length_tour, tour);
        if (best_tour == -1 || tours[a].first < best_length) {
            best_tour = a;
            best_length = tours[a].first;
        }
        if (print)
            printf("Best length = %ld, best_tour = %d\n", best_length, best_tour);
    }
    return tours[best_tour];
}

// Greedy tour
pair<long int, int *> greedy(int N, double** points, long int** distances) {
    if (print)
        printf("greedy_tour\n");
    clock_t start_greedy = clock();
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
    int *tour =neighbors_to_tour(N, neighbor);

    if (print_length)
        printf("Length tour for greedy = %ld\n", length_tour);
    if (print_time) {
        clock_t end_greedy = clock();
        long elapsed = timediff(start_greedy, end_greedy);
        if (print_time)
            printf("Greedy took %ld microseconds\n", elapsed);
    }
    if(version == 3) return enhance3(N, distances, length_tour, tour);
    else if(version == 2) return enhance2(N, distances, length_tour, tour);
    else if (version == 1) return enhance(N, distances, length_tour, tour);
    else return make_pair(length_tour, tour);
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
        printf("Length tour for CW = %ld\n", length_tour);
    if (print_time) {
        clock_t end_greedy = clock();
        long elapsed = timediff(start_greedy, end_greedy);
        if (print_time)
            printf("Clarke-Wright took %ld microseconds\n", elapsed);
    }
    if(version == 3) return enhance3(N, distances, length_tour, tour);
    else if(version == 2) return enhance2(N, distances, length_tour, tour);
    else if (version == 1) return enhance(N, distances, length_tour, tour);
    else return make_pair(length_tour, tour);
}

int count_reachable(int v, int *visited, vector <int> neighbor[1000]){
  // Mark the current node as visited
  visited[v] = 1;
  int count = 1;

  // Recur for all vertices adjacent to this vertex
  for (int i = 0; i < neighbor[v].size(); i++){
      if (neighbor[v][i] != -1 && visited[neighbor[v][i]]==0)
          count += count_reachable(neighbor[v][i], visited, neighbor);
  }
  return count;
}

void remove_edge(vector <int> neighbor[1000], int current, int next){
    bool not_found = true;
    for (int i = 0; i < neighbor[current].size() && not_found; i++){
        if (neighbor[current][i]==next) {
            neighbor[current][i] = -1;
            not_found = false;
        }
    }
    not_found = true;
    for (int i = 0; i < neighbor[next].size() && not_found; i++){
        if (neighbor[next][i]==current) {
            neighbor[next][i] = -1;
            not_found = false;
        }
    }
}

// The function to check if edge current-next can be considered as next edge in the Euler tour
bool valid_edge(int current, int next, vector <int> neighbor[1000]) {
  // The edge current-next is valid in one of the following two cases:
  // 1) If next is the only adjacent vertex of current
  int count = 0;
  for (int i = 0; i < neighbor[current].size(); i++){
     if (neighbor[current][i] != -1)
        count++;
  }
  if (count == 1)
    return true;

  // 2) If there are multiple adjacents, then current-next is not a bridge
  // count of vertices reachable from current
  int *visited = (int *)calloc(1000, sizeof(int));
  //int visited[1000];
  //memset(visited, false, 1000);
  int count1 = count_reachable(current, visited, neighbor);
  // Remove edge (current, next) and after removing the edge, count
  // vertices reachable from u
  remove_edge(neighbor, current, next);
  //memset(visited, false, 1000);
  int *visited2 = (int *)calloc(1000, sizeof(int));
  int count2 = count_reachable(current, visited2, neighbor);
  //Add the edge back to the graph
  neighbor[current].push_back(next);
  neighbor[next].push_back(current);
  //If count1 is greater, then edge (current, next) is a bridge
  if (count1 > count2)
    return false;
  else
    return true;
}

//Christofides
pair<long int,int *> christofides(int N, double** points, long int** distances){
	if(print)
		printf("christofides\n");
	clock_t start = clock();
	int *tour = (int *)malloc(N*sizeof(int));
	int *degrees = (int *)calloc(N, sizeof(int));
	long int length_tour = 0;
    int i,j,k,r;

    // Find the minimum spanning tree
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
            if (print)
                printf("edge taken : (%d, %d)\n", i, j);
        }
    }

    // Find vertices with odd degree
    vector <int> odd_vertices;
    for (i=0; i < N; i++){
        if ((degrees[i] % 2) == 1) {
            odd_vertices.push_back(i);
            if (print)
                printf("odd vertice : %d\n", i);
        }
    }

    // Minimal matching
    vector<pair<int, pair<int, int> > > v_odd;
    L = (odd_vertices.size()*(odd_vertices.size()-1))/2;
    v_odd.resize(L);
    if (print){
        printf("number of odd vertices : %lu\n", odd_vertices.size());
    }
    k = 0;
    for (i=0; i < odd_vertices.size(); i++){
        for (j=i+1; j < odd_vertices.size(); j++){
            v_odd[k].first = distances[odd_vertices[i]][odd_vertices[j]];
            v_odd[k].second.first = odd_vertices[i];
            v_odd[k].second.second = odd_vertices[j];
            if (print)
                printf("add to v_odd : (%d, %d, %d)\n", v_odd[k].first, v_odd[k].second.first, v_odd[k].second.second);
            k++;
        }
    }
    sort(v_odd.begin(),v_odd.end());
    int *odd_used = (int *)calloc(N, sizeof(int));
    for (i=0; i < L; i++) {
        if (odd_used[v_odd[i].second.first] == 0 && odd_used[v_odd[i].second.second] == 0) {
            // We take the edge
            odd_used[v_odd[i].second.first] = 1;
            odd_used[v_odd[i].second.second] = 1;
            neighbor[v_odd[i].second.first].push_back(v_odd[i].second.second);
            neighbor[v_odd[i].second.second].push_back(v_odd[i].second.first);
            if (print)
                printf("matching : (%d, %d)\n", v_odd[i].second.first, v_odd[i].second.second);
        }
    }

    // Euler tour
    if (print)
        printf("Euler Tour :\n");
    vector <int> euler_tour;
    bool cont = true;
    int current = 0;
    euler_tour.push_back(0);
    if (print)
        printf("%d\n", 0);
    while(cont){
        cont = false;
        for (i = 0; i != neighbor[current].size() && !cont; i++){
              int next = neighbor[current][i];
              // If edge current-next is not removed and it is a valid next edge
              if (next != -1 && valid_edge(current, next, neighbor)){
                  // remove edge
                  remove_edge(neighbor, current, next);
                  // add edge to euler tour
                  euler_tour.push_back(next);
                  if (print)
                    printf("%d\n", next);
                  cont = true;
                  current = next;
              }
          }
    }

    // Create the final tour
    if (print)
        printf("Final tour :\n");
    k = 0;
    int *seen = (int *)calloc(1000, sizeof(int));
    for (i = 0; i < euler_tour.size(); i++){
        if (seen[euler_tour[i]] == 0){
            tour[k] = euler_tour[i];
            if (print)
                printf("%d\n", tour[k]);
            seen[euler_tour[i]] = 1;
            if (k > 0)
                length_tour += distances[tour[k-1]][tour[k]];
            k++;
        }
    }
    length_tour += distances[tour[N-1]][tour[0]];
    if (print_length)
        printf("Length tour for Christofides = %ld\n", length_tour);

    clock_t end = clock();
    long elapsed = timediff(start, end);
    if (print_time)
        printf("Christophides took %ld microseconds\n", elapsed);

    if(version == 3) return enhance3(N, distances, length_tour, tour);
    else if(version == 2) return enhance2(N, distances, length_tour, tour);
    else if (version == 1) return enhance(N, distances, length_tour, tour);
    else return make_pair(length_tour, tour);
}

int main(int argc, char *argv[]) {
    /* Get arguments */
    int N, i, j, k,r;
    scanf("%d", &N);
    N = 1000;
    double **points = init_matrix_double(N, 2);
    int K = min(20, N-1);
    int count1=0,count2=0,count3=0,count4=0;
for(r=0;r<2;r++){

    int x = time(NULL) % 1000;
    srand(x);
    printf("Seed = %d\n", x);
    for(i = 0; i < N; i++) {
        points[i][0] = rand() % 1000;
        points[i][1] = rand() % 1000;
    }



    /* srand(time(NULL));

    for (i=0; i<N; i++) {
        scanf("%lf", points[i]);
        scanf("%lf", points[i]+1);
        if (print)
            printf("Node %d: (%f, %f)\n", i, points[i][0], points[i][1]);
    }
	*/

    /* Compute distances between all the points */
    long **distances = compute_distances(N, points);

    clock_t start = clock();

    for(i = 0; i < N; i++){
    	closeto[i].resize(N-1);
    	k = 0;
    	for(j = 0; j < N; j++){
    		if(i != j){
    			closeto[i][k] = make_pair(distances[i][j], j);
    			k++;
    		}
    	}

    	nth_element(closeto[i].begin(), closeto[i].begin() + K, closeto[i].end());
    	closeto[i].resize(K);
    	sort(closeto[i].begin(), closeto[i].end());
    }

    if (print_time) {
        clock_t end = clock();
        long elapsed = timediff(start, end);
        printf("Initialization of enhance 2 took %ld microseconds\n", elapsed);
    }

    /* Find tours */
    pair <long int, int*> tour_nn = nearest_neighbor(N, points, distances, min(6,N));

    pair <long int, int*> tour_greedy = greedy(N, points, distances);

    /*long int test_length = 1000000000;
    pair <long int, int*> tour_cw;
    for(i = 0; i < 1; i++){
        pair <long int, int*> tour_cw_test = clarke_wright(N, points, distances);
        if(tour_cw_test.first < test_length) {
            tour_cw = tour_cw_test;
            test_length = tour_cw_test.first;
        }
    }*/

    pair <long int, int*> tour_christofides = christofides(N, points, distances);

    /* Find the best one */
    //long int best_length = min(tour_nn.first, min(tour_greedy.first, min(tour_cw.first, tour_christofides.first)));
    long int best_length = min(tour_nn.first, min(tour_greedy.first, tour_christofides.first));
    int *best_tour;
    if (best_length == tour_nn.first){
        best_tour = tour_nn.second;
    	count1++;
        }
    else if (best_length == tour_greedy.first){
        best_tour = tour_greedy.second;
    	count2++;
        }
    //else if (best_length == tour_cw.first){
      //  best_tour = tour_cw.second;
      //  count3++;}
    else{
        best_tour = tour_christofides.second;
    	count4++;
        }

    /* Print the best tour */
    long int verify_length = 0;
    for (i=0; i<N; i++){
        printf("%d\n", best_tour[i]);
        verify_length += distances[best_tour[i]][best_tour[(i+1)%N]];
    }
    if (print_length)
        printf("Best length = %ld = %ld\n", best_length, verify_length);

    printf("count1=%d count2=%d count3=%d count4=%d count4=%d\n");
    return EXIT_SUCCESS;
}

}
