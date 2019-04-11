/// COMPILE: mpicxx mpi.cpp -o mpi 
/// RUN: mpiexec -np 4 ./mpi <num_of_points>
#include <stdio.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <time.h>
#include <unistd.h>
#include "mpi.h"
using namespace std;

int POINTS;

typedef struct
{
    float x;
    float y;
} POINT;
vector<POINT> p;
int n=0;

void GenPoints();
void SortPoint();
int CreateConvexHull(int,int);

timespec diff(timespec start, timespec end);

main(int argc, char* argv[])
{
    int number, NUMBER;
    int myid, numprocs;
    int lower_bound, interval;
    timespec time1, time2;
    POINTS = atoi(argv[1]);
    p.resize(POINTS);

    MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    if(myid==0) {
        
        GenPoints();
  
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
        SortPoint();
    }

    interval = POINTS/numprocs;
	lower_bound = myid * interval;
	
	//number = CreateConvexHull(lower_bound,interval);
    int j,Num,k,t;
    Num=n;
    for(k=3+lower_bound; k<lower_bound+interval || k<=Num; k++)
    {
        j=2;
        for(j=2; j<=k-1; j++)
        {
            if(((p[k-j+1].x-p[k-j].x)*(p[k-j].y-p[0].y)-(p[k-j].x-p[0].x)*(p[k-j+1].y-p[k-j].y))*((p[k-j+1].x-p[k-j].x)*(p[k].y-p[k-j+1].y)-(p[k].x-p[k-j+1].x)*(p[k-j+1].y-p[k-j].y))>0)
            {
                for(t=k-1; t<Num; t++)
                    p[t]=p[t+1];
                Num--;
                k=k-1;
                j=j-1;
            }
        }
    }

    MPI_Reduce(&Num, &NUMBER, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    if(myid==0) {
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2);

        printf("Points: %d \t Boundary Points: %d \t ",POINTS,NUMBER);
        printf("Time: %ld s + %.9ld ns",diff(time1,time2).tv_sec,diff(time1,time2).tv_nsec);

    }


    MPI_Finalize();
}


void GenPoints()
{
    const int maxX = 6400-100, maxY = 4800-300;
    vector<int> xx(maxX);
    for(int i=0; i<maxX; i++)
        xx[i] = i+50;

    vector<int> yy(maxY);
    for(int j=0; j<maxY; j++)
        yy[j] = j+50;

    vector<POINT> seeds(maxX*maxY);
    int s = 0;
    for(int i=0; i<maxX; i++)
    {
        for(int j=0; j<maxY; j++)
        {
            seeds[s].x = xx[i];
            seeds[s].y = yy[j];
            s++;
        }
    }

    random_shuffle ( seeds.begin(), seeds.end() );

    for(int i=0; i<POINTS; i++)
    {
        p[i].x = seeds[i].x;
        p[i].y = seeds[i].y;
    }
    n = POINTS-1;
}



void SortPoint()
{
    POINT *min,temp;
    int i,j;
    min=&p[0];
    for(i=1; i<n; i++)
        if(p[i].y<min->y)
            min=&p[i];
    temp=*min;
    *min=p[0];
    p[0]=temp;
    for(i=2; i<n; i++)
    {
        j=i-1;
        temp=p[i];
        while(j>0&&((temp.x-p[0].x)*(p[j].y-p[0].y)-(p[j].x-p[0].x)*(temp.y-p[0].y))>0)
        {
            p[j+1]=p[j];
            j--;
        }
        p[j+1]=temp;
    }
    p[n]=p[0];
}



int CreateConvexHull(int lower_bound, int interval)
{
    int j,Num,k,t;
    Num=n;
    for(k=3+lower_bound; k<lower_bound+interval || k<=Num; k++)
    {
        j=2;
        for(j=2; j<=k-1; j++)
        {
            if(( (p[k-j+1].x-p[k-j].x) * (p[k-j].y-p[0].y) - (p[k-j].x-p[0].x)*(p[k-j+1].y-p[k-j].y))*((p[k-j+1].x-p[k-j].x)*(p[k].y-p[k-j+1].y)-(p[k].x-p[k-j+1].x)*(p[k-j+1].y-p[k-j].y))>0)
            {
                for(t=k-1; t<Num; t++)
                    p[t]=p[t+1];
                Num--;
                k=k-1;
                j=j-1;
            }
        }
    }
    return Num;
}


timespec diff(timespec start, timespec end)
{
	timespec temp;
	if ((end.tv_nsec-start.tv_nsec)<0) {
		temp.tv_sec = end.tv_sec-start.tv_sec-1;
		temp.tv_nsec = 1000000000+end.tv_nsec-start.tv_nsec;
	} else {
		temp.tv_sec = end.tv_sec-start.tv_sec;
		temp.tv_nsec = end.tv_nsec-start.tv_nsec;
	}
	return temp;
}
