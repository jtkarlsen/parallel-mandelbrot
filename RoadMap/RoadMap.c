#include "graphicsScreen.h"
#include "StopWatch.h"
#include <stdio.h>
#include <mpi.h>

#define WIDTH 1500
#define HEIGHT 1500
#define SIZE WIDTH*HEIGHT

#define DATA_TAG  1
#define TERMINATE_TAG  2
#define RESULT_TAG  3

//The number of iterations decides how long we will compute before giving up
#define ITERATIONS 255

double deltaxmin, deltaymin, deltaxmax, deltaymax;
double box_x_min, box_x_max, box_y_min, box_y_max;
int zooms=10;
int crc = 0;
int rank;
double data[5];


inline double translate_x(int x) {
  return (((box_x_max-box_x_min)/WIDTH)*x)+box_x_min;
}

inline double translate_y(int y) {
  return (((box_y_max-box_y_min)/HEIGHT)*y)+box_y_min;
}


//Simple Mandelbrot divergence test
//DO NOT MODIFY THIS
inline int solve(double x, double y)
{
  double r=0.0,s=0.0;

  double next_r,next_s;
  int itt=0;

  while((r*r+s*s)<=4.0) {
    next_r=r*r-s*s+x;
    next_s=2*r*s+y;
    r=next_r; s=next_s;
    if(++itt==ITERATIONS)break;
  }

  return itt;
}


// Our "main" function
void CreateMap() {
  int x,y,color, rowcolor;

  //Loops over pixels
  for(y=0;y<HEIGHT;y++) {
	rowcolor = 0;
    for(x=0;x<WIDTH;x++) {
      //The color is determined by the number of iterations
      //More iterations means brighter color
      color = solve(translate_x(x),translate_y(y));
      crc+=color;
      rowcolor+=color;
      //printf("Slave return: %d\n", color);
    }
    //printf("Row %d is %d\n", y, rowcolor);
  }
}


int
RoadMap ()
{
  int i;

  //Sets the bounding box,
  //Note that the x_min is -1.5 instead of -2.5, since we are using a square window
  box_x_min=-1.5; box_x_max=0.5;
  box_y_min=-1.0; box_y_max=1.0;

  deltaxmin=(-0.9-box_x_min)/zooms;
  deltaxmax=(-0.65-box_x_max)/zooms;
  deltaymin=(-0.4-box_y_min)/zooms;
  deltaymax=(-0.1-box_y_max)/zooms;

  //Updates the map for every zoom level
  CreateMap();
  for(i=0;i<zooms;i++){
    box_x_min+=deltaxmin;
    box_x_max+=deltaxmax;
    box_y_min+=deltaymin;
    box_y_max+=deltaymax;
    CreateMap();
  }

  return 0;
}


void
PrepareSend(int y)
{
	data[4] = (double)y;
}


void
ChangeZoom()
{
	box_x_min+=deltaxmin;
	box_x_max+=deltaxmax;
	box_y_min+=deltaymin;
	box_y_max+=deltaymax;
	data[0] = box_x_min;
	data[1] = box_x_max;
	data[2] = box_y_min;
	data[3] = box_y_max;
}


int
PrepareSlavework()
{
	box_x_min = data[0];
	box_x_max = data[1];
	box_y_min = data[2];
	box_y_max = data[3];
	return (int)data[4];
}


int
DoSlavework(int y)
{
	int x = 0;
	int color = 0;
	int temp = 0;

	for(x=0;x<WIDTH;x++){
		temp=solve(translate_x(x),translate_y(y));
		color+=temp;
		//printf("Slavework result: %d\n", temp);
	}
	return color;
}


void
master()
{
	char buf[256];
	sw_init();
	sw_start();

	int number_processors, i, y;
	int color = 0;
	int work = HEIGHT * (zooms+1);
	int work_left = work;
	int work_returned = 0;

	MPI_Status status;
	MPI_Comm_size(MPI_COMM_WORLD, &number_processors);

	//Sets the bounding box
	box_x_min=-1.5; box_x_max=0.5;
	box_y_min=-1.0; box_y_max=1.0;

	data[0] = box_x_min;
	data[1] = box_x_max;
	data[2] = box_y_min;
	data[3] = box_y_max;

	deltaxmin=(-0.9-box_x_min)/zooms;
	deltaxmax=(-0.65-box_x_max)/zooms;
	deltaymin=(-0.4-box_y_min)/zooms;
	deltaymax=(-0.1-box_y_max)/zooms;

	for(i=1; i<number_processors; i++)
	{
		PrepareSend(i-1);
		MPI_Send(&data, 5, MPI_DOUBLE, i, DATA_TAG, MPI_COMM_WORLD);
		work_left--;
	}

	while(work_returned < work)
	{
		if (work_left > 0)
		{
			for(y=i-1;y<HEIGHT;y++)
			{
				MPI_Recv(&color, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
				work_returned++;
				crc+=color;
				//printf("Slave return: %d\n", color);
				if (work_left > 0)
				{
					PrepareSend(y);
					MPI_Send(&data, 5, MPI_DOUBLE, status.MPI_SOURCE, DATA_TAG, MPI_COMM_WORLD);
					work_left--;
				}
			}
			ChangeZoom();
			i=1;
		} else {
			MPI_Recv(&color, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			crc+=color;
			work_returned++;
		}
	}
	for(i=1; i<number_processors; i++)
	{
		MPI_Send(&data, 5, MPI_DOUBLE, i, TERMINATE_TAG, MPI_COMM_WORLD);
	}

	sw_stop();

	sw_timeString(buf);

	printf("Time taken: %s\n",buf);
	printf("CRC is %x - %d\n",crc, crc);
}


void
slave()
{
	int y = 0;
	int color = 0; //color contains the values for the row colors
	MPI_Status status;

	char buf[256];

	MPI_Recv(data, 5, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

	while(status.MPI_TAG==DATA_TAG)
	{
		y = PrepareSlavework();
		//printf("%f, %f, %f, %f\n", box_x_min, box_x_max, box_y_min, box_y_max);
		color = DoSlavework(y);
		//printf("For %d, color is %d\n", y, color);
		//printf("Row %d is %d\n", y, color);
		if (rank == 50)
		{
			sw_init();
			sw_start();
		}
		MPI_Send(&color, 1, MPI_INT, 0, RESULT_TAG, MPI_COMM_WORLD);
		MPI_Recv(data, 5, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		if (rank == 50)
		{
			sw_stop();
			sw_timeString(buf);
			printf("Time taken: %s\n",buf);
		}
		if(status.MPI_TAG==TERMINATE_TAG)
		{
			return;
		}
	}

}


int
main (int argc, char *argv[]){

    int i, number_processors;
    char buf[256];

  	MPI_Init(&argc, &argv );
  	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  	MPI_Comm_size(MPI_COMM_WORLD, &number_processors);

  	for(i = 0; i<5; i++)
  	{
  		data[i] = 0;
  	}

  	if(number_processors == 1)
  	{
  		sw_init();
  		sw_start();
  		RoadMap();
  		sw_stop();
  		sw_timeString(buf);
  		printf("Time taken: %s\n",buf);
  		printf("CRC is %x - %d\n",crc, crc);

  	} else {
  		if(rank==0) //MASTER
  		{
  			master();
  		}
  		else //FILTY SLAVE
  		{
  			slave();
  		}
  	}

  MPI_Finalize();
  return 0;

}
