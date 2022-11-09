
#include <cstdlib>
#include <iostream>
#include <fstream>
#include "mpi.h"
#include <cstring>
#include <math.h>
#include <algorithm>
#include <vector>
using namespace std;

//convert 2d array into 1d array
void array_convert(int** sub_panel,int* transfer, int row, int col){
	for(int i=0;i<row;i++){
		for(int j=0;j<col;j++){
			transfer[i*col+j]=sub_panel[i][j];
		}
	}
}
//convert 1d array back to 2d array
void convert_array(int** sub_panel,int* receiver,int row,int col){
	int counter=0;
	for(int i=0;i<row;i++){
		for(int j=0;j<col;j++){
			sub_panel[i][j]=receiver[counter];
			counter++;
		}
	}
	
}
//initialize a empty 2d array with size N
int** make_panel(int m)
{
	
    int **arr = (int**)malloc(m * sizeof(int*));  
    for (int i = 0; i < m; ++i) // arr[1] +n; arr[2]+n;
        arr[i] = (int*)malloc(m * sizeof(int));
    return arr;
}
//read data from input file
void read_panel(int** panel, string filename, int size)
{
    string line;
    ifstream infile (filename);
    int i = 0;
    if (infile.is_open())
    {
		
        while (getline (infile,line))
        {	
            for (int j = 0; j < size; j++)
            {
                if (line[j] == '1'){
                    panel[i][j] = 1;
                    // cout << "ij "<< i << " " << j << endl;
                }
            }
            i++;
        }
        infile.close();
    } else {
        exit(1);
    }
    return;
}
//split bigger map into smaller ones
void sync_sub_panel_data(int** all_panel, int** sub_panel, int rows, int rowe, int cols, int cole,int size)
{
	//for blocks that the coordinate that is bigger than the side length of each block
	if(rows >= size){
		for (int ri = rows; ri <rowe; ri++){
			for (int ci = cols; ci < cole; ci++){
				sub_panel[ri-size][ci] = all_panel[ri][ci];
			}
		}
	}
	else if(cols >= size){
		for (int ri = rows; ri <rowe; ri++){
			for (int ci = cols; ci < cole; ci++){
				sub_panel[ri][ci-size] = all_panel[ri][ci];
			}
		}
		
	}else{
		for (int ri = rows; ri <rowe; ri++){
			for (int ci = cols; ci < cole; ci++){
				sub_panel[ri][ci] = all_panel[ri][ci];
			}
		}
	}
}
//read the data after each generation
void read_sub_panel_data(int** all_panel, int** sub_panel, int rows, int rowe, int cols, int cole,int size)
{
	//for blocks that the coordinate that is bigger than the side length of each block
	if(rows >= size){
		for (int ri = rows; ri <rowe; ri++){
			for (int ci = cols; ci < cole; ci++){
				all_panel[ri][ci]=sub_panel[ri-size][ci];
			}
		}
	}
	else if(cols >= size){
		for (int ri = rows; ri <rowe; ri++){
			for (int ci = cols; ci < cole; ci++){
				all_panel[ri][ci]=sub_panel[ri][ci-size];
			}
		}
		
	}else{
		for (int ri = rows; ri <rowe; ri++){
			for (int ci = cols; ci < cole; ci++){
				all_panel[ri][ci]=sub_panel[ri][ci];
			}
		}
	}  
}


int main(int argc, char *argv[])
{

    int N = atoi(argv[1]);
    int k = atoi(argv[2]);
    int m = atoi(argv[3]);
    
    int p;
    int id;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &p);  //  Get the number of processes.
    MPI_Comm_rank(MPI_COMM_WORLD, &id); //  Get the individual process ID.
    
	
    double wtime;
    if (id == 0)
    {
        wtime = MPI_Wtime();
    }
	int splitcore=sqrt(p);
    int chunk = N /  splitcore;
    if (N > chunk*p){
        chunk += 1;
    }
	//calculate the correspond x y coordinate for sub blocks
    int *start = new int[splitcore];
    int *end = new int[splitcore];
    for(int i=0;i< splitcore;i++){
		start[i]=chunk*i;
		end[i]=start[i]+chunk;
	}
	end[ splitcore-1]=N;
    
    int mysize = chunk; // border

    int** all_panel = NULL;
    vector<int**> sub_panels;
    vector<int> sub_sizes;
    if (id == 0)
    {
        string filename = argv[4];
        all_panel = make_panel(N);
        read_panel(all_panel, filename, N);
        // split to child.
        for (int i=0; i < p; i++){
            int ** sub_panel = make_panel(chunk);
            sub_panels.push_back(sub_panel);
            sub_sizes.push_back(chunk);
        }
    } 
    
    int **sub_panel = make_panel(mysize);
    int **new_sub_panel = make_panel(mysize);

    int *rows=new int[p];
    int *rowe=new int[p];
    int *cols=new int[p];
    int *cole=new int[p];
    int counter=0;
    //store coordinates into array
	for(int i=0;i<splitcore;i++){
		for(int j=0;j<splitcore;j++){
			rows[counter]=start[i];
			rowe[counter]=end[i];
			cols[counter]=start[j];
			cole[counter]=end[j];
			counter++;
			}
		}

	int transfer[mysize*mysize];
	
	int subsize=mysize*mysize;
	int **sub_pa = make_panel(mysize);
    for (int ki = 0; ki <= k; ki++)
    {
        // sync and send the sub_panel.
        if (id == 0)
        {
            // sync panel
            	
            for (int index=0; index < p; index++)
            {                
                sync_sub_panel_data(all_panel, sub_panels[index], rows[index],rowe[index],cols[index],cole[index],sub_sizes[index]);
                
                if (index != 0)
                {

					array_convert(sub_panels[index],transfer,rowe[index]-rows[index],cole[index]-cols[index]);
                    MPI_Send(transfer, mysize*mysize, MPI_INT, index, 0, MPI_COMM_WORLD);
                } else {
					
                    sub_pa = sub_panels[0];
                }
            }
        } else {
			int receiver[mysize*mysize];
            MPI_Recv(receiver, mysize*mysize, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);            
            convert_array(sub_pa,receiver,rowe[id]-rows[id],cole[id]-cols[id]);

        }
        
        // clac the sub_panel;
        
        for (int ri = 0; ri < mysize; ri++)
        {
            for (int ci = 0; ci < mysize; ci++)
            {
                int up = ci -1;
                int down = ci +1;
                int left = ri -1;
                int right = ri + 1;
                int count =0;
                //if the block is at top left corner 
                if(left <0 && up<0 ){
					 count= sub_pa[right][ci] + sub_pa[ri][down] + sub_pa[right][down] ;
				//if the block is at top right corner
				}else if(right >= mysize && up<0){
					 count= sub_pa[left][ci] + sub_pa[left][down] + sub_pa[ri][down] ;
				//if the block is at bottom right corner
				}else if(right >= mysize && down >=N){
					 count= sub_pa[left][up] + sub_pa[ri][up] + sub_pa[left][ci] ; 
				//if the block is at bottom left corner
				}else if(left <0 && down >=N){
					 count= sub_pa[ri][up] + sub_pa[right][up] + sub_pa[right][ci] ;	
			    // for blocks at left border 
				}else if(left <0){
					count= sub_pa[ri][up]  + sub_pa[ri][down] + sub_pa[right][up] +  sub_pa[right][ci] +sub_pa[right][down];
				// for blocks at top border
				}else if(up<0){
					 count= sub_pa[left][ci]  + sub_pa[right][ci] + sub_pa[left][down] +  sub_pa[ri][down] +sub_pa[right][down];
				//for blocks at right border
				}else if(right>=mysize){
					 count= sub_pa[left][up] + sub_pa[left][ci] + sub_pa[left][down] + sub_pa[ri][up] +sub_pa[ri][down];
				//for blocks at bottom border
				}else if(down>=N){
					 count= sub_pa[left][up] + sub_pa[ri][up] + sub_pa[right][up] + sub_pa[left][ci]  +sub_pa[right][ci];
				}else{
					 count = sub_pa[left][up] + sub_pa[ri][up] + sub_pa[right][up] + sub_pa[left][ci]+sub_pa[right][ci]+sub_pa[left][down]+sub_pa[ri][down]+sub_pa[right][down];
				} 
                int survive = 0;
                if (sub_pa[ri][ci])
                {
                    if  (count ==2 || count == 3)
                    {
                        survive = 1;
                    }
                } else if (count == 3)
                {
                    survive = 1;
                }
                new_sub_panel[ri][ci] = survive;
            }

            
        }
        int **new_pa = make_panel(mysize);
        int **new_map = make_panel(N);
        // send to the master.
        if (id == 0)
        {
            // recive panel
            for (int index=0; index < p; index++)
            {                
                if (index != 0)
                { 
					int new_receive[mysize*mysize];            
                    MPI_Recv(new_receive, mysize*mysize, MPI_INT, index, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);				
					convert_array(new_pa,new_receive,rowe[index]-rows[index],cole[index]-cols[index]);
                    read_sub_panel_data(all_panel, new_pa, rows[index],rowe[index], cols[index], cole[index],mysize);
                } else {
                    read_sub_panel_data(all_panel, new_sub_panel, rows[0],rowe[0], cols[0], cole[0],mysize);
                }
            }
        } else {
			int new_sub[mysize*mysize];
			array_convert(new_sub_panel,new_sub,rowe[id]-rows[id],cole[id]-cols[id]);
            MPI_Send(new_sub, mysize*mysize, MPI_INT, 0, 0, MPI_COMM_WORLD);
          //  cout << "send " << id << " to 0" << endl;

        }

        MPI_Barrier(MPI_COMM_WORLD);
        if (id == 0)
        {
			m--;
            if (m==0)
            {
				m=atoi(argv[3]);
                string filename = "test.output.N." + to_string(N) + ".k." + to_string(ki+1) + ".txt";
                ofstream myfile;
                myfile.open (filename);
                for (int ri = 0; ri < N; ri++)
                {
                    for (int ci = 0; ci < N; ci++)
                       myfile << all_panel[ri][ci];
                    myfile << endl;
                }   
                myfile.close();
            }
        }
    }

    if (id == 0)
    {
        double wtime2 = MPI_Wtime();
        cout << "  Elapsed wall clock time = " << wtime2 - wtime << " seconds.\n";                
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();    
    return 0;
}
