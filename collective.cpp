#include "collective.h"
#include "utils.h"
#include<vector>
#include<iostream>
#include<set>
#include <math.h>

using namespace std;


const int COLLECTIVE_DEBUG = 0;

/*************************** DECLARE YOUR HELPER FUNCTIONS HERE ************************/



/*************************** collective.h functions ************************/


 void HPC_Bcast(void* buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm) {
    // TODO: Implement this function using only sends and receives for communication instead of MPI_Bcast.

    int world_rank;
    MPI_Comm_rank(comm, &world_rank);

    int world_size;
    MPI_Comm_size(comm, &world_size);

    MPI_Status status;
    int d = log2(world_size);
    int flip = 1 << (d-1);
    int mask = flip - 1;

    for(int j=d-1; j>=0; j--)
    {
        if (((world_rank ^ root) & mask) == 0)
        {
            int pair = world_rank ^ flip;
            if (((world_rank ^ root) & flip) == 0)
            {
                MPI_Send(buffer, count, datatype, pair, 111, comm);
            }
            else
            {
                MPI_Recv(buffer, count, datatype, pair, 111, comm, &status);
            }
        }
        mask = mask >> 1;
        flip = flip >> 1;
    }

    //MPI_Bcast(buffer, count, datatype, root, comm);
}
    

void HPC_Prefix(const HPC_Prefix_func* prefix_func, const void *sendbuf, void *recvbuf, int count,
                MPI_Datatype datatype, MPI_Comm comm, void* wb1, void* wb2, void* wb3) {
    if (count <= 0) return;

    /* Step 1. Run user function on local data with a NULL previous prefix. */
    const void* local_last_prefix = prefix_func(NULL, sendbuf, recvbuf, count);

    // TODO: Implement the rest of this function using sends and receives for communication.

    // Initialize MPI parameters
    int world_rank;
    MPI_Comm_rank(comm, &world_rank);

    int world_size;
    MPI_Comm_size(comm, &world_size);

    MPI_Status status;

    // To handle the case of p not being a power of 2
    int d = ceil(log2(world_size));

    // Initialize tota sum to wb1
    prefix_func(NULL, local_last_prefix, wb1, 1);

    // Counter to keep track if we are updating prefix sum for the first time or not
    int c = 0;

    /* Step 2. Perform parallel prefix sum using the last local prefix sum on each processor */

    for(int j = 0; j<= d-1; j++)
    {
        // Find corresponding pair to send and recive from
        int new_rank = world_rank ^ (1<<j);

        // Do the communication only if it's a real processor 
        if (world_rank < world_size && new_rank < world_size)
        {
            MPI_Send(wb1, 1, datatype, new_rank, 111, comm);
            MPI_Recv(wb2, 1, datatype, new_rank, 111, comm, &status);

            // Handle non-commutativity 
            if (world_rank > new_rank)
            {
                prefix_func(wb2, wb1, wb1, 1);
            }
            else
            {
                prefix_func(wb1, wb2, wb1, 1);
            }

            // Initialize prefix to first element of the recieved data I*x = x 
            if (c == 0 && world_rank > new_rank)
            {
                prefix_func(NULL, wb2, wb3, 1);
                c += 1;
            }
            
            // Once prefix has been initialized keep updating it with the recieved data
            else if (c > 0 && world_rank > new_rank)
            {
                prefix_func(wb2, wb3, wb3, 1);
                c += 1;
            }

        }

    }   

    /* Step 3: Add the reult of parallel prefix sum on a processor to each of its local sum*/
    if (world_rank > 0)
        prefix_func(wb3, sendbuf, recvbuf, count);
    
}


/*************************** DEFINE YOUR HELPER FUNCTIONS HERE ************************/

