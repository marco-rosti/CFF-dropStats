!stab randomly in order to find drops, if exhaustive search takes too long
!#define RAND_SEARCH

!drops are contiguous even if they only touch at a cell corner 
!#define DIAGONALS_CONNECTED

!share fields between processors on the same node
#define MPI_WINDOW

!make a joint pdf of drop curvature vs diameter, slow due to allreduce
!#define CURVATURE_PDF
