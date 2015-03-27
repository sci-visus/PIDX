/*****************************************************
 **  PIDX Parallel I/O Library                      **
 **  Copyright (c) 2010-2014 University of Utah     **
 **  Scientific Computing and Imaging Institute     **
 **  72 S Central Campus Drive, Room 3750           **
 **  Salt Lake City, UT 84112                       **
 **                                                 **
 **  PIDX is licensed under the Creative Commons    **
 **  Attribution-NonCommercial-NoDerivatives 4.0    **
 **  International License. See LICENSE.md.         **
 **                                                 **
 **  For information about this project see:        **
 **  http://www.cedmav.com/pidx                     **
 **  or contact: pascucci@sci.utah.edu              **
 **  For support: PIDX-support@visus.net            **
 **                                                 **
 *****************************************************/ 

#if defined(BGL) || defined(BGP) || defined(BGQ)

//void identity(MPI_Comm comm, int *iotask)
void identity(MPI_Comm comm)
{
  int rank;
  int np;
  int my_name_len;
  char my_name[255];
    
  MPI_Comm_rank(comm,&rank);
  MPI_Comm_size(comm,&np);
  MPI_Get_processor_name(my_name, &my_name_len);


  MPIX_Hardware_t hw;
  MPIX_Hardware(&hw);

  /*  Get the personality  */
  Personality pers;
  char message[100];

  /* Number of MPI tasks per Pset */
  int coreId;
  int *TasksPerPset;
  int *tmp;
  int i,ierr;

  Personality personality;
  Kernel_GetPersonality(&pers, sizeof(pers));

  int numIONodes,numPsets,numNodesInPset,rankInPset;


  int myrank, numpsets, psetID, psetsize, psetrank;

  MPI_Comm_rank(comm,&myrank);
  bgq_pset_info (comm, &numpsets, &psetID, &psetsize, &psetrank);

  numIONodes = numpsets; 
  numNodesInPset = psetsize; 
  rankInPset = myrank; 


  numPsets = numpsets; 
    
  if(rank == 0) 
  {
    printf("number of IO nodes in block: %i \n",numIONodes);
    printf("number of Psets in block : %i \n",numPsets);
    printf("number of compute nodes in Pset: %i \n",numNodesInPset);
  }

  int psetNum;
  psetNum = psetID;

  /*
  if((*iotask)>0) 
  {
    printf( "%04i (%-50s %s) %i yes\n", rank, my_name, message, psetNum );
  } 
  else 
  {
  */
  printf( "%04i (%-50s %s) %i --\n", rank, my_name, message, psetNum);
  //}
  printf("MPI task %6i is rank %3i in Pset: %3i \n",rank, rankInPset,psetNum);

  
/* Determine which core on node....  I don't want to put more than one io-task per node */
  coreId = get_processor_id ();

  TasksPerPset = malloc(numPsets*sizeof(int));
  tmp = malloc(numPsets*sizeof(int));
  for(i=0;i<numPsets;i++) 
    tmp[i]=0;
  if(coreId == 0) 
    tmp[psetNum]=1;
  
  ierr = MPI_Allreduce(tmp,TasksPerPset,numPsets,MPI_INT,MPI_SUM,comm);
  
  if(rank == 0) 
  for(i=0;i<numPsets;i++) 
    printf("Pset: %3i has %3i nodes \n",i,TasksPerPset[i]);
  
  free(tmp);
  free(TasksPerPset);

}

#endif
