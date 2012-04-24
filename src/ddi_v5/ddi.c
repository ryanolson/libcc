



void DDI_Init(int argc, char **argv)
{
        int status;
        MPI_Initialized(&status);
        if(!status) {
           MPI_Init(&argc, &argv);
        }

        status = ARMCI_Init();
        if(status) {
           
