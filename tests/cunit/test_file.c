/*
 * Tests for the file functions PIOc_create, PIOc_open, and
 * PIOc_close.
 *
 * Ed Hartnett
 */
#include <pio.h>
#include <pio_tests.h>

/* The number of tasks this test should run on. */
#define TARGET_NTASKS 4

/* The name of this test. */
#define TEST_NAME "test_file"

#define NUM_NETCDF_FLAVORS 4
#define NDIM 3
#define X_DIM_LEN 400
#define Y_DIM_LEN 400
#define VAR_NAME "foo"
#define ATT_NAME "bar"

/* The dimension names. */
char dim_name[NDIM][NC_MAX_NAME + 1] = {"timestep", "x", "y"};

/* Length of the dimensions in the sample data. */
int dim_len[NDIM] = {NC_UNLIMITED, X_DIM_LEN, Y_DIM_LEN};

/* Define metadata for the test file. */
int
define_metadata(int ncid, int my_rank)
{
    int dimids[NDIM]; /* The dimension IDs. */
    int varid; /* The variable ID. */
    int ret;

    for (int d = 0; d < NDIM; d++)
        if ((ret = PIOc_def_dim(ncid, dim_name[d], (PIO_Offset)dim_len[d], &dimids[d])))
            ERR(ret);

    if ((ret = PIOc_def_var(ncid, VAR_NAME, NC_INT, NDIM, dimids, &varid)))
        ERR(ret);

    return PIO_NOERR;
}

/* Check the metadata in the test file. */
int
check_metadata(int ncid, int my_rank)
{
    int ndims, nvars, ngatts, unlimdimid, natts, dimid[NDIM];
    PIO_Offset len_in;
    char name_in[NC_MAX_NAME + 1];
    nc_type xtype_in;
    int ret;

    /* Check how many dims, vars, global atts there are, and the id of
     * the unlimited dimension. */
    if ((ret = PIOc_inq(ncid, &ndims, &nvars, &ngatts, &unlimdimid)))
        return ret;
    if (ndims != NDIM || nvars != 1 || ngatts != 0 || unlimdimid != 0)
        return ERR_AWFUL;

    /* Check the dimensions. */
    for (int d = 0; d < NDIM; d++)
    {
        if ((ret = PIOc_inq_dim(ncid, d, name_in, &len_in)))
            ERR(ret);
        if (len_in != dim_len[d] || strcmp(name_in, dim_name[d]))
            return ERR_AWFUL;
    }

    /* Check the variable. */
    if ((ret = PIOc_inq_var(ncid, 0, name_in, &xtype_in, &ndims, dimid, &natts)))
        ERR(ret);
    if (strcmp(name_in, VAR_NAME) || xtype_in != NC_INT || ndims != NDIM ||
        dimid[0] != 0 || dimid[1] != 1 || dimid[2] != 2 || natts != 0)
        return ERR_AWFUL;

    return PIO_NOERR;
}

/* Run Tests for PIO file operations. */
int main(int argc, char **argv)
{
    int my_rank;    /* Zero-based rank of processor. */
    int ntasks;     /* Number of processors involved in current execution. */
    int iotype;     /* Specifies the flavor of netCDF output format. */
    int niotasks;   /* Number of processors that will do IO. */
    int ioproc_stride = 1;    /* Stride in the mpi rank between io tasks. */
    int numAggregator = 0;    /* Number of the aggregator? Always 0 in this test. */
    int ioproc_start = 0;     /* Zero based rank of first processor to be used for I/O. */
    int dimids[NDIM];         /* The dimension IDs. */
    PIO_Offset elements_per_pe;  /* Array index per processing unit. */
    int iosysid;    /* The ID for the parallel I/O system. */
    int ncid;       /* The ncid of the netCDF file. */
    int varid;      /* The ID of the netCDF varable. */
    int ioid;       /* The I/O description ID. */
    PIO_Offset *compdof;  /* The decomposition mapping. */
    int ret;        /* Return code. */
    int fmt, d, d1, i;    /* Index for loops. */
    int num_flavors;      /* Number of PIO netCDF flavors in this build. */
    int flavor[NUM_FLAVORS]; /* iotypes for the supported netCDF IO flavors. */
    MPI_Comm test_comm;   /* A communicator for this test. */

    /* Initialize test. */
    if ((ret = pio_test_init(argc, argv, &my_rank, &ntasks, TARGET_NTASKS,
                             &test_comm)))
        ERR(ERR_INIT);

    /* Only do something on TARGET_NTASKS tasks. */
    if (my_rank < TARGET_NTASKS)
    {
        /* Figure out iotypes. */
        if ((ret = get_iotypes(&num_flavors, flavor)))
            ERR(ret);

        /* keep things simple - 1 iotask per MPI process */
        niotasks = TARGET_NTASKS;

        /* Initialize the PIO IO system. This specifies how
         * many and which processors are involved in I/O. */
        if ((ret = PIOc_Init_Intracomm(test_comm, niotasks, ioproc_stride,
                                       ioproc_start, PIO_REARR_SUBSET, &iosysid)))
            ERR(ret);

        /* Describe the decomposition. This is a 1-based array, so add 1! */
        elements_per_pe = X_DIM_LEN * Y_DIM_LEN / TARGET_NTASKS;
        if (!(compdof = malloc(elements_per_pe * sizeof(PIO_Offset))))
            return PIO_ENOMEM;
        for (i = 0; i < elements_per_pe; i++)
            compdof[i] = my_rank * elements_per_pe + i + 1;

        /* Create the PIO decomposition for this test. */
        printf("rank: %d Creating decomposition...\n", my_rank);
        if ((ret = PIOc_InitDecomp(iosysid, PIO_FLOAT, 2, &dim_len[1], (PIO_Offset)elements_per_pe,
                                   compdof, &ioid, NULL, NULL, NULL)))
            ERR(ret);
        free(compdof);

        /* Use PIO to create the example file in each of the four
         * available ways. */
        for (fmt = 0; fmt < num_flavors; fmt++)
        {
            char filename[NC_MAX_NAME + 1]; /* Test filename. */
            char iotype_name[NC_MAX_NAME + 1];

            /* Overwrite existing test file. */
            int mode = PIO_CLOBBER;

            /* If this is netCDF-4, add the netCDF4 flag. */
            if (flavor[fmt] == PIO_IOTYPE_NETCDF4C || flavor[fmt] == PIO_IOTYPE_NETCDF4P)
            {
                printf("%d adding NC_NETCDF4 flag\n", my_rank);
                mode |= NC_NETCDF4;
            }

            /* If this is pnetcdf or netCDF-4 parallel, add the MPIIO flag. */
            if (flavor[fmt] == PIO_IOTYPE_PNETCDF || flavor[fmt] == PIO_IOTYPE_NETCDF4P)
            {
                printf("%d adding NC_MPIIO flag\n", my_rank);
                mode |= NC_MPIIO;
            }

            /* Create a filename. */
            if ((ret = get_iotype_name(flavor[fmt], iotype_name)))
                return ret;
            sprintf(filename, "%s_%s.nc", TEST_NAME, iotype_name);

            /* Create the netCDF output file. */
            printf("rank: %d Creating sample file %s with format %d...\n",
                   my_rank, filename, flavor[fmt]);
            if ((ret = PIOc_create(iosysid, filename, mode, &ncid)))
                ERR(ret);

            /* Define the test file metadata. */
            if ((ret = define_metadata(ncid, my_rank)))
                ERR(ret);

            /* End define mode. */
            if ((ret = PIOc_enddef(ncid)))
                ERR(ret);

            /* Close the netCDF file. */
            printf("rank: %d Closing the sample data file...\n", my_rank);
            if ((ret = PIOc_closefile(ncid)))
                ERR(ret);

            /* Reopen the test file. */
            printf("rank: %d Re-opening sample file %s with format %d...\n",
                   my_rank, filename, flavor[fmt]);
            if ((ret = PIOc_open(iosysid, filename, mode, &ncid)))
                ERR(ret);

            /* Check the test file metadata. */
            if ((ret = check_metadata(ncid, my_rank)))
                ERR(ret);

            /* Close the netCDF file. */
            printf("rank: %d Closing the sample data file...\n", my_rank);
            if ((ret = PIOc_closefile(ncid)))
                ERR(ret);

            /* Put a barrier here to make output look better. */
            if ((ret = MPI_Barrier(test_comm)))
                MPIERR(ret);

        }

        /* Free the PIO decomposition. */
        printf("rank: %d Freeing PIO decomposition...\n", my_rank);
        if ((ret = PIOc_freedecomp(iosysid, ioid)))
            ERR(ret);
    } /* endif my_rank < TARGET_NTASKS */

    /* Wait for everyone to catch up. */
    printf("%d %s waiting for all processes!\n", my_rank, TEST_NAME);
    MPI_Barrier(test_comm);

    /* Finalize the MPI library. */
    printf("%d %s Finalizing...\n", my_rank, TEST_NAME);
    if ((ret = pio_test_finalize(&test_comm)))
        return ret;

    printf("%d %s SUCCESS!!\n", my_rank, TEST_NAME);

    return 0;
}