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

#include <stdio.h>
#include <getopt.h>
#include <string.h>
#include <assert.h>
#include <errno.h>
#include <stdlib.h>
#include <string.h>
#include <stdlib.h>

#include "PIDX_Visus.h"

#define PIDX_MAX_DIMENSIONS 5
#define MAX_TEMPLATE_DEPTH 6

static int parse_args(int argc, char **argv);
static void usage(void);

/* global dimensions of 3D volume */
static int* extents = 0;

/* output file */
static char output_file[512] = {0};

static int resolution = 0;
static int variable_count = 0;
static int bits_per_block, b_per_file;
char filename_template[1024];

double swap(double d)
{
   double a;
   unsigned char *dst = (unsigned char *)&a;
   unsigned char *src = (unsigned char *)&d;

   dst[0] = src[7];
   dst[1] = src[6];
   dst[2] = src[5];
   dst[3] = src[4];
   dst[4] = src[3];
   dst[5] = src[2];
   dst[6] = src[1];
   dst[7] = src[0];

   return a;
}

int main(int argc, char **argv) {
    char bitSequence[512], bitPattern[512];
    ;
    int i, j, nfiles;
    FILE* idx_file;
    int var = 0, ret = 0, file_no = 0;
    char variable_name[1024][1024];
    int *sample_per_variable;
    int bounding_box[2][5];
    block_bitmap* bitmap;
    int* file_status;
    int maxh;

   
    ret = parse_args(argc, argv);
    if (ret < 0) {
        usage();
    }

    extents = (int*) malloc(sizeof (int) * PIDX_MAX_DIMENSIONS);
    assert(extents);
    memset(extents, 0, sizeof (int) * PIDX_MAX_DIMENSIONS);

    /*
    printf("Variable Count %d\n", variable_count);
    variable_name = (char**)malloc(sizeof(char*) * variable_count);
    assert(variable_name);
    memset(variable_name, 0, sizeof(char*) * variable_count);
    for(i = 0 ; i < variable_count ; i++)
    {
        variable_name[i] =  (char*)malloc(sizeof(char) * 1024);
        assert(variable_name[i]);
        memset(variable_name[i], 0, sizeof(char) * 1024);
    }
     */

    sample_per_variable = malloc(sizeof (*sample_per_variable) * variable_count);


    idx_file = fopen(output_file, "r");
    if (!idx_file) {
        ret = -1;
    }
    fscanf(idx_file, "(version)\n6\n");
    fscanf(idx_file, "(box)\n0 %d 0 %d 0 %d 0 %d 0 %d\n", &extents[0], &extents[1], &extents[2], &extents[3], &extents[4]);
    fscanf(idx_file, "(fields)\n");
    for (var = 0; var < variable_count; var++) {
        fscanf(idx_file, "%s %d*", &variable_name[var], &sample_per_variable[var]);
        fscanf(idx_file, "float64 ");
        if (var != variable_count - 1)
            fscanf(idx_file, " + \n");
    }
    fscanf(idx_file, "\n(bits)\n%s\n", &bitSequence);
    fscanf(idx_file, "(bitsperblock)\n%d\n(blocksperfile)\n%d\n", &bits_per_block, &b_per_file);
    fscanf(idx_file, "(filename_template)\n./%s", &filename_template);


    fclose(idx_file);
    extents[0] = extents[0] + 1;
    extents[1] = extents[1] + 1;
    extents[2] = extents[2] + 1;
    extents[3] = extents[3] + 1;
    extents[4] = extents[4] + 1;
    int samples_per_block = pow(2, bits_per_block);
    
    {
        printf("(box)\n0 %d 0 %d 0 %d 0 %d 0 %d\n", extents[0], extents[1], extents[2], extents[3], extents[4]);
        printf("(fields)\n");
        for(var = 0 ; var < variable_count ; var++)
        {
            printf("%s %d*", variable_name[var], sample_per_variable[var]);
            printf("float64 ");
            if(var != variable_count - 1)
              printf(" + \n");
        }
        printf("\n(bits)\n%s\n", bitSequence);
        printf("(bitsperblock)\n%d\n(blocksperfile)\n%d\n", bits_per_block, b_per_file);
        printf("(filename_template)\n./%s", filename_template);
    }
     

    for (i = 0; i < PIDX_MAX_DIMS; i++) {
        bounding_box[0][i] = 0;
        bounding_box[1][i] = extents[i];
    }
    bitmap = (block_bitmap*) malloc(sizeof (block_bitmap));
    if (!bitmap) {
        fprintf(stderr, "[File : %s] [Line : %d] bitmap\n", __FILE__, __LINE__);
    }
    memset(bitmap, 0, sizeof (block_bitmap));

    PointND extentsP;
    extentsP.x = extents[0];
    extentsP.y = extents[1];
    extentsP.z = extents[2];
    extentsP.u = extents[3];
    extentsP.v = extents[4];
    GuessBitmaskPattern(bitSequence, extentsP);
    maxh = strlen(bitSequence);


    printf("\nTotal Levels [%d] Bit Sequence [%s]\n", maxh, bitSequence);

    for (i = 0; i <= maxh; i++)
        bitPattern[i] = RegExBitmaskBit(bitSequence, i);

    createBlockBitmap(bounding_box, b_per_file, bits_per_block, maxh, resolution, bitPattern, bitmap);

    //maximum number of files possible
    nfiles = (getPowerOf2(extents[0]) * getPowerOf2(extents[1]) * getPowerOf2(extents[2]) * getPowerOf2(extents[3]) * getPowerOf2(extents[4])) / ((long long) pow(2, bits_per_block) * (long long) b_per_file);
    if ((getPowerOf2(extents[0]) * getPowerOf2(extents[1]) * getPowerOf2(extents[2]) * getPowerOf2(extents[3]) * getPowerOf2(extents[4])) % ((long long) pow(2, bits_per_block) * (long long) b_per_file))
        nfiles++;
    assert(nfiles != 0);


    file_status = (int*) malloc(nfiles * sizeof (int));
    if (!file_status) {
        fprintf(stderr, "[File : %s] [Line : %d] file_status\n", __FILE__, __LINE__);
    }
    memset(file_status, 0, nfiles * sizeof (int));
    file_status[0] = 1;
    for (i = 1; i < bitmap->levels; i++) {
        for (j = 0; j < bitmap->hz_block_count_array[i]; j++) {
            file_no = bitmap->hz_block_number_array[i][j] / b_per_file;
            file_status[file_no] = 1;
        }
    }
    long long miss_count = 0;
    int fd;
    uint32_t* binheader;
    int binheader_count;
    long long *ZYX, lost_element_count = 0, element_count = 0;
    binheader_count = 10 + 10 * b_per_file * variable_count;
    binheader = (uint32_t*) malloc(sizeof (*binheader)*(binheader_count));
    if (!binheader) {
        fprintf(stderr, "[File : %s] [Line : %d] binheader\n", __FILE__, __LINE__);
    }

    
    //nfiles = 1;
    
    int nExistingFiles = 0;

    if (bits_per_block + (unsigned int) floor(log2(b_per_file)) + 1 >= maxh)
	nExistingFiles = 1;
    else 
    {
	for (i = 0; i < nfiles; i++)
	    if (file_status[i] == 1)
		nExistingFiles++;
    }
    printf("Number of Files = %d (%d)\n", nExistingFiles, nfiles);
    for (i = 0; i < /*nfiles*/nExistingFiles; i++) {
        if (file_status[i] == 1) {
            char bin_file[PATH_MAX];
            ret = generate_file_name(i, bin_file, PATH_MAX);
            if (ret == -1) {
                fprintf(stderr, "[File : %s] [Line : %d] generate_file_name\n", __FILE__, __LINE__);
            }
            printf("File Name : %s\n", bin_file);
            fd = open(bin_file, O_RDONLY);
            if (fd < 0) {
                fprintf(stderr, "[File : %s] [Line : %d] open\n", __FILE__, __LINE__);
            }

            ret = read(fd, binheader, (sizeof (*binheader) * binheader_count));
            if (ret < 0) {
                fprintf(stderr, "[File : %s] [Line : %d] read\n", __FILE__, __LINE__);
            }
            if (ret < (int) (sizeof (*binheader) * binheader_count)) {
                fprintf(stderr, "[File : %s] [Line : %d] read : size\n", __FILE__, __LINE__);
            }

            int bpf = 0, bpf_offset = 0, offsetted_bpf = 0;
            double* data_buffer = NULL;
            size_t data_size;
            off_t data_offset;


            for (var = 0; var < variable_count; var++) {
                printf("\n");
                for (bpf = 0; bpf < b_per_file; bpf++) {
                    if (is_block_present((bpf + (b_per_file * i)), bitmap)) {
                        //printf("X");

                        //bpf_offset = find_block_negative_offset(bits_per_block, b_per_file,  (bpf+(b_per_file * i)), getLevelFromBlock((bpf+(b_per_file * i)), bits_per_block), bitmap);
                        //offsetted_bpf = bpf + bpf_offset;


                        data_offset = /*32768 + bpf * 32768 * 8 * sample_per_variable[var];*/ ntohl(binheader[(bpf + var * b_per_file)*10 + 12]);
                        data_size = /*32768*8*sample_per_variable[var];*/ ntohl(binheader[(bpf + var * b_per_file)*10 + 14]);

                        //printf("[%d : %d] [%d] Present Block : %d[%d] : %d %d\n", var, sample_per_variable[var], sizeof (long long), bpf, b_per_file, data_offset, data_size);

                        data_buffer = (double*) malloc(data_size);
                        assert(data_buffer != NULL);
                        memset(data_buffer, 0, data_size);

                        ret = pread(fd, data_buffer, data_size, data_offset);
                        printf("[%d] :: [%d]\n", ret, data_size);
			 assert(ret == data_size);
                        long long hz_index, hz_val;
                        int s;
                        //printf("SIZE : %d\n", data_size/(sizeof(double) * sample_per_variable[var]));
                        for (hz_val = 0; hz_val < data_size/(sizeof(double) * sample_per_variable[var]); hz_val++) {
                            hz_index = (b_per_file * i * samples_per_block) + (bpf * samples_per_block) + hz_val;
                            ZYX = Hz_to_xyz(bitPattern, maxh - 1, hz_index);

                            
                            if ((int) data_buffer[hz_val * sample_per_variable[var] + 0] == 0) {
                                miss_count++;
                                //printf("ALL ZEROS\n");
                                //continue;
                            }
                            //printf("XXXXXXXX\n");
                            int check_bit = 1, s = 0;

                            for (s = 0; s < sample_per_variable[var]; s++)
                                check_bit = check_bit && ((int) data_buffer[hz_val * sample_per_variable[var] + s] == s + var + 100 + (extents[0] * extents[1] * ZYX[2])+(extents[0]*(ZYX[1])) + ZYX[0]);
                                //check_bit = check_bit && ((int) swap(data_buffer[hz_val * sample_per_variable[var] + s]) == s + 100 + (extents[0] * extents[1] * ZYX[2])+(extents[0]*(ZYX[1])) + ZYX[0]);

                            if (check_bit == 0) {
                                lost_element_count++;
//                                 printf("[%d] X [%d] : [%lld] : [%lld] : %lld %lld %lld : S%d : %d : %lld\n", var, bpf, hz_val, hz_index, ZYX[0], ZYX[1], ZYX[2], sample_per_variable[var], (int) data_buffer[hz_val * sample_per_variable[var] + 0], (long long) 100 + (extents[0] * extents[1] * ZYX[2])+(extents[0]*(ZYX[1])) + ZYX[0]);
                                //printf("[%d] X [%d] : [%lld] : [%lld] : %lld %lld %lld : S%d : %d : %lld\n", var, bpf, hz_val, hz_index, ZYX[0], ZYX[1], ZYX[2], sample_per_variable[var], (int) swap(data_buffer[hz_val * sample_per_variable[var] + 0]), (long long) 100 + (extents[0] * extents[1] * ZYX[2])+(extents[0]*(ZYX[1])) + ZYX[0]);
                            } else {
//                                printf("Y [%d] : [%lld] : [%lld] : %lld %lld %lld : S%d : %lld : %lld\n", bpf, hz_val, hz_index, ZYX[0], ZYX[1], ZYX[2], sample_per_variable[var], (long long)(data_buffer[hz_val * sample_per_variable[var] + 0]), (long long)100 + (extents[0] * extents[1] * ZYX[2])+(extents[0]*(ZYX[1])) + ZYX[0]);
                                element_count++;
                            }

                            free(ZYX);
                        }
                        free(data_buffer);
                        data_buffer = 0;
                    }
                }
            }

        }
    }
    int total_samples = 0;
    for (var = 0; var < variable_count; var++)
        total_samples = total_samples + sample_per_variable[var];


    //  printf("Element Count : Lost Element Count :: %d : %d : %lld\n", element_count, lost_element_count, miss_count);
    printf("%lld : %lld\n", (long long) (element_count), (long long) extents[0] * extents[1] * extents[2] * extents[3] * extents[4] * variable_count);
    assert(element_count == (long long) extents[0] * extents[1] * extents[2] * extents[3] * extents[4] * variable_count);
    free(sample_per_variable);
    sample_per_variable = 0;


}

static int parse_args(int argc, char **argv) {
    char flags[] = "v:f:r:";
    int one_opt = 0, i = 0;

    while ((one_opt = getopt(argc, argv, flags)) != EOF) {
        /* postpone error checking for after while loop */
        switch (one_opt) {
            case('v'):
                sscanf(optarg, "%d", &variable_count);
                break;
            case('f'):
                sprintf(output_file, "%s", optarg);
                break;
	    case('r'):
                sscanf(optarg, "%d", &resolution);
                break;
            case('?'):
                return (-1);
        }
    }
    /* need positive dimensions */

    return (0);
}

/* prints usage instructions */
static void usage(void) {
    printf("Usage: idx-verify -f Filename.idx\n");
    printf("  -f: IDX Filename\n");
    printf("\n");
    return;
}

int generate_file_name(int file_number, char* filename, int maxlen) {
    long long address = 0;
    unsigned int segs[MAX_TEMPLATE_DEPTH] = {0};
    int seg_count = 0;
    char* pos;
    int digit;
    int ret;

    // determine the first HZ address for the file in question 
    address = file_number * b_per_file;

    // walk backwards through the file name template to find places where we need to substitute strings


    for (pos = &filename_template[strlen(filename_template) - 1];
            pos != filename_template;
            pos--) {
        // be careful not to look past the end of the array 
        if (pos - filename_template > (strlen(filename_template) - 3))
            continue;

        if (pos[0] == '%' && pos[1] == '0' && pos[3] == 'x') {
            // TODO: for now we have a hard coded max depth 
            if (seg_count >= MAX_TEMPLATE_DEPTH) {
                fprintf(stderr, "Error: generate_filename() function can't handle this template yet: %s\n", filename_template);
                return (-1);
            }

            // found an occurance of %0 in the template; check the next character to see how many bits to use here 

            switch (pos[2]) {
                case '1':
                    segs[seg_count] += address & 0xf;
                    address = address >> 4;
                    break;
                case '2':
                    segs[seg_count] += address & 0xff;
                    address = address >> 8;
                    break;
                case '3':
                    segs[seg_count] += address & 0xfff;
                    address = address >> 12;
                    break;
                case '4':
                    segs[seg_count] += address & 0xffff;
                    address = address >> 16;
                    break;
                case '5':
                    segs[seg_count] += address & 0xfffff;
                    address = address >> 20;
                    break;
                default:
                    // TODO: generalize this to work for any value 
                    fprintf(stderr, "Error: generate_filename() function can't handle this template yet: %s\n", filename_template);
                    return (-1);
            }
            seg_count++;
        }
    }
    switch (seg_count) {
        case 0:
            ret = strlen(filename_template);
            if (ret < maxlen) {
                strcpy(filename, filename_template);
            }
            break;
        case 1:
            ret = snprintf(filename, maxlen, filename_template, segs[0]);
            break;
        case 2:
            ret = snprintf(filename, maxlen, filename_template,
                    segs[1], segs[0]);
            break;
        case 3:
            ret = snprintf(filename, maxlen, filename_template,
                    segs[2], segs[1], segs[0]);
            break;
        case 4:
            ret = snprintf(filename, maxlen, filename_template,
                    segs[3], segs[2], segs[1], segs[0]);
            break;
        case 5:
            ret = snprintf(filename, maxlen, filename_template,
                    segs[4], segs[3], segs[2], segs[1], segs[0]);
            break;
        case 6:
            ret = snprintf(filename, maxlen, filename_template,
                    segs[5], segs[4], segs[3], segs[2], segs[1], segs[0]);
            break;
        default:
            // TODO: generalize this 
            fprintf(stderr, "Error: generate_filename() function can't handle this template yet: %s\n", filename_template);
            return (-1);
            break;
    }
    // make sure that the resulting string fit into the buffer ok 
    if (ret >= maxlen - 1) {
        fprintf(stderr, "Error: filename too short in generate_filename()\n");
        return (-1);
    }
    //printf("filename: %s\n", filename);
    return (0);
}
