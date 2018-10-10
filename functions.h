#ifndef FUNCTIONS_H
#define FUNCTIONS_H


#endif // FUNCTIONS_H

CELL **open_raster_C_variable(char *input_name){

    struct Cell_head cellhd;	/* it stores region information,*/
    char    *raster_mapset;		/* mapset name */
    void    *inrast = NULL;            /* input buffer */
    int     row, col;
    int     nrows, ncols;
    int     infd;
    RASTER_MAP_TYPE data_type;	/* type of the map (CELL/DCELL/...) */
    CELL   **output;

    /* returns NULL if the map was not found in any mapset,
     * mapset name otherwise */
    raster_mapset = (char *) G_find_raster2(input_name, "");
    if (raster_mapset == NULL)
        G_fatal_error(_("Raster map <%s> not found"), input_name);

    /* determine the inputmap type (CELL/FCELL/DCELL) */
    data_type = Rast_map_type(input_name, raster_mapset);

    /* Allocate input buffer */
    inrast = Rast_allocate_buf(data_type);

    /* Allocate output buffer, use input map data_type */
    nrows = Rast_window_rows();
    ncols = Rast_window_cols();

    /* allocating memory for the output variable */
    output = (CELL **)G_malloc (sizeof(CELL *)*nrows);

    /* Rast_open_old - returns file destriptor (>0) */
    infd = Rast_open_old(input_name, raster_mapset);

    /* controlling, if we can open input raster */
    Rast_get_cellhd(input_name, raster_mapset, &cellhd);

    G_debug(3, "number of rows %d", cellhd.rows);

    /* for each row */
    for (row = 0; row < nrows; row++) {

        G_percent(row, nrows, 2);

        /* allocating memory for the output variable */
        output[row] = (CELL *)G_malloc (sizeof(CELL)*ncols);

    /* read input map */
    Rast_get_row(infd, inrast, row, data_type);

    // process the data 
       for (col = 0; col < ncols; col++) {
        
        output[row][col] = ((CELL *) inrast)[col];
       }
    }

    /* memory cleanup */
    G_free(inrast);

    /* closing raster maps */
    Rast_close(infd);

    return (output);

}

DCELL **open_raster_D_variable(char *input_name){

    struct Cell_head cellhd;	/* it stores region information,*/
    char    *raster_mapset;		/* mapset name */
    void    *inrast = NULL;            /* input buffer */
    int     row, col;
    int     nrows, ncols;
    int     infd;
    RASTER_MAP_TYPE data_type;	/* type of the map (CELL/DCELL/...) */
    DCELL   **output;

    /* returns NULL if the map was not found in any mapset,
     * mapset name otherwise */
    raster_mapset = (char *) G_find_raster2(input_name, "");
    if (raster_mapset == NULL)
        G_fatal_error(_("Raster map <%s> not found"), input_name);

    /* determine the inputmap type (CELL/FCELL/DCELL) */
    data_type = Rast_map_type(input_name, raster_mapset);

    /* Allocate input buffer */
    inrast = Rast_allocate_buf(data_type);

    /* Allocate output buffer, use input map data_type */
    nrows = Rast_window_rows();
    ncols = Rast_window_cols();

    /* allocating memory for the output variable */
    output = (DCELL **)G_malloc (sizeof(DCELL *)*nrows);

    /* Rast_open_old - returns file destriptor (>0) */
    infd = Rast_open_old(input_name, raster_mapset);

    /* controlling, if we can open input raster */
    Rast_get_cellhd(input_name, raster_mapset, &cellhd);

    G_debug(3, "number of rows %d", cellhd.rows);

    /* for each row */
    for (row = 0; row < nrows; row++) {

        G_percent(row, nrows, 2);

        /* allocating memory for the output variable */
        output[row] = (DCELL *)G_malloc (sizeof(DCELL)*ncols);

    /* read input map */
    Rast_get_row(infd, inrast, row, data_type);

    /* process the data */
    for (col = 0; col < ncols; col++) {
                    output[row][col] = ((DCELL *) inrast)[col];
              
    }
    }

    /* memory cleanup */
    G_free(inrast);

    /* closing raster maps */
    Rast_close(infd);

    return (output);

}

void save_raster_C_variable(CELL **input,  char *input_name) {

    int     row;
    int     nrows;
    int     outfd;
    struct  History history;	/* holds meta-data (title, comments,..) */

    G_message(_("Saving <%s> raster"), input_name);
	
    /* Allocate output buffer, use input map data_type */
    nrows = Rast_window_rows();

    /* controlling, if we can write the raster */
    outfd = Rast_open_new(input_name, CELL_TYPE);

    /* for each row */
    for (row = 0; row < nrows; row++) {

        G_percent(row, nrows, 2);

    /* write raster row to output raster map */
    Rast_put_row(outfd, input[row], CELL_TYPE);
    }

    /* closing raster maps */
    Rast_close(outfd);

    /* memory cleanup */
    //close_raster_C_variable(input);

    /* add command line incantation to history file */
    Rast_short_history(input_name, "raster", &history);
    Rast_command_history(&history);
    Rast_write_history(input_name, &history);

    return;

}

void save_raster_D_variable(DCELL **input,  char *input_name) {

    int     row;
    int     nrow;
    int     outfd1;
    struct  History history;	/* holds meta-data (title, comments,..) */
//    RASTER_MAP_TYPE data_type;	/* type of the map (CELL/DCELL/...) */

    G_message(_("Saving <%s> raster"), input_name);

    /* Allocate output buffer, use input map data_type */
    nrow = Rast_window_rows();

    /* controlling, if we can write the raster */
    outfd1 = Rast_open_new(input_name, DCELL_TYPE);

    /* for each row */
    for (row = 0; row < nrow; row++) {

        G_percent(row, nrow, 2);

    /* write raster row to output raster map */
		Rast_put_row(outfd1, input[row], DCELL_TYPE);
    }

    /* closing raster maps */
    Rast_close(outfd1);

    /* memory cleanup */
    //close_raster_D_variable(input);

    /* add command line incantation to history file */
    Rast_short_history(input_name, "raster", &history);
    Rast_command_history(&history);
    Rast_write_history(input_name, &history);

    //return;

}

void close_raster_C_variable(CELL **input) {

    int	row, nrows;

    nrows = Rast_window_rows();

    for(row = 0; row < nrows; row++){
        G_free(input[row]);
    }
    G_free(input);

}

void close_raster_D_variable(DCELL **input) {

    int	row, nrows;

    nrows = Rast_window_rows();

    for(row = 0; row < nrows; row++){
        G_free(input[row]);
    }
    G_free(input);

}
