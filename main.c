
/****************************************************************************
 *
 * MODULE:       r.recovery
 * 
 * AUTHOR(S):    Luiz Augusto Richit ------------ luizaugustorichit@gmail.com
 *               Roberto Valmir da Silva -------- roberto.silva@uffs.edu.br
 *               Tomas Carlotto ----------------- thomas.carl@hotmail.com
 *               José Mario vicensi Grzybowski -- jose.grzybowski@uffs.edu.br
 *               
 *               
 * PURPOSE:      Simulates forest recovey/growth based on vegetation indeces
 * 
 * 
 * COPYRIGHT:    (C) 2017 by the GRASS Development Team
 *
 *               This program is free software under the GNU General Public
 *   	    	 License (>=v2). Read the file COPYING that comes with GRASS
 *   	    	 for details.
 *
 ****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <grass/gis.h>
#include <grass/raster.h>
#include <grass/glocale.h>
#include <grass/gmath.h>

#include "Simulator.h"
#include "functions.h"
#include "Calib_BiCGStab.h"

typedef double ValueType;

//  Feature
typedef struct{
    
//-------------INPUTS FOR SIMULATIATION ('Simulation Inputs guisection')----------------
    
    struct Option *DensityMap;
    struct Option *SoilUseMap;
    struct Option *Growth_rate;
    struct Option *D_Forest;
    struct Option *nyears;
    struct Option *Carrying_Capacity;
    struct Option *Res;
    struct Option *ClassWater;
    struct Option *ClassElem;
    struct Option *OutputMap;
    struct Flag *Simulator_flag;
    
//------------INPUTS FOR CALIBRATION ('Parameters Calibration guisection')---------------

    struct Option *Dens_Initial_Map;
    struct Option *Dens_End_Map;
    struct Option *SoilUse_Calib;
    struct Option *Time_Growth;
    struct Option *Cal_Carrying_Capacity;
    struct Option *Res_calib;
    struct Option *ClassWater_Calib;
    struct Option *ClassElem_Calib;
    struct Option *DensOutput;
    struct Option *ErrorOutput;
    struct Flag *Calibrator_flag;

    struct Option *input, *outpute;	// options 

} 

paramType;
paramType param;		//Parameters 


void set_parameters(){
    
//---------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------
//----INPUTS OF SIMULATION MODULE -------------
//---------------------------------------------------------------------------------------------
    
    param.Simulator_flag = G_define_flag();
    param.Simulator_flag->key = 's';
    param.Simulator_flag->label = _("Perform ONLY the simulation of forest growth for the reported inputs.");
    param.Simulator_flag->description = _("By marking only this option, the user limits the application of the module to the simulation of forest growth of a system of interest, but the values of 'Growth Rate - r' and 'Forest Diffusion - D' must be known or previous obtained");
    param.Simulator_flag->guisection = _("Simulation inputs");
    
    param.DensityMap = G_define_standard_option(G_OPT_R_INPUT);
    param.DensityMap->key = "density";
    param.DensityMap->label = _("Name of input density raster map");
    param.DensityMap->description = _("The values of the density map must have vegetation density in the range [0,1] and raster map must be greater than 3x3 pixels.");
    param.DensityMap->required = NO;
    param.DensityMap->guisection = _("Simulation inputs");

    param.SoilUseMap = G_define_standard_option(G_OPT_R_INPUT);
    param.SoilUseMap->key = "soil_use";
    param.SoilUseMap->label = _("Name of input soil use raster map");
    param.SoilUseMap->description = _("Map of soil use  with integer values of classes classification");
    param.SoilUseMap->required = NO;
    param.SoilUseMap->guisection = _("Simulation inputs");

    param.Growth_rate = G_define_option();
    param.Growth_rate->key = "ru";
    param.Growth_rate->label = _("Intrinsic forest growth rate (0.0<ru<=1.0) [1/year]");
    param.Growth_rate->required = NO;
    param.Growth_rate->type = TYPE_DOUBLE;
    param.Growth_rate->guisection = _("Simulation inputs");
    
    param.D_Forest = G_define_option();
    param.D_Forest->key = "du";
    param.D_Forest->label = _("Forest diffusion or dispersal coefficient (du>0.0) [m^2/year]");
    param.D_Forest->type = TYPE_DOUBLE;
    param.D_Forest->guisection = _("Simulation inputs");

    param.Carrying_Capacity = G_define_option();
    param.Carrying_Capacity->key = "ku";
    param.Carrying_Capacity->label = _("Environmental carrying capacity (0.0<ku<=1.0)[-]");
    param.Carrying_Capacity->required = NO;
    param.Carrying_Capacity->type = TYPE_DOUBLE;
    param.Carrying_Capacity->answer = "1.0";
    param.Carrying_Capacity->description = _("Carrying Capacity is a dimensionless coefficient. For dens forests, the recommended value is 1.0.");
    param.Carrying_Capacity->guisection = _("Simulation inputs");

    param.nyears = G_define_option();
    param.nyears->key = "nyears";
    param.nyears->label = _("Simulation time of forest growth [years]");
    param.nyears->required = NO;
    param.nyears->answer = "50";
    param.nyears->type = TYPE_INTEGER;
    param.nyears->guisection = _("Simulation inputs");
    
    param.ClassElem= G_define_option();
    param.ClassElem->key = "elem_class";
    param.ClassElem->label = _("Soil use map classes that represent elements of the system");
    param.ClassElem->type = TYPE_INTEGER;
    param.ClassElem->required = NO;
    param.ClassElem->multiple = YES;
    param.ClassElem->guisection = _("Simulation inputs");

    param.ClassWater= G_define_option();
    param.ClassWater->key = "water_class";
    param.ClassWater->label = _("Soil use map classes that represent permeable elements of the system");
    param.ClassWater->description = _("The system elements corresponding of this class(es) are not elements of system but allows interactivity between separated elements (i.e. rivers, lakes, streets or roads, etc.)");
    param.ClassWater->type = TYPE_INTEGER;
    param.ClassWater->required = NO;
    param.ClassWater->multiple = YES;
    param.ClassWater->guisection = _("Simulation inputs");

    param.OutputMap = G_define_option();
    param.OutputMap->key = "output";
    param.OutputMap->label = _("Name of output density raster map");
    param.OutputMap->required = NO;
    param.OutputMap->description = _("Saves the final condition of forest evolution to the value of reported growth time (nyears). This raster density result can be used as a new initial condition to simulate the evolution for the same or any other increment of time.As output, annual growth maps are saved with the addition of an indicative suffix to the output name entered by the user (i.e. for the 'Park' name the outputs are 'Park_Forest_1_yrs_Growth, Park_Forest_2_yrs_Growth, ... Park_Forest_x_yrs_Growth (x = nyears informed)).For the density maps resulting from the simulation, the initial density values for the non-system elements will be copied.");
    param.OutputMap->guisection = _("Simulation inputs");
    
//---------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------
//------------INPUTS OF CALIBRATION MODULE-----------------------------------------------------
//---------------------------------------------------------------------------------------------
    
    param.Calibrator_flag = G_define_flag();
    param.Calibrator_flag->key = 'c';
    param.Calibrator_flag->label = _("Perform ONLY the calibration of forest growth.");
    param.Calibrator_flag->description = _("By marking only this option the module do not run any simulation growth, the user only obtains as output the calibrated parameters (ru, du) for the evaluated forest system.");
    param.Calibrator_flag->guisection= _("Parameters calibration");
    
    param.Dens_Initial_Map = G_define_standard_option(G_OPT_R_INPUT);
    param.Dens_Initial_Map->key = "initial_density_calib";
    param.Dens_Initial_Map->label = _("Name of input density raster map for initial condition");
    param.Dens_Initial_Map->description = _("In case of calibration this option is required and density map must be in the range [0,1]");
    param.Dens_Initial_Map->required = NO;
    param.Dens_Initial_Map->guisection= _("Parameters calibration");

    param.Dens_End_Map= G_define_standard_option(G_OPT_R_INPUT);
    param.Dens_End_Map->key = "end_density_calib";
    param.Dens_End_Map->label = _("Name of input density raster map for the end condition");
    param.Dens_End_Map->description = _("In case of calibration this option is required and density map must be in the range [0,1]");
    param.Dens_End_Map->required = NO;
    param.Dens_End_Map->guisection= _("Parameters calibration");
    
    param.SoilUse_Calib = G_define_standard_option(G_OPT_R_INPUT);
    param.SoilUse_Calib->key = "soil_use_calib";
    param.SoilUse_Calib->label = _("Name of input soil use raster map for calibration");
    param.SoilUse_Calib->description = _("In case of calibration this option is required.");
    param.SoilUse_Calib->required = NO;
    param.SoilUse_Calib->guisection= _("Parameters calibration");
    
    param.Time_Growth = G_define_option();
    param.Time_Growth->key = "nyears_calib";
    param.Time_Growth->label = _("Time interval between the two system conditions [years]");
    param.Time_Growth->answer = "15";
    param.Time_Growth->type = TYPE_INTEGER;
    param.Time_Growth->description = _(" - In case of calibration this option is required; - The time growth interval of few years are not recommended.");
    param.Time_Growth->guisection= _("Parameters calibration");
    
    param.Cal_Carrying_Capacity = G_define_option();
    param.Cal_Carrying_Capacity->key = "k_calib";
    param.Cal_Carrying_Capacity->label = _("Environmental carrying capacity (0.0<k_calib<=1.0)[-]");
    param.Cal_Carrying_Capacity->type = TYPE_DOUBLE;
    param.Cal_Carrying_Capacity->answer = "1.0";
    param.Cal_Carrying_Capacity->description = _("In case of calibration this option is required.");
    param.Cal_Carrying_Capacity->guisection= _("Parameters calibration");
    
    param.ClassElem_Calib= G_define_option();
    param.ClassElem_Calib->key = "elem_class_calib";
    param.ClassElem_Calib->label = _("Soil use map classes that represent elements of the system");
    param.ClassElem_Calib->type = TYPE_INTEGER;
    param.ClassElem_Calib->required = NO;
    param.ClassElem_Calib->multiple = YES;
    param.ClassElem_Calib->guisection = _("Parameters calibration");

    param.ClassWater_Calib= G_define_option();
    param.ClassWater_Calib->key = "water_class_calib";
    param.ClassWater_Calib->label = _("Soil use map classes that represent permeable elements of the system");
    param.ClassWater_Calib->description = _("The system elements corresponding of this class(es) are not elements of system but allows interactivity between separate elements (i.e. rivers, lakes, etc.)");
    param.ClassWater_Calib->type = TYPE_INTEGER;
    param.ClassWater_Calib->required = NO;
    param.ClassWater_Calib->multiple = YES;
    param.ClassWater_Calib->guisection = _("Parameters calibration");
    
    param.DensOutput= G_define_standard_option(G_OPT_R_OUTPUT);
    param.DensOutput->key = "dens_output_calib";
    param.DensOutput->label = _("Name of raster map simulated with parameters from calibration");
    param.DensOutput->description = _("Stores the simulation result for forest growth for growth time and others inputs reported.For the density map resulting from the calibration, the initial density values for the non-system elements will be copied.");
    param.DensOutput->required = NO;
    param.DensOutput->answer= "Density_cal";
    param.DensOutput->guisection = _("Parameters calibration");
    
    param.ErrorOutput= G_define_standard_option(G_OPT_R_OUTPUT);
    param.ErrorOutput->key = "error_output_calib";
    param.ErrorOutput->label = _("Name of raster map of calibration errors");
    param.ErrorOutput->description = _("Return the error between observed and predict end conditions ('O(i,j)-P(i,j)') for calibrated parameters. Negative Values means underestimation and positive means overestimation od density");
    param.ErrorOutput->required = NO;
    param.ErrorOutput->answer= "Error_cal";
    param.ErrorOutput->guisection = _("Parameters calibration");
    
}


int main(int argc, char *argv[])
{
    
    struct Cell_head window; 
    int nrows, ncols;
    int Temp_value;
   
    struct History history;	/* holds meta-data (title, comments,..) */

    

    struct GModule *module;	 // GRASS module for parsing arguments

    // Initialize GIS environment
    G_gisinit(argv[0]);		 // reads grass env, stores program name to G_program_name()

    // Header of window module
    module = G_define_module();
    G_add_keyword(_("raster"));
    G_add_keyword(_("forest"));
    G_add_keyword(_("growth"));
    G_add_keyword(_("DLG-model"));
    module->description = _("Simulates forest growth from the Diffusive-Logistic Growth model. It enables calibration of model parameters if they are not known.");
    
    // set parameters for user (user window)
    set_parameters();

    // options and flags parser 
    if (G_parser(argc, argv))
	exit(EXIT_FAILURE);

    
//---------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------
// -----MANAGEMENT OF SELECT FLAGS-------------------------------------------------------------
//---------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------

//---IF Calibration and Simulation flags is both select----
 
    if (param.Calibrator_flag->answer && param.Simulator_flag->answer){
        G_fatal_error(_("\nYou  have marked the option to  perform  ONLY the calibration\n" 
                        "and ONLY the simulation concomitantly. To run the  module you\n"
                        "must select an  option for each  execution.  If  you  want to\n"
                        "perform the calibration check the option and run the  module.\n"
                        "Keep the  calibrated  values  and then run the  simulation by\n"
                        "entering these calibrated parameters or inform the parameters\n"
                        "directly (if they are already known) and mark and execute ONLY\n"
                        "the simulation option."));
    }


//---------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------
// -----VERIFY THE INPUTS FOR CALIBRATION (MESSAGE ERRORS)-------------------------------------
//---------------------------------------------------------------------------------------------
 if (param.Calibrator_flag->answer){  
 
    double k_cal; // Carrying Capacity employed for Calibration
    int dt;       // Time interval between two conditions of density evolution
    
    
    if (!param.Dens_Initial_Map->answer){
        G_fatal_error(_("For the calibration option is required the initial density rater\n map input (<%s>) for interest calibration area."),param.Dens_Initial_Map->key);
    }
    if (!param.Dens_End_Map->answer){
        G_fatal_error(_("For the calibration option is required the end density rater map input\n (<%s>) for interest calibration area."),param.Dens_End_Map->key);
    }
    if (!param.SoilUse_Calib->answer){
        G_fatal_error(_("For the calibration option is required the soil use rater map input\n (<%s>) for interest calibration area."),param.SoilUse_Calib->key);
    }
    if(!param.Time_Growth->answer){
       G_fatal_error(_("You must enter the <%s> value."),param.Time_Growth->key);
    }
    //Load the values from paramete to variable (dt)
    sscanf(param.Time_Growth->answer, "%d", &dt);
    
    if (dt==15){
       G_warning(_("You are running the module using the default value for <%s>.\nThe interval between two conditions is the default value?"), param.Time_Growth->key);
    }
    
    if(!param.Cal_Carrying_Capacity->answer){
        G_fatal_error(_("You must enter the <%s> value."),param.Cal_Carrying_Capacity->key);
    }
    //Load the values from paramete to variable (k_cal)
    sscanf(param.Cal_Carrying_Capacity->answer, "%lf", &k_cal);
    
    if (k_cal==1.0){
        G_message(_("You are running the module using the default value for <%s>."), param.Cal_Carrying_Capacity->key);
    }
    if (k_cal>1.0){
        G_fatal_error(_("Attention: <%s> must be less than or equal to 1.0."),param.Cal_Carrying_Capacity->key);
    }
    if (!param.ClassWater_Calib->answer){
        //G_message(_("You did not inform any Element at <%s>, are you sure?"), param.ClassWater_Calib->key);
    }
    if (!param.ClassElem_Calib->answer){
        G_fatal_error(_("You need inform some integer value to indicate the \nelements of the system at option <%s>."),param.ClassElem_Calib->key);
    }
    if(param.Cal_Carrying_Capacity->answer && param.Time_Growth->answer){
         if (k_cal<=0.0){
             G_fatal_error(_("The <%s> input have negative or null value. Check the reported value. "), param.Cal_Carrying_Capacity->key);
         }
         if (dt<=0){
             G_fatal_error(_("The <%s> input have negative or null value. Check the reported value."), param.Time_Growth->key);
         }
    } 
    //if (!param.DensOutput->answer){
      //  G_fatal_error(_("For the Calibration option is required a name for density \nraster map output (<%s>) simulated for optimazed parameters"),param.DensOutput->key);
    //}
    //if (!param.ErrorOutput->answer){
     //   G_fatal_error(_("For the Calibration option is required a name for density \nraster map output (<%s>) simulated for optimazed parameters"),param.DensOutput->key);
    //}

 } // End of 'if(C==1)'

//---------------------------------------------------------------------------------------------
// -----VERIFY THE INPUTS FOR SIMULATION (MESSAGE ERRORS)--------------------------------------
//---------------------------------------------------------------------------------------------
 if (param.Simulator_flag->answer){
    
    int n_years;
    double ru;
    double Du;
    double ku;
    
    //Load the values entered into the variables:
    if (!param.DensityMap->answer){
        G_fatal_error(_("For the simulation option is required the density rater map input (<%s>)."),param.DensityMap->key);
    }
    if (!param.SoilUseMap->answer){
        G_fatal_error(_("For the simulation option is required the soil use rater map input (<%s>)."),param.SoilUseMap->key);
    }
    if (!param.Growth_rate->answer){
        G_fatal_error(_("The parameter <%s> is required for simulation and therefore needs to be specified to perform simulation of forest growth ('<%s> section ')."),param.Growth_rate->key,param.Calibrator_flag->guisection);
    }
    if (!param.D_Forest->answer){ 
        G_fatal_error(_("The parameter <%s> is required for simulation and therefore needs to be specified to perform simulation of forest growth ('<%s> section ')."), param.D_Forest->key,param.Calibrator_flag->guisection);
    }
    if (!param.Carrying_Capacity->answer){
        G_fatal_error(_("For the simulation option is required the carrying capacity value (<%s>)."),param.Carrying_Capacity->answer);
    }
        sscanf(param.Carrying_Capacity->answer, "%lf", &ku); // copy value from k_calib
    
    if (!param.nyears->answer){
        G_fatal_error(_("For the simulation is required the time simulation at (<%s>)."),param.nyears->key);
    }
        sscanf(param.nyears->answer, "%d", &n_years);    // copy value from nyears
    
    if (param.Growth_rate->answer){
        sscanf(param.Growth_rate->answer, "%lf", &ru);   // copy value from ru
        if (ru<=0.0){
            G_fatal_error(_("The <%s> input have negative or null value. Check the reported value."),param.Growth_rate->key);
        }
        if (ru>1.0){
            G_fatal_error(_("Attention: <%s> must be less than or equal to 1.0."),param.Growth_rate->key);
        }
    }
    if (param.D_Forest->answer){
        sscanf(param.D_Forest->answer, "%lf", &Du);      // copy value from Du
        if (Du<=0.0){
            G_fatal_error(_("The <%s> input have negative or null value. Check the reported value."), param.D_Forest->key);
        }
    }
    if (n_years==50){
        G_message(_("You are running the module using the default value for <%s>."), param.nyears->key);
    }
    if (ku==1.0){
        G_message(_("You are running the module using the default value for <%s>."), param.Carrying_Capacity->key);
    }
    if (!param.ClassWater->answer){
        G_message(_("You did not inform any Element at <%s>, are you sure?"), param.ClassWater->key);
    }
    if (!param.ClassElem->answer){
        G_fatal_error(_("You need inform some integer value to indicate \n the elements of the system at option <%s>."),param.ClassElem->key);
    }
    if(param.Carrying_Capacity->answer && param.nyears->answer){
         if (ku<=0.0){
             G_fatal_error(_("The <%s> input have negative or null value. Check the reported value."), param.Carrying_Capacity->key);
         }
         if (n_years<=0){
             G_fatal_error(_("The <%s> input have negative or null value. Check the reported value."), param.nyears->key);
         }
    } 
    if (ku>1.0){
        G_fatal_error(_("Attention: <%s> must be less than or equal to 1.0."),param.Carrying_Capacity->key);
    }
    if (!param.OutputMap->answer){
        G_fatal_error(_("For the simulation option is required a name\n for density rater map output (<%s>)."),param.OutputMap->key);
    }
 }// end of 'if (S==1)'


//---------------------------------------------------------------------------------------------
// GET MAP INFORMATIONS ( nº of rows and cols and map resolution)------------------------------
//---------------------------------------------------------------------------------------------

    nrows = Rast_window_rows();
    ncols = Rast_window_cols();
     
    
    G_get_window(&window);
    double yres = (double)window.ns_res;
    double xres = (double)window.ew_res;
    double res = (xres+yres)/2.0;
    
    if (xres!=yres){
        G_warning(_("\nNumerical  formalization was  designed for square  resolution    maps,\n"
                    "but  the  map  resolution  differs  between the x  and y   directions:\n"
                    "xresolution = %.3f and  yresolution=%.3f. The  module apply the mean \n"
                    "value as standard resolution, but  this can only be  tolerated in  the \n"
                    "case of near-value resolution in  both directions. If  this is not the \n"
                    "case, the simulation results becomes invalid."),xres,yres);  
    }
    
//---------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------
// RUN THE CALIBRATOR FUNCTION (calibrator.h) - return D and r coefficients
//---------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------
if (param.Calibrator_flag->answer){// run calibrator.h
    
    double k_cal;        // Carrying Capacity employed for Calibration;
    int dt;              // Time interval between two conditions of density evolution;
    int nClassElem_cal;  // number of elements of classes vector 'ClassElem';
    int nClassWater_cal; // number of elements of classes vector 'ClassWater';
    int i, j;            // counter of rows and columns;
    DCELL** U0;
    DCELL** Uf;
    CELL** Soil_Use;
    char* Error_name;   // 2D-array for error raster map;
    char* Density_name; // 2D-array for name density raster map;

    //Load the values entered into the variables
    sscanf(param.Time_Growth->answer, "%d", &dt);
    sscanf(param.Cal_Carrying_Capacity->answer, "%lf", &k_cal);
    
    
//--------------------Copy parameters with MULTIPLE values------------------------------- 
    
    for (i = 0; param.ClassElem_Calib->answers[i]; i++) {
    }
    nClassElem_cal=i;
    int *ClassElem_vec_cal = (int*)G_malloc(nClassElem_cal* sizeof(int));

    for (i=0; i<nClassElem_cal; i++){
        sscanf(param.ClassElem_Calib->answers[i], "%d", &Temp_value);
        ClassElem_vec_cal[i]=Temp_value;;
    }
    
    
    nClassWater_cal=0;
    i=0;
    
    if (param.ClassWater_Calib->answer){
        for (i = 0; param.ClassWater_Calib->answers[i]; i++) {
        }
    }
    nClassWater_cal=i;
    
    if (nClassWater_cal<=0){
        nClassWater_cal=1;
    }
    int *ClassWater_vec_cal = (int*)G_malloc(nClassWater_cal* sizeof(int));
    
    if (param.ClassWater_Calib->answer){
    for (i=0; i<nClassWater_cal; i++){
         sscanf(param.ClassWater_Calib->answers[i], "%d", &Temp_value);
         ClassWater_vec_cal[i]=Temp_value;
    }
    }
    if (!param.ClassWater_Calib->answer){
        ClassWater_vec_cal[0]=-1;
        nClassWater_cal=-1;
    }
    
    //----------------------------------------------------------------
    // Verify if some value are the same
    for (i=0; i<nClassElem_cal; i++){
        for(j=0; j<nClassWater_cal; j++){
            if(ClassElem_vec_cal[i]==ClassWater_vec_cal[j]){
                G_fatal_error(_("One or more values ​​reported in <%s> and <%s>\n"
                                "are the same (see <%s> section)."),param.ClassWater_Calib->key, param.ClassElem_Calib->key,param.Calibrator_flag->guisection);
            }
                
        }
    }
    
    // stores options and flags to variables

    char *name_U0 = param.Dens_Initial_Map->answer;        // name of initial density raster map (Calibration)
    U0 = open_raster_D_variable(name_U0);
    
    char *name_Uf = param.Dens_End_Map->answer;            // name of end density raster map (Calibration)
    Uf = open_raster_D_variable(name_Uf);
    
    char *name_soil_use = param.SoilUse_Calib->answer;     // name of soil use raster map (Calibration)
    Soil_Use = open_raster_C_variable(name_soil_use);
    
    DCELL** Error = (DCELL**) G_malloc( nrows * sizeof(DCELL*));
    DCELL** Density = (DCELL**) G_malloc( nrows * sizeof(DCELL*));
    
        for (i = 0; i < nrows; i++){
            
            Error[i] = (DCELL*) G_malloc( ncols * sizeof(DCELL));
            Density[i] = (DCELL*) G_malloc( ncols * sizeof(DCELL));
            
               for (j = 0; j < ncols ; j++){
                   Error[i][j] = U0[i][j];
                 Density[i][j] = U0[i][j];
               }
        }

    
    //run calibrator.h
    Calibrator (Soil_Use, U0, Uf, k_cal, dt, nrows, ncols, ClassWater_vec_cal, ClassElem_vec_cal, res, nClassElem_cal, nClassWater_cal, Error, Density);
    
    
    if(param.ErrorOutput->answer){
        Error_name = param.ErrorOutput->answer; 
        save_raster_D_variable(Error, Error_name);
    }
    if(param.DensOutput->answer){
        Density_name = param.DensOutput->answer;
        save_raster_D_variable(Density, Density_name);
    }

    
}

//---------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------
// RUN THE SIMULATOR FUNCTION (simulator.h) - employ D and r coefficients declared by user 
//---------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------
if (param.Simulator_flag->answer){ // run simulator.h
    
    
    int n_years, i, j;
    double ru, Du, ku;
    int nClassElem, nClassWater;
    DCELL** Dens_0;
    CELL** SoilUseRaster;
    char* name_output;
    
    //Load the values entered into the variables
    sscanf(param.nyears->answer, "%d", &n_years);
    sscanf(param.Carrying_Capacity->answer, "%lf", &ku);
    
    if(param.D_Forest->answer){
       sscanf(param.D_Forest->answer,"%lf", &Du);
    }
    if(param.Growth_rate->answer){
       sscanf(param.Growth_rate->answer,"%lf", &ru);
    }
    
    //--------------------Copy parameters with MULTIPLE values------------------------------- 
    for (i = 0; param.ClassElem->answers[i]; i++) {
        // count number of elements of ClassElem vec
    }
    nClassElem=i;
    
    
    int *ClassElem_vec = (int*) G_malloc(nClassElem * sizeof(int));

    for (i=0; i<nClassElem; i++){
        sscanf(param.ClassElem->answers[i], "%d", &Temp_value);
        ClassElem_vec[i]=Temp_value;;
    }
    
    
    nClassWater=0;
    i=0;
    
    if (param.ClassWater->answer){
        for (i = 0; param.ClassWater->answers[i]; i++) {
        }
    }
    nClassWater=i;
    
    if (nClassWater<=0){
        nClassWater=1;
    }

    int *ClassWater_vec = (int*) G_malloc(nClassWater * sizeof(int));
    
    if (param.ClassWater->answer){
        for (i=0; i<nClassWater; i++){
         sscanf(param.ClassWater->answers[i], "%d", &Temp_value);
         ClassWater_vec[i]=Temp_value;
        }
    }
    if (!param.ClassWater->answer){
        ClassWater_vec[0]=-1;
        nClassWater=-1;
    }

    //----------------------------------------------------------------
    // Verify if some value are the same
    for (i=0; i<nClassElem; i++){
        for(j=0; j<nClassWater; j++){
            if(ClassElem_vec[i]==ClassWater_vec[j]){
                G_fatal_error(_("One or more values ​​reported in <%s> and <%s>\n"
                                "are the same (see <%s> section)."),param.ClassWater->key,param.ClassElem->key,param.Simulator_flag->guisection);
            }
                
        }
    }

    // stores options and flags to variables

    char *name_U0_simul = param.DensityMap->answer;        // name of initial density raster map (Calibration)
    Dens_0 = open_raster_D_variable(name_U0_simul);
    
    char *name_soilusesimul = param.SoilUseMap->answer;     // name of soil use raster map (Calibration)
    SoilUseRaster = open_raster_C_variable(name_soilusesimul);
    
    name_output=param.OutputMap->answer;
    
    // run simulator.h.
    
    Simulator(n_years, Du, ru,  ku, res, nrows, ncols, SoilUseRaster, Dens_0, ClassWater_vec, ClassElem_vec, nClassElem, nClassWater, name_output);
    G_message(_("\n The Model Exit Success!\n"));
    
}

//---------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------


    exit(EXIT_SUCCESS);
}

