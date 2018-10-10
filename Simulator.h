#ifndef SIMULATOR_H_INCLUDED
#define SIMULATOR_H_INCLUDED

typedef double ValueType;

//--------------------------------------------------------------------------------------------------
//begin {Author:  Roberto Valmir da Silva - roberto.silva@uffs.edu.br}

void save_raster_DCELL(DCELL **input,  char *input_name) {

    int     row;
    int     nrows;
    int     outfd;
    struct  History history;	// holds meta-data (title, comments,..) 
//    RASTER_MAP_TYPE data_type;	// type of the map (CELL/DCELL/...) 

    G_message(_("Saving <%s> raster"), input_name);

    // Allocate output buffer, use input map data_type 
    nrows = Rast_window_rows();

    // controlling, if we can write the raster 
    outfd = Rast_open_new(input_name, DCELL_TYPE);

    // for each row 
    for (row = 0; row < nrows; row++) {

        G_percent(row, nrows, 2);

    // write raster row to output raster map 
		Rast_put_row(outfd, input[row], DCELL_TYPE);
    }

    // closing raster maps 
    Rast_close(outfd);

    // memory cleanup 
    //close_raster_D_variable(input);

    // add command line incantation to history file 
    Rast_short_history(input_name, "raster", &history);
    Rast_command_history(&history);
    Rast_write_history(input_name, &history);

    //return;

}
//end {Author:  Roberto Valmir da Silva - roberto.silva@uffs.edu.br}
//--------------------------------------------------------------------------------------------------


//---------------------------------------------------------------------------------------------------
//begin {Author:  Luiz Augusto Richit - luizaugustorichit@gmail.com}

// Reformulate Aux to count Variables 
void Reformulate_Aux(int ncols, int nrows, int** Aux){

    int i, j;
    int k=0;

     for (i=0 ; i<nrows ; i++){
        for (j=0 ; j<ncols ; j++){

            if (Aux[i][j]==1){
                k++;
                Aux[i][j]=k; // Storage the corresponding  element number of system;
            }
        }
     }
}

// Fill 'ElemDiag' and 'Str' vectors 
void FuncCounter_0(int* Counter, int* ElemDiag, int* Str, int iter){

    int i;
    int y=0;
    int u=0;

    for (i=0; i<5; i++){
      if(Counter[i]!=0 && i!=2){
        y++;
      }
    }
    y+=1; //For central point;

    Str[iter+1]=Str[iter]+y;

    if (iter>0){

       for (i=0; i<2; i++){
         if(Counter[i]!=0){
           u++;
         }
       }
       ElemDiag[iter]=Str[iter]+u;
    }

}

// Fill 'Value_aux' and 'Col' vectors 
void FuncCounter_1(int* Counter, int* Str, int* Value_aux, int* Col, int* TempCol, int iter){

int i;
    // CC for Distant Diagonal
    if (Counter[0]==1 && Counter[4]==0){
        Counter[0]=2;
    }
    if (Counter[0]==0 && Counter[4]==1){
        Counter[4]=2;
    }
    if (Counter[0]==0 && Counter[4]==0){
        Counter[2]+=2;
    }
    // CC for Adjacent Diagonal
    if (Counter[1]==1 && Counter[3]==0){
        Counter[1]=2;
    }
    if (Counter[1]==0 && Counter[3]==1){
        Counter[3]=2;
    }
    if (Counter[1]==0 && Counter[3]==0){
        Counter[2]+=2;
    }

  int g=0;

    for (i=0; i<5; i++){
        if(Counter[i]!=0 || i==2){
           g++;
           Value_aux[Str[iter]+g-1] = Counter[i];
                 Col[Str[iter]+g-1] = TempCol[i];
        }
    }
}

// Fill 'TempCol' and 'Counter' vectors 
void Create_CRS_Vecs_0(int** Aux, int nrows, int ncols, int* ElemDiag, int* Str){

    // Used to count the number of non-null elements of characterized Matrix  

int *Counter = (int*) G_malloc(5* sizeof(int));

int i, j, k, r1, a, g1;
int iter=-1;

     for (i=0 ; i<nrows ; i++){
        for (j=0 ; j<ncols ; j++){

        // clear values of Counter for begin next verification 
        for (k=0; k<5; k++){
          Counter[k]=0;
        }
        // Begin verification of adjacencies 
           if (Aux[i][j]==1){
                iter++;

if (i==0){

        if (j==0){

              if (Aux[i][j+1]==1){
                  Counter[3]=1;
              }
              if (Aux[i+1][j]==1){
                  Counter[4]=1;
              }

              if (Aux[i][j+1]==-1){
                    for (r1 = j+2; r1 < ncols; r1++){
                            if (Aux[i][r1]==1){
                                Counter[3]=1;
                                r1=ncols;
                            }
                    }
              }
              if (Aux[i+1][j]==-1){
                  for (r1=i+2; r1<nrows; r1++){
                      if (Aux[r1][j]==1){
                          Counter[4]=1;
                          r1=nrows;
                      }
                  }
              }

        }//----End of j==0 within if-Statement (i==0)-----
        if (j>0 && j<ncols-1){

              if (Aux[i][j+1]==1){
                 Counter[3]=1;
              }
              if (Aux[i][j-1]==1){
                 Counter[1]=1;
              }
              if (Aux[i+1][j]==1){
                 Counter[4]=1;
              }
              if (Aux[i][j-1]==-1 && j>1){
                      for (r1 = j-2; r1 >= 0; r1--){
                           if (Aux[i][r1]==1){
                                Counter[1]=1;
                                r1=-1;
                           }
                     }
              }
              if (Aux[i][j+1]==-1 && j<=ncols-3){// if (j+2<=ncols-1)
                   for (r1 = j+2; r1 <ncols; r1++){
                            if (Aux[i][r1]==1){
                                Counter[3]=1;
                                r1=ncols;
                            }
                    }
              }
              if (Aux[i+1][j]==-1 && i<=nrows-3){
                  for (r1=i+2; r1<nrows; r1++){
                        if (Aux[r1][j]>0){
                                Counter[4]=1;
                                r1=nrows;
                        }
                  }
              }

        }//----End of  0<j<ncols  if-Statement-----
        if (j==ncols-1){

            if (Aux[i][j-1]==1){
                  Counter[1]=1;
            }
            if (Aux[i+1][j]==1){
                  Counter[4]=1;
            }
            if (Aux[i][j-1]==-1){
                for (r1 = j-2; r1 >= 0; r1--){
                    if (j-2>=0){
                        if (Aux[i][r1]==1){
                            Counter[1]=1;
                            r1=-1;
                        }
                    }
                }
            }
            if (Aux[i+1][j]==-1){
                for (r1=i+2; r1<nrows; r1++){
                        if (Aux[r1][j]==1){
                            Counter[4]=1;
                            r1=nrows;
                        }
                }
            }
        }//----End of  j==(ncols-1)  if-Statement-----

}//----End of  i==0  if-Statement-----

if (i>0 && i<nrows-1){

        if (j==0){

                  if (Aux[i][j+1]==1){
                      Counter[3]=1;
                  }
                  if (Aux[i+1][j]==1){
                      Counter[4]=1;
                  }
                  if (Aux[i-1][j]==1){
                      Counter[0]=1;
                  }
                  if (Aux[i][j+1]==-1){
                        for (r1=j+2; r1<ncols; r1++){
                          if (Aux[i][r1]>0){
                             Counter[3]=1;
                             r1=ncols;
                          }
                        }
                  }
                  if (Aux[i+1][j]==-1){
                        for (r1=i+2; r1<nrows; r1++){
                          if (Aux[r1][j]==1){
                             Counter[4]=1;
                             r1=nrows;
                          }
                        }
                  }
                  if (Aux[i-1][j]==-1 && i>=2){
                        for (r1=i-2; r1>=0; r1--){
                          if (Aux[r1][j]==1){
                             Counter[0]=1;
                             r1=-1;
                          }
                        }
                  }

        }//----End of  j==0  if-Statement------
        if (j>0 && j<ncols-1){

                  if (Aux[i][j+1]==1){
                     Counter[3]=1;
                  }
                  if (Aux[i][j-1]==1){
                     Counter[1]=1;
                  }
                  if (Aux[i+1][j]==1){
                     Counter[4]=1;
                  }
                  if (Aux[i-1][j]==1){
                     Counter[0]=1;
                  }
                  if (Aux[i][j+1]==-1){
                        for (r1=j+2; r1<ncols; r1++){
                          if (Aux[i][r1]==1){
                             Counter[3]=1;
                             r1=ncols;
                          }
                        }

                  }
                  if (Aux[i][j-1]==-1){
                        for (r1=j-1; r1>=0; r1--){
                          if (Aux[i][r1]==1){
                             Counter[1]=1;
                             r1=-1;
                          }
                        }
                  }
                  if (Aux[i-1][j]==-1 && i>=2){
                        for (r1=i-2; r1>=0; r1--){
                          if (Aux[r1][j]==1){
                             Counter[0]=1;
                             r1=-1;
                          }
                        }
                  }
                  if (Aux[i+1][j]==-1){
                        for (r1=i+2; r1<nrows; r1++){
                            if (Aux[r1][j]==1){
                                Counter[4]=1;
                                r1=nrows;
                            }
                        }
                  }

       }//----End of  0<j<ncols  if-Statement-----
        if (j==ncols-1){

                    if (Aux[i][j-1]==1){
                          Counter[1]=1;
                    }
                    if (Aux[i+1][j]==1){
                          Counter[4]=1;
                    }
                    if (Aux[i-1][j]==1){
                          Counter[0]=1;
                    }
                    if (Aux[i][j-1]==-1){
                        for (r1=j-2; r1>=0; r1--){
                          if (Aux[i][r1]==1){
                             Counter[1]=1;
                             r1=-1;
                          }
                        }

                    }
                    if (Aux[i-1][j]==-1 && i>=2){
                        for (r1=i-2; r1>=0; r1--){
                          if (Aux[r1][j]>0){
                              Counter[0]=1;
                              r1=-1;
                          }
                        }
                    }
                    if (Aux[i+1][j]==-1){
                       if(i<=nrows-3){
                            for (r1=i+2; r1<nrows; r1++){
                                if (Aux[r1][j]==1){
                                    Counter[4]=1;
                                    r1=nrows;
                                }
                            }
                       }
                    }
        }//----End of  j=ncols-1  if-Statement-----

}//----End of  0<i<nrows  if-Statement-----

if (i==nrows-1){

        if (j==0){

                  if (Aux[i][j+1]==1){
                      Counter[3]=1;
                  }
                  if (Aux[i-1][j]==1){
                      Counter[0]=1;
                  }
                  if (Aux[i][j+1]==-1){
                        for (r1=j+2; r1<ncols; r1++){
                          if (Aux[i][r1]==1){
                                Counter[3]=1;
                                r1=ncols;
                          }
                        }
                  }
                  if (Aux[i-1][j]==-1 && i>1){
                        for (r1=i-2; r1>-1; r1--){
                          if (Aux[r1][j]>0){
                               Counter[0]=1;
                               r1=-1;
                          }
                        }
                  }

        } //----End of  j==0  if-Statement------
        if (j>0 && j<ncols-1){

                  if (Aux[i][j+1]==1){
                     Counter[3]=1;
                  }
                  if (Aux[i][j-1]==1){
                     Counter[1]=1;
                  }
                  if (Aux[i-1][j]==1){
                     Counter[0]=1;
                  }
                  if (Aux[i][j+1]==-1 && j<=ncols-3){
                      for (r1=j+2; r1<ncols; r1++){
                          if (Aux[i][r1]>0){
                              Counter[3]=1;
                              r1=ncols;
                          }
                      }
                  }
                  if (Aux[i][j-1]==-1 && j>1){
                       for (r1=j-2; r1>=0; r1--){
                          if (Aux[i][r1]>0){
                              Counter[1]=1;
                              r1=-1;
                          }
                      }
                  }
                  if (Aux[i-1][j]==-1 && i>1){
                        for (r1=i-2; r1>=0; r1--){
                          if (Aux[r1][j]>0){
                              Counter[0]=1;
                              r1=-1;
                          }
                        }
                  }

       }//----End of  0<j<ncols  if-Statement-----
        if (j==ncols-1){

                  if (Aux[i][j-1]==1){
                          Counter[1]=1;
                  }
                  if (Aux[i-1][j]==1){
                          Counter[0]=1;
                  }
                  if (Aux[i][j-1]==-1){
                        for (r1=j-2; r1>=0; r1--){
                          if (Aux[i][r1]>0){
                              Counter[1]=1;
                              r1=-1;
                          }
                        }
                  }
                  if (Aux[i-1][j]==-1 && i>1){
                        for (r1=i-2; r1>=0; r1--){
                          if (Aux[r1][j]>0){
                             Counter[0]=1;
                             r1=-1;
                          }
                        }
                  }

        }//----End of  j=ncols-1  if-Statement-----

}//----End of  i=nrows-1  if-Statement-----

FuncCounter_0(Counter,ElemDiag,Str,iter);

           }//----End of  (Aux[i][j]==1)  if-Statement-----
        }//----End of the for loop (columns)-----
      } //----End of the for loop (rows)----
}

// Fill ColTemp and Counter (temp) to Fill Value_ and 'Col'   
void Create_CRS_Vecs_1(int** Aux, int nrows, int ncols, int* Str, int* Value_aux, int* Col){

int *Counter = (int*) G_malloc(5* sizeof(int));
int *TempCol = (int*) G_malloc(5* sizeof(int));


int i, j, k, r1;
int iter=-1;

// Uses to complete 'Value' and 'Col' Vectors 

iter=-1;

    for (i=0 ; i<nrows ; i++){
        for (j=0 ; j<ncols ; j++){


        // Clear values of Counter for begin next verification
        for (k=0; k<5; k++){
          Counter[k] = 0;
          TempCol[k] = 0;
        }

        // Begin verification of adjacencies 

           if (Aux[i][j]>0){
                iter++;
                TempCol[2]=Aux[i][j]-1;

if (i==0){

        if (j==0){

              if (Aux[i][j+1]>0 ){
                  Counter[3]=1;
                  TempCol[3]=Aux[i][j+1]-1;
              }
              if (Aux[i+1][j]>0){
                  Counter[4]=1;
                  TempCol[4]=Aux[i+1][j]-1;
              }
              if (Aux[i][j+1]==-1){
                    for (r1 = j+2; r1 < ncols; r1++){
                            if (Aux[i][r1]>0){
                                Counter[3]=1;
                                TempCol[3]=Aux[i][r1]-1;
                                r1=ncols;
                            }
                    }
              }
              if (Aux[i+1][j]==-1){
                        for (r1=i+2; r1<nrows; r1++){      //r1=i+2:n see that if i==0 so 0+2=2
                          if (Aux[r1][j]>0){
                                Counter[4]=1;
                                TempCol[4]=Aux[r1][j]-1;
                                r1=nrows;
                          }
                        }
              }
        }//----End of j==0 within if-Statement (i==0)-----
        if (j>0 && j<ncols-1){

              if (Aux[i][j+1]>0){
                 Counter[3]=1;
                 TempCol[3]=Aux[i][j+1]-1;
              }
              if (Aux[i][j-1]>0){
                 Counter[1]=1;
                 TempCol[1]=Aux[i][j-1]-1;
              }
              if (Aux[i+1][j]>0){
                 Counter[4]=1;
                 TempCol[4]=Aux[i+1][j]-1;
              }

              if (Aux[i][j-1]==-1 && j>1){

                   for (r1 = j-2; r1 >= 0; r1--){
                            if (Aux[i][r1]>0){
                                Counter[1]=1;
                                TempCol[1]=Aux[i][r1]-1;
                                r1=-1;
                            }
                    }
              }
              if (Aux[i][j+1]==-1 && j<=ncols-3){

                   for (r1 = j+2; r1 <ncols; r1++){
                            if (Aux[i][r1]>0){
                                Counter[3]=1;
                                TempCol[3]=Aux[i][r1]-1;
                                r1=ncols;
                            }
                    }
              }
              if (Aux[i+1][j]==-1 && i<=nrows-3){
                        for (r1=i+2; r1<nrows; r1++){ // for r1=i+1:n   i=0 ... i+1=1
                          if (Aux[r1][j]>0){
                                Counter[4]=1;
                                TempCol[4]=Aux[r1][j]-1;
                          }
                        }

              }

        }//----End of  0<j<ncols  if-Statement-----
        if (j==ncols-1){

            if (Aux[i][j-1]>0){
                  Counter[1]=1;
                  TempCol[1]=Aux[i][j-1]-1;
            }
            if (Aux[i+1][j]>0){
                  Counter[4]=1;
                  TempCol[4]=Aux[i+1][j]-1;
            }

            if (Aux[i][j-1]==-1){
                for (r1 = j-2; r1 >= 0; r1--){
                        if (Aux[i][r1]>0){
                            Counter[1]=1;
                            TempCol[1]=Aux[i][r1]-1;
                            r1=-1;
                        }
                }
            }
            if (Aux[i+1][j]==-1){

                    for (r1=i+2; r1<nrows; r1++){
                          if (Aux[r1][j]>0){
                            Counter[4]=1;
                            TempCol[4]=Aux[r1][j]-1;
                            r1=nrows;
                          }
                        }

            }
        }//----End of  j==(ncols-1)  if-Statement-----

}//----End of  i==0  if-Statement-----

if (i>0 && i<nrows-1){

        if (j==0){

                  if (Aux[i][j+1]>0){
                      Counter[3]=1;
                      TempCol[3]=Aux[i][j+1]-1;
                  }
                  if (Aux[i+1][j]>0){
                      Counter[4]=1;
                      TempCol[4]=Aux[i+1][j]-1;
                  }
                  if (Aux[i-1][j]>0){
                      Counter[0]=1;
                      TempCol[0]=Aux[i-1][j]-1;
                  }
                  if (Aux[i][j+1]==-1){

                        for (r1=j+2; r1<ncols; r1++){

                          if (Aux[i][r1]>0){
                             Counter[3]=1;
                             TempCol[3]=Aux[i][r1]-1;
                             r1=ncols;
                          }
                        }
                  }
                  if (Aux[i+1][j]==-1){

                        for (r1=i+2; r1<nrows; r1++){

                          if (Aux[r1][j]>0){
                             Counter[4]=1;
                             TempCol[4]=Aux[r1][j]-1;
                             r1=nrows;
                          }
                        }
                  }
                  if (Aux[i-1][j]==-1 && i>=2){
                        for (r1=i-2; r1>=0; r1--){

                          if (Aux[r1][j]>0){
                             Counter[0]=1;
                             TempCol[0]=Aux[r1][j]-1;
                             r1=-1;
                          }
                        }
                  }

        }//----End of  j==0  if-Statement------
        if (j>0 && j<ncols-1){

                  if (Aux[i][j+1]>0){
                     Counter[3]=1;
                     TempCol[3]=Aux[i][j+1]-1;
                  }
                  if (Aux[i][j-1]>0){
                     Counter[1]=1;
                     TempCol[1]=Aux[i][j-1]-1;
                  }
                  if (Aux[i+1][j]>0){
                     Counter[4]=1;
                     TempCol[4]=Aux[i+1][j]-1;
                  }
                  if (Aux[i-1][j]>0){
                     Counter[0]=1;
                     TempCol[0]=Aux[i-1][j]-1;
                  }
                  if (Aux[i][j+1]==-1){

                        for (r1=j+2; r1<ncols; r1++){

                          if (Aux[i][r1]>0){
                             Counter[3]=1;
                             TempCol[3]=Aux[i][r1]-1;
                             r1=ncols;
                          }
                        }

                  }
                  if (Aux[i][j-1]==-1){

                        for (r1=j-2; r1>=0; r1--){

                          if (Aux[i][r1]>0){
                             Counter[1]=1;
                             TempCol[1]=Aux[i][r1]-1;
                             r1=-1;
                          }
                        }
                  }
                  if (Aux[i-1][j]==-1 && i>=2){

                        for (r1=i-2; r1>=0; r1--){

                          if (Aux[r1][j]>0){
                             Counter[0]=1;
                             TempCol[0]=Aux[r1][j]-1;
                             r1=-1;
                          }
                        }
                  }
                  if (Aux[i+1][j]==-1){

                        for (r1=i+2; r1<nrows; r1++){
                            if (Aux[r1][j]>0){
                                Counter[4]=1;
                                TempCol[4]=Aux[r1][j]-1;
                                r1=nrows;
                            }
                        }

                  }

       }//----End of  0<j<ncols  if-Statement-----
        if (j==ncols-1){

                    if (Aux[i][j-1]>0){
                          Counter[1]=1;
                          TempCol[1]=Aux[i][j-1]-1;
                    }
                    if (Aux[i+1][j]>0){
                          Counter[4]=1;
                          TempCol[4]=Aux[i+1][j]-1;
                    }
                    if (Aux[i-1][j]>0){
                          Counter[0]=1;
                          TempCol[0]=Aux[i-1][j]-1;
                    }
                    if (Aux[i][j-1]==-1){
                        for (r1=j-2; r1>=0; r1--){
                          if (Aux[i][r1]>0){
                             Counter[1]=1;
                             TempCol[1]=Aux[i][r1]-1;
                             r1=-1;
                          }
                        }

                    }
                    if (Aux[i-1][j]==-1 && i>=2){
                        for (r1=i-2; r1>=0; r1--){
                          if (Aux[r1][j]>0){
                              Counter[0]=1;
                              TempCol[0]=Aux[r1][j]-1;
                              r1=-1;
                          }
                        }
                    }
                    if (Aux[i+1][j]==-1){
                        for (r1=i+2; r1<nrows; r1++){
                           if (Aux[r1][j]>0){
                                      Counter[4]=1;
                                      TempCol[4]=Aux[r1][j]-1;
                           }
                        }

                    }
        }//----End of  j=ncols-1  if-Statement-----

}//----End of  0<i<nrows  if-Statement-----

if (i==nrows-1){

        if (j==0){

                  if (Aux[i][j+1]>0){
                      Counter[3]=1;
                      TempCol[3]=Aux[i][j+1]-1;
                  }
                  if (Aux[i-1][j]>0){
                      Counter[0]=1;
                      TempCol[0]=Aux[i-1][j]-1;
                  }
                  if (Aux[i][j+1]==-1){
                            for (r1=j+2; r1<ncols; r1++){
                                if (Aux[i][r1]>0){
                                    Counter[3]=1;
                                    TempCol[3]=Aux[i][r1]-1;
                                    r1=ncols;
                                }
                            }
                  }
                  if (Aux[i-1][j]==-1 && i>1){

                        for (r1=i-2; r1>=0; r1--){

                          if (Aux[r1][j]>0){
                               Counter[0]=1;
                               TempCol[0]=Aux[r1][j]-1;
                               r1=-1;
                          }
                        }
                  }

        } //----End of  j==0  if-Statement------
        if (j>0 && j<ncols-1){

                  if (Aux[i][j+1]>0){
                     Counter[3]=1;
                     TempCol[3]=Aux[i][j+1]-1;
                  }
                  if (Aux[i][j-1]>0){
                     Counter[1]=1;
                     TempCol[1]=Aux[i][j-1]-1;
                  }
                  if (Aux[i-1][j]>0){
                     Counter[0]=1;
                     TempCol[0]=Aux[i-1][j]-1;
                  }
                  if (Aux[i][j+1]==-1 && j<=ncols-3){

                      for (r1=j+2; r1<ncols; r1++){

                          if (Aux[i][r1]>0){
                              Counter[3]=1;
                              TempCol[3]=Aux[i][r1]-1;
                              r1=ncols;
                          }
                      }
                  }
                  if (Aux[i][j-1]==-1 && j>1){

                       for (r1=j-2; r1>=0; r1--){

                          if (Aux[i][r1]>0){
                              Counter[1]=1;
                              TempCol[1]=Aux[i][r1]-1;
                              r1=ncols;
                          }
                      }
                  }
                  if (Aux[i-1][j]==-1 && i>1){

                        for (r1=i-2; r1>=0; r1--){

                          if (Aux[r1][j]>0){
                              Counter[0]=1;
                              TempCol[0]=Aux[r1][j]-1;
                              r1=-1;
                          }
                        }
                  }



       }//----End of  0<j<ncols  if-Statement-----
        if (j==ncols-1){

                  if (Aux[i][j-1]>0){
                          Counter[1]=1;
                          TempCol[1]=Aux[i][j-1]-1;
                  }
                  if (Aux[i-1][j]>0){
                          Counter[0]=1;
                          TempCol[0]=Aux[i-1][j]-1;
                  }
                  if (Aux[i][j-1]==-1){

                        for (r1=j-2; r1>=0; r1--){

                          if (Aux[i][r1]>0){
                              Counter[1]=1;
                              TempCol[1]=Aux[i][r1]-1;
                              r1=-1;
                          }
                        }
                  }
                  if (Aux[i-1][j]==-1 && i>1){

                        for (r1=i-2; r1>=0; r1--){

                          if (Aux[r1][j]>0){
                             Counter[0]=1;
                             TempCol[0]=Aux[r1][j]-1;
                             r1=-1;
                          }
                        }
                  }

        }//----End of  j=ncols-1  if-Statement-----

}//----End of  i=nrows-1  if-Statement-----

FuncCounter_1(Counter,Str,Value_aux,Col,TempCol,iter);
//--------------------------------------------------------

           }//----End of  (Aux[i][j]>0)  if-Statement-----
        }//----End of the for loop (columns)-----
    }//----End of the for loop (rows)----


}

// Calculate the Number of Variables of System within a Rectangle
int Number_of_Variables(int rows_0,int rows_f, int col_0, int col_f,int **Aux){
    
      int row, col;
      int Sum=0;
      // Add the number of elements of input data that are elements of the System within a interval
      
      for (row=rows_0 ; row<rows_f ; row++){
        for (col=col_0 ; col<col_f ; col++){
           if (Aux[row][col]==1){
            Sum++;
           }
        }
      }
   return (Sum);
   
}// End Number_of_Variables() Function

// Concatenate a Matrix to 1D-Vector 
void Matrix2Vec(int var,int nrows, int ncols, ValueType **Dens, int **Aux, ValueType *U0){

    int i,j;   // i and j indexes
    int k=-1; //counter of variables
    ValueType u;

// Concatenates an array of scattered variables 
    for (i=0; i<nrows; i++){
        for (j=0; j<ncols; j++){
            if (Aux[i][j]==1){     // Condition to be a Variable of system
                k++;
                U0[k]= Dens[i][j]; // receive the Density value
            }
        }
    }

}//End of Matrix2Vec() Function 

// Turn 1D-Vector of original Matrix
void Vec2Matrix(int ncols, int nrows, int var, ValueType *U0, ValueType **Dens, int **Aux){

    // this function return a output matrix map with simulated recovery
    //We use the same matrix Dens to return the output matrix to economy memory
    //U0 is the vector of values of forest in end condition, because it receive the x-values at tf
    

    int i, j, q;

    q=0;
    for (i=0; i<nrows; i++){
        for (j=0; j<ncols; j++){
            if (Aux[i][j]>0){
                Dens[i][j]= U0[q];
                q++;
            }
        }
    }
}

// Make corrections of conflicting elements 
void Corrections(int ncols, int nrows, ValueType **Dens, int **Aux){
     int i, j;

     for ( i = 0 ; i < nrows; i++){
        for ( j = 0 ; j < ncols ; j++){
           if(Dens[i][j]<0.0){
                Aux[i][j]=0;
           }
           if(Dens[i][j]>1.0){
                Aux[i][j]=0;
           }
        }
     }

}

// Generate Aux Matrix 
void Generate_Aux_Matrix(int** Aux, int** SoilUse, int nrows, int ncols, int* ClassWater, int* ClassElem, int nClassElem, int nClassWater){

    int z, i, j;
if(nClassWater!=-1){
    for ( z = 0 ; z < nClassWater ; z++){
        for (i=0; i<nrows; i++){
            for (j=0; j<ncols; j++){
                if( SoilUse[i][j]==ClassWater[z] ){ // if the class if considered Water Surface
                    Aux[i][j]=-1;
                }
            }
        }
    }
}
    for ( z = 0 ; z < nClassElem ; z++){
        for (i=0; i<nrows; i++){
            for (j=0; j<ncols; j++){
                if(SoilUse[i][j]==ClassElem[z]){ // if the class if considered Element of System
                    Aux[i][j]=1;
                }
            }
        }
    }
}

// Complete definitively Value_B 
void CRS_Value_B(int var, ValueType coef, int sizeLines, int* Value_aux, ValueType *Value_B,int* ElemDiag){

    int i,j;
    ValueType m = 0.5*coef;
    ValueType n = 1.0-(2.0*coef);

  for (i = 0; i<sizeLines ; i++){
       Value_B[i]=(ValueType) Value_aux[i]*m;
  }

  for (i=0 ; i<var ; i++){
       j=ElemDiag[i];
       Value_B[j]+=n;
  }
}

// Complete definitively Value_A 
void CRS_Value_A(int var, ValueType coef, int sizeLines, int* Value_aux, ValueType *Value_A,int *ElemDiag){

    int i,j;
    ValueType m = -0.5*coef;
    ValueType n = 1.0+(2.0*coef);

  for (i = 0; i<sizeLines ; i++){
       Value_A[i]=(ValueType) Value_aux[i]*m;
  }

  for (i=0 ; i<var ; i++){
       j=ElemDiag[i];
       Value_A[j]+=n;
  }
}

// Formulate b auxiliary vector based in Model Equation of Phenomenon
void Fomulate_b (int var, ValueType r, ValueType ku, ValueType dt, ValueType *U0x, ValueType *b, ValueType *Value_B, int *Str, int *Col){

    int i,j,q;

     for ( i = 0; i < var ; i++){
            for (j=Str[i]; j<Str[i+1]; j++){
              q=Col[j];
              b[i]+=Value_B[j]*U0x[q];;
            }
        b[i]+=r*dt*U0x[i]*(1-(U0x[i]/ku));
     }
}
//end {Author:  Luiz Augusto Richit - luizaugustorichit@gmail.com}





//---------------------------------------------------------------------------------------------------
// begin {Author:  José Mario V. Grzybowski - jose.grzybowski@uffs.edu.br}

// Sub-routines to BICGStab
ValueType dot_p(ValueType *x, ValueType *y, int n){

	ValueType res = 0.0;
    int i;
    for (i = 0; i < n; i++){
        res += x[i] * y[i];
    }
    return res;
}

void add(ValueType *a, ValueType *b, int rows, ValueType *c){
    int i;
		for (i = 0; i<rows;i++){
			c[i] = a[i] + b[i];
	}
}

void subtr(ValueType *a, ValueType *b, int rows, ValueType *c){

    int i;
		for (i = 0; i<rows; i++){
			c[i] = a[i] - b[i];
	}
}

void atr(ValueType *a, int rows, ValueType *c){
    int i;
		for ( i = 0; i<rows; i++){
			c[i] = a[i];
	}
}

void sc_vec(ValueType a, ValueType *b, int rows, ValueType *c){
    int i;
		for (i = 0; i<rows; i++){
			c[i] = a*b[i];
	}
}

void mv_mult(ValueType *Value_A, ValueType *x, int var, int* Str, int* Col, ValueType *result){ // in matrix form: result = mat * vec;

    int i,j;
    int q;

    for ( i = 0; i < var ; i++){
        result[i]=0.0;
    }

    for ( i = 0; i < var ; i++){
            for (j=Str[i]; j<Str[i+1]; j++){
              q=Col[j];
              result[i]+=Value_A[j]*x[q];
            }
     }
}

// BiCGStab function
void BiCGStab(ValueType *Value_A, int* Str, int* Col, ValueType *b, ValueType *x0, ValueType tol, ValueType *x1, int n ){

int i;

	ValueType *r0 = (ValueType*) G_malloc(n* sizeof(ValueType));
	ValueType *r1 = (ValueType*) G_malloc(n* sizeof(ValueType));
	ValueType *p0 = (ValueType*) G_malloc(n* sizeof(ValueType));

	ValueType *p1 = (ValueType*) G_malloc(n* sizeof(ValueType));
	ValueType *rb = (ValueType*) G_malloc(n* sizeof(ValueType));
	ValueType *v0 = (ValueType*) G_malloc(n* sizeof(ValueType));
	ValueType *v1 = (ValueType*) G_malloc(n* sizeof(ValueType));

	 ValueType *s = (ValueType*) G_malloc(n* sizeof(ValueType));
	 ValueType *z = (ValueType*) G_malloc(n* sizeof(ValueType));
	 ValueType *t = (ValueType*) G_malloc(n* sizeof(ValueType));
	 ValueType *y = (ValueType*) G_malloc(n* sizeof(ValueType));

   ValueType alpha = 1.0;
   ValueType ro0 = 1.0;
   ValueType omega0 = 1.0;
   ValueType ro1;
   ValueType omega1;
   ValueType beta;
   ValueType norma;

	for (i = 0; i < n; i++){
	    p0[i] = 0;
		v0[i] = 0;
		x0[i] = 0;

	}

    ValueType * buff1 = (ValueType*) G_malloc(n* sizeof(ValueType));
    ValueType * buff2 = (ValueType*) G_malloc(n* sizeof(ValueType));

	int flag = 1;
	int j;
	i = 0;

	mv_mult(Value_A,x0,n,Str,Col,r0);

	subtr(b,r0,n,r0);

	atr(r0,n,rb); 


	while (flag==1){

		i++;
		ro1 = dot_p(rb,r0,n);
		beta = (ro1/ro0)*(alpha/omega0);

		sc_vec(omega0,v0,n,p1);

		subtr(p0,p1,n,p1);

		sc_vec(beta,p1,n,p1);
                
		add(r0,p1,n,p1);

		atr(p1,n,y);

		mv_mult(Value_A,y,n,Str,Col,v1);

		alpha = ro1/(dot_p(rb,v1,n));
                
		sc_vec(alpha,v1,n,s);
                
		subtr(r0,s,n,s);
                
		atr(s,n,z);

		mv_mult(Value_A,z,n,Str,Col,t);
                
		omega1 = (dot_p(t,s,n))/(dot_p(t,t,n));

		for ( j = 0; j < n; j++){
			buff1[j] = 0.0;
		}

		for (j = 0; j < n; j++){
			buff2[j] = 0.0;
		}

		sc_vec(omega1,z,n,buff1);
                
		sc_vec(alpha,y,n,buff2);

		add(buff1,buff2,n,x1);
                
		add(x0,x1,n,x1);
                
		mv_mult(Value_A,x1,n,Str,Col,buff1);

		subtr(b,buff1,n,buff1);
                
		norma = dot_p(buff1,buff1,n);
                
		norma = sqrt(norma);
                
		if (norma < tol || i>20 ){
			flag = 0;
		}

		else {
			sc_vec(omega1,t,n,r1);
			subtr(s,r1,n,r1);
			ro0 = ro1;
			atr(p1,n,p0);
			atr(v1,n,v0);
			omega0 = omega1;
			atr(x1,n,x0);
			atr(r1,n,r0);
			}
	}
	
}
// end {Author:  José Mario V. Grzybowski - jose.grzybowski@uffs.edu.br}
//--------------------------------------------------------------------------------------------------------------------------



//--------------------------------------------------------------------------------------------------------------------------
// begin {Author:  Luiz Augusto Richit - luizaugustorichit@gmail.com}

// Begin the Simulator() Function
void Simulator(ValueType nyears, ValueType Du, ValueType ru,ValueType ku, ValueType Res, int nrows, int ncols, int** SoilUse, ValueType** Dens, int* ClassWater, int* ClassElem, int nClassElem, int nClassWater, char* name_output){

    ValueType dx = Res;   // dx in meters
    ValueType dt = 1.0;   // dt in year
    ValueType TolBiCGStab = 0.000000001;

    ValueType coef;
    coef=((dt)*Du)/(dx*dx);
    
    int i,j;

    // Create Aux Matrix
    int   **Aux = (int**) G_malloc( nrows * sizeof(int*));
        for (i = 0; i < nrows; i++){
            Aux[i] = (int*) G_malloc( ncols * sizeof(int));
               for (j = 0; j<ncols ; j++){
                 Aux[i][j] = 0;
               }
        }

    Generate_Aux_Matrix(Aux,SoilUse,nrows,ncols,ClassWater,ClassElem,nClassElem,nClassWater);

    Corrections(ncols,nrows,Dens,Aux);

    // Number of Variables of the System 
     int var;
     var=Number_of_Variables(0,nrows,0,ncols,Aux);
     printf("Number of variables = %d",var);
     

    // Preallocating a Concatenated Vector with Initial Density Values - U0, b of Ax=b  and solution vec 'x_sol' 

     ValueType *x_sol = (ValueType *) G_malloc(var* sizeof(ValueType));
        ValueType *U0 = (ValueType *) G_malloc(var* sizeof(ValueType));
        ValueType *x0 = (ValueType *) G_malloc(var* sizeof(ValueType));
         ValueType *b = (ValueType *) G_malloc(var* sizeof(ValueType));
             for (i = 0; i < var; i++){
                  U0[i] = 0.00;
                   b[i] = 0.00;
               x_sol[i] = 0.00;
                  x0[i] = 0.00;
             }

   // Generate one 1D-array U0 to initial Density Values 
     Matrix2Vec(var,nrows,ncols,Dens,Aux,U0);

   // Pre-allocating ElenDiag and Str 
     int *ElemDiag = (int*) G_malloc( var* sizeof(int));
          int *Str = (int*) G_malloc ((var+1)* sizeof(int));
                  ElemDiag[0] = 0;
                       Str[0] = 0;

   // Create 'a fortiori' ElemDiag and Str 
    Create_CRS_Vecs_0(Aux, nrows, ncols, ElemDiag, Str);


    int sizeLines=Str[var]+1;
        int *Value_aux = (int*) G_malloc(sizeLines* sizeof(int));
            for (i = 0; i < sizeLines; i++){
                Value_aux[i]=0;

            }

    ValueType* Value_B = (ValueType *) G_malloc(sizeLines* sizeof(ValueType));
    ValueType* Value_A = (ValueType *) G_malloc(sizeLines* sizeof(ValueType));
              int* Col = (int*) G_malloc(sizeLines* sizeof(int));


   Reformulate_Aux(ncols,nrows,Aux);

   Create_CRS_Vecs_1(Aux, nrows, ncols, Str, Value_aux, Col);


      ValueType time = 0.0;    // initial
          int ntimes = 0;  // initialize and pre-allocating with zero

   // Count the time steps Number = ntimes 
      while (time < (ValueType) nyears){
        ntimes++;
        time+=dt;
       }
       
   // Solve the system recursive times (ntimes) 

   //  Fill CRS Value_A-Matrix
     CRS_Value_A(var, coef, sizeLines, Value_aux, Value_A, ElemDiag);

   //  Fill CRS Value_B-Matrix
     CRS_Value_B(var, coef, sizeLines, Value_aux, Value_B, ElemDiag);
     
     
     
     int nt; //Counter of time steps
     char* nameTemp;
     char* buf[1024];
     char* time_middleffix[4];
     //char* Rast_color[2014];
     
        for (nt=0; nt<ntimes; nt++){

                 // Formulate the b vector B*U0+r*dt*U0(1-U0/ku);
                   Fomulate_b(var,ru,ku,dt,U0,b,Value_B,Str,Col);

                 // Solve the Equation Ax=b for step 'nt'
                   BiCGStab(Value_A, Str, Col, b, x0, TolBiCGStab, x_sol, var);

                 // Loading the solution x to U0 to solve the next step
                    for (i = 0; i < var ; i++){
                         U0[i] = x_sol[i];
                         b[i] = 0.0; // clear the b vector
                    }
                  // name of outputs rasters
                  
                  sprintf(buf,"");
                  sprintf(buf,"%s",name_output);
                  strcat(buf,"_Forest_");
                  sprintf(time_middleffix, "%i", nt+1);
		  strcat(buf, time_middleffix); 
                  strcat(buf,"_yrs_Growth");
                  nameTemp=buf;
                  
                  // run manage colors
                  
                  //sprintf(Rast_color,"r.resample --overwrite input=%s output=%s", "raster_temp_tmax", "raster_temp_tmax2");
		//		if(run(Rast_color))
                   //                exit(1);
                  
                  // fill the 2d array with its density values  
                  Vec2Matrix(ncols,nrows,var,U0,Dens,Aux);

                  // save each raster
                  save_raster_DCELL(Dens, nameTemp);
                  
                  
       }
   // Use the 'Dens' Map to  receive density values. Elements that are no variables are keep with initial density values 
     
   //Vec2Matrix(ncols,nrows,var,U0,Dens,Aux);
   //save_raster_DCELL(Dens, nameTemp);

}
//end {Author:  Luiz Augusto Richit - luizaugustorichit@gmail.com}

#endif // SIMULATOR_H_INCLUDED
