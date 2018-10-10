#ifndef CALIB_BICGSTAB_H_INCLUDED
#define CALIB_BICGSTAB_H_INCLUDED

typedef double ValueType;

/* ***************************************************************************
	 *
	 * MODULE:       header file (.h) with calibration functions for r.richit module
         * 
         * 
         * AUTHOR(S):    Luiz Augusto Richit ------------ luizaugustorichit@gmail.com
         *               Roberto Valmir da Silva -------- roberto.silva@uffs.edu.br
         *               Tomas Carlotto ----------------- thomas.carl@hotmail.com
         *               José Mario vicensi Grzybowski -- jose.grzybowski@uffs.edu.br
         * 
	 *
	 * PURPOSE:      funtions to Calibrate the parameters r and D of Diffusive-Logistic Growth Model to simulate
	 *               forest recovery with r.richit module
	 *
	 * COPYRIGHT:    (C) Richit, L. A. et al., 2016.
	 *
	 *
 *************************************************************************** */

 /* Reformulate Aux to count Variables */
void Reformulate_Aux_cal(int ncols, int nrows, int** Aux){

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

/* Fill 'ElemDiag' and 'Str' vectors */
void FuncCounter_0_cal(int* Counter, int* ElemDiag, int* Str, int iter){

    int i;
    int y=0;
    int u=0;

    for (i=0; i<5; i++){
            if(Counter[i]!=0 && i!=2){
                y++;
            }
    }
    y+=1; // for central point

    Str[iter+1]=Str[iter]+y;

    if (iter>0){

       for (i=0; i<2; i++){
          if(Counter[i]!=0){
                u++;
          }
       }
       ElemDiag[iter]=Str[iter]+u; //
    }
}

/* Fill 'Value_aux' and 'Col' vectors */
void FuncCounter_1_cal(int* Counter, int* Str, int* Value_aux, int* Col, int* TempCol, int iter){

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

/* Fill 'TempCol' and 'Counter' vectors */
void Create_CRS_Vecs_0_cal(int** Aux, int nrows, int ncols, int* ElemDiag, int* Str){
/* Used to count the number of non-null elements of characterized Matr_calix  */

int *Counter = (int*) G_malloc(5* sizeof(int));

int i, j, k, r1;
int iter=-1;

     for (i=0 ; i<nrows ; i++){
        for (j=0 ; j<ncols ; j++){


        /* clear values of Counter for begin next verification */
        for (k=0; k<5; k++){
          Counter[k]=0;
        }

        /* Begin verification of adjacencies */
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

        }/*----End of j==0 within if-Statement (i==0)-----*/
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

        }/*----End of  0<j<ncols  if-Statement-----*/
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
        }/*----End of  j==(ncols-1)  if-Statement-----*/

}/*----End of  i==0  if-Statement-----*/

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

        }/*----End of  j==0  if-Statement------*/
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

       }/*----End of  0<j<ncols  if-Statement-----*/
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
        }/*----End of  j=ncols-1  if-Statement-----*/

}/*----End of  0<i<nrows  if-Statement-----*/

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

}/*----End of  i=nrows-1  if-Statement-----*/

FuncCounter_0_cal(Counter,ElemDiag,Str,iter);

           }/*----End of  (Aux[i][j]==1)  if-Statement-----*/
        }/*----End of the for loop (columns)-----*/
      } /*----End of the for loop (rows)----*/

//G_free(Counter);
//G_free(TempCol);

}

void Create_CRS_Vecs_1_cal(int** Aux, int nrows, int ncols, int* Str, int* Value_aux, int* Col){

int *Counter = (int*) G_malloc(5* sizeof(int));
int *TempCol = (int*) G_malloc(5* sizeof(int));

int i, j, k, r1;
int iter=-1;

/* Uses to complete 'Value' and 'Col' Vectors */
iter=-1;

    for (i=0 ; i<nrows ; i++){
        for (j=0 ; j<ncols ; j++){

        /* Clear values of Counter for begin next verification*/
        for (k=0; k<5; k++){
          Counter[k]=0;
          TempCol[k]=0;
        }
        /* Begin verification of adjacencies */

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
        }/*----End of j==0 within if-Statement (i==0)-----*/
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

        }/*----End of  0<j<ncols  if-Statement-----*/
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
        }/*----End of  j==(ncols-1)  if-Statement-----*/

}/*----End of  i==0  if-Statement-----*/

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

        }/*----End of  j==0  if-Statement------*/
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

       }/*----End of  0<j<ncols  if-Statement-----*/
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
        }/*----End of  j=ncols-1  if-Statement-----*/

}/*----End of  0<i<nrows  if-Statement-----*/

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

        } /*----End of  j==0  if-Statement------*/
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



       }/*----End of  0<j<ncols  if-Statement-----*/
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

        }/*----End of  j=ncols-1  if-Statement-----*/

}/*----End of  i=nrows-1  if-Statement-----*/

FuncCounter_1_cal(Counter,Str,Value_aux,Col,TempCol,iter);
/*--------------------------------------------------------*/

           }/*----End of  (Aux[i][j]>0)  if-Statement-----*/
        }/*----End of the for loop (columns)-----*/
    }/*----End of the for loop (rows)----*/

}

/* Equations for error */
void Hist_Distrib(int subdiv, int var, int* NumberHist, ValueType* U_vec){

    int i,j;
    ValueType q1,q2,w,L;
    L = (ValueType) subdiv;
    w=1.0/L;

    /*Clear the pre-values */
    for (i=0; i<subdiv; i++){
        NumberHist[i]=0;
    }

    /*Divide into classes */
    for (i=0; i<var; i++){
       for (j=0; j<subdiv; j++){
           q1=(ValueType)j*w;
           q2=(ValueType)(j+1)*w;

         if(U_vec[i]>=q1 && U_vec[i]<q2 ){
            NumberHist[j]+=1;
         }
       }
    }
}

ValueType HistNumberAnalysis(int subdiv,int var,int* NumberHistNat, int* NumberHistCalc){

   int res=0;
   ValueType o;
   o = (ValueType) subdiv;
   ValueType Area;

   int i;

    for (i=0; i<subdiv; i++){
        if (NumberHistNat[i]<NumberHistCalc[i]){
            res+=NumberHistCalc[i]-NumberHistNat[i];
        }
        if (NumberHistNat[i]>=NumberHistCalc[i]){
            res+=NumberHistNat[i]-NumberHistCalc[i];
        }
    }

    Area=(ValueType) res; //S=b*h,  h=res or in percentage h=res/var;  b=1/subdiv;
    Area/=(ValueType) var;
    Area/=o;

    return Area;
}

/* Calculate a Root Mean Square Error of real and simulated values*/
ValueType RMSE_Error(int var, ValueType  *Uf_vec, ValueType *U0x){

    int i;
    ValueType ErrorRMSE=0.0;
    ValueType a;

    for (i=0; i<var; i++){
            a=(Uf_vec[i]-U0x[i])*(Uf_vec[i]-U0x[i]);
            ErrorRMSE+=a;
    }
    ErrorRMSE/=(ValueType) var;
    ErrorRMSE=sqrt(ErrorRMSE);

  return ErrorRMSE;
}

ValueType MAPE_Error (int var, ValueType* Uf_vec, ValueType* U0x){

    int i;
    ValueType ErrorMAPE;
    ValueType Temp=0.0;
    ValueType x=0.0;
    ValueType y=0.0;

    for (i=0; i<var; i++){

         Temp=(Uf_vec[i]-U0x[i]);
         y+=Uf_vec[i];

         if(Temp<0.0){
            Temp*=-1.0;
         }
         x+=Temp;
    }
     ErrorMAPE=x/y;
  return ErrorMAPE;
}

ValueType Bias_Error (int var, ValueType* Uf_vec, ValueType* U0x){

    int i;
    ValueType ErrorBias;
    ValueType Temp=0.0;
    ValueType x=0.0;
    ValueType y=(ValueType) var;

        for (i=0; i<var; i++){

         Temp=(Uf_vec[i]-U0x[i]);

               if(Temp<0.0){
                  Temp*=-1.0;
               }
         x+=(Temp/Uf_vec[i]);
        }

    ErrorBias=x/y;

  return ErrorBias;
}

ValueType MSE_Error (int var, ValueType* Uf_vec, ValueType* U0x){

    int i;
    ValueType ErrorMSE=0.0;

    for (i=0; i<var; i++){
            ErrorMSE+=((Uf_vec[i]-U0x[i])*(Uf_vec[i]-U0x[i]));
    }
    ErrorMSE/=(ValueType) var;

  return ErrorMSE;
}

/* Calculate the Number of Variables of System within a Rectangle*/
int Number_of_Variables_cal(int rows_0,int rows_f, int col_0, int col_f,int **Aux){

      int row, col;
      int Sum=0;
      /* Add the number of elements of input data that are elements of the System within a interval*/
      for (row=rows_0 ; row<rows_f ; row++){
        for (col=col_0 ; col<col_f ; col++){
           if (Aux[row][col]==1){
            Sum++;
           }
        }
      }
   return (Sum);
}/* End Number_of_Variables_cal() Function*/

/* Make corrections of conflicting elements */
void Corrections_cal(int ncols, int nrows, ValueType **U0map, ValueType **Ufmap, int **Aux){

     int i, j;

     for ( i = 0 ; i < nrows; i++){
        for ( j = 0 ; j < ncols ; j++){
           if(U0map[i][j]<0.0 || Ufmap[i][j]<0.0){ // Values out off range (value<0)
                Aux[i][j]=0;
           }
           if(U0map[i][j]>1.0){ // Values out off range (value>1)
                Aux[i][j]=0;
           }
           if(Ufmap[i][j]>1.0){ // Values out off range (value>1)
                Aux[i][j]=0;
           }
        }
     }
}

/* Concatenate a Matr_calix to 1D-Vector */
void Matrix2Vec_cal(int var,int nrows, int ncols, ValueType **Dens, int **Aux, ValueType *U0){

    int i,j;   // i and j indexes
    ValueType u;
    int k=-1; //counter of variables

/* Concatenates an array of scattered variables */
    for (i=0; i<nrows; i++){
        for (j=0; j<ncols; j++){
                u=(ValueType)Dens[i][j];
            if (Aux[i][j]==1){     // Condition to be a Variable of system
                k++;
                U0[k]=u; // receive the Density value
            }
        }
    }

}/*End of Matrix2Vec_cal() Function */

/* Create the outputs maps */
void CreateOutputsMaps(int ncols, int nrows, ValueType **UfSimul,ValueType **ErrDif, int **Aux, ValueType* U0x_tot, ValueType* Uf_tot){

/* This function return a output matr_calix map with simulated recovery &
              distortion between  calculated and real density values */

    int i, j, q;

    q=0;
    for (i=0; i<nrows; i++){
        for (j=0; j<ncols; j++){
            if (Aux[i][j]>0){

               UfSimul[i][j] = U0x_tot[q];
                ErrDif[i][j] = U0x_tot[q] - Uf_tot[q];
                q++;
            }
        }
    }
}

/* Generate Aux Matr_calix */
void Generate_Aux_Matrix_cal(int** Aux, int** SoilUse, int nrows, int ncols, int* ClassWater, int* ClassElem, int nClassElem, int nClassWater){

    int z, i, j;
    
if (nClassWater!=-1){
    for ( z = 0 ; z < nClassWater ; z++){
        for (i=0; i<nrows; i++){
            for (j=0; j<ncols; j++){
                if( SoilUse[i][j]==ClassWater[z] ){ // if the class if considered water surface
                        Aux[i][j]=-1;
                }
            }
        }
    }
}
    for ( z = 0 ; z < nClassElem ; z++){

        for (i=0; i<nrows; i++){
            for (j=0; j<ncols; j++){
                if(SoilUse[i][j]==ClassElem[z]){ // if the class if considered element of System
                        Aux[i][j]=1;
                }
            }
        }
    }

}

/* Complete definitively Vec_B */
void CRS_Value_B_cal(int var, ValueType coef, int sizeLines, int* Value_aux, ValueType *Value_B,int* ElemDiag){

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

/* Complete definitively Vec_B */
void CRS_Value_A_cal(int var, ValueType coef, int sizeLines, int* Value_aux, ValueType *Value_A,int *ElemDiag){

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

/* Formulate b auxiliary vector based in Model Equation of Phenomenon*/
void Fomulate_b_cal (int var, ValueType r, ValueType ku, ValueType dt, ValueType *U0x, ValueType *b, ValueType *Value_B, int *Str, int *Col){

    int i,j,q;

     for ( i = 0; i < var ; i++){

            for (j=Str[i]; j<Str[i+1]; j++){
              q=Col[j];
              b[i]+=Value_B[j]*U0x[q];;
            }
        b[i]+=r*dt*U0x[i]*(1-(U0x[i]/ku));
     }
}

/* Redefine the value for test*/
void Parameters_range(int nDu, int nru, ValueType* ru_s, ValueType* Du_s, ValueType Du0, ValueType Duf, ValueType ru0, ValueType ruf){

    ValueType dr=1.0;
    int i;

/* create new rectangle to test values of Du and ru*/
dr = (Duf-Du0)/(nDu-1);

    for (i=0; i<nDu; i++){
        Du_s[i]=Du0+(dr*i);// similar to linspace(Du0,Duf,nDu);
    }

dr = (ruf-ru0)/(nru-1);

    for (i=0; i<nru; i++){
        ru_s[i]=ru0+(dr*i);// similar to linspace(ru0,ruf,nDu);
    }
}

/* subroutines to BICGStab*/
ValueType dot_p_cal(ValueType *x, ValueType *y, int n){

	ValueType res = 0.0;
    int i;
    for (i = 0; i < n; i++){
        res += x[i] * y[i];
    }
    return res;
}

void add_cal(ValueType *a, ValueType *b, int rows, ValueType *c){
    int i;
		for (i = 0; i<rows;i++){
			c[i] = a[i] + b[i];
	}
}

void subtr_cal(ValueType *a, ValueType *b, int rows, ValueType *c){

    int i;
		for (i = 0; i<rows; i++){
			c[i] = a[i] - b[i];
	}
}

void atr_cal(ValueType *a, int rows, ValueType *c){
    int i;
		for ( i = 0; i<rows; i++){
			c[i] = a[i];
	}
}

void sc_vec_cal(ValueType a, ValueType *b, int rows, ValueType *c){
    int i;
		for (i = 0; i<rows; i++){
			c[i] = a*b[i];
	}
}

void mv_mult_cal(ValueType *Value_A, ValueType *x, int var, int* Str, int* Col, ValueType *result){ // in matr_calix form: result = mat * vec;

    int i,j,q;

     for ( i = 0; i < var ; i++){
          result[i]=0.00;
     }
     for ( i = 0; i < var ; i++){
        for (j=Str[i]; j<Str[i+1]; j++){
            q=Col[j];
            result[i]+=Value_A[j]*x[q];
        }
     }
}

//void BiCGStab_cal
void BiCGStab_cal(ValueType *Value_A, int* Str, int* Col, ValueType *b, ValueType *x0, ValueType tol, ValueType *x1, int n ){

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
		r0[i] = 0;

	}

    ValueType * buff1 = (ValueType*) G_malloc(n* sizeof(ValueType));
    ValueType * buff2 = (ValueType*) G_malloc(n* sizeof(ValueType));

	mv_mult_cal(Value_A,x0,n,Str,Col,r0);

	subtr_cal(b,r0,n,r0);

	atr_cal(r0,n,rb); // rb receive r0

	int flag = 1;
	int j;
	i = 0;

	while (flag==1){

		i++;
		ro1 = dot_p_cal(rb,r0,n);
		beta = (ro1/ro0)*(alpha/omega0);

		sc_vec_cal(omega0,v0,n,p1);

		subtr_cal(p0,p1,n,p1);

		sc_vec_cal(beta,p1,n,p1);
		add_cal(r0,p1,n,p1);

		//p1 = r0 + beta*(p0 - omega0*v0);
		atr_cal(p1,n,y);

		mv_mult_cal(Value_A,y,n,Str,Col,v1);
		//mv_mult_cal(Value_A, col, str, y, sizeLines, v1);

		alpha = ro1/(dot_p_cal(rb,v1,n));

		//printf("%f\n",alpha);
		sc_vec_cal(alpha,v1,n,s);
		subtr_cal(r0,s,n,s);
		//s = r0 - alpha*v1;
		atr_cal(s,n,z);
		//z = s;

		mv_mult_cal(Value_A,z,n,Str,Col,t);
		//mv_mult_cal(Value_A, col, str, z, sizeLines, t);
		omega1 = (dot_p_cal(t,s,n))/(dot_p_cal(t,t,n));

		for ( j = 0; j < n; j++){
			buff1[j] = 0.0;
		}
		//ValueType buff1[n][1] = {0};

		for (j = 0; j < n; j++){
			buff2[j] = 0.0;
		}
		//float buff2[n][1] = {0};

		sc_vec_cal(omega1,z,n,buff1);
		sc_vec_cal(alpha,y,n,buff2);

		add_cal(buff1,buff2,n,x1);
		add_cal(x0,x1,n,x1);
		//x1 = x0 + alpha*y + omega1*z;

		mv_mult_cal(Value_A,x1,n,Str,Col,buff1);

		//mv_mult_cal(Value_A, col, str, x1, sizeLines, buff1);

		subtr_cal(b,buff1,n,buff1);

		//norma = b - mv_mult_cal(A,x1,n,n); // buffer

		norma = dot_p_cal(buff1,buff1,n);
		norma = sqrt(norma);

		//norma = (norma*norma)^(0.5);
		if (norma < tol || i>20 ){
			flag = 0;
		}
		else {
			sc_vec_cal(omega1,t,n,r1);
			subtr_cal(s,r1,n,r1);
			ro0 = ro1;
			atr_cal(p1,n,p0);
			atr_cal(v1,n,v0);
			omega0 = omega1;
			atr_cal(x1,n,x0);
			atr_cal(r1,n,r0);
			}
	}

  // G_free(p0); G_free(r0); G_free(x0); G_free(v0); G_free(p1);G_free(v1); G_free(t); G_free(r1); G_free(x1); G_free(s);
  // G_free(buff1); G_free(buff2); G_free(rb); G_free(z);  G_free(y);

}


/* Obtain the D and r parameters of System*/
void Calibrator (int **SoilUse, ValueType** U0map, ValueType** Ufmap, ValueType ku, ValueType nyears, int nrows, int ncols, int *ClassWater, int *ClassElem, ValueType Res,int nClassElem, int nClassWater,ValueType** ErrDif, ValueType** UfSimul){
     
     ValueType TolBiCGStab_cal=0.000000001;
     ValueType Tol=0.01;
     int i,j;
     int subdiv=25;
     

    ValueType Du0 = 0.0001;   // small value
    ValueType Duf = 10.0;    // big value
    ValueType ru0 = 0.001;  // small value
    ValueType ruf = 0.10;  // big value
     

    ValueType dx=Res;   // in meters (for UTM porojections)
    ValueType dt=1.0;   // don't change

    ValueType ErCal;
    ValueType ErVal;
    ValueType ru;
    ValueType D;

// Create a vector of Du values with a interval based in input range
int pDu1, pDun, nDu;

pDu1=(int)(log10(Du0));
pDun=(int)(log10(Duf));

nDu=pDun-pDu1;

nDu*=2;

if (nDu<5){
    nDu=5;
}

ValueType dr=1.0;
dr*=(pDun-pDu1);
dr/=(nDu-1);

// Use the initial and final Du-values to create a vector with subdivisions to test
ValueType *Du_s = (ValueType*)G_malloc( nDu*sizeof(ValueType));
        for (i = 0; i < nDu; i++){
            Du_s[i]=pow(10,pDu1+i*dr);
        }

// Force Du0 and Duf as extremes values (Due to rounding, discrepancies may occur)
 Du_s[0]=Du0;
 Du_s[nDu-1]=Duf;

// Create a vector of ru values with a interval based in input range
int pru1, prun, nru;

pru1=(int)(log10(ru0));
prun=(int)(log10(ruf));

nru=prun-pru1;

if (nru<5){
    nru=5; // não viciar dentro de um intervalo
}

dr=1.0;
dr*=(prun-pru1);
dr/=(nru-1);

// Use the initial and final ru-values to create a vector with subdivisions to test
ValueType *ru_s = (ValueType*)G_malloc( nru*sizeof(ValueType));//  divisions for each decimal order
        for (i = 0; i < (nru); i++){
            ru_s[i]=pow(10,pru1+i*dr);
        }
// Force ru0 and ruf as extremes values (Due to rounding, discrepancies may occur)
 ru_s[0]=ru0;
 ru_s[nru-1]=ruf;

/*Create Aux Matr_calix*/
int   **Aux = (int**) G_malloc( nrows * sizeof(int*));
    for (i = 0; i < nrows; i++){
        Aux[i] = (int*) G_malloc( ncols * sizeof(int));
           for (j=0; j < ncols ; j++){
                 Aux[i][j] = 0;
           }
    }

/* Create Auxiliar matrix 'Aux' to A and B CRS formulation  */
Generate_Aux_Matrix_cal(Aux,SoilUse,nrows,ncols,ClassWater,ClassElem,nClassElem,nClassWater);

/*Corrections_cal for conflicting elements of raster maps*/
Corrections_cal(ncols,nrows,U0map,Ufmap,Aux);

/*Break the density maps in two means to calibrate and validate*/
int col,row;

if(ncols>=nrows){

        if(ncols%2==0){
          col=(ncols/2);
          row=nrows;
        }
        else{
          col=(ncols-1)/2;
          row=nrows;
        }
}

if(ncols<nrows){

        if(nrows%2==0){
          row=(nrows/2);
          col=ncols;
        }
        else{
          row=(nrows-1)/2;
          col=ncols;
        }
}

/* preallocating maps*/

ValueType **U0mapcal = (ValueType **) G_malloc( row * sizeof(ValueType*));
ValueType **Ufmapcal = (ValueType **) G_malloc( row * sizeof(ValueType*));
      int   **Auxcal = (int**) G_malloc( row * sizeof(int*));
ValueType **U0mapval = (ValueType **) G_malloc( row * sizeof(ValueType*));
ValueType **Ufmapval = (ValueType **) G_malloc( row * sizeof(ValueType*));
      int   **Auxval = (int**) G_malloc( row * sizeof(int*));

    for (i = 0; i < row; i++){
      U0mapcal[i] = (ValueType *) G_malloc( col * sizeof(ValueType));
      Ufmapcal[i] = (ValueType *) G_malloc( col * sizeof(ValueType));
        Auxcal[i] = (int*) G_malloc( col * sizeof(int));
      U0mapval[i] = (ValueType *) G_malloc( col * sizeof(ValueType));
      Ufmapval[i] = (ValueType *) G_malloc( col * sizeof(ValueType));
        Auxval[i] = (int*) G_malloc( col * sizeof(int));

           for (j=0; j<col ; j++){
               U0mapcal[i][j]=0.0;
               Ufmapcal[i][j]=0.0;
                 Auxcal[i][j]=0;
               U0mapval[i][j]=0.0;
               Ufmapval[i][j]=0.0;
                 Auxval[i][j]=0;
           }
    }

int d,e;
d=ncols%2;
e=nrows%2;

if (row==nrows){

    for (i = 0; i < row; i++){
           for (j=0; j<col; j++){

               if (d==0){
               U0mapcal[i][j]=U0map[i][j];
               Ufmapcal[i][j]=Ufmap[i][j];
                 Auxcal[i][j]=Aux[i][j];
               U0mapval[i][j]=U0map[i][j+col];
               Ufmapval[i][j]=Ufmap[i][j+col];
                 Auxval[i][j]=Aux[i][j+col];
               }
               if (d!=0){
               U0mapcal[i][j]=U0map[i][j];
               Ufmapcal[i][j]=Ufmap[i][j];
                 Auxcal[i][j]=Aux[i][j];
               U0mapval[i][j]=U0map[i][j+col+1];
               Ufmapval[i][j]=Ufmap[i][j+col+1];
                 Auxval[i][j]=Aux[i][j+col+1];
               }
           }
        }
}

if (col==ncols){

         for (i = 0; i < row; i++){
           for (j=0; j<col; j++){

             if (e==0){
               U0mapcal[i][j]=U0map[i][j];
               Ufmapcal[i][j]=Ufmap[i][j];
                 Auxcal[i][j]=Aux[i][j];
               U0mapval[i][j]=U0map[i+row][j];
               Ufmapval[i][j]=Ufmap[i+row][j];
                 Auxval[i][j]=Aux[i+row][j];
             }
             if (e!=0){
               U0mapcal[i][j]=U0map[i][j];
               Ufmapcal[i][j]=Ufmap[i][j];
                 Auxcal[i][j]=Aux[i][j];
               U0mapval[i][j]=U0map[i+row+1][j];
               Ufmapval[i][j]=Ufmap[i+row+1][j];
                 Auxval[i][j]=Aux[i+row+1][j];
             }
         }
        }
}

/* Number of variables within calibration and validation maps */
int var_cal,var_val,var_tot;
var_cal = Number_of_Variables_cal(0,row,0,col,Auxcal);
var_val = Number_of_Variables_cal(0,row,0,col,Auxval);
var_tot = Number_of_Variables_cal(0,nrows,0,ncols,Aux);


/* -C-A-L-I-B-R-A-T-I-O-N-/*

/* Preallocating vector to concatenate maps*/
ValueType *U0_vec_cal = (ValueType*) G_malloc( var_cal* sizeof(ValueType));
   ValueType *U0x_cal = (ValueType*) G_malloc( var_cal* sizeof(ValueType));
ValueType *Uf_vec_cal = (ValueType*) G_malloc( var_cal* sizeof(ValueType));
     ValueType *b_cal = (ValueType*) G_malloc( var_cal* sizeof(ValueType));


        for (i = 0; i < var_cal; i++){
            U0_vec_cal[i]=0.0;
            Uf_vec_cal[i]=0.0;
               U0x_cal[i]=0.0;
                 b_cal[i]=0.0;
        }

/* Concatenate U0map and Ufmap to U0vec and Ufvec */
Matrix2Vec_cal(var_cal, row, col, U0mapcal, Auxcal, U0_vec_cal);
Matrix2Vec_cal(var_cal, row, col, Ufmapcal, Auxcal, Uf_vec_cal);

/* Apply input Raster maps to formulate the System Equations and Matr_calices */
int *ElemDiag_cal = (int*) G_malloc( var_cal* sizeof(int));
     int *Str_cal = (int*) G_malloc ((var_cal+1)* sizeof(int));

Str_cal[0]=0;
ElemDiag_cal[0]=0;


// Create 'a fortiori' ElemDiag_x and Str_x:
Create_CRS_Vecs_0_cal(Auxcal, row, col, ElemDiag_cal, Str_cal);


int sizeLines_cal=Str_cal[var_cal]+1;

    int *Value_aux_cal = (int*) G_malloc(sizeLines_cal* sizeof(int));
        for (i = 0; i < sizeLines_cal; i++){
                Value_aux_cal[i]=0;
        }
ValueType* Value_B_cal = (ValueType*) G_malloc(sizeLines_cal* sizeof(ValueType));
ValueType* Value_A_cal = (ValueType*) G_malloc(sizeLines_cal* sizeof(ValueType));
         int * Col_cal = (int*) G_malloc(sizeLines_cal* sizeof(int));


Reformulate_Aux_cal(col,row,Auxcal);

Create_CRS_Vecs_1_cal(Auxcal, row,  col,  Str_cal, Value_aux_cal, Col_cal);

 /* Count the time steps Number = ntimes */
 ValueType time=0.0;
 int ntimes=0;
      while (time<(ValueType)nyears){
        ntimes++;
        time+=dt;
       }

/* Initialize solution x-vector  */
    //CALIBRATION
       ValueType *x_cal = (ValueType*) G_malloc(var_cal* sizeof(ValueType));
      ValueType *x0_cal = (ValueType*) G_malloc(var_cal* sizeof(ValueType));
           for (i = 0; i < var_cal; i++){
                 x0_cal[i] = 0.0;
           }
    //VALIDATION
       ValueType *x_val = (ValueType*) G_malloc(var_val* sizeof(ValueType));
      ValueType *x0_val = (ValueType*) G_malloc(var_val* sizeof(ValueType));
           for (i = 0; i < var_val; i++){
                 x0_val[i] = 0.0;
           }
    // FULLY MAP
       ValueType *x_tot = (ValueType*) G_malloc(var_tot* sizeof(ValueType));
      ValueType *x0_tot = (ValueType*) G_malloc(var_tot* sizeof(ValueType));
           for (i = 0; i < var_tot; i++){
                 x0_tot[i] = 0.0;
           }

/* Test, successively the rectangle of Du and ru-values*/
int *NumberHistNat = (int*) G_malloc(subdiv*sizeof(int));
int *NumberHistCalc = (int*) G_malloc(subdiv*sizeof(int));

Hist_Distrib(subdiv,var_cal,NumberHistNat,Uf_vec_cal);


int line,column;

ValueType coef,r,dif_r,dif_D,Ecal,Error;
Error=1.0;


int Ermin=1;  //flag for control calibration steps
int iter=0;  //Counter of calibration steps
int nt;     //Counter of time integration steps
int u,k;   //Counter of evaluate the ru and Du within their vectors

while (Ermin==1){

    iter++;

/*---Test inferior limit of range test----*/
for (u=0;  u < nDu; u++){
  for (k=0; k < nru; k++){

            for (i=0; i<var_cal; i++){
                U0x_cal[i]=U0_vec_cal[i];
            }

               coef = (Du_s[u]*dt)/(dx*dx);
               r = ru_s[k];

        /* Fill CRS A-Matr_calix*/
        CRS_Value_A_cal(var_cal, coef, sizeLines_cal, Value_aux_cal, Value_A_cal, ElemDiag_cal);

        /* Complete Value_B with its values based in Madj_cal for B_cal-Matr_calix*/
        CRS_Value_B_cal(var_cal, coef, sizeLines_cal, Value_aux_cal, Value_B_cal, ElemDiag_cal);

        /* Solve the system recursive times (ntimes) */
               for ( nt=0; nt<ntimes; nt++){

                 // Formulate the b vector B*U0+r*dt*U0(1-U0/ku);
                   Fomulate_b_cal(var_cal,r,ku,dt,U0x_cal,b_cal,Value_B_cal,Str_cal,Col_cal);

                 // Solve the Equation Ax=b for step 'nt'
                   BiCGStab_cal(Value_A_cal, Str_cal, Col_cal, b_cal, x0_cal, TolBiCGStab_cal, x_cal, var_cal);

                 // Loading the solution x to U0 to solve the next step
                   for (i=0; i<var_cal ; i++){
                      U0x_cal[i]=x_cal[i];
                           b_cal[i]=0.0;
                           x_cal[i]=0.0;
                   }
               }

Hist_Distrib(subdiv, var_cal, NumberHistCalc, U0x_cal);
Ecal=HistNumberAnalysis(subdiv, var_cal, NumberHistNat, NumberHistCalc);
//Ecal=RMSE_Error(var_cal,Uf_vec_cal, U0x_cal);
//Ecal=MAPE_Error(var_cal, Uf_vec_cal, U0x_cal);
//Ecal=MSE_Error(var_cal, Uf_vec_cal, U0x_cal);
//Ecal=Bias_Error(var_cal, Uf_vec_cal, U0x_cal);

          if (Ecal<Error){//
               line=u; // D
               column=k; // r
               Error=Ecal; // Error[0]=Error[1]; Error[1]=Ecal; Error_step[1]=Error[1];
               ErCal=Ecal;
          }

     } // end of k-for statement for iter==1;
    } // end of u-for statement for iter==1;

/* Redefine ru and Du for calculate news vectors to ru_s and Du_s */
    if (column==0){
        ru0=ru_s[column];
        ruf=ru_s[column+1];
    }
    if (column==(nru-1)){
        ru0=ru_s[column-1];
        ruf=ru_s[column];
    }
    if (column>0 && column<(nru-1)){
        ru0=ru_s[column-1];
        ruf=ru_s[column+1];
    }
    if (line==0){
        Du0=Du_s[line];
        Duf=Du_s[line+1];
    }
    if (line==(nDu-1)){
        Du0=Du_s[line-1];
        Duf=Du_s[line];
    }
    if (line>0 && line<(nDu-1)){
        Du0=Du_s[line-1];
        Duf=Du_s[line+1];
    }

/* Calculate news vectors to ru_s and Du_s */
Parameters_range(nDu,nru,ru_s,Du_s,Du0,Duf,ru0,ruf);

Error=1.0;

/*  Evaluate the stop criteria*/
if (iter>1){

dif_r=(ruf-ru0)/ruf;
dif_D=(Duf-Du0)/Duf;

             if (dif_r<Tol && dif_D<Tol){ //Stores the values of r and D that meet the lowest tolerable error condition.
                Ermin=0;  // Ermin receive a new value to End the test within the while-statement
                D=Du_s[line];
                ru=ru_s[column];
             }
}

}// end of While-statement

coef=(D*dt)/(dx*dx);

/* -C-A-L-I-B-R-A-T-I-O-N-RMSE */

/* Create and Fill CRS A-Matr_calix*/
CRS_Value_A_cal(var_cal, coef, sizeLines_cal, Value_aux_cal, Value_A_cal, ElemDiag_cal);

/* Complete Value_B with its values based in Madj_cal for B_cal-Matr_calix, and too Col_cal and Lin_cal  */
CRS_Value_B_cal(var_cal, coef, sizeLines_cal, Value_aux_cal, Value_B_cal, ElemDiag_cal);


/* Solve the system recursive times (ntimes) */

               for (nt=0; nt<ntimes; nt++){

                  if (nt==0){
                    for (i=0; i<var_cal; i++){
                        U0x_cal[i]=U0_vec_cal[i];
                    }
                  }

                 // Formulate the b vector B*U0+r*dt*U0(1-U0/ku);
                   Fomulate_b_cal(var_cal,ru,ku,dt,U0x_cal,b_cal,Value_B_cal,Str_cal,Col_cal);


                 // Solve the Equation Ax=b for step 'nt'
                   BiCGStab_cal(Value_A_cal, Str_cal, Col_cal, b_cal, x0_cal, TolBiCGStab_cal, x_cal, var_cal);


                 // Loading the solution x to U0 to solve the next step
                   for (i=0; i<var_cal ; i++){
                      U0x_cal[i]=x_cal[i];
                           b_cal[i]=0.0;
                           x_cal[i]=0.0;
                   }
               }

ErCal=RMSE_Error(var_cal,Uf_vec_cal, U0x_cal);

printf("\n Calibration is Finished, wait for display calibrated parameters and errors!\n");


/* -V-A-L-I-D-A-T-I-O-N- */

ValueType *U0_vec_val = (ValueType *) G_malloc( var_val* sizeof(ValueType));
   ValueType *U0x_val = (ValueType *) G_malloc( var_val* sizeof(ValueType));
ValueType *Uf_vec_val = (ValueType *) G_malloc( var_val* sizeof(ValueType));
     ValueType *b_val = (ValueType *) G_malloc( var_val* sizeof(ValueType));

        for (i = 0; i < var_val; i++){
            U0_vec_val[i]=0.0;
            Uf_vec_val[i]=0.0;
               U0x_val[i]=0.0;
                 b_val[i]=0.0;
        }

/* Concatenate Initial and End Density maps Condition */
Matrix2Vec_cal(var_val, row, col, U0mapval, Auxval, U0_vec_val);
Matrix2Vec_cal(var_val, row, col, Ufmapval, Auxval, Uf_vec_val);

/* Create ElemDiag and Str (Compressed row storage vector Str)*/
int *ElemDiag_val = (int*) G_malloc( var_val* sizeof(int));
     int *Str_val = (int*) G_malloc ((var_val+1)* sizeof(int));

Str_val[0]=0;
ElemDiag_val[0]=0;

/* Complete ElemDiag and Str*/
Create_CRS_Vecs_0_cal(Auxval, row, col, ElemDiag_val, Str_val);

int sizeLines_val=Str_val[var_val]+1;

    int* Value_aux_val = (int*) G_malloc(sizeLines_val* sizeof(int));
            for (i = 0; i < sizeLines_val; i++){
                Value_aux_val[i]=0;
            }
ValueType* Value_B_val = (ValueType*) G_malloc(sizeLines_val* sizeof(ValueType));
ValueType* Value_A_val = (ValueType*) G_malloc(sizeLines_val* sizeof(ValueType));
         int * Col_val = (int*) G_malloc(sizeLines_val* sizeof(int));

/* Reformulate Aux Matr_calix*/
Reformulate_Aux_cal(col,row,Auxval);

/*Complete 'Value_aux_x' and 'Col_x' based on 'Str_x' */
Create_CRS_Vecs_1_cal(Auxval, row,  col,  Str_val, Value_aux_val, Col_val);

/* Create and Fill CRS A-Matr_calix*/
CRS_Value_A_cal(var_val, coef, sizeLines_val, Value_aux_val, Value_A_val, ElemDiag_val);

/* Complete Value_B with its values based in Madj_cal for B_cal-Matr_calix, and too Col_cal and Lin_cal  */
CRS_Value_B_cal(var_val, coef, sizeLines_val, Value_aux_val, Value_B_val, ElemDiag_val);

/* Solve the system recursive times (ntimes) */

               for (nt=0; nt<ntimes; nt++){

                  if (nt==0){
                    for (i=0; i<var_val; i++){
                        U0x_val[i]=U0_vec_val[i];
                    }
                  }
                 // Formulate the b vector B*U0+r*dt*U0(1-U0/ku);
                   Fomulate_b_cal(var_val,ru,ku,dt,U0x_val,b_val,Value_B_val,Str_val,Col_val);

                 // Solve the Equation Ax=b for step 'nt'
                   BiCGStab_cal(Value_A_val, Str_val, Col_val, b_val, x0_val, TolBiCGStab_cal, x_val, var_val);

                 // Loading the solution x to U0 to solve the next step
                   for (i=0; i<var_val ; i++){
                      U0x_val[i]=x_val[i];
                           b_val[i]=0.0;
                           x_val[i]=0.0;
                   }
               }

ErVal=RMSE_Error(var_val,Uf_vec_val, U0x_val);


/* Re-S-I-M-U-L-A-T-I-O-N- */

ValueType *U0_vec_tot = (ValueType*) G_malloc( var_tot* sizeof(ValueType));
ValueType *Uf_vec_tot = (ValueType*) G_malloc( var_tot* sizeof(ValueType));
   ValueType *U0x_tot = (ValueType*) G_malloc( var_tot* sizeof(ValueType));
     ValueType *b_tot = (ValueType*) G_malloc( var_tot* sizeof(ValueType));

        for (i = 0; i < var_tot; i++){
            U0_vec_tot[i]=0.0;
            Uf_vec_tot[i]=0.0;
               U0x_tot[i]=0.0;
                 b_tot[i]=0.0;
        }

Matrix2Vec_cal(var_tot, nrows, ncols, U0map, Aux, U0_vec_tot);
Matrix2Vec_cal(var_tot, nrows, ncols, Ufmap, Aux, Uf_vec_tot);

int *ElemDiag_tot = (int*) G_malloc( var_tot* sizeof(int));
     int *Str_tot = (int*) G_malloc ((var_tot+1)* sizeof(int));

Str_tot[0]=0;
ElemDiag_tot[0]=0;

Create_CRS_Vecs_0_cal(Aux, nrows,ncols, ElemDiag_tot, Str_tot);


int sizeLines_tot=Str_tot[var_tot]+1;

    int *Value_aux_tot = (int*) G_malloc(sizeLines_tot* sizeof(int));
        for (i = 0; i < sizeLines_tot; i++){
                Value_aux_tot[i]=0;
        }
ValueType* Value_B_tot = (ValueType*) G_malloc(sizeLines_tot* sizeof(ValueType));
ValueType* Value_A_tot = (ValueType*) G_malloc(sizeLines_tot* sizeof(ValueType));
         int * Col_tot = (int*) G_malloc(sizeLines_tot* sizeof(int));

Reformulate_Aux_cal(ncols,nrows,Aux);

Create_CRS_Vecs_1_cal(Aux,  nrows, ncols, Str_tot, Value_aux_tot, Col_tot);

/*  Fill CRS Value_A-Matr_calix*/
CRS_Value_A_cal(var_tot, coef, sizeLines_tot, Value_aux_tot, Value_A_tot, ElemDiag_tot);

/*  Fill CRS Value_B-Matr_calix*/
CRS_Value_B_cal(var_tot, coef, sizeLines_tot, Value_aux_tot, Value_B_tot, ElemDiag_tot);

/* Solve the system recursive times (ntimes) */

               for (nt=0; nt<ntimes; nt++){

                  if (nt==0){
                    for (i=0; i<var_tot; i++){
                        U0x_tot[i]=U0_vec_tot[i];
                    }
                  }
                  // Formulate the b vector B*U0+r*dt*U0(1-U0/ku);
                   Fomulate_b_cal(var_tot,ru,ku,dt,U0x_tot,b_tot,Value_B_tot,Str_tot,Col_tot);

                 // Solve the Equation Ax=b for step 'nt'
                   BiCGStab_cal(Value_A_tot, Str_tot, Col_tot, b_tot, x0_tot, TolBiCGStab_cal, x_tot, var_tot);


                 // Loading the solution x to U0 to solve the next step
                   for ( i=0; i<var_tot ; i++){
                      U0x_tot[i]=x_tot[i];
                           b_tot[i]=0.0;
                           x_tot[i]=0.0;
                   }
               }

/* Memory cleanup */

//G_free(U0x_tot); G_free(ElemDiag_tot); G_free(U0_vec_tot); G_free(Uf_vec_tot); G_free(x_tot); G_free(x0_tot); G_free(b_tot); G_free(Value_B_tot); G_free(Value_A_tot); G_free(Col_tot); G_free(Str_tot); G_free(Value_aux_tot);

//G_free(U0x_cal); G_free(ElemDiag_cal); G_free(U0_vec_cal); G_free(Uf_vec_cal); G_free(x_cal); G_free(x0_cal); G_free(b_cal); G_free(Value_B_cal); G_free(Value_A_cal); G_free(Col_cal); G_free(Str_cal); G_free(Value_aux_cal);

//G_free(U0x_val); G_free(ElemDiag_val); G_free(U0_vec_val); G_free(Uf_vec_val); G_free(x_val); G_free(x0_val); G_free(b_val); G_free(Value_B_val); G_free(Value_A_val); G_free(Col_val); G_free(Str_val); G_free(Value_aux_val);

//

printf("\n Creating  Outputs maps, Wait a few!\n");

CreateOutputsMaps(ncols,nrows,UfSimul,ErrDif,Aux,U0x_tot,Uf_vec_tot);


 printf("\n D =  %f\n",D);
 printf("\n r =  %f\n",ru);
 printf("\n Error (RMSE) Calibration =  %f\n",ErCal);
 printf("\n Error (RMSE) Validation =  %f\n",ErVal);

   
}

#endif // CALIB_BICGSTAB_H_INCLUDED
