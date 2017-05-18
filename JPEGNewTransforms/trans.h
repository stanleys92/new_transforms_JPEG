/*   jidctint.h
Added 2016-2017 by Stanislav Svoboda
  - added functions for allocation of whole image
Copyright (C) 2016, 2017 Stanislav Svoboda
For conditions of distribution and use, see the accompanying LICENSE.md file.
*/

#include <stdlib.h>
#include <stdio.h>

#ifndef IMAGE_ARRAY_H
#define IMAGE_ARRAY_H 1
static double **mallocMatrix(int w,int h){
  double **matrix=malloc(sizeof(double*)*h);
  for(int i=0;i<h;i++){
    matrix[i]=malloc(sizeof(double)*w);
  }
  return matrix;
}

static int isEven(int x){
  if((x%2)==0) return 1;
  else return 0;
}
static int isOdd(int x){
  if((x%2)==0) return 0;
  else return 1;
}

double **matImg;
double **matCHBImg;
double **matCHRImg;
int widthImg;
int heightImg;
int widthCHImg;
int heightCHImg;


int transV;

int transIV;



int load_matrix;
#endif
