/*   jfdctint.h
Added 2016-2017 by Stanislav Svoboda
  - added functions for new transforms (FORWARD)
Copyright (C) 2016, 2017 Stanislav Svoboda
For conditions of distribution and use, see the accompanying LICENSE.md file.
*/

#include <stdio.h>
#include <stdlib.h>
#include "trans.h"
#include <math.h>
#define PI 3.14159 

int nstate=0;
int wpos=0;
int hpos=0;
int wCHpos=0;
int hCHpos=0;
int enable=1;
int enableWAV=1;
double fastMatrix[8][8];
int shifting=128;

//beta constant
double betaVec[8]={0.00000012795,0.00074839045,0.03338096005,0.2718346173,0.7281653827,0.96661904,0.9992516096,0.9999998721};
double LDCT_scale=0.8;

/*separable CDF 5/3 and CDF 9/7 constants*/
double alfa=-0.5;
double beta=0.25;
const double alfa97=-1.586134342;
const double beta97=-0.05298011854;
const double gama97=0.8829110762;
const double delta97=0.4435068522;
const double scale97=1.149604398;
/*nonseparable CDF 5/3 and CDF 9/7 constants*/
double alfaRB=-0.25;
double betaRB=0.125;
const double alfa97RB=-1.586134342/2;
const double beta97RB=-0.05298011854/2;
const double gama97RB=0.8829110762/2;
const double delta97RB=0.44350685222/2;
const double scale97RB=1.149604398;

/*Scale matrices*/
double dctScale[8][8]={{1024.00,924.25,942.36,924.25,1020.00,924.25,942.35,924.24}, 
                       {924.25,837.49,853.89,837.49,924.25,837.49,853.89,837.48}, 
                       {942.36,853.89,870.62,853.90,942.36,853.89,870.62,853.89}, 
                       {924.25,837.49,853.90,837.49,924.25,837.49,853.89,837.48}, 
                       {1020.00,924.25,942.36,924.25,1020.00,924.25,942.35,924.24}, 
                       {924.25,837.49,853.89,837.49,924.25,837.49,853.89,837.48}, 
                       {942.35,853.89,870.62,853.89,942.35,853.89,870.62,853.89}, 
                       {924.24,837.48,853.89,837.48,924.24,837.48,853.89,837.48}};


double dstScale[8][8]={{870.03,872.34,872.82,873.04,873.17,873.26,873.33,873.39}, 
                       {872.34,873.09,873.25,873.32,873.36,873.39,873.42,873.43}, 
                       {872.82,873.25,873.34,873.38,873.40,873.42,873.43,873.44}, 
                       {873.04,873.32,873.38,873.41,873.42,873.43,873.44,873.45}, 
                       {873.17,873.36,873.40,873.42,873.43,873.44,873.44,873.45}, 
                       {873.26,873.39,873.42,873.43,873.44,873.45,873.45,873.45}, 
                       {873.33,873.42,873.43,873.44,873.44,873.45,873.45,873.45}, 
                       {873.39,873.43,873.44,873.45,873.45,873.45,873.45,873.45}};

double whtScale[8][8]={{1016.0,1020.0,1020.0,1020.0,1020.0,1020.0,1020.0,1020.0}, 
                       {1020.0,1020.0,1020.0,1020.0,1020.0,1020.0,1020.0,1020.0},
                       {1020.0,1020.0,1020.0,1020.0,1020.0,1020.0,1020.0,1020.0},
                       {1020.0,1020.0,1020.0,1020.0,1020.0,1020.0,1020.0,1020.0},
                       {1020.0,1020.0,1020.0,1020.0,1020.0,1020.0,1020.0,1020.0},
                       {1020.0,1020.0,1020.0,1020.0,1020.0,1020.0,1020.0,1020.0},
                       {1020.0,1020.0,1020.0,1020.0,1020.0,1020.0,1020.0,1020.0},
                       {1020.0,1020.0,1020.0,1020.0,1020.0,1020.0,1020.0,1020.0}};


double DCHTMatrix[8][8]={{0.35355,0.35355,0.35355,0.35355,0.35355,0.35355,0.35355,0.35355},
                         {-0.54006,-0.38576,-0.23146,-0.07715,0.07715,0.23146,0.38576,0.54006},
                         {0.54006,0.07715,-0.23146,-0.38576,-0.38576,-0.23146,0.07715,0.54006},
                         {-0.43082,0.30773,0.43082,0.18464,-0.18464,-0.43082,-0.30773,0.43082},
                         {0.28204,-0.52378,-0.12087,0.36262,0.36262,-0.12087,-0.52378,0.28204},
                         {-0.14979,0.49215,-0.36377,-0.32097,0.32097,0.36377,-0.49215,0.14979},
                         {0.06155,-0.30773,0.55391,-0.30773,-0.30773,0.55391,-0.30773,0.06155},
                         {-0.01707,0.11949,-0.35846,0.59744,-0.59744,0.35846,-0.11949,0.01707}};


//modified CDF 5/3 and 9/7
static int wavIndex[21][3]={{0,0,1},
                            {0,0,2},
                            {1,1,2},
                            {0,1,2},
                            {1,0,2},
                            {0,0,4},
                            {2,2,4},
                            {0,2,4},
                            {2,0,4},
                            {1,1,4},
                            {3,3,4},
                            {1,3,4},
                            {3,1,4},
                            {0,1,4},
                            {2,3,4},
                            {0,3,4},
                            {2,1,4},
                            {1,0,4},
                            {3,2,4},
                            {1,2,4},
                            {3,0,4}};

//Walsh matrix
static int walsh[8][8] = {
        {1,  1,  1,  1,  1,  1,  1,  1},
        {1,  1,  1,  1, -1, -1, -1, -1},
        {1,  1, -1, -1, -1, -1,  1,  1},
        {1,  1, -1, -1,  1,  1, -1, -1},
        {1, -1, -1,  1,  1, -1, -1,  1},
        {1, -1, -1,  1, -1,  1,  1, -1},
        {1, -1,  1, -1, -1,  1, -1,  1},
        {1, -1,  1, -1,  1, -1,  1, -1}
    };

//approximated DCT
static double APPROX[8][8]={{0.35355339,0.35355339,0.35355339,0.35355339,0.35355339,0.35355339,0.35355339,0.35355339},
                       {1.41421356,0,0,0,0,0,0,-1.41421356},
                       {0.5,0,0,-0.5,-0.5,0,0,0.5},
                       {0,0,-0.70710678,0,0,0.70710678,0,0},
                       {0.35355339,-0.35355339,-0.35355339,0.35355339,0.35355339,-0.35355339,-0.35355339,0.35355339},
                       {0,-1.41421356,0,0,0,0,1.41421356,0},
                       {0,-0.5,0.5,0,0,0.5,-0.5,0},
                       {0,0,0,-0.70710678,0.70710678,0,0,0}};


//modified SMRT X index
int SMRT_X[8][8]={{0, 0, 0, 0, 0, 0, 0, 0},
{1, 1, 3, 5, 7, 1, 3, 1},
{1, 1, 3, 5, 7, 1, 3, 1},
{1, 1, 3, 5, 7, 1, 3, 1},
{1, 1, 3, 5, 7, 1, 3, 1},
{2, 2, 2, 2, 2, 2, 6, 2},
{2, 6, 6, 6, 6, 2, 6, 2},
{4, 4, 4, 4, 4, 4, 4, 4}};

//modified SMRT Y index
int SMRT_Y[8][8]={{0, 1, 1, 1, 1, 2, 2, 4},
{0, 1, 1, 1, 1, 2, 2, 4},
{0, 1, 1, 1, 1, 2, 2, 4},
{0, 1, 1, 1, 1, 2, 2, 4},
{0, 1, 1, 1, 1, 2, 2, 4},
{0, 1, 1, 1, 1, 2, 2, 4},
{0, 1, 1, 1, 1, 2, 2, 4},
{0, 1, 1, 1, 1, 2, 2, 4}};

//modified SMRT matrix
int SMRT_p[8][8]={{0, 0, 1, 2, 3, 0, 2, 0},
{0, 0, 0, 0, 0, 0, 0, 0},
{1, 1, 1, 1, 1, 1, 1, 1},
{2, 2, 2, 2, 2, 2, 2, 2},
{3, 3, 3, 3, 3, 3, 3, 3},
{0, 0, 1, 2, 3, 0, 0, 0},
{2, 0, 1, 2, 3, 2, 2, 2},
{0, 0, 1, 2, 3, 0, 2, 0}};


/*whole image structures*/
typedef struct listMem{
  int type;
  int wpos;
  int hpos;
  struct listMem *next;
}*lMem;


typedef struct headList{
  struct listMem *first;
  struct listMem *act;
}*hLst;


struct headList *initHead(){
  struct headList *newHead;
  newHead=malloc(sizeof(struct headList));
  newHead->first=NULL;
  newHead->act=NULL;
  return newHead;
}

struct headList *HList=NULL;

void insertIntoList(int type,int w,int h){
  struct listMem *newMem;
  newMem=malloc(sizeof(struct listMem));
  newMem->type=type;
  newMem->wpos=w;
  newMem->hpos=h;
  newMem->next=NULL;
  if(HList->first==NULL){
    HList->first=newMem;
    HList->act=newMem;
  }
  else{
    HList->act->next=newMem;
    HList->act=newMem;
  }
}











/* computation of nonseparable CDF 5/3
 * @param1: whole image
 * @param2: image width
 * @param3: image height
 * @param4: step size
 * return
 */
void wavPassRB(double **matrix,int w,int h,int step){
  double scaleRB=sqrt(2);
//horizontal/vertical lifting
  int predX;
  int succX;
  int predY;
  int succY;
  
  for(int i=0,indH=0;i<h;i=i+step,indH++){
    for(int j=0,indW=0;j<w;j=j+step,indW++){
      if((indH % 2)!=(indW % 2)){
        if((i-step)<0) predY=i+step;
        else predY=i-step;
        if((i+step)>=h) succY=i-step;
        else succY=i+step;

        if((j-step)<0) predX=j+step;
        else predX=j-step;
        if((j+step)>=w) succX=j-step;
        else succX=j+step;
        matrix[i][j]=matrix[i][j]+alfaRB*(matrix[predY][j]+matrix[succY][j]+matrix[i][predX]+matrix[i][succX]);
      }
    }
  }
  for(int i=0,indH=0;i<h;i=i+step,indH++){
    for(int j=0,indW=0;j<w;j=j+step,indW++){
      if((indH % 2)==(indW % 2)){
        if((i-step)<0) predY=i+step;
        else predY=i-step;
        if((i+step)>=h) succY=i-step;
        else succY=i+step;

        if((j-step)<0) predX=j+step;
        else predX=j-step;
        if((j+step)>=w) succX=j-step;
        else succX=j+step;
        matrix[i][j]=matrix[i][j]+betaRB*(matrix[predY][j]+matrix[succY][j]+matrix[i][predX]+matrix[i][succX]);
      }
    }
  }
  //scaling
  for(int i=0,indH=0;i<h;i=i+step,indH++){
    for(int j=0,indW=0;j<w;j=j+step,indW++){
      if((indH % 2)==(indW % 2)){
        matrix[i][j]=matrix[i][j]*scaleRB;
      }
      else if((indH % 2)!=(indW % 2)){
        matrix[i][j]=matrix[i][j]/(-scaleRB);
      } 
    }
  }

  //diagonal lifting


  for(int i=0,indH=0;i<h;i=i+step,indH++){
    for(int j=0,indW=0;j<w;j=j+step,indW++){
      if(((indH % 2)==1)&&((indW % 2)==1)){
        if((i-step)<0) predY=i+step;
        else predY=i-step;
        if((i+step)>=h) succY=i-step;
        else succY=i+step;

        if((j-step)<0) predX=j+step;
        else predX=j-step;
        if((j+step)>=w) succX=j-step;
        else succX=j+step;
        matrix[i][j]=matrix[i][j]+alfaRB*(matrix[predY][predX]+matrix[predY][succX]+matrix[succY][predX]+matrix[succY][succX]);
      }
    }
  }
  for(int i=0,indH=0;i<h;i=i+step,indH++){
    for(int j=0,indW=0;j<w;j=j+step,indW++){
      if(((indH % 2)==0)&&((indW % 2)==0)){
        if((i-step)<0) predY=i+step;
        else predY=i-step;
        if((i+step)>=h) succY=i-step;
        else succY=i+step;

        if((j-step)<0) predX=j+step;
        else predX=j-step;
        if((j+step)>=w) succX=j-step;
        else succX=j+step;
        matrix[i][j]=matrix[i][j]+betaRB*(matrix[predY][predX]+matrix[predY][succX]+matrix[succY][predX]+matrix[succY][succX]);
      }
    }
  }

  //scaling
  for(int i=0,indH=0;i<h;i=i+step,indH++){
    for(int j=0,indW=0;j<w;j=j+step,indW++){
      if(((indH % 2)==0)&&((indW % 2)==0)){
        matrix[i][j]=matrix[i][j]*scaleRB;
      }
      else if(((indH % 2)==1)&&((indW % 2)==1)){
        matrix[i][j]=matrix[i][j]/(-scaleRB);
      }
    }
  }
  return;
}



/* computation of nonseparable CDF 9/7
 * @param1: whole image
 * @param2: image width
 * @param3: image height
 * @param4: step size
 * return
 */
void wavPass97RB(double **matrix,int w,int h,int step){
//horizontal/vertical lifting
  int predX;
  int succX;
  int predY;
  int succY;
  
  for(int i=0,indH=0;i<h;i=i+step,indH++){
    for(int j=0,indW=0;j<w;j=j+step,indW++){
      if((indH % 2)!=(indW % 2)){
        if((i-step)<0) predY=i+step;
        else predY=i-step;
        if((i+step)>=h) succY=i-step;
        else succY=i+step;

        if((j-step)<0) predX=j+step;
        else predX=j-step;
        if((j+step)>=w) succX=j-step;
        else succX=j+step;
        matrix[i][j]=matrix[i][j]+alfa97RB*(matrix[predY][j]+matrix[succY][j]+matrix[i][predX]+matrix[i][succX]);
      }
    }
  }
  for(int i=0,indH=0;i<h;i=i+step,indH++){
    for(int j=0,indW=0;j<w;j=j+step,indW++){
      if((indH % 2)==(indW % 2)){
        if((i-step)<0) predY=i+step;
        else predY=i-step;
        if((i+step)>=h) succY=i-step;
        else succY=i+step;

        if((j-step)<0) predX=j+step;
        else predX=j-step;
        if((j+step)>=w) succX=j-step;
        else succX=j+step;
        matrix[i][j]=matrix[i][j]+beta97RB*(matrix[predY][j]+matrix[succY][j]+matrix[i][predX]+matrix[i][succX]);
      }
    }
  }

  for(int i=0,indH=0;i<h;i=i+step,indH++){
    for(int j=0,indW=0;j<w;j=j+step,indW++){
      if((indH % 2)!=(indW % 2)){
        if((i-step)<0) predY=i+step;
        else predY=i-step;
        if((i+step)>=h) succY=i-step;
        else succY=i+step;

        if((j-step)<0) predX=j+step;
        else predX=j-step;
        if((j+step)>=w) succX=j-step;
        else succX=j+step;
        matrix[i][j]=matrix[i][j]+gama97RB*(matrix[predY][j]+matrix[succY][j]+matrix[i][predX]+matrix[i][succX]);
      }
    }
  }
  for(int i=0,indH=0;i<h;i=i+step,indH++){
    for(int j=0,indW=0;j<w;j=j+step,indW++){
      if((indH % 2)==(indW % 2)){
        if((i-step)<0) predY=i+step;
        else predY=i-step;
        if((i+step)>=h) succY=i-step;
        else succY=i+step;

        if((j-step)<0) predX=j+step;
        else predX=j-step;
        if((j+step)>=w) succX=j-step;
        else succX=j+step;
        matrix[i][j]=matrix[i][j]+delta97RB*(matrix[predY][j]+matrix[succY][j]+matrix[i][predX]+matrix[i][succX]);
      }
    }
  }
  //scaling
  for(int i=0,indH=0;i<h;i=i+step,indH++){
    for(int j=0,indW=0;j<w;j=j+step,indW++){
      if((indH % 2)==(indW % 2)){
        matrix[i][j]=matrix[i][j]*scale97RB;
      }
      else if((indH % 2)!=(indW % 2)){
        matrix[i][j]=matrix[i][j]/(-scale97RB);
      } 
    }
  }
  
  //diagonal lifting


  for(int i=0,indH=0;i<h;i=i+step,indH++){
    for(int j=0,indW=0;j<w;j=j+step,indW++){
      if(((indH % 2)==1)&&((indW % 2)==1)){
        if((i-step)<0) predY=i+step;
        else predY=i-step;
        if((i+step)>=h) succY=i-step;
        else succY=i+step;

        if((j-step)<0) predX=j+step;
        else predX=j-step;
        if((j+step)>=w) succX=j-step;
        else succX=j+step;
        matrix[i][j]=matrix[i][j]+alfa97RB*(matrix[predY][predX]+matrix[predY][succX]+matrix[succY][predX]+matrix[succY][succX]);
      }
    }
  }
  for(int i=0,indH=0;i<h;i=i+step,indH++){
    for(int j=0,indW=0;j<w;j=j+step,indW++){
      if(((indH % 2)==0)&&((indW % 2)==0)){
        if((i-step)<0) predY=i+step;
        else predY=i-step;
        if((i+step)>=h) succY=i-step;
        else succY=i+step;

        if((j-step)<0) predX=j+step;
        else predX=j-step;
        if((j+step)>=w) succX=j-step;
        else succX=j+step;
        matrix[i][j]=matrix[i][j]+beta97RB*(matrix[predY][predX]+matrix[predY][succX]+matrix[succY][predX]+matrix[succY][succX]);
      }
    }
  }
  
  for(int i=0,indH=0;i<h;i=i+step,indH++){
    for(int j=0,indW=0;j<w;j=j+step,indW++){
      if(((indH % 2)==1)&&((indW % 2)==1)){
        if((i-step)<0) predY=i+step;
        else predY=i-step;
        if((i+step)>=h) succY=i-step;
        else succY=i+step;

        if((j-step)<0) predX=j+step;
        else predX=j-step;
        if((j+step)>=w) succX=j-step;
        else succX=j+step;
        matrix[i][j]=matrix[i][j]+gama97RB*(matrix[predY][predX]+matrix[predY][succX]+matrix[succY][predX]+matrix[succY][succX]);
      }
    }
  }
  for(int i=0,indH=0;i<h;i=i+step,indH++){
    for(int j=0,indW=0;j<w;j=j+step,indW++){
      if(((indH % 2)==0)&&((indW % 2)==0)){
        if((i-step)<0) predY=i+step;
        else predY=i-step;
        if((i+step)>=h) succY=i-step;
        else succY=i+step;

        if((j-step)<0) predX=j+step;
        else predX=j-step;
        if((j+step)>=w) succX=j-step;
        else succX=j+step;
        matrix[i][j]=matrix[i][j]+delta97RB*(matrix[predY][predX]+matrix[predY][succX]+matrix[succY][predX]+matrix[succY][succX]);
      }
    }
  }
  //scaling
  for(int i=0,indH=0;i<h;i=i+step,indH++){
    for(int j=0,indW=0;j<w;j=j+step,indW++){
      if(((indH % 2)==0)&&((indW % 2)==0)){
        matrix[i][j]=matrix[i][j]*scale97RB;
      }
      else if(((indH % 2)==1)&&((indW % 2)==1)){
        matrix[i][j]=matrix[i][j]/(-scale97RB);
      }
    }
  }
  return;
}



/* computation of separable and modified CDF 5/3
 * @param1: whole image
 * @param2: start row
 * @param3: start col
 * @param4: image width
 * @param5: image height
 * @param6: step size
 * return
 */
void wavPass(double **matrix,int r,int c,int w,int h,int step){
double scale=sqrt(2);
scale=1.23;
//process rows
  for(int i=r;i<h;i=i+step){
    for(int j=c,ind=0;j<w;j=j+step,ind++){
      if(isOdd(ind)){
        int pred;
        int succ;
        if((j+step)>=w)succ=j-step;
        else succ=j+step;
        pred=j-step;
        matrix[i][j]=matrix[i][j]+matrix[i][pred]*alfa+matrix[i][succ]*alfa;
      }
    }
  }
  for(int i=r;i<h;i=i+step){
    for(int j=c,ind=0;j<w;j=j+step,ind++){
      if(isEven(ind)){
        int pred;
        int succ;
        if((j-step)<0)pred=j+step;
        else pred=j-step;
        succ=j+step;
        matrix[i][j]=matrix[i][j]+matrix[i][pred]*beta+matrix[i][succ]*beta;
      }
    }
  }
  //scaling
  for(int i=r;i<h;i=i+step){
    for(int j=c,ind=0;j<w;j=j+step,ind++){
      if(isEven(ind)){
        matrix[i][j]=matrix[i][j]*scale;
      }
      else{
        matrix[i][j]=matrix[i][j]/(-scale);
      }
    }
  }
//process columns
  for(int i=c;i<w;i=i+step){
    for(int j=r,ind=0;j<h;j=j+step,ind++){
      if(isOdd(ind)){
        int pred;
        int succ;
        if((j+step)>=h)succ=j-step;
        else succ=j+step;
        pred=j-step;
        matrix[j][i]=matrix[j][i]+matrix[pred][i]*alfa+matrix[succ][i]*alfa;
      }
    }
  }
  
  for(int i=c;i<w;i=i+step){
    for(int j=r,ind=0;j<h;j=j+step,ind++){
      if(isEven(ind)){
        int pred;
        int succ;
        if((j-step)<0)pred=j+step;
        else pred=j-step;
        succ=j+step;
        matrix[j][i]=matrix[j][i]+matrix[pred][i]*beta+matrix[succ][i]*beta;
      }
    }
  }
  //scaling
  for(int i=c;i<w;i=i+step){
    for(int j=r,ind=0;j<h;j=j+step,ind++){
      if(isEven(ind)){
        matrix[j][i]=matrix[j][i]*scale;
      }
      else{
        matrix[j][i]=matrix[j][i]/(-scale);
      }
    }
  }
  return;
}


/* computation of separable and modified CDF 9/7
 * @param1: whole image
 * @param2: start row
 * @param3: start col
 * @param4: image width
 * @param5: image height
 * @param6: step size
 * return
 */
void wavPass97(double **matrix,int r,int c,int w,int h,int step){
//process rows
  int pred;
  int succ;
  for(int i=r;i<h;i=i+step){
    for(int j=c,ind=0;j<w;j=j+step,ind++){
      if(isOdd(ind)){
        if((j+step)>=w)succ=j-step;
        else succ=j+step;
        pred=j-step;
        matrix[i][j]=matrix[i][j]+matrix[i][pred]*alfa97+matrix[i][succ]*alfa97;
      }
    }
  }
  for(int i=r;i<h;i=i+step){
    for(int j=c,ind=0;j<w;j=j+step,ind++){
      if(isEven(ind)){
        if((j-step)<0)pred=j+step;
        else pred=j-step;
        succ=j+step;
        matrix[i][j]=matrix[i][j]+matrix[i][pred]*beta97+matrix[i][succ]*beta97;
      }
    }
  }
  for(int i=r;i<h;i=i+step){
    for(int j=c,ind=0;j<w;j=j+step,ind++){
      if(isOdd(ind)){
        if((j+step)>=w)succ=j-step;
        else succ=j+step;
        pred=j-step;
        matrix[i][j]=matrix[i][j]+matrix[i][pred]*gama97+matrix[i][succ]*gama97;
      }
    }
  }
  for(int i=r;i<h;i=i+step){
    for(int j=c,ind=0;j<w;j=j+step,ind++){
      if(isEven(ind)){
        if((j-step)<0)pred=j+step;
        else pred=j-step;
        succ=j+step;
        matrix[i][j]=matrix[i][j]+matrix[i][pred]*delta97+matrix[i][succ]*delta97;
      }
    }
  }
  //scaling
  for(int i=r;i<h;i=i+step){
    for(int j=c,ind=0;j<w;j=j+step,ind++){
      if(isEven(ind)){
        matrix[i][j]=matrix[i][j]*scale97;
      }
      else{
        matrix[i][j]=matrix[i][j]/(-scale97);
      }
    }
  }
//process columns
  for(int i=c;i<w;i=i+step){
    for(int j=r,ind=0;j<h;j=j+step,ind++){
      if(isOdd(ind)){
        if((j+step)>=h)succ=j-step;
        else succ=j+step;
        pred=j-step;
        matrix[j][i]=matrix[j][i]+matrix[pred][i]*alfa97+matrix[succ][i]*alfa97;
      }
    }
  }
  
  for(int i=c;i<w;i=i+step){
    for(int j=r,ind=0;j<h;j=j+step,ind++){
      if(isEven(ind)){
        if((j-step)<0)pred=j+step;
        else pred=j-step;
        succ=j+step;
        matrix[j][i]=matrix[j][i]+matrix[pred][i]*beta97+matrix[succ][i]*beta97;
      }
    }
  }
  for(int i=c;i<w;i=i+step){
    for(int j=r,ind=0;j<h;j=j+step,ind++){
      if(isOdd(ind)){
        if((j+step)>=h)succ=j-step;
        else succ=j+step;
        pred=j-step;
        matrix[j][i]=matrix[j][i]+matrix[pred][i]*gama97+matrix[succ][i]*gama97;
      }
    }
  }
  
  for(int i=c;i<w;i=i+step){
    for(int j=r,ind=0;j<h;j=j+step,ind++){
      if(isEven(ind)){
        if((j-step)<0)pred=j+step;
        else pred=j-step;
        succ=j+step;
        matrix[j][i]=matrix[j][i]+matrix[pred][i]*delta97+matrix[succ][i]*delta97;
      }
    }
  }
  //scaling
  for(int i=c;i<w;i=i+step){
    for(int j=r,ind=0;j<h;j=j+step,ind++){
      if(isEven(ind)){
        matrix[j][i]=matrix[j][i]*scale97;
      }
      else{
        matrix[j][i]=matrix[j][i]/(-scale97);
      }
    }
  }
  return;
}



/* get new input value, neighbor block(for LDCT precomputation)
 * @param1: value from actual block
 * @param2: value from neighbor block
 * @param3: index offset
 * return new input value
 */
double getNegFX(double fst,double sec,int pos){
  double res=0.0;
  switch(pos){
    case -1:
      res=(betaVec[3]*sec-betaVec[4]*fst)/(betaVec[3]-betaVec[4]);
      break;
    case -2:
      res=(betaVec[2]*sec-betaVec[5]*fst)/(betaVec[2]-betaVec[5]);
      break;
    case -3:
      res=(betaVec[1]*sec-betaVec[6]*fst)/(betaVec[1]-betaVec[6]);
      break;
    case -4:
      res=(betaVec[0]*sec-betaVec[7]*fst)/(betaVec[0]-betaVec[7]);
      break;

  }
  return res;
}

/* get new input value, actual block(for LDCT precomputation)
 * @param1: value from actual block
 * @param2: value from neighbor block
 * @param3: index offset
 * return new input value
 */
double getPosFX(double fst,double sec,int pos){
  double res=0.0;
  switch(pos){
    case 1:
      res=(betaVec[4]*fst-betaVec[3]*sec)/(betaVec[4]-betaVec[3]);
      break;
    case 2:
      res=(betaVec[5]*fst-betaVec[2]*sec)/(betaVec[5]-betaVec[2]);
      break;
    case 3:
      res=(betaVec[6]*fst-betaVec[1]*sec)/(betaVec[6]-betaVec[1]);
      break;
    case 4:
      res=(betaVec[7]*fst-betaVec[0]*sec)/(betaVec[7]-betaVec[0]);
      break;

  }
  return res;
}

/* LDCT precomputation (whole image)
 * @param1: whole image matrix
 * @param3: input matrix width
 * @param2: input matrix height
 */
void localdct(double **matrix,int w,int h){
  double tempMat[8];
  //rows
  for(int i=0;i<h;i=i+8){
    for(int j=0;j<w;j=j+8){
      for(int k=0;k<8;k++){
        if((j-4)>=0){
          tempMat[0]=getNegFX(matrix[i+k][j-4],matrix[i+k][j+3],-4);
          tempMat[1]=getNegFX(matrix[i+k][j-3],matrix[i+k][j+2],-3);
          tempMat[2]=getNegFX(matrix[i+k][j-2],matrix[i+k][j+1],-2);
          tempMat[3]=getNegFX(matrix[i+k][j-1],matrix[i+k][j],-1);
          tempMat[4]=getPosFX(matrix[i+k][j],matrix[i+k][j-1],1);
          tempMat[5]=getPosFX(matrix[i+k][j+1],matrix[i+k][j-2],2);
          tempMat[6]=getPosFX(matrix[i+k][j+2],matrix[i+k][j-3],3);
          tempMat[7]=getPosFX(matrix[i+k][j+3],matrix[i+k][j-4],4);
        
          matrix[i+k][j-4]=tempMat[0];
          matrix[i+k][j-3]=tempMat[1];
          matrix[i+k][j-2]=tempMat[2];
          matrix[i+k][j-1]=tempMat[3];
          matrix[i+k][j]=tempMat[4];
          matrix[i+k][j+1]=tempMat[5];
          matrix[i+k][j+2]=tempMat[6];
          matrix[i+k][j+3]=tempMat[7];
        }
      } 
    }
  }
  //columns
  for(int i=0;i<h;i=i+8){
    for(int j=0;j<w;j=j+8){
      for(int k=0;k<8;k++){
        if((i-4)>=0){
          tempMat[0]=getNegFX(matrix[i-4][j+k],matrix[i+3][j+k],-4);
          tempMat[1]=getNegFX(matrix[i-3][j+k],matrix[i+2][j+k],-3);
          tempMat[2]=getNegFX(matrix[i-2][j+k],matrix[i+1][j+k],-2);
          tempMat[3]=getNegFX(matrix[i-1][j+k],matrix[i][j+k],-1);
          tempMat[4]=getPosFX(matrix[i][j+k],matrix[i-1][j+k],1);
          tempMat[5]=getPosFX(matrix[i+1][j+k],matrix[i-2][j+k],2);
          tempMat[6]=getPosFX(matrix[i+2][j+k],matrix[i-3][j+k],3);
          tempMat[7]=getPosFX(matrix[i+3][j+k],matrix[i-4][j+k],4);
 
          matrix[i-4][j+k]=tempMat[0];
          matrix[i-3][j+k]=tempMat[1];
          matrix[i-2][j+k]=tempMat[2];
          matrix[i-1][j+k]=tempMat[3];
          matrix[i][j+k]=tempMat[4];
          matrix[i+1][j+k]=tempMat[5];
          matrix[i+2][j+k]=tempMat[6];
          matrix[i+3][j+k]=tempMat[7];
        }
      }   
    }
  }
}


/* DCT computation (whole image)
 * @param1: whole image matrix
 * @param3: input matrix width
 * @param2: input matrix height
 */
void DCT_II(double **matrix,int w,int h){

  int length=8;
  for(int i=0;i<8;i++){
    for(int j=0;j<8;j++){

      fastMatrix[i][j]=cos((i*PI*(2*j+1))/((double)length*2));
    }
  }

  double sum=0.0;  
  double li=0;
  double lj=0;
  double tempMat[8][8];

  for(int r=0;r<h;r=r+8){
    for(int c=0;c<w;c=c+8){
      //------------DCT--------------
      for(int i=0;i<8;i++){
        for(int j=0;j<8;j++){
          if(i==0) li=1/(sqrt(2));
          else li=1;
          if(j==0) lj=1/(sqrt(2));
          else lj=1;
          for(int n1=0;n1<length;n1++){
            for(int n2=0;n2<length;n2++){
              sum=sum+(matrix[n1+r][n2+c]-128)*0.25*li*lj*fastMatrix[i][n1]*fastMatrix[j][n2];
            }
          }
          //sum=sum*8;
          tempMat[i][j]=sum;
          sum=0;
        }
      }
      for(int i=0;i<8;i++){
        for(int j=0;j<8;j++){
          matrix[i+r][j+c]=tempMat[i][j];
        }
      }
      //------------DCT END--------------
    }
  }
}



#define EIGHT 8
/* Modification DHT coefficient position
 * @param1: matrix for modification
 */
void prepMatHartley(double tempMat[EIGHT][EIGHT]){
  double pom[8];
  for(int i=0;i<8;i++){
    for(int j=0;j<8;j++){
      pom[j]=tempMat[i][j]; 
    }
    tempMat[i][2]=pom[7];
    tempMat[i][3]=pom[2];
    tempMat[i][4]=pom[6];
    tempMat[i][6]=pom[3];
    tempMat[i][7]=pom[4];
  }
  for(int i=0;i<8;i++){
    for(int j=0;j<8;j++){
      pom[j]=tempMat[j][i]; 
    }
    tempMat[2][i]=pom[7];
    tempMat[3][i]=pom[2];
    tempMat[4][i]=pom[6];
    tempMat[6][i]=pom[3];
    tempMat[7][i]=pom[4];
  }
}

/* Wavelet index bit reversal
 * @param1: matrix for modification
 */
void prepMat(double tempMat[EIGHT][EIGHT]){
  double pom=0.0;
  for(int i=0;i<8;i++){
    pom=tempMat[i][1];
    tempMat[i][1]=tempMat[i][4];
    tempMat[i][4]=pom;
    pom=tempMat[i][3];
    tempMat[i][3]=tempMat[i][6];
    tempMat[i][6]=pom;
  }
  for(int i=0;i<8;i++){
    pom=tempMat[1][i];
    tempMat[1][i]=tempMat[4][i];
    tempMat[4][i]=pom;
    pom=tempMat[3][i];
    tempMat[3][i]=tempMat[6][i];
    tempMat[6][i]=pom;
  }
}

/* Get MRT constant
 * @param1: input matrix - row index 
 * @param2: input matrix - column index
 * @param3: base matrix - row index 
 * @param4: base matrix - column index
 * @param5: MRT matrix identificator
 * @return: MRT constant
 */
double getA(int k1,int k2,int n1,int n2,int p){
  int A=0;
  int pom=0;
  pom=(n1*k1+n2*k2)%8;
  if(pom==p){
    return 1.0;
  }
  else if(pom==(p+4)){
    return -1.0;
  }
  else return 0.0;
}


/* Modification SMRT coefficient position
 * @param1: matrix for modification
 */
void mySMRT(double M[8][8]){
  double p[8][8];
  for(int i=0;i<8;i++){
    for(int j=0;j<8;j++){
      p[i][j]=M[i][j];
    }
  }
  M[5][1]=p[3][2];
  M[6][1]=p[4][2];
  M[1][2]=p[2][4];
  M[2][2]=p[1][4];
  M[3][2]=p[3][3];
  M[4][2]=p[4][3];
  M[1][3]=p[3][4];
  M[2][3]=p[2][2];
  M[3][3]=p[5][1];
  M[4][3]=p[1][5];
  M[1][4]=p[4][4];
  M[2][4]=p[2][3];
  M[3][4]=p[6][1];
  M[4][4]=p[1][6];
  M[1][5]=p[1][2];
  M[1][6]=p[1][3];



}


/* Get SMRT constant
 * @param1: input matrix - row index 
 * @param2: input matrix - column index
 * @param3: MRT matrix identificator
 * @param4: unique coefficient - row index 
 * @param5: unique coefficient - column index
 * @return: is unique coefficient
 */
int isUnique(int r, int c, int p, int *rt, int *ct){
int res=0;
for(int i=0;i<8;i++){
  for(int j=0;j<8;j++){
    if((SMRT_X[i][j]==r)&&(SMRT_Y[i][j]==c)&&(SMRT_p[i][j]==p)){
      *rt=i;
      *ct=j;
      res=1;
      return res;
    }
  }
}
return res;
}



