/*   jidctint.h
Added 2016-2017 by Stanislav Svoboda
  - added functions for new transforms (INVERSE)
Copyright (C) 2016, 2017 Stanislav Svoboda
For conditions of distribution and use, see the accompanying LICENSE.md file.
*/

#include <stdio.h>
#include <stdlib.h>
#include "trans.h"
#include <math.h>
#define PI 3.14159 


int Instate=0;
int Iwpos=0;
int Ihpos=0;
int IwCHpos=0;
int IhCHpos=0;

/*nonseparable CDF 5/3 and CDF 9/7 constants*/
double IalfaRB=-0.25;
double IbetaRB=0.125;
const double Ialfa97RB=-1.586134342/2;
const double Ibeta97RB=-0.05298011854/2;
const double Igama97RB=0.8829110762/2;
const double Idelta97RB=0.4435068522/2;
const double Iscale97RB=1.149604398;

/*separable CDF 5/3 and CDF 9/7 constants*/
double Ialfa=-0.5;
double Ibeta=0.25;
const double Ialfa97=-1.586134342;
const double Ibeta97=-0.05298011854;
const double Igama97=0.8829110762;
const double Idelta97=0.4435068522;
const double Iscale97=1.149604398;

int ienable=1;
double ifastMatrix[8][8];
int IenableWAV=1;

double ILDCT_scale=0.8;//LDCT scale constant

//beta constant
double IbetaVec[8]={0.00000012795,0.00074839045,0.03338096005,0.2718346173,0.7281653827,0.96661904,0.9992516096,0.9999998721};



/*Scale matrices*/
double IdctScale[8][8]={{1024.00,924.25,942.36,924.25,1020.00,924.25,942.35,924.24}, 
                       {924.25,837.49,853.89,837.49,924.25,837.49,853.89,837.48}, 
                       {942.36,853.89,870.62,853.90,942.36,853.89,870.62,853.89}, 
                       {924.25,837.49,853.90,837.49,924.25,837.49,853.89,837.48}, 
                       {1020.00,924.25,942.36,924.25,1020.00,924.25,942.35,924.24}, 
                       {924.25,837.49,853.89,837.49,924.25,837.49,853.89,837.48}, 
                       {942.35,853.89,870.62,853.89,942.35,853.89,870.62,853.89}, 
                       {924.24,837.48,853.89,837.48,924.24,837.48,853.89,837.48}};

double IdstScale[8][8]={{870.03,872.34,872.82,873.04,873.17,873.26,873.33,873.39}, 
                       {872.34,873.09,873.25,873.32,873.36,873.39,873.42,873.43}, 
                       {872.82,873.25,873.34,873.38,873.40,873.42,873.43,873.44}, 
                       {873.04,873.32,873.38,873.41,873.42,873.43,873.44,873.45}, 
                       {873.17,873.36,873.40,873.42,873.43,873.44,873.44,873.45}, 
                       {873.26,873.39,873.42,873.43,873.44,873.45,873.45,873.45}, 
                       {873.33,873.42,873.43,873.44,873.44,873.45,873.45,873.45}, 
                       {873.39,873.43,873.44,873.45,873.45,873.45,873.45,873.45}};

double IwhtScale[8][8]={{1016.0,1020.0,1020.0,1020.0,1020.0,1020.0,1020.0,1020.0}, 
                       {1020.0,1020.0,1020.0,1020.0,1020.0,1020.0,1020.0,1020.0},
                       {1020.0,1020.0,1020.0,1020.0,1020.0,1020.0,1020.0,1020.0},
                       {1020.0,1020.0,1020.0,1020.0,1020.0,1020.0,1020.0,1020.0},
                       {1020.0,1020.0,1020.0,1020.0,1020.0,1020.0,1020.0,1020.0},
                       {1020.0,1020.0,1020.0,1020.0,1020.0,1020.0,1020.0,1020.0},
                       {1020.0,1020.0,1020.0,1020.0,1020.0,1020.0,1020.0,1020.0},
                       {1020.0,1020.0,1020.0,1020.0,1020.0,1020.0,1020.0,1020.0}};

double IDCHTMatrix[8][8]={{0.35355,0.35355,0.35355,0.35355,0.35355,0.35355,0.35355,0.35355},
                         {-0.54006,-0.38576,-0.23146,-0.07715,0.07715,0.23146,0.38576,0.54006},
                         {0.54006,0.07715,-0.23146,-0.38576,-0.38576,-0.23146,0.07715,0.54006},
                         {-0.43082,0.30773,0.43082,0.18464,-0.18464,-0.43082,-0.30773,0.43082},
                         {0.28204,-0.52378,-0.12087,0.36262,0.36262,-0.12087,-0.52378,0.28204},
                         {-0.14979,0.49215,-0.36377,-0.32097,0.32097,0.36377,-0.49215,0.14979},
                         {0.06155,-0.30773,0.55391,-0.30773,-0.30773,0.55391,-0.30773,0.06155},
                         {-0.01707,0.11949,-0.35846,0.59744,-0.59744,0.35846,-0.11949,0.01707}};



//modified CDF 5/3 and 9/7
static int IwavIndex[21][3]={{0,0,1},
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

//approximated DCT matrix
static double IAPPROX[8][8]={{0.35355339,0.35355339,0.35355339,0.35355339,0.35355339,0.35355339,0.35355339,0.35355339},
  {0.40824829,0.40824829,0.40824829,0.0,0.0,-0.40824829,-0.40824829,-0.40824829},
  {0.44721359,0.22360679,-0.22360679,-0.44721359,-0.44721359,-0.22360679,0.22360679,0.44721359},
  {0.40824829,0.0,-0.40824829,-0.40824829,0.40824829,0.40824829,0.0,-0.40824829},
  {0.35355339,-0.35355339,-0.35355339,0.35355339,0.35355339,-0.35355339,-0.35355339,0.35355339},
  {0.40824829,-0.40824829,0.0,0.40824829,-0.40824829,0.0,0.40824829,-0.40824829},
  {0.22360679,-0.44721359,0.44721359,-0.22360679,-0.22360679,0.44721359,-0.44721359,0.22360679},
  {0.0,-0.40824829,0.40824829,-0.40824829,0.40824829,-0.40824829,0.40824829,0.0}};

static double IAPPROXT[8][8]={{0.35355339,0.40824829,0.44721359,0.40824829,0.35355339,0.40824829,0.22360679,0.0},
  {0.35355339,0.40824829,0.22360679,0.0,-0.35355339,-0.40824829,-0.44721359,-0.40824829,},
  {0.35355339,0.40824829,-0.22360679,-0.40824829,-0.35355339,0.0,0.44721359,0.40824829,},
  {0.35355339,0.0,-0.44721359,-0.40824829,0.35355339,0.40824829,-0.22360679,-0.40824829,},
  {0.35355339,0.0,-0.44721359,0.40824829,0.35355339,-0.40824829,-0.22360679,0.40824829,},
  {0.35355339,-0.40824829,-0.22360679,0.40824829,-0.35355339,0.0,0.44721359,-0.40824829,},
  {0.35355339,-0.40824829,0.22360679,0.0,-0.35355339,0.40824829,-0.44721359,0.40824829,},
  {0.35355339,-0.40824829,0.44721359,-0.40824829,0.35355339,-0.40824829,0.22360679,0.0}};




/*whole image structures*/
typedef struct IlistMem{
  int type;
  int wpos;
  int hpos;
  struct IlistMem *next;
}*IlMem;


typedef struct IheadList{
  struct IlistMem *first;
  struct IlistMem *act;
}*IhLst;


struct IheadList *IinitHead(){
  struct IheadList *newHead;
  newHead=malloc(sizeof(struct IheadList));
  newHead->first=NULL;
  newHead->act=NULL;
  return newHead;
}

struct IheadList *IHList=NULL;

void IinsertIntoList(int type,int w,int h){
  struct IlistMem *newMem;
  newMem=malloc(sizeof(struct IlistMem));
  newMem->type=type;
  newMem->wpos=w;
  newMem->hpos=h;
  newMem->next=NULL;
  if(IHList->first==NULL){
    IHList->first=newMem;
    IHList->act=newMem;
  }
  else{
    IHList->act->next=newMem;
    IHList->act=newMem;
  }
}






#define IEIGHT 8
/* Wavelet index bit reversal
 * @param1: matrix for modification
 */
void IprepMat(double tempMat[IEIGHT][IEIGHT]){
  double pom=0.0;
  for(int i=0;i<8;i++){
    pom=tempMat[1][i];
    tempMat[1][i]=tempMat[4][i];
    tempMat[4][i]=pom;
    pom=tempMat[3][i];
    tempMat[3][i]=tempMat[6][i];
    tempMat[6][i]=pom;
  }
  for(int i=0;i<8;i++){
    pom=tempMat[i][1];
    tempMat[i][1]=tempMat[i][4];
    tempMat[i][4]=pom;
    pom=tempMat[i][3];
    tempMat[i][3]=tempMat[i][6];
    tempMat[i][6]=pom;
  }
}



/* computation of nonseparable CDF 5/3
 * @param1: whole image
 * @param2: image width
 * @param3: image height
 * @param4: step size
 * return
 */
void wavPassInvRB(double **matrix,int w,int h,int step){
  double scaleRB=sqrt(2);

  int predX;
  int succX;
  int predY;
  int succY;
  //scaling  
  for(int i=0,indH=0;i<h;i=i+step,indH++){
    for(int j=0,indW=0;j<w;j=j+step,indW++){
      if(((indH % 2)==0)&&((indW % 2)==0)){
        matrix[i][j]=matrix[i][j]/scaleRB;
      }
      else if(((indH % 2)==1)&&((indW % 2)==1)){
        matrix[i][j]=matrix[i][j]*(-scaleRB);
      }
    }
  }
  //diagonal lifting
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
        matrix[i][j]=matrix[i][j]-IbetaRB*(matrix[predY][predX]+matrix[predY][succX]+matrix[succY][predX]+matrix[succY][succX]);
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
        matrix[i][j]=matrix[i][j]-IalfaRB*(matrix[predY][predX]+matrix[predY][succX]+matrix[succY][predX]+matrix[succY][succX]);
      }
    }
  }

    //scaling
  for(int i=0,indH=0;i<h;i=i+step,indH++){
    for(int j=0,indW=0;j<w;j=j+step,indW++){
      if((indH % 2)==(indW % 2)){
        matrix[i][j]=matrix[i][j]/scaleRB;
      }
      else if((indH % 2)!=(indW % 2)){
        matrix[i][j]=matrix[i][j]*(-scaleRB);
      } 
    }
  }
  //horizontal/vertical lifting
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
        matrix[i][j]=matrix[i][j]-IbetaRB*(matrix[predY][j]+matrix[succY][j]+matrix[i][predX]+matrix[i][succX]);
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
        matrix[i][j]=matrix[i][j]-IalfaRB*(matrix[predY][j]+matrix[succY][j]+matrix[i][predX]+matrix[i][succX]);
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
void wavPassInv97RB(double **matrix,int w,int h,int step){

  int predX;
  int succX;
  int predY;
  int succY;
  
//scaling
  for(int i=0,indH=0;i<h;i=i+step,indH++){
    for(int j=0,indW=0;j<w;j=j+step,indW++){
      if(((indH % 2)==0)&&((indW % 2)==0)){
        matrix[i][j]=matrix[i][j]/Iscale97RB;
      }
      else if(((indH % 2)==1)&&((indW % 2)==1)){
        matrix[i][j]=matrix[i][j]*(-Iscale97RB);
      }
    }
  }
  //diagonal lifting
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
        matrix[i][j]=matrix[i][j]-Idelta97RB*(matrix[predY][predX]+matrix[predY][succX]+matrix[succY][predX]+matrix[succY][succX]);
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
        matrix[i][j]=matrix[i][j]-Igama97RB*(matrix[predY][predX]+matrix[predY][succX]+matrix[succY][predX]+matrix[succY][succX]);
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
        matrix[i][j]=matrix[i][j]-Ibeta97RB*(matrix[predY][predX]+matrix[predY][succX]+matrix[succY][predX]+matrix[succY][succX]);
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
        matrix[i][j]=matrix[i][j]-Ialfa97RB*(matrix[predY][predX]+matrix[predY][succX]+matrix[succY][predX]+matrix[succY][succX]);
      }
    }
  }


  //scaling
  for(int i=0,indH=0;i<h;i=i+step,indH++){
    for(int j=0,indW=0;j<w;j=j+step,indW++){
      if((indH % 2)==(indW % 2)){
        matrix[i][j]=matrix[i][j]/Iscale97RB;
      }
      else if((indH % 2)!=(indW % 2)){
        matrix[i][j]=matrix[i][j]*(-Iscale97RB);
      } 
    }
  }
  //horizontal/vertical
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
        matrix[i][j]=matrix[i][j]-Idelta97RB*(matrix[predY][j]+matrix[succY][j]+matrix[i][predX]+matrix[i][succX]);
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
        matrix[i][j]=matrix[i][j]-Igama97RB*(matrix[predY][j]+matrix[succY][j]+matrix[i][predX]+matrix[i][succX]);
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
        matrix[i][j]=matrix[i][j]-Ibeta97RB*(matrix[predY][j]+matrix[succY][j]+matrix[i][predX]+matrix[i][succX]);
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
        matrix[i][j]=matrix[i][j]-Ialfa97RB*(matrix[predY][j]+matrix[succY][j]+matrix[i][predX]+matrix[i][succX]);
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
void wavPassInv(double **matrix,int r,int c,int w,int h,int step){
  //scaling
  double scale=sqrt(2);
  scale=1.23;
  int pred;
  int succ;
  for(int i=c;i<w;i=i+step){
    for(int j=r,ind=0;j<h;j=j+step,ind++){
      if(isEven(ind)){
        matrix[j][i]=matrix[j][i]/scale;
      }
      else{
        matrix[j][i]=matrix[j][i]*(-scale);
      }
    }
  }
  //process columns
  for(int i=c;i<w;i=i+step){
    for(int j=r,ind=0;j<h;j=j+step,ind++){
      if(isEven(ind)){
        if((j-step)<0)pred=j+step;
        else pred=j-step;
        succ=j+step;
        matrix[j][i]=matrix[j][i]-matrix[pred][i]*Ibeta-matrix[succ][i]*Ibeta;
      }
    }
  }
  for(int i=c;i<w;i=i+step){
    for(int j=r,ind=0;j<h;j=j+step,ind++){
      if(isOdd(ind)){
        if((j+step)>=h)succ=j-step;
        else succ=j+step;
        pred=j-step;
        matrix[j][i]=matrix[j][i]-matrix[pred][i]*Ialfa-matrix[succ][i]*Ialfa;
      }
    }
  }


  //scaling
  for(int i=r;i<h;i=i+step){
    for(int j=c,ind=0;j<w;j=j+step,ind++){
      if(isEven(ind)){
        matrix[i][j]=matrix[i][j]/scale;
      }
      else{
        matrix[i][j]=matrix[i][j]*(-scale);
      }
    }
  }
  //process rows
  for(int i=r;i<h;i=i+step){
    for(int j=c,ind=0;j<w;j=j+step,ind++){
      if(isEven(ind)){
        if((j-step)<0)pred=j+step;
        else pred=j-step;
        succ=j+step;
        matrix[i][j]=matrix[i][j]-matrix[i][pred]*Ibeta-matrix[i][succ]*Ibeta;
      }
    }
  }
  for(int i=r;i<h;i=i+step){
    for(int j=c,ind=0;j<w;j=j+step,ind++){
      if(isOdd(ind)){
        if((j+step)>=w)succ=j-step;
        else succ=j+step;
        pred=j-step;
        matrix[i][j]=matrix[i][j]-matrix[i][pred]*Ialfa-matrix[i][succ]*Ialfa;
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
void wavPassInv97(double **matrix,int r,int c,int w,int h,int step){
  int pred;
  int succ;
  //scaling
  for(int i=c;i<w;i=i+step){
    for(int j=r,ind=0;j<h;j=j+step,ind++){
      if(isEven(ind)){
        matrix[j][i]=matrix[j][i]/Iscale97;
      }
      else{
        matrix[j][i]=matrix[j][i]*(-Iscale97);
      }
    }
  }
  //process columns
  for(int i=c;i<w;i=i+step){
    for(int j=r,ind=0;j<h;j=j+step,ind++){
      if(isEven(ind)){
        if((j-step)<0)pred=j+step;
        else pred=j-step;
        succ=j+step;
        matrix[j][i]=matrix[j][i]-matrix[pred][i]*Idelta97-matrix[succ][i]*Idelta97;
      }
    }
  }
  for(int i=c;i<w;i=i+step){
    for(int j=r,ind=0;j<h;j=j+step,ind++){
      if(isOdd(ind)){
        if((j+step)>=h)succ=j-step;
        else succ=j+step;
        pred=j-step;
        matrix[j][i]=matrix[j][i]-matrix[pred][i]*Igama97-matrix[succ][i]*Igama97;
      }
    }
  }
  for(int i=c;i<w;i=i+step){
    for(int j=r,ind=0;j<h;j=j+step,ind++){
      if(isEven(ind)){
        if((j-step)<0)pred=j+step;
        else pred=j-step;
        succ=j+step;
        matrix[j][i]=matrix[j][i]-matrix[pred][i]*Ibeta97-matrix[succ][i]*Ibeta97;
      }
    }
  }
  for(int i=c;i<w;i=i+step){
    for(int j=r,ind=0;j<h;j=j+step,ind++){
      if(isOdd(ind)){
        if((j+step)>=h)succ=j-step;
        else succ=j+step;
        pred=j-step;
        matrix[j][i]=matrix[j][i]-matrix[pred][i]*Ialfa97-matrix[succ][i]*Ialfa97;
      }
    }
  }
  //scaling
  for(int i=r;i<h;i=i+step){
    for(int j=c,ind=0;j<w;j=j+step,ind++){
      if(isEven(ind)){
        matrix[i][j]=matrix[i][j]/Iscale97;
      }
      else{
        matrix[i][j]=matrix[i][j]*(-Iscale97);
      }
    }
  }
  //process rows
  for(int i=r;i<h;i=i+step){
    for(int j=c,ind=0;j<w;j=j+step,ind++){
      if(isEven(ind)){
        if((j-step)<0)pred=j+step;
        else pred=j-step;
        succ=j+step;
        matrix[i][j]=matrix[i][j]-matrix[i][pred]*Idelta97-matrix[i][succ]*Idelta97;
      }
    }
  }
  for(int i=r;i<h;i=i+step){
    for(int j=c,ind=0;j<w;j=j+step,ind++){
      if(isOdd(ind)){
        if((j+step)>=w)succ=j-step;
        else succ=j+step;
        pred=j-step;
        matrix[i][j]=matrix[i][j]-matrix[i][pred]*Igama97-matrix[i][succ]*Igama97;
      }
    }
  }
  for(int i=r;i<h;i=i+step){
    for(int j=c,ind=0;j<w;j=j+step,ind++){
      if(isEven(ind)){
        if((j-step)<0)pred=j+step;
        else pred=j-step;
        succ=j+step;
        matrix[i][j]=matrix[i][j]-matrix[i][pred]*Ibeta97-matrix[i][succ]*Ibeta97;
      }
    }
  }
  for(int i=r;i<h;i=i+step){
    for(int j=c,ind=0;j<w;j=j+step,ind++){
      if(isOdd(ind)){
        if((j+step)>=w)succ=j-step;
        else succ=j+step;
        pred=j-step;
        matrix[i][j]=matrix[i][j]-matrix[i][pred]*Ialfa97-matrix[i][succ]*Ialfa97;
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
double IgetNegFX(double fst,double sec,int pos){
  double res=0.0;
  switch(pos){
    case -1:
      res=IbetaVec[3]*sec+IbetaVec[4]*fst;
      break;
    case -2:
      res=IbetaVec[2]*sec+IbetaVec[5]*fst;
      break;
    case -3:
      res=IbetaVec[1]*sec+IbetaVec[6]*fst;
      break;
    case -4:
      res=IbetaVec[0]*sec+IbetaVec[7]*fst;
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
double IgetPosFX(double fst,double sec,int pos){
  double res=0.0;
  switch(pos){
    case 1:
      res=IbetaVec[4]*fst+IbetaVec[3]*sec;
      break;
    case 2:
      res=IbetaVec[5]*fst+IbetaVec[2]*sec;
      break;
    case 3:
      res=IbetaVec[6]*fst+IbetaVec[1]*sec;
      break;
    case 4:
      res=IbetaVec[7]*fst+IbetaVec[0]*sec;
      break;

  }
  return res;
}

/* LDCT precomputation (whole image)
 * @param1: whole image matrix
 * @param3: input matrix width
 * @param2: input matrix height
 */
void Ilocaldct(double **matrix,int w,int h){
 double tempMat[8];
  //rows
  for(int i=0;i<h;i=i+8){
    for(int j=0;j<w;j=j+8){
      for(int k=0;k<8;k++){
        if((i-4)>=0){
          tempMat[0]=IgetNegFX(matrix[i-4][j+k],matrix[i+3][j+k],-4);
          tempMat[1]=IgetNegFX(matrix[i-3][j+k],matrix[i+2][j+k],-3);
          tempMat[2]=IgetNegFX(matrix[i-2][j+k],matrix[i+1][j+k],-2);
          tempMat[3]=IgetNegFX(matrix[i-1][j+k],matrix[i][j+k],-1);
          tempMat[4]=IgetPosFX(matrix[i][j+k],matrix[i-1][j+k],1);
          tempMat[5]=IgetPosFX(matrix[i+1][j+k],matrix[i-2][j+k],2);
          tempMat[6]=IgetPosFX(matrix[i+2][j+k],matrix[i-3][j+k],3);
          tempMat[7]=IgetPosFX(matrix[i+3][j+k],matrix[i-4][j+k],4);

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
//columns
  for(int i=0;i<h;i=i+8){
    for(int j=0;j<w;j=j+8){
      for(int k=0;k<8;k++){
        if((j-4)>=0){
          tempMat[0]=IgetNegFX(matrix[i+k][j-4],matrix[i+k][j+3],-4);
          tempMat[1]=IgetNegFX(matrix[i+k][j-3],matrix[i+k][j+2],-3);
          tempMat[2]=IgetNegFX(matrix[i+k][j-2],matrix[i+k][j+1],-2);
          tempMat[3]=IgetNegFX(matrix[i+k][j-1],matrix[i+k][j],-1);
          tempMat[4]=IgetPosFX(matrix[i+k][j],matrix[i+k][j-1],1);
          tempMat[5]=IgetPosFX(matrix[i+k][j+1],matrix[i+k][j-2],2);
          tempMat[6]=IgetPosFX(matrix[i+k][j+2],matrix[i+k][j-3],3);
          tempMat[7]=IgetPosFX(matrix[i+k][j+3],matrix[i+k][j-4],4);

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

  
}




/* Get MRT coefficients
 * @param1: UMRT matrix
 * @param2: MRT 1. matrix
 * @param3: MRT 2. matrix
 * @param4: MRT 3. matrix
 * @param5: MRT 4. matrix
 */
void getMRT(double M[8][8],double M0[8][8],double M1[8][8],double M2[8][8],double M3[8][8]){
  double pom0[8][8]={{M[0][0],M[0][1],M[0][5],M[0][1],M[0][7],M[0][1],M[0][5],M[0][1]},
                     {M[1][0],M[1][1],M[1][5],M[1][2],M[1][7],M[1][3],M[1][6],M[1][4]},
                     {M[5][0],M[5][1],M[5][5],M[6][1],M[5][7],M[5][1],M[5][6],M[6][1]},
                     {M[1][0],M[1][2],M[1][6],M[1][1],M[1][7],M[1][4],M[1][5],M[1][3]},
                     {M[7][0],M[7][1],M[7][5],M[7][1],M[7][7],M[7][1],M[7][5],M[7][1]},
                     {M[1][0],M[1][3],M[1][5],M[1][4],M[1][7],M[1][1],M[1][6],M[1][2]},
                     {M[5][0],M[6][1],M[5][6],M[5][1],M[5][7],M[6][1],M[5][5],M[5][1]},
                     {M[1][0],M[1][4],M[1][6],M[1][3],M[1][7],M[1][2],M[1][5],M[1][1]}};

  double pom1[8][8]={{0,M[0][2],0,M[0][4],0,-M[0][2],0,-M[0][4]},
                     {M[2][0],M[2][1],M[2][5],M[4][2],M[2][7],-M[2][3],M[4][6],-M[4][4]},
                     {0,M[5][2],0,M[6][4],0,-M[5][2],0,-M[6][4]},
                     {M[4][0],M[2][2],M[2][6],M[4][1],M[4][7],-M[2][4],M[4][5],-M[4][3]},
                     {0,M[7][2],0,M[7][4],0,-M[7][2],0,-M[7][4]},
                     {-M[2][0],M[2][3],-M[2][5],M[4][4],-M[2][7],-M[2][1],-M[4][6],-M[4][2]},
                     {0,M[6][2],0,M[5][4],0,-M[6][2],0,-M[5][4]},
                     {-M[4][0],M[2][4],-M[2][6],M[4][3],-M[4][7],-M[2][2],-M[4][5],-M[4][1]}};

  double pom2[8][8]={{0,M[0][3],M[0][6],-M[0][3],0,M[0][3],-M[0][6],-M[0][3]},
                     {M[3][0],M[3][1],M[3][5],-M[3][2],M[3][7],M[3][3],-M[3][6],-M[3][4]},
                     {M[6][0],M[5][3],M[6][5],-M[6][3],M[6][7],M[5][3],-M[6][6],-M[6][3]},
                     {-M[3][0],M[3][2],M[3][6],-M[3][1],-M[3][7],M[3][4],-M[3][5],-M[3][3]},
                     {0,M[7][3],M[7][6],-M[7][3],0,M[7][3],-M[7][6],-M[7][3]},
                     {M[3][0],M[3][3],M[3][5],-M[3][4],M[3][7],M[3][1],-M[3][6],-M[3][2]},
                     {-M[6][0],M[6][3],M[6][6],-M[5][3],-M[6][7],M[6][3],-M[6][5],-M[5][3]},
                     {-M[3][0],M[3][4],M[3][6],-M[3][3],-M[3][7],M[3][2],-M[3][5],-M[3][1]}};

  double pom3[8][8]={{0,M[0][4],0,M[0][2],0,-M[0][4],0,-M[0][2]},
                     {M[4][0],M[4][1],M[4][5],M[2][2],M[4][7],-M[4][3],M[2][6],-M[2][4]},
                     {0,M[5][4],0,M[6][2],0,-M[5][4],0,-M[6][2]},
                     {M[2][0],M[4][2],M[4][6],M[2][1],M[2][7],-M[4][4],M[2][5],-M[2][3]},
                     {0,M[7][4],0,M[7][2],0,-M[7][4],0,-M[7][2]},
                     {-M[4][0],M[4][3],-M[4][5],M[2][4],-M[4][7],-M[4][1],-M[2][6],-M[2][2]},
                     {0,M[6][4],0,M[5][2],0,-M[6][4],0,-M[5][2]},
                     {-M[2][0],M[4][4],-M[4][6],M[2][3],-M[2][7],-M[4][2],-M[2][5],-M[2][1]}};

  for(int i=0;i<8;i++){
    for(int j=0;j<8;j++){
      M0[i][j]=pom0[i][j];
      M1[i][j]=pom1[i][j];
      M2[i][j]=pom2[i][j];
      M3[i][j]=pom3[i][j];
    }
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
double IgetA(int k1,int k2,int n1,int n2,int p){
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
void ImySMRT(double M[8][8]){
  double p[8][8];
  for(int i=0;i<8;i++){
    for(int j=0;j<8;j++){
      p[i][j]=M[i][j];
    }
  }
  M[3][2]=p[5][1];
  M[4][2]=p[6][1];
  M[2][4]=p[1][2];
  M[1][4]=p[2][2];
  M[3][3]=p[3][2];
  M[4][3]=p[4][2];
  M[3][4]=p[1][3];
  M[2][2]=p[2][3];
  M[5][1]=p[3][3];
  M[1][5]=p[4][3];
  M[4][4]=p[1][4];
  M[2][3]=p[2][4];
  M[6][1]=p[3][4];
  M[1][6]=p[4][4];
  M[1][2]=p[1][5];
  M[1][3]=p[1][6];



}
