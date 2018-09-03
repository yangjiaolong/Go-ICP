/****************************************************************
3D Euclidean Distance Transform Class
Last modified: Feb 13, 2013

Functions computing DT are derived from "dt.c" written by 
Alexander Vasilevskiy.

Jiaolong Yang <yangjiaolong@gmail.com>
****************************************************************/

#include <math.h>
#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include <memory.h>

#include "jly_3ddt.h"

#define FALSE 0
#define TRUE 1

#define sqrt1(x) sqrt((double)x)

#define ROUND(x) (int((x)+0.5))
#define FLOOR(x) (int((x)))
#define CEIL(x) (int((x+0.99999)))

// ************************************
//  
//       D-Euclidean BEGIN
//
// ************************************


void initDE(Array3dDEucl3D & inDE){
  int z,y,x;
  for (z=0;z<inDE.Zdim;z++) 
    for(y=0;y<inDE.Ydim;y++) 
      for(x=0;x<inDE.Xdim;x++) {
	inDE.data[z][y][x].v=infty;
	inDE.data[z][y][x].h=infty;
	inDE.data[z][y][x].d=infty;
	inDE.data[z][y][x].distance=infty;
      }
  inDE.data[2][2][3].v=0;
  inDE.data[2][2][3].h=0;
  inDE.data[2][2][3].d=0;
  inDE.data[2][2][3].distance=0.0;
}

DEucl3D MINforwardDE3(Array3dDEucl3D& A, int z,int y,int x)
{
  DEucl3D*** inDE = A.data;
  int Xdim = A.Xdim;
  int Ydim = A.Ydim;
  int Zdim = A.Zdim;

  //              MASK:
  //
  //             (0,0,0)(0,0,1)
  //      (0,1,1)(0,1,0)(0,1,1)
  //
  
  DEucl3D mask[5];
  DEucl3D min;
  int looper;
  min.distance=infty;
 
  if (z<Zdim-1) {
    mask[0].v=inDE[z+1][y][x].v; 
    mask[0].h=inDE[z+1][y][x].h;
    mask[0].d=inDE[z+1][y][x].d+1;
    mask[0].distance=sqrt1(mask[0].v*mask[0].v+mask[0].h*mask[0].h+mask[0].d*mask[0].d);
  }
  else {
    mask[0].v=infty;
    mask[0].h=infty;
    mask[0].d=infty;
    mask[0].distance=infty;
  }

  if ((y<Ydim-1)&&(z<Zdim-1)) {
    mask[1].v=inDE[z+1][y][x].v; 
    mask[1].h=inDE[z+1][y][x].h+1;
    mask[1].d=inDE[z+1][y][x].d+1;
    mask[1].distance=sqrt1(mask[1].v*mask[1].v+mask[1].h*mask[1].h+mask[1].d*mask[1].d);
  }
  else {
    mask[1].v=infty;
    mask[1].h=infty;
    mask[1].d=infty;
    mask[1].distance=infty;

  }
  
  if (y<Ydim-1) {
    mask[2].v=inDE[z][y+1][x].v; 
    mask[2].h=inDE[z][y+1][x].h+1; 
    mask[2].d=inDE[z][y+1][x].d; 
    mask[2].distance=sqrt1(mask[2].v*mask[2].v+mask[2].h*mask[2].h+mask[2].d*mask[2].d);
  }
  else {
    mask[2].v=infty;
    mask[2].h=infty;
    mask[2].d=infty;
    mask[2].distance=infty;
  }
  mask[3].v=inDE[z][y][x].v;
  mask[3].h=inDE[z][y][x].h;
  mask[3].d=inDE[z][y][x].d;
  mask[3].distance=sqrt1(mask[3].v*mask[3].v+mask[3].h*mask[3].h+mask[3].d*mask[3].d);
  if ((z>0)&&(y<Ydim-1)) {
    mask[4].v=inDE[z-1][y+1][x].v;
    mask[4].h=inDE[z-1][y+1][x].h+1;
    mask[4].d=inDE[z-1][y+1][x].d+1;
    mask[4].distance=sqrt1(mask[4].v*mask[4].v+mask[4].h*mask[4].h+mask[4].d*mask[4].d);
  }
 else {
   mask[4].v=infty; 
   mask[4].h=infty;
   mask[4].d=infty;
   mask[4].distance=infty;
 }

  for(looper=0;looper<5;looper++) {
    if (mask[looper].distance<min.distance) min=mask[looper];
  }
 
  return min;
}


DEucl3D MINforwardDE4(Array3dDEucl3D& A, int z,int y,int x)
{
  DEucl3D*** inDE = A.data;
  int Xdim = A.Xdim;
  int Ydim = A.Ydim;
  int Zdim = A.Zdim;
  
  //
  //   MASK:   (0,0,1)(0,0,0)
  //
  DEucl3D mask[2];
  DEucl3D min;
  int looper;
  min.distance=infty;
  if (z>0) {
    mask[0].v=inDE[z-1][y][x].v; 
    mask[0].h=inDE[z-1][y][x].h;
    mask[0].d=inDE[z-1][y][x].d+1;
    mask[0].distance=sqrt1(mask[0].v*mask[0].v+mask[0].h*mask[0].h+mask[0].d*mask[0].d);         
  }
  else {
    mask[0].v=infty;
    mask[0].h=infty;
    mask[0].d=infty;
    mask[0].distance=infty;
  }
  mask[1].v=inDE[z][y][x].v;
  mask[1].h=inDE[z][y][x].h;
  mask[1].d=inDE[z][y][x].d; 
  mask[1].distance=sqrt1(mask[1].v*mask[1].v+mask[1].h*mask[1].h+mask[1].d*mask[1].d);
  
  for(looper=0;looper<2;looper++) {
    if (mask[looper].distance<min.distance) min=mask[looper];
  }
 
  return min;
}

DEucl3D MINbackwardDE3(Array3dDEucl3D& A, int z,int y,int x)
{
  DEucl3D*** inDE = A.data;
  int Xdim = A.Xdim;
  int Ydim = A.Ydim;
  int Zdim = A.Zdim;

  //
  //                     MASK:
  //
  //             (0,1,1)(0,1,0)(0,1,1)
  //             (0,0,1)(0,0,0)
  //
  //
  
  
  DEucl3D mask[5];
  DEucl3D min;
  float temp;
  int looper;
  min.distance=infty;
  if ((z>0)&&(y>0)) {
    mask[0].v=inDE[z-1][y-1][x].v; 
    mask[0].h=inDE[z-1][y-1][x].h+1;
    mask[0].d=inDE[z-1][y-1][x].d+1;
    mask[0].distance=sqrt1(mask[0].v*mask[0].v+mask[0].h*mask[0].h+mask[0].d*mask[0].d);
  }
  else {
    mask[0].v=infty;
    mask[0].h=infty;
    mask[0].d=infty;
    mask[0].distance=infty;
  }
  if (y>0){
    mask[1].v=inDE[z][y-1][x].v; 
    mask[1].h=inDE[z][y-1][x].h+1;
    mask[1].d=inDE[z][y-1][x].d;
    mask[1].distance=sqrt1(mask[1].v*mask[1].v+mask[1].h*mask[1].h+mask[1].d*mask[1].d);
  }
  else {
    mask[1].v=infty;
    mask[1].h=infty;
    mask[1].d=infty;
    mask[1].distance=infty;

  }
  if ((z<Zdim-1)&&(y>0)) {
    mask[2].v=inDE[z+1][y-1][x].v; 
    mask[2].h=inDE[z+1][y-1][x].h+1; 
    mask[2].d=inDE[z+1][y-1][x].d+1; 
    mask[2].distance=sqrt1(mask[2].v*mask[2].v+mask[2].h*mask[2].h+mask[2].d*mask[2].d);
  }
  else {
    mask[2].v=infty;
    mask[2].h=infty;
    mask[2].d=infty;
    mask[2].distance=infty;
  }
  mask[3].v=inDE[z][y][x].v;
  mask[3].h=inDE[z][y][x].h;
  mask[3].d=inDE[z][y][x].d;
  mask[3].distance=sqrt1(mask[3].v*mask[3].v+mask[3].h*mask[3].h+mask[3].d*mask[3].d);
  if (z>0) {
    mask[4].v=inDE[z-1][y][x].v;
    mask[4].h=inDE[z-1][y][x].h;
    mask[4].d=inDE[z-1][y][x].d+1;
    mask[4].distance=sqrt1(mask[4].v*mask[4].v+mask[4].h*mask[4].h+mask[4].d*mask[4].d);
  }
 else {
   mask[4].v=infty; 
   mask[4].h=infty;
   mask[4].d=infty;
   mask[4].distance=infty;
 }

  for(looper=0;looper<5;looper++) {
    if (mask[looper].distance<min.distance) min=mask[looper];
  }
 
  return min;
}


DEucl3D MINforwardDE2(Array3dDEucl3D& A, int z,int y,int x)
{
  DEucl3D*** inDE = A.data;
  int Xdim = A.Xdim;
  int Ydim = A.Ydim;
  int Zdim = A.Zdim;

  //
  // MASK:   (0,0,0)(0,0,1)
  //

  DEucl3D mask[2];
  DEucl3D min;
  int looper;
  min.distance=infty;
  if (z<Zdim-1) {
    mask[0].v=inDE[z+1][y][x].v; 
    mask[0].h=inDE[z+1][y][x].h;
    mask[0].d=inDE[z+1][y][x].d+1;
    mask[0].distance=sqrt1(mask[0].v*mask[0].v+mask[0].h*mask[0].h+mask[0].d*mask[0].d);         
  }
  else {
    mask[0].v=infty;
    mask[0].h=infty;
    mask[0].d=infty;
    mask[0].distance=infty;
  }
  mask[1].v=inDE[z][y][x].v;
  mask[1].h=inDE[z][y][x].h;
  mask[1].d=inDE[z][y][x].d;
  mask[1].distance=sqrt1(mask[1].v*mask[1].v+mask[1].h*mask[1].h+mask[1].d*mask[1].d);
  
  for(looper=0;looper<2;looper++) {
    if (mask[looper].distance<min.distance) min=mask[looper];
  }
 
 
  return min;
}

 
DEucl3D MINbackwardDE1(Array3dDEucl3D& A, int z,int y,int x)
{
  DEucl3D*** inDE = A.data;
  int Xdim = A.Xdim;
  int Ydim = A.Ydim;
  int Zdim = A.Zdim;

  //
  //                     MASK:
  //
  //             (1,1,1)(1,1,0)(1,1,1)
  //             (1,0,1)(1,0,0)(1,0,1)
  //             (1,1,1)(1,1,0)(1,1,1)
  //
  
  
  DEucl3D mask[14];
  DEucl3D min;
  float temp;
  int looper;
  min.distance=infty;
 
 if ((z>0)&&(y>0)&&(x<Xdim-1)) {
    mask[0].v=inDE[z-1][y-1][x+1].v+1; 
    mask[0].h=inDE[z-1][y-1][x+1].h+1;
    mask[0].d=inDE[z-1][y-1][x+1].d+1;
    mask[0].distance=sqrt1(mask[0].v*mask[0].v+mask[0].h*mask[0].h+mask[0].d*mask[0].d);
  }
  else {
    mask[0].v=infty;
    mask[0].h=infty;
    mask[0].d=infty;
    mask[0].distance=infty;
  }

  if ((y>0)&&(x<Xdim-1)){
    mask[1].v=inDE[z][y-1][x+1].v+1; 
    mask[1].h=inDE[z][y-1][x+1].h+1;
    mask[1].d=inDE[z][y-1][x+1].d;
    mask[1].distance=sqrt1(mask[1].v*mask[1].v+mask[1].h*mask[1].h+mask[1].d*mask[1].d);
  }
  else {
    mask[1].v=infty;
    mask[1].h=infty;
    mask[1].d=infty;
    mask[1].distance=infty;

  }
   if ((z<Zdim-1)&&(y>0)&&(x<Xdim-1)) {
    mask[2].v=inDE[z+1][y-1][x+1].v+1; 
    mask[2].h=inDE[z+1][y-1][x+1].h+1; 
    mask[2].d=inDE[z+1][y-1][x+1].d+1; 
    mask[2].distance=sqrt1(mask[2].v*mask[2].v+mask[2].h*mask[2].h+mask[2].d*mask[2].d);
  }
  else {
    mask[2].v=infty;
    mask[2].h=infty;
    mask[2].d=infty;
    mask[2].distance=infty;
  }
  
  if ((z>0)&&(x<Xdim-1)) {
    mask[3].v=inDE[z-1][y][x+1].v+1;
    mask[3].h=inDE[z-1][y][x+1].h;
    mask[3].d=inDE[z-1][y][x+1].d+1;
    mask[3].distance=sqrt1(mask[3].v*mask[3].v+mask[3].h*mask[3].h+mask[3].d*mask[3].d);
  }
 else {
   mask[3].v=infty; 
   mask[3].h=infty;
   mask[3].d=infty;
   mask[3].distance=infty;
 }
 
 if (x<Xdim-1) {
    mask[4].v=inDE[z][y][x+1].v+1;
    mask[4].h=inDE[z][y][x+1].h;
    mask[4].d=inDE[z][y][x+1].d;
    mask[4].distance=sqrt1(mask[4].v*mask[4].v+mask[4].h*mask[4].h+mask[4].d*mask[4].d);
  }
 else {
   mask[4].v=infty; 
   mask[4].h=infty;
   mask[4].d=infty;
   mask[4].distance=infty;
 }

 if ((x<Xdim-1)&&(z<Zdim-1)) {
    mask[5].v=inDE[z+1][y][x+1].v+1;
    mask[5].h=inDE[z+1][y][x+1].h;
    mask[5].d=inDE[z+1][y][x+1].d+1;
    mask[5].distance=sqrt1(mask[5].v*mask[5].v+mask[5].h*mask[5].h+mask[5].d*mask[5].d);
  }
 else {
   mask[5].v=infty; 
   mask[5].h=infty;
   mask[5].d=infty;
   mask[5].distance=infty;
 }

if ((x<Xdim-1)&&(z>0)&&(y<Ydim-1)) {
    mask[6].v=inDE[z-1][y+1][x+1].v+1;
    mask[6].h=inDE[z-1][y+1][x+1].h+1;
    mask[6].d=inDE[z-1][y+1][x+1].d+1;
    mask[6].distance=sqrt1(mask[6].v*mask[6].v+mask[6].h*mask[6].h+mask[6].d*mask[6].d);
  }
 else {
   mask[6].v=infty; 
   mask[6].h=infty;
   mask[6].d=infty;
   mask[6].distance=infty;
 }

if ((x<Xdim-1)&&(y<Ydim-1)) {
    mask[7].v=inDE[z][y+1][x+1].v+1;
    mask[7].h=inDE[z][y+1][x+1].h+1;
    mask[7].d=inDE[z][y+1][x+1].d;
    mask[7].distance=sqrt1(mask[7].v*mask[7].v+mask[7].h*mask[7].h+mask[7].d*mask[7].d);
  }
 else {
   mask[7].v=infty; 
   mask[7].h=infty;
   mask[7].d=infty;
   mask[7].distance=infty;
 }

if ((x<Xdim-1)&&(y<Ydim-1)&&(z<Zdim-1)) {
    mask[8].v=inDE[z+1][y+1][x+1].v+1;
    mask[8].h=inDE[z+1][y+1][x+1].h+1;
    mask[8].d=inDE[z+1][y+1][x+1].d+1;
    mask[8].distance=sqrt1(mask[8].v*mask[8].v+mask[8].h*mask[8].h+mask[8].d*mask[8].d);
  }
 else {
   mask[8].v=infty; 
   mask[8].h=infty;
   mask[8].d=infty;
   mask[8].distance=infty;
 }

  //              MASK :
  //
  //             (0,0,0)(0,0,1)
  //      (0,1,1)(0,1,0)(0,1,1)
  //
  
  if (z<Zdim-1) {
    mask[9].v=inDE[z+1][y][x].v; 
    mask[9].h=inDE[z+1][y][x].h;
    mask[9].d=inDE[z+1][y][x].d+1;
    mask[9].distance=sqrt1(mask[9].v*mask[9].v+mask[9].h*mask[9].h+mask[9].d*mask[9].d);
  }
  else {
    mask[9].v=infty;
    mask[9].h=infty;
    mask[9].d=infty;
    mask[9].distance=infty;
  }
  
  if ((y<Ydim-1)&&(z<Zdim-1)) {
    mask[10].v=inDE[z+1][y][x].v; 
    mask[10].h=inDE[z+1][y][x].h+1;
    mask[10].d=inDE[z+1][y][x].d+1;
    mask[10].distance=sqrt1(mask[10].v*mask[10].v+mask[10].h*mask[10].h+mask[10].d*mask[10].d);
  }
  else {
    mask[10].v=infty;
    mask[10].h=infty;
    mask[10].d=infty;
    mask[10].distance=infty;

  }
  if (y<Ydim-1) {
    mask[11].v=inDE[z][y+1][x].v; 
    mask[11].h=inDE[z][y+1][x].h+1; 
    mask[11].d=inDE[z][y+1][x].d; 
    mask[11].distance=sqrt1(mask[11].v*mask[11].v+mask[11].h*mask[11].h+mask[11].d*mask[11].d);
  }
  else {
    mask[11].v=infty;
    mask[11].h=infty;
    mask[11].d=infty;
    mask[11].distance=infty;
  }
  mask[12].v=inDE[z][y][x].v;
  mask[12].h=inDE[z][y][x].h;
  mask[12].d=inDE[z][y][x].d;
  mask[12].distance=sqrt1(mask[12].v*mask[12].v+mask[12].h*mask[12].h+mask[12].d*mask[12].d);
 
  if ((z>0)&&(y<Ydim-1)) {
    mask[13].v=inDE[z-1][y+1][x].v;
    mask[13].h=inDE[z-1][y+1][x].h+1;
    mask[13].d=inDE[z-1][y+1][x].d+1;
    mask[13].distance=sqrt1(mask[13].v*mask[13].v+mask[13].h*mask[13].h+mask[13].d*mask[13].d);
  }
 else {
   mask[13].v=infty; 
   mask[13].h=infty;
   mask[13].d=infty;
   mask[13].distance=infty;
 }

  for(looper=0;looper<14;looper++) {
    if (mask[looper].distance<min.distance) min=mask[looper];
  }
  return min;
}

DEucl3D MINforwardDE1(Array3dDEucl3D& A, int z,int y,int x)
{
  DEucl3D*** inDE = A.data;
  int Xdim = A.Xdim;
  int Ydim = A.Ydim;
  int Zdim = A.Zdim;

  DEucl3D mask[14];
  DEucl3D min;
  float temp;
  int looper;
  min.distance=infty;
 
  
  //                     MASK:
  //
  //             (1,1,1)(1,1,0)(1,1,1)
  //             (1,0,1)(1,0,0)(1,0,1)
  //             (1,1,1)(1,1,0)(1,1,1)
  //
 if ((z>0)&&(y>0)&&(x>0)) {
    mask[0].v=inDE[z-1][y-1][x-1].v+1; 
    mask[0].h=inDE[z-1][y-1][x-1].h+1;
    mask[0].d=inDE[z-1][y-1][x-1].d+1;
    mask[0].distance=sqrt1(mask[0].v*mask[0].v+mask[0].h*mask[0].h+mask[0].d*mask[0].d);
  }
  else {
    mask[0].v=infty;
    mask[0].h=infty;
    mask[0].d=infty;
    mask[0].distance=infty;
  }
  if ((y>0)&&(x>0)){
    mask[1].v=inDE[z][y-1][x-1].v+1; 
    mask[1].h=inDE[z][y-1][x-1].h+1;
    mask[1].d=inDE[z][y-1][x-1].d;
    mask[1].distance=sqrt1(mask[1].v*mask[1].v+mask[1].h*mask[1].h+mask[1].d*mask[1].d);
  }
  else {
    mask[1].v=infty;
    mask[1].h=infty;
    mask[1].d=infty;
    mask[1].distance=infty;

  }
  if ((z<Zdim-1)&&(y>0)&&(x>0)) {
    mask[2].v=inDE[z+1][y-1][x-1].v+1; 
    mask[2].h=inDE[z+1][y-1][x-1].h+1; 
    mask[2].d=inDE[z+1][y-1][x-1].d+1; 
    mask[2].distance=sqrt1(mask[2].v*mask[2].v+mask[2].h*mask[2].h+mask[2].d*mask[2].d);
  }
  else {
    mask[2].v=infty;
    mask[2].h=infty;
    mask[2].d=infty;
    mask[2].distance=infty;
  }
  if ((z>0)&&(x>0)) {
    mask[3].v=inDE[z-1][y][x-1].v+1;
    mask[3].h=inDE[z-1][y][x-1].h;
    mask[3].d=inDE[z-1][y][x-1].d+1;
    mask[3].distance=sqrt1(mask[3].v*mask[3].v+mask[3].h*mask[3].h+mask[3].d*mask[3].d);
  }
 else {
   mask[3].v=infty; 
   mask[3].h=infty;
   mask[3].d=infty;
   mask[3].distance=infty;
 }
 if (x>0) {
    mask[4].v=inDE[z][y][x-1].v+1;
    mask[4].h=inDE[z][y][x-1].h;
    mask[4].d=inDE[z][y][x-1].d;
    mask[4].distance=sqrt1(mask[4].v*mask[4].v+mask[4].h*mask[4].h+mask[4].d*mask[4].d);
  }
 else {
   mask[4].v=infty; 
   mask[4].h=infty;
   mask[4].d=infty;
   mask[4].distance=infty;
 }
 if ((x>0)&&(z<Zdim-1)) {
    mask[5].v=inDE[z+1][y][x-1].v+1;
    mask[5].h=inDE[z+1][y][x-1].h;
    mask[5].d=inDE[z+1][y][x-1].d+1;
    mask[5].distance=sqrt1(mask[5].v*mask[5].v+mask[5].h*mask[5].h+mask[5].d*mask[5].d);
  }
 else {
   mask[5].v=infty; 
   mask[5].h=infty;
   mask[5].d=infty;
   mask[5].distance=infty;
 }
if ((x>0)&&(z>0)&&(y<Ydim-1)) {
    mask[6].v=inDE[z-1][y+1][x-1].v+1;
    mask[6].h=inDE[z-1][y+1][x-1].h+1;
    mask[6].d=inDE[z-1][y+1][x-1].d+1;
    mask[6].distance=sqrt1(mask[6].v*mask[6].v+mask[6].h*mask[6].h+mask[6].d*mask[6].d);
  }
 else {
   mask[6].v=infty; 
   mask[6].h=infty;
   mask[6].d=infty;
   mask[6].distance=infty;
 }
if ((x>0)&&(y<Ydim-1)) {
    mask[7].v=inDE[z][y+1][x-1].v+1;
    mask[7].h=inDE[z][y+1][x-1].h+1;
    mask[7].d=inDE[z][y+1][x-1].d;
    mask[7].distance=sqrt1(mask[7].v*mask[7].v+mask[7].h*mask[7].h+mask[7].d*mask[7].d);
  }
 else {
   mask[7].v=infty; 
   mask[7].h=infty;
   mask[7].d=infty;
   mask[7].distance=infty;
 }
if ((x>0)&&(y<Ydim-1)&&(z<Zdim-1)) {
    mask[8].v=inDE[z+1][y+1][x-1].v+1;
    mask[8].h=inDE[z+1][y+1][x-1].h+1;
    mask[8].d=inDE[z+1][y+1][x-1].d+1;
    mask[8].distance=sqrt1(mask[8].v*mask[8].v+mask[8].h*mask[8].h+mask[8].d*mask[8].d);
  }
 else {
   mask[8].v=infty; 
   mask[8].h=infty;
   mask[8].d=infty;
   mask[8].distance=infty;
 }

  //
  //                     MASK:
  //
  //             (0,1,1)(0,1,0)(0,1,1)
  //             (0,0,1)(0,0,0)
  //
  //
  
  
  if ((z>0)&&(y>0)) {
    mask[9].v=inDE[z-1][y-1][x].v; 
    mask[9].h=inDE[z-1][y-1][x].h+1;
    mask[9].d=inDE[z-1][y-1][x].d+1;
    mask[9].distance=sqrt1(mask[9].v*mask[9].v+mask[9].h*mask[9].h+mask[9].d*mask[9].d);
  }
  else {
    mask[9].v=infty;
    mask[9].h=infty;
    mask[9].d=infty;
    mask[9].distance=infty;
  }
 

  if (y>0){
    mask[10].v=inDE[z][y-1][x].v; 
    mask[10].h=inDE[z][y-1][x].h+1;
    mask[10].d=inDE[z][y-1][x].d;
    mask[10].distance=sqrt1(mask[10].v*mask[10].v+mask[10].h*mask[10].h+mask[10].d*mask[10].d);
  }
  else {
    mask[10].v=infty;
    mask[10].h=infty;
    mask[10].d=infty;
    mask[10].distance=infty;

  }

 
  if ((z<Zdim-1)&&(y>0)) {
    mask[11].v=inDE[z+1][y-1][x].v; 
    mask[11].h=inDE[z+1][y-1][x].h+1; 
    mask[11].d=inDE[z+1][y-1][x].d+1; 
    mask[11].distance=sqrt1(mask[11].v*mask[11].v+mask[11].h*mask[11].h+mask[11].d*mask[11].d);
  }
  else {
    mask[11].v=infty;
    mask[11].h=infty;
    mask[11].d=infty;
    mask[11].distance=infty;
  }
  
  mask[12].v=inDE[z][y][x].v;
  mask[12].h=inDE[z][y][x].h;
  mask[12].d=inDE[z][y][x].d;
  mask[12].distance=sqrt1(mask[12].v*mask[12].v+mask[12].h*mask[12].h+mask[12].d*mask[12].d);
  
  if (z>0) {
    mask[13].v=inDE[z-1][y][x].v;
    mask[13].h=inDE[z-1][y][x].h;
    mask[13].d=inDE[z-1][y][x].d+1;
    mask[13].distance=sqrt1(mask[13].v*mask[13].v+mask[13].h*mask[13].h+mask[13].d*mask[13].d);
  }
 else {
   mask[13].v=infty; 
   mask[13].h=infty;
   mask[13].d=infty;
   mask[13].distance=infty;
 }

  for(looper=0;looper<14;looper++) {
    if (mask[looper].distance<min.distance) min=mask[looper];
  }
 
  return min;
}

const bool Debug= false;

void DEuclidean(Array3dDEucl3D& A)
{
  DEucl3D*** inDE = A.data;
  int Xdim = A.Xdim;
  int Ydim = A.Ydim;
  int Zdim = A.Zdim;

  int y,x,z;
  if (Debug) printf("DEuclidean started ...\n");
  for (x=0;x<Xdim;x++) {
    for(y=0;y<Ydim;y++){ 
      for(z=0;z<Zdim;z++) inDE[z][y][x]=MINforwardDE1(A, z,y,x);
      for(z=Zdim-1;z>-1;z--) inDE[z][y][x]=MINforwardDE2(A, z,y,x);
    }
    for(y=Ydim-1;y>-1;y--) { 
      for(z=Zdim-1;z>-1;z--) inDE[z][y][x]=MINforwardDE3(A, z,y,x);
      for(z=0;z<Zdim;z++) inDE[z][y][x]=MINforwardDE4(A, z,y,x);
    }
  }
  for (x=Xdim-1;x>-1;x--) {
    for(y=Ydim-1;y>-1;y--){ 
      for(z=Zdim-1;z>-1;z--) inDE[z][y][x]=MINbackwardDE1(A, z,y,x);
      for(z=0;z<Zdim;z++) inDE[z][y][x]=MINforwardDE4(A, z,y,x);
    }
    for(y=0;y<Ydim;y++) { 
      for(z=0;z<Zdim;z++) inDE[z][y][x]=MINbackwardDE3(A, z,y,x);
      for(z=Zdim-1;z>-1;z--) inDE[z][y][x]=MINforwardDE2(A, z,y,x);
      
    }
  }
  if (Debug) printf("DEuclidean finished\n");
 
}

// ************************************
//  
//       D-Euclidean END
//
// ***********************************

void printArray3D(const Array3dDEucl3D& A, int x, char type)
{
  DEucl3D*** inDE = A.data;
  int Xdim = A.Xdim;
  int Ydim = A.Ydim;
  int Zdim = A.Zdim;

  int z,y;
  printf("\nslice x= %d\n",x);
  for (y=0;y<Ydim;y++) {
    printf("\n");
    for(z=0;z<Zdim;z++) 
      if (type=='f') {
	if (inDE[z][y][x].distance==infty) 
	  printf("* "); 
	else
	  printf("%f ",inDE[z][y][x].distance);
      } 
  }
  printf("\n");
}

// ****************************
//
//  float functions BEGIN
//
// ****************************
float MINforward3Dfloat(Array3dDEucl3D& A, float d1,float d2,float d3,int z,int y,int x){
  DEucl3D*** inDE = A.data;
  int Xdim = A.Xdim;
  int Ydim = A.Ydim;
  int Zdim = A.Zdim;

  float mask[14];
  float min=infty;
  int looper;

  // MASK
  //
  //   (+d3)(+d2)(+d3)   (+d2)(+d1)(+d2)
  //   (+d2)(+d1)(+d2)   (+d1)( 0 )
  //   (+d3)(+d2)(+d3)


  // apply the mask
  if ((z>0)&&(y>0)) mask[0]=inDE[z-1][y-1][x].distance+d2; else mask[0]=infty;
  if (y>0) mask[1]=inDE[z][y-1][x].distance+d1; else mask[1]=infty;
  if (( z<Zdim-1)&&(y>0)) mask[2]=inDE[z+1][y-1][x].distance+d2; else mask[2]=infty;
  mask[3]=inDE[z][y][x].distance;
  if (z>0) mask[4]=inDE[z-1][y][x].distance+d1; else mask[4]=infty;

  if ((x>0)&&(z>0)&&(y>0)) mask[5]=inDE[z-1][y-1][x-1].distance+d3; else mask[5]=infty;
  if ((x>0)&&(y>0)) mask[6]=inDE[z][y-1][x-1].distance+d2; else mask[6]=infty;
  if ((x>0)&&(z<Zdim-1)&&(y>0)) mask[7]=inDE[z+1][y-1][x-1].distance+d3; else mask[7]=infty;
  if ((x>0)&&(z>0)) mask[8]=inDE[z-1][y][x-1].distance+d2; else mask[8]=infty;
  if (x>0) mask[9]=inDE[z][y][x-1].distance+d1; else mask[9]=infty;
  if ((x>0)&&(z<Zdim-1)) mask[10]=inDE[z+1][y][x-1].distance+d2; else mask[10]=infty;
  if ((x>0)&&(z>0)&&(y<Ydim-1)) mask[11]=inDE[z-1][y+1][x-1].distance+d3; else mask[11]=infty;
  if ((x>0)&&(y<Ydim-1)) mask[12]=inDE[z][y+1][x-1].distance+d2; else mask[12]=infty;
  if ((x>0)&&(z<Zdim-1)&&(y<Ydim-1)) mask[13]=inDE[z+1][y+1][x-1].distance+d3; else mask[13]=infty;
  // find minimum
  for(looper=0;looper<14;looper++) if (mask[looper]<min) min=mask[looper];
 
  
  return min;
}

float MINbackward3Dfloat(Array3dDEucl3D& A, float d1,float d2,float d3,int z,int y,int x){
  DEucl3D*** inDE = A.data;
  int Xdim = A.Xdim;
  int Ydim = A.Ydim;
  int Zdim = A.Zdim;

  float mask[14];
  float min=infty;
  int looper;
  
  //                MASK:
  //
  //                       (+d3)(+d2)(+d3)   
  //          ( 0 )(+d1)   (+d2)(+d1)(+d2)   
  //     (+d2)(+d1)(+d2)   (+d3)(+d2)(+d3)

  // apply the mask
  
  mask[0]=inDE[z][y][x].distance;
  if (z<Zdim-1) mask[1]=inDE[z+1][y][x].distance+d1; else mask[1]=infty;
  if ((y<Ydim-1)&&(z<Zdim-1)) mask[2]=inDE[z+1][y+1][x].distance+d2; else mask[2]=infty;
  if (y<Ydim-1) mask[3]=inDE[z][y+1][x].distance+d1; else mask[3]=infty;
  if ((z>0)&&(y<Ydim-1)) mask[4]=inDE[z-1][y+1][x].distance+d2; else mask[4]=infty;

  if ((x<Xdim-1)&&(z>0)&&(y>0)) mask[5]=inDE[z-1][y-1][x+1].distance+d3; else mask[5]=infty;
  if ((x<Xdim-1)&&(y>0)) mask[6]=inDE[z][y-1][x+1].distance+d2; else mask[6]=infty;
  if ((x<Xdim-1)&&(z<Zdim-1)&&(y>0)) mask[7]=inDE[z+1][y-1][x+1].distance+d3; else mask[7]=infty;
  if ((x<Xdim-1)&&(z>0)) mask[8]=inDE[z-1][y][x+1].distance+d2; else mask[8]=infty;
  if (x<Xdim-1) mask[9]=inDE[z][y][x+1].distance+d1; else mask[9]=infty;
  if ((x<Xdim-1)&&(z<Zdim-1)) mask[10]=inDE[z+1][y][x+1].distance+d2; else mask[10]=infty;
  if ((x<Xdim-1)&&(z>0)&&(y<Ydim-1)) mask[11]=inDE[z-1][y+1][x+1].distance+d3; else mask[11]=infty;
  if ((x<Xdim-1)&&(y<Ydim-1)) mask[12]=inDE[z][y+1][x+1].distance+d2; else mask[12]=infty;
  if ((x<Xdim-1)&&(z<Zdim-1)&&(y<Ydim-1)) mask[13]=inDE[z+1][y+1][x+1].distance+d3; else mask[13]=infty;
  // find minimum
  for(looper=0;looper<14;looper++) if (mask[looper]<min) min=mask[looper];
 
 
  return min;
}

void DistanceTransform3Dfloat(Array3dDEucl3D& A, float d1,float d2,float d3){
  DEucl3D*** inDE = A.data;
  int Xdim = A.Xdim;
  int Ydim = A.Ydim;
  int Zdim = A.Zdim;

  int z,y,x;
  if (Debug) printf("chamfer float started ...\n");
 
 for (x=0;x<Xdim;x++) 
   for (y=0;y<Ydim;y++) 
     for(z=0;z<Zdim;z++) inDE[z][y][x].distance=MINforward3Dfloat(A, d1,d2,d3,z,y,x);

 
 for (x=Xdim-1;x>-1;x--)
   for (y=Ydim-1;y>-1;y--) 
     for(z=Zdim-1;z>-1;z--) inDE[z][y][x].distance=MINbackward3Dfloat(A, d1,d2,d3,z,y,x);
  if (Debug) printf("chamfer float finished\n");
 
}

// ***************************
//
//  float functions END
//
// ***************************

DT3D::DT3D()
{
	A.data = NULL;
}

void DT3D::Build(double* _x, double* _y, double* _z, int num)
{
	xMin=_x[0]; xMax=_x[0]; yMin=_y[0]; yMax=_y[0]; zMin=_z[0]; zMax=_z[0];
	int i;
	for(i = 1; i < num; i ++)
	{
		if(xMin > _x[i]) xMin = _x[i];
		if(xMax < _x[i]) xMax = _x[i];
		if(yMin > _y[i]) yMin = _y[i];
		if(yMax < _y[i]) yMax = _y[i];
		if(zMin > _z[i]) zMin = _z[i];
		if(zMax < _z[i]) zMax = _z[i];
	}

	double xCenter = (xMin+xMax)/2;
	double yCenter = (yMin+yMax)/2;
	double zCenter = (zMin+zMax)/2;
	xMin = xCenter - expandFactor*(xMax-xCenter);
	xMax = xCenter + expandFactor*(xMax-xCenter);
	yMin = yCenter - expandFactor*(yMax-yCenter);
	yMax = yCenter + expandFactor*(yMax-yCenter);
	zMin = zCenter - expandFactor*(zMax-zCenter);
	zMax = zCenter + expandFactor*(zMax-zCenter);

	double max = xMax-xMin > yMax-yMin ? xMax-xMin : yMax-yMin;
	max = max > zMax-zMin ? max : zMax-zMin;

	xMin = xCenter - max/2;
	xMax = xCenter + max/2;
	yMin = yCenter - max/2;
	yMax = yCenter + max/2;
	zMin = zCenter - max/2;
	zMax = zCenter + max/2;

	scale = SIZE/max;

	//printf("DTaccu:%lf\n",sqrt(3.0)/2/scale);

	if(A.data == NULL)
	{
		A.Init(SIZE, SIZE, SIZE);
	}

	int x,y,z;

	DEucl3D*** inDE = A.data;
	int Xdim = A.Xdim;
	int Ydim = A.Ydim;
	int Zdim = A.Zdim;
	
	for(z=0; z<Zdim; z++) {
		//	printf("x=%d y=%d z=%d\n",x,y,z);
		for(y=0; y<Ydim; y++) {
			//	printf("x=%d y=%d z=%d\n",x,y,z);
			for(x=0; x<Xdim; x++) {
				//	printf("x=%d y=%d z=%d\n",x,y,z);
					inDE[z][y][x].distance=infty;
					inDE[z][y][x].h=infty;
					inDE[z][y][x].v=infty;
					inDE[z][y][x].d=infty;
			}
		}
	}
	for(i = 0; i < num; i++)
	{
		x = ROUND((_x[i]-xMin)*scale);
		y = ROUND((_y[i]-yMin)*scale);
		z = ROUND((_z[i]-zMin)*scale);

		if(x<0 || x>=Xdim || y<0 || y>=Ydim || z<0 || z>=Zdim)
			continue;

		inDE[z][y][x].distance=0;
		inDE[z][y][x].h=0;
		inDE[z][y][x].v=0;
		inDE[z][y][x].d=0;
	}

	DEuclidean(A);

	for(z=0; z<Zdim; z++) {
		for(y=0; y<Ydim; y++) {
			for(x=0; x<Xdim; x++) {
				//inDE[z][y][x].distance = (inDE[z][y][x].distance  - sqrt(3.0))/scale;
				inDE[z][y][x].distance = (inDE[z][y][x].distance)/scale;
				if(inDE[z][y][x].distance < 0)
					inDE[z][y][x].distance = 0;
			}
		}
	}
}

float DT3D::Distance(double _x, double _y, double _z)
{
	int x, y, z;
	x = ROUND((_x-xMin)*scale);
	y = ROUND((_y-yMin)*scale);
	z = ROUND((_z-zMin)*scale);

	if(x > -1 && x < SIZE && y > -1 && y < SIZE && z > -1 && z < SIZE)
		return A.data[z][y][x].distance;

	float a = 0, b = 0, c = 0;
	if(x < 0)
	{
		a = x;
		x = 0;
	}
	else if(x >= SIZE)
	{
		a = x-SIZE+1;
		x = SIZE-1;
	}
		
	if(y < 0)
	{
		b = y;
		y = 0;
	}
	else if(y >= SIZE)
	{
		b = y-SIZE+1;
		y = SIZE-1;
	}
		
	if(z < 0)
	{
		c = z;
		z = 0;
	}
	else if(z >= SIZE)
	{
		c = z-SIZE+1;
		z = SIZE-1;
	}
		
	return sqrt(a*a+b*b+c*c)/scale + A.data[z][y][x].distance;
}
