#include "poe.h"
#include<stdio.h>
#include<stdlib.h>
#include"mex.h"
void error(char * msg)
{
  mexErrMsgTxt(msg);
  /*fprintf(stderr, "%s\n", msg);*/
  /*exit(1);*/
}

image_double_ptr new_image_double(int X, int Y)
{
  image_double_ptr img;
  img=(image_double_ptr ) malloc(sizeof(struct image_double));
  img->X=X;
  img->Y=Y;
  img->data=(double *) malloc(X*Y*sizeof(double));
  if (img==NULL || img->data==NULL)
    error("not enough memory");
  return img;
}

void image_double_initialize(image_double_ptr  img, double value_ini)
{
  if (img==NULL || img->data==NULL || img->X<0 || img->Y<0)
    error("image_double_initialize: invalid input");
  int i,j;
  for (i=0;i<img->X;i++)
    for (j=0;j<img->Y;j++)
      img->data[i+j*img->X]=value_ini;
}

void free_image_double(image_double_ptr  img)
{
  if (img==NULL || img->data==NULL)
    error("free_image_double: invalid input");
  free((void *)img->data);
  free((void *)img);
}

void InitializeList(rect_list_ptr plist)
{
  if (plist==NULL)
    {
      error("InitializeList: invalid input");
    }
  plist->start=NULL;
  plist->end=NULL;
  plist->count=0;
}

int ListIsEmpty(rect_list_ptr plist)
{
  if (plist->count==0)
    return 1;
  else
    return 0;
}

static void CopyToNode(rect_item_ptr pitem, rect_node_ptr pnode)
{
  if (pnode==NULL||pitem==NULL)
    error("CopyToNode: invalid input");
  pnode->item=*pitem;
  pnode->next=NULL;
}


void AddNode(rect_list_ptr plist, rect_item_ptr pitem)
{
  rect_node_ptr pnode;
  pnode=(rect_node_ptr) malloc(sizeof(struct rect_node));
  CopyToNode(pitem, pnode);
  if (plist==NULL)
    {
      error("AddNode: invalid input");
    }
  if (ListIsEmpty(plist))
    {
      plist->start=pnode;
    }
  else
    plist->end->next=pnode;
  plist->end=pnode;
  plist->count++;
}

void EmptyTheList(rect_list_ptr plist)
{
  rect_node_ptr pnode;
  while(plist->count>0)
    {
      pnode=plist->start;
      plist->start=plist->start->next;
      free(pnode);
      plist->count--;
      /*printf("%d\n", plist->count);*/
    }
  plist->end=NULL;
}


void InitializePlist(point_list_ptr plist)
{
  if (plist==NULL)
    error("not enough memory");
  plist->start=NULL;
  plist->end=NULL;
  plist->count=0;
}

int PlistIsEmpty(const point_list_ptr plist)
{
  if (plist==NULL)
    error("PlistIsEmpty: invalid input");

  if (plist->count==0)
    return 1;
  else
    return 0;
}

static void CopyPointToNode(point_node_ptr pnode, point_item_ptr pitem)
{
  if (pnode==NULL || pitem==NULL)
    error("CopyPointToNode: invalid input");
  pnode->item=*pitem;
  pnode->next=NULL;
}
  
void AddPointNode(point_list_ptr plist, point_item_ptr pitem)
{
  if (plist==NULL || pitem==NULL)
    error("AddPointNode: invalid input");
  point_node_ptr pnode;
  pnode=(point_node_ptr) malloc(sizeof(struct point_node));
  if (pnode==NULL)
    error("not enough memory");
  CopyPointToNode(pnode,pitem);
  if (PlistIsEmpty(plist))
    {
      plist->start=pnode;
    }
  else
    plist->end->next=pnode;
  plist->end=pnode;
  plist->count++;
}

void EmptyPlist(point_list_ptr plist)
{
  point_node_ptr pnode;
  while(plist->count>0)
    {
      pnode=plist->start;
      plist->start=plist->start->next;
      free(pnode);
      plist->count--;
    }
  plist->end=NULL;
}
