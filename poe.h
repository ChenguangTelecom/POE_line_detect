#ifndef POE_H_
#define POE_H_
#include<stdio.h>

struct image_double{
  double * data;
  int X;
  int Y;
};

typedef struct image_double * image_double_ptr;

struct rect{
  double x1,y1,x2,y2;
  double cx;
  double cy;
  double theta;
};

typedef struct rect  rect_item;
typedef struct rect * rect_item_ptr;

struct rect_node{
  rect_item item;
  struct rect_node * next;
};

typedef struct  rect_node * rect_node_ptr; 
struct rect_list{
  rect_node_ptr start;
  rect_node_ptr end;
  int count;
};

  typedef struct rect_list * rect_list_ptr;

struct point{
  int x;
  int y;
};

typedef struct point point_item;
typedef struct point * point_item_ptr;

struct point_node
{
  point_item item;
  struct point_node * next;
};

typedef struct point_node * point_node_ptr;

struct point_list
{
  point_node_ptr start;
  point_node_ptr end;
  int count;
};

typedef struct point_list * point_list_ptr;


void error(char * msg);
image_double_ptr new_image_double(int X, int Y);
void image_double_initialize(image_double_ptr  img, double value_ini);
void free_image_double(image_double_ptr  img);
image_double_ptr * kernel_generation(int W, int P);
void free_kernel(image_double_ptr *kernel, int P);
static int max_response(double * p_response, int P);
  void poe(image_double_ptr o_matrix,image_double_ptr BW, image_double_ptr * kernel, int P,image_double_ptr flag_matrix);

void InitializeList(rect_list_ptr plist);
int ListIsEmpty(rect_list_ptr plist);
void AddNode(rect_list_ptr plist, rect_item_ptr pitem);
void EmptyTheList(rect_list_ptr plist);

void InitializePlist(point_list_ptr plist);
int PlistIsEmpty(const point_list_ptr plist);
void AddPointNode(point_list_ptr plist, point_item_ptr pitem);
void EmptyPlist(point_list_ptr plist);
rect_list_ptr line_segment_detection(image_double_ptr o_matrix, image_double_ptr flag_matrix, double tau);
#endif
