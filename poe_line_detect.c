#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include "mex.h"
#include "poe.h"
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#define EDGEMAP prhs[0]
#define W_VALUE prhs[1]
#define P_VALUE prhs[2]
#define TAU prhs[3]
#define LINES plhs[0]

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray * prhs[])
{
  double * edgemap, *out;
  double value_p, value_w, tau;
  if (nrhs!=4)
    mexErrMsgTxt("must have four inputs: edgemap, W, P. tau\n");
  if (nlhs!=1)
    mexErrMsgTxt("must have one output\n");
  if (mxIsComplex(EDGEMAP) || mxGetNumberOfDimensions(EDGEMAP)!=2 || !mxIsDouble(EDGEMAP) )
    mexErrMsgTxt("the edge map data must be a 2D double matrix.\n");
  if (!mxIsDouble(P_VALUE) || !mxIsDouble(W_VALUE) || !mxIsDouble(TAU) || mxGetNumberOfElements(P_VALUE)!=1 || mxGetNumberOfElements(W_VALUE)!=1 || mxGetNumberOfElements(TAU)!=1)
    mexErrMsgTxt("P_VALUE, W_VALUE and TAU must be double scalars.\n");
  else
    {
      value_w=mxGetScalar(W_VALUE);
      value_p=mxGetScalar(P_VALUE);
      tau=mxGetScalar(TAU);
    }
  int X, Y,i,j;
  X=mxGetM(EDGEMAP);
  Y=mxGetN(EDGEMAP);
  edgemap=mxGetPr(EDGEMAP);
  image_double_ptr BW;
  image_double_ptr flag_matrix;
  BW=new_image_double((int) (X+2*value_w), (int) (Y+2*value_w));
  flag_matrix=new_image_double(X,Y);
  image_double_initialize(BW,0.0);
  for (i=0;i<X;i++)
    for (j=0;j<Y;j++)
      {
	BW->data[i+(int) value_w+(j+(int) value_w)*(X+2* (int) value_w)]=edgemap[i+j*X];
	flag_matrix->data[i+j*X]=edgemap[i+j*X]>0.5?1.0:0.0;
      }
  
  image_double_ptr *kernel;
  kernel=kernel_generation( (int) value_w, (int) value_p);
  image_double_ptr o_matrix;
  o_matrix=new_image_double(X,Y);
  image_double_initialize(o_matrix,-100.0);
  poe(o_matrix,BW, kernel, (int) value_p, flag_matrix);
  rect_list_ptr lines_list;
  lines_list=line_segment_detection(o_matrix, flag_matrix, tau);
  LINES=mxCreateDoubleMatrix(lines_list->count, 4, mxREAL);
  out=mxGetPr(LINES);
  rect_node_ptr pnode;
  pnode=lines_list->start;
  int num_count;
  num_count=lines_list->count;
  i=0;
  while(pnode!=NULL)
    {
      out[i+0*num_count]=(pnode->item).x1+1.0;
      out[i+1*num_count]=(pnode->item).y1+1.0;
      out[i+2*num_count]=(pnode->item).x2+1.0;
      out[i+3*num_count]=(pnode->item).y2+1.0;
      i++;
      pnode=pnode->next;
    }
  EmptyTheList(lines_list);
  free(lines_list);
  
  
  
  free_image_double(BW);
  free_image_double(flag_matrix);
  free_image_double(o_matrix);
  free_kernel(kernel, (int) value_p);  
}


image_double_ptr * kernel_generation(int W, int P)
{
  image_double_ptr * kernel;
  kernel=(image_double_ptr *) malloc(P*sizeof(image_double_ptr));
  
  if (kernel==NULL)
    error("not enough memory");
  int i,j,k, curr_x, curr_y;
  for (k=0;k<P;k++)
    kernel[k]=new_image_double(2*W+1,2*W+1);
  double theta;
  double dist;
  for (k=0;k<P;k++)
    {
      theta=(double) k*M_PI/(double) P;
      image_double_initialize(kernel[k],0.0);
      for (i=0;i<2*W+1;i++)
	for (j=0;j<2*W+1;j++)
	  {
	    curr_x=i-W;
	    curr_y=j-W;
	    dist=-(double) curr_x*sin(theta)+curr_y*cos(theta);
	    if (fabs(dist)<0.5 && pow((double) curr_x,2)+pow((double) curr_y,2)<= pow((double) W,2)+0.000001)
	      kernel[k]->data[i+j*(2*W+1)]=1.0;
	  }
    }
  return kernel;
}

void free_kernel(image_double_ptr *kernel, int P)
{
  if (kernel==NULL)
    error("free_kernel: invalid input kernel");
  int k;
  for (k=0;k<P;k++)
    {
      free_image_double(kernel[k]);
    }
  free((void *) kernel);
}

static int max_response(double * p_response, int P)
{
  if (p_response==NULL)
    error("max_response: invalid p_response");
  int i, max_index=0;
  double temp=p_response[0];
  for (i=1;i<P;i++)
    {
      if(temp<p_response[i])
	{
	  max_index=i;
	  temp=p_response[i];
	}
    }
  return max_index;
}

void poe(image_double_ptr o_matrix,image_double_ptr BW, image_double_ptr * kernel, int P,image_double_ptr flag_matrix)
{
  int i,j,k,m,n,curr_i,curr_j,test_i, test_j;
  int W, X,Y,X_W,Y_W;
  int max_index;
  X=flag_matrix->X;
  Y=flag_matrix->Y;
  X_W=BW->X;
  W=(X_W-X)/2;
  Y_W=BW->Y;
  double *theta;
  double * p_response;
  p_response=(double *) malloc(P*sizeof(double));
  theta=(double *) malloc(P*sizeof(double));
  for (i=0;i<P;i++)
    {
      p_response[i]=0.0;
      theta[i]=(double) i*M_PI/(double) P;
    }
  for (i=0;i<X;i++)
    for (j=0;j<Y;j++)
      {
	if (flag_matrix->data[i+j*X]>0.5)
	  {
	    curr_i=i+W;
	    curr_j=j+W;
	    for (k=0;k<P;k++)
	      p_response[k]=0.0;
	    for (k=0;k<P;k++)
	      {
		for (m=-W;m<W+1;m++)
		  for(n=-W;n<W+1;n++)
		    {
		      test_i=curr_i+m;
		      test_j=curr_j+n;
		      p_response[k]+=BW->data[test_i+test_j*X_W]*kernel[k]->data[m+W+(n+W)*(2*W+1)];     
		    }
	      }
	    max_index=max_response(p_response,P);
	    o_matrix->data[i+j*X]=theta[max_index];
	  }
      }
  free(p_response);
  free(theta);
}
static int max(int x, int y)
{
  return x<y?y:x;
}
static int min(int x, int y)
{
  return x<y?x:y;
}

static int isaligned(double reg_angle, double theta, double tau)
{
  theta-=reg_angle;
  if (theta<0.0)
    theta=-theta;
  if(theta>M_PI/2.0)
    theta=M_PI-theta;
  if (theta<=tau+0.01)
    return 1;
  else
    return 0;
}

static int region_grow(point_list_ptr preg, double reg_angle, image_double_ptr o_matrix, image_double_ptr flag_matrix, double tau)
{
  int i, j,X,Y, curr_x,curr_y,test_x,test_y;
  X=flag_matrix->X;
  Y=flag_matrix->Y;
  double theta;
  point_item pix_coor;
  point_node_ptr pnode;
  pnode=preg->start;
  while(pnode!=NULL)
    {
      curr_x=(pnode->item).x;
      curr_y=(pnode->item).y;
      for (i=-1;i<2;i++)
	for (j=-1;j<2;j++)
	  {
	    test_x=min(max(curr_x+i,0),X-1);
	    test_y=min(max(curr_y+j,0),Y-1);
	    theta=o_matrix->data[test_x+test_y*X];
	    if (flag_matrix->data[test_x+test_y*X]>0.5 && isaligned(reg_angle, theta, tau) )
	      {
		pix_coor.x=test_x;
		pix_coor.y=test_y;
		flag_matrix->data[test_x+test_y*X]=0.0;
		AddPointNode(preg,&pix_coor);
	      }
	  }
      pnode=pnode->next;
    }
  double reg_threshold;
  reg_threshold=-log(sqrt(pow((double) X* (double) Y,5)))/log(3.0/16.0);
  if ((double) preg->count>=reg_threshold)
    return 1;
  else
    {
      pnode=preg->start;
      while(pnode!=NULL)
	{
	  curr_x=(pnode->item).x;
	  curr_y=(pnode->item).y;
	  flag_matrix->data[curr_x+curr_y*X]=1.0;
	  pnode=pnode->next;
	}
      EmptyPlist(preg);
      return 0;
    }
}

static double get_theta(point_list_ptr preg, double cx, double cy)
{
  if (preg==NULL)
    error("get_theta: invalid input");
  double theta;
  double mxx=0.0;
  double mxy=0.0;
  double myy=0.0;
  int curr_x, curr_y;
  point_node_ptr pnode;
  pnode=preg->start;
  while(pnode!=NULL)
    {
      curr_x=(pnode->item).x;
      curr_y=(pnode->item).y;
      mxx+=pow((double) curr_y-cy,2);
      myy+=pow((double) curr_x-cx,2);
      mxy-=((double) curr_x-cx)*((double) curr_y-cy);
      pnode=pnode->next;
    }
  mxx/=(double) preg->count;
  mxy/=(double) preg->count;
  myy/=(double) preg->count;
  double lambda;
  lambda=0.5*(mxx-myy-sqrt(pow(mxx-myy,2)+4.0*pow(mxy,2)));
  theta=fabs(mxx)>fabs(myy)?atan2(lambda-mxx,mxy):atan2(mxy,lambda-myy);
  return theta;
}

static void region_to_rectangle(rect_list_ptr lines_list, point_list_ptr preg)
{
  if (lines_list==NULL || preg==NULL)
    error("region_to_rectangle: invalid input");
  double cx, cy;
  double theta;
  cx=cy=0.0;
  point_node_ptr pnode;
  pnode=preg->start;
  while(pnode!=NULL)
    {
      cx+=(double) (pnode->item).x;
      cy+=(double) (pnode->item).y;
      pnode=pnode->next;
    }
  cx/=(double) preg->count;
  cy/=(double) preg->count;
  theta=get_theta(preg, cx, cy);
  double dx;
  dx=cos(theta);
  double dy;
  dy=sin(theta);
  double l_max, l_min, w_max, w_min,l,w;
  l_max=l_min=w_max=w_min=0.0;
  pnode=preg->start;
  int curr_x, curr_y;
  while(pnode!=NULL)
    {
      curr_x=(pnode->item).x;
      curr_y=(pnode->item).y;
      l=((double) curr_x-cx)*dx+((double) curr_y-cy)*dy;
      w=-((double) curr_x-cx)*dy+((double) curr_y-cy)*dx;
      if (l>l_max)
	l_max=l;
      if (l<l_min)
	l_min=l;
      if (w>w_max)
	w_max=w;
      if (w<w_min)
	w_min=w;
      pnode=pnode->next;
    }
  double x1,y1,x2,y2;
  x1=cx+l_min*dx;
  y1=cy+l_min*dy;
  x2=cx+l_max*dx;
  y2=cy+l_max*dy;
  rect_item line;
  line.x1=x1;
  line.y1=y1;
  line.x2=x2;
  line.y2=y2;
  line.cx=cx;
  line.cy=cy;
  line.theta=theta;
  AddNode(lines_list,&line);
}

rect_list_ptr line_segment_detection(image_double_ptr o_matrix, image_double_ptr flag_matrix, double tau)
{
  rect_list_ptr lines_list;
  lines_list=(rect_list_ptr) malloc(sizeof(struct rect_list));
  InitializeList(lines_list);
  int i,j, X, Y;
  X=flag_matrix->X;
  Y=flag_matrix->Y;
 
  point_list_ptr preg;
  preg=(point_list_ptr) malloc(sizeof(struct point_list));
  InitializePlist(preg);
  double reg_theta;
  point_item pix_coor;
  for (i=0;i<X;i++)
    for (j=0;j<Y;j++)
      {
	if(flag_matrix->data[i+j*X]>0.5)
	  {
	    reg_theta=o_matrix->data[i+j*X];
	    pix_coor.x=i;
	    pix_coor.y=j;
	    flag_matrix->data[i+j*X]=0.0;
	    InitializePlist(preg);
	    AddPointNode(preg,&pix_coor);
	    if (region_grow(preg, reg_theta, o_matrix, flag_matrix, tau))
	      {
		region_to_rectangle(lines_list, preg);
		EmptyPlist(preg);
	      }
	    
	  }
      }
  free(preg);
  return lines_list;
}  
