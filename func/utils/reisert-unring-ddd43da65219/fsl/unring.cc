//     unring.cc  unrings 4D volumes
//     Marco Reisert, Elias Kellner, Medical Physics, University Medical Center Freiburg

/*  
    LICENCE
        
    The Software remains the property of the University Medical Center 
    Freiburg ("the University").
    
    The Software is distributed "AS IS" under this Licence solely for
    non-commercial use in the hope that it will be useful, but in order
    that the University as a charitable foundation protects its assets for
    the benefit of its educational and research purposes, the University
    makes clear that no condition is made or to be implied, nor is any
    warranty given or to be implied, as to the accuracy of the Software,
    or that it will be suitable for any particular purpose or for use
    under any specific conditions. Furthermore, the University disclaims
    all responsibility for the use which is made of the Software. It
    further disclaims any liability for the outcomes arising from using
    the Software.
    
    The Licensee agrees to indemnify the University and hold the
    University harmless from and against any and all claims, damages and
    liabilities asserted by third parties (including claims for
    negligence) which arise directly or indirectly from the use of the
    Software or the sale of any products based on the Software.
    
    No part of the Software may be reproduced, modified, transmitted or
    transferred in any form or by any means, electronic or mechanical,
    without the express permission of the University. The permission of
    the University is not required if the said reproduction, modification,
    transmission or transference is done without financial return, the
    conditions of this Licence are imposed upon the receiver of the
    product, and all original and amended source code is included in any
    transmitted product. You may be held legally responsible for any
    copyright infringement that is caused or encouraged by your failure to
    abide by these terms and conditions.
    
    You are not permitted under this Licence to use this Software
    commercially. Use for which any financial return is received shall be
    defined as commercial use, and includes (1) integration of all or part
    of the source code or the Software into a product for sale or license
    by or on behalf of Licensee to third parties or (2) use of the
    Software or any derivative of it for research with the final aim of
    developing software products for sale or license to a third party or
    (3) use of the Software or any derivative of it for research with the
    final aim of developing non-software products for sale or license to a
    third party, or (4) use of the Software to provide any service to an
    external organisation for which payment is received. */

#include "newimage/newimageall.h"
#include "newimage/fmribmain.h"
#include <fftw3.h>
#include <stdio.h>
#define PI  3.1416

using namespace NEWIMAGE;


void unring_1D(fftw_complex *data,int n, int numlines,int nsh,int minW, int maxW)
{
    
   
    fftw_complex *in, *out;
    fftw_plan p,pinv;
    
    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);    
    p = fftw_plan_dft_1d(n, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    pinv = fftw_plan_dft_1d(n, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
    
    fftw_complex *sh = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n *(2*nsh+1));
    fftw_complex *sh2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n *(2*nsh+1));
    
    
    double nfac = 1/double(n);
    
    int *shifts = (int*) malloc(sizeof(int)*(2*nsh+1));
    shifts[0] = 0;
    for (int j = 0; j < nsh;j++)
    {
        shifts[j+1] = j+1;
        shifts[1+nsh+j] = -(j+1);
    }
    
    double *TV1arr = new double[2*nsh+1]; 
    double *TV2arr = new double[2*nsh+1]; 
    
    for (int k = 0; k < numlines; k++)
    {
     
        
        fftw_execute_dft(p,&(data[n*k]),sh);
        
        int maxn = (n%2 == 1)? (n-1)/2 : n/2 -1;
        
        for (int j = 1; j < 2*nsh+1; j++)
        {
            double phi = PI/double(n) * double(shifts[j])/double(nsh);
            fftw_complex u = {cos(phi),sin(phi)};
            fftw_complex e = {1,0};

            sh[j*n ][0] = sh[0][0];
            sh[j*n ][1] = sh[0][1];

            if (n%2 == 0)
            {
                sh[j*n + n/2][0] = 0;
                sh[j*n + n/2][1] = 0;
            }
            
            for (int l = 0; l < maxn; l++)
            {            
                
                double tmp = e[0];
                e[0] = u[0]*e[0] - u[1]*e[1];
                e[1] = tmp*u[1] + u[0]*e[1];
                
                int L ;
                L = l+1;
                sh[j*n +L][0] = (e[0]*sh[L][0] - e[1]*sh[L][1]);
                sh[j*n +L][1] = (e[0]*sh[L][1] + e[1]*sh[L][0]);
                L = n-1-l;
                sh[j*n +L][0] = (e[0]*sh[L][0] + e[1]*sh[L][1]);
                sh[j*n +L][1] = (e[0]*sh[L][1] - e[1]*sh[L][0]);                
                                
            }
        }                
        
                
        for (int j = 0; j < 2*nsh+1; j++)
        {
            fftw_execute_dft(pinv,&(sh[j*n]),&sh2[j*n]);
        }
  
        for (int j=0;j < 2*nsh+1;j++)
        {
            TV1arr[j] = 0;
            TV2arr[j] = 0;           
            const int l = 0;
            for (int t = minW; t <= maxW;t++)                
            {
                TV1arr[j] += fabs(sh2[j*n + (l-t+n)%n ][0] - sh2[j*n + (l-(t+1)+n)%n ][0]);
                TV1arr[j] += fabs(sh2[j*n + (l-t+n)%n ][1] - sh2[j*n + (l-(t+1)+n)%n ][1]);
                TV2arr[j] += fabs(sh2[j*n + (l+t+n)%n ][0] - sh2[j*n + (l+(t+1)+n)%n ][0]);
                TV2arr[j] += fabs(sh2[j*n + (l+t+n)%n ][1] - sh2[j*n + (l+(t+1)+n)%n ][1]);
            }
        }
                  

        
        
        for(int l=0; l < n; l++)
        {
            double minTV = 999999999999;
            int minidx= 0;
            for (int j=0;j < 2*nsh+1;j++)
            {                                                                    
                
//                 double TV1 = 0;
//                 double TV2 = 0;
//                 for (int t = minW; t <= maxW;t++)                
//                 {
//                     TV1 += fabs(sh2[j*n + (l-t)%n ][0] - sh2[j*n + (l-(t+1))%n ][0]);
//                     TV1 += fabs(sh2[j*n + (l-t)%n ][1] - sh2[j*n + (l-(t+1))%n ][1]);
//                     TV2 += fabs(sh2[j*n + (l+t)%n ][0] - sh2[j*n + (l+(t+1))%n ][0]);
//                     TV2 += fabs(sh2[j*n + (l+t)%n ][1] - sh2[j*n + (l+(t+1))%n ][1]);
// 
//                 }
//                   
//                 
//                 if (TV1 < minTV)
//                 {
//                     minTV = TV1;
//                     minidx = j;
//                 }
//                 if (TV2 < minTV)
//                 {
//                     minTV = TV2;
//                     minidx = j;
//                 }
                
                
                if (TV1arr[j] < minTV)
                {
                    minTV = TV1arr[j];
                    minidx = j;
                }
                if (TV2arr[j] < minTV)
                {
                    minTV = TV2arr[j];
                    minidx = j;
                }
                
                TV1arr[j] += fabs(sh2[j*n + (l-minW+1+n)%n ][0] - sh2[j*n + (l-(minW)+n)%n ][0]);
                TV1arr[j] -= fabs(sh2[j*n + (l-maxW+n)%n ][0] - sh2[j*n + (l-(maxW+1)+n)%n ][0]);
                TV2arr[j] += fabs(sh2[j*n + (l+maxW+1+n)%n ][0] - sh2[j*n + (l+(maxW+2)+n)%n ][0]);
                TV2arr[j] -= fabs(sh2[j*n + (l+minW+n)%n ][0] - sh2[j*n + (l+(minW+1)+n)%n ][0]);
                
                TV1arr[j] += fabs(sh2[j*n + (l-minW+1+n)%n ][1] - sh2[j*n + (l-(minW)+n)%n ][1]);
                TV1arr[j] -= fabs(sh2[j*n + (l-maxW+n)%n ][1] - sh2[j*n + (l-(maxW+1)+n)%n ][1]);
                TV2arr[j] += fabs(sh2[j*n + (l+maxW+1+n)%n ][1] - sh2[j*n + (l+(maxW+2)+n)%n ][1]);
                TV2arr[j] -= fabs(sh2[j*n + (l+minW+n)%n ][1] - sh2[j*n + (l+(minW+1)+n)%n ][1]);
            
            }
             
           
            double a0r = sh2[minidx*n + (l-1+n)%n ][0];
            double a1r = sh2[minidx*n + l][0];
            double a2r = sh2[minidx*n + (l+1+n)%n ][0];
            double a0i = sh2[minidx*n + (l-1+n)%n ][1];
            double a1i = sh2[minidx*n + l][1];
            double a2i = sh2[minidx*n + (l+1+n)%n ][1];
            double s = double(shifts[minidx])/nsh/2;
            
            //data[k*n + l][0] =  (a1r - 0.5*(a2r-a0r)*s + (0.5*(a2r+a0r) - a1r)*s*s)*nfac;
            //data[k*n + l][1] =  (a1i - 0.5*(a2i-a0i)*s + (0.5*(a2i+a0i) - a1i)*s*s)*nfac;

            
            if (s>0)
            {
                data[k*n + l][0] =  (a1r*(1-s) + a0r*s)*nfac;
                data[k*n + l][1] =  (a1i*(1-s) + a0i*s)*nfac;
            }
            else
            {
                s = -s;
                data[k*n + l][0] =  (a1r*(1-s) + a2r*s)*nfac;
                data[k*n + l][1] =  (a1i*(1-s) + a2i*s)*nfac;
            }
            
        }
        
        
       
        
    }
    
    
     delete TV1arr;
     delete TV2arr;
     free(shifts);
     fftw_destroy_plan(p);
     fftw_destroy_plan(pinv);
     fftw_free(in); 
     fftw_free(out);
     fftw_free(sh);
     fftw_free(sh2);
    
    
    
    
}


void unring_2d(fftw_complex *data1,fftw_complex *tmp2, const int *dim_sz, int nsh, int minW, int maxW)
{

    
        double eps = 0;
        fftw_complex *tmp1 =  (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * dim_sz[0]*dim_sz[1]);        
        fftw_complex *data2 =  (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * dim_sz[0]*dim_sz[1]);    
        
        fftw_plan p,pinv,p_tr,pinv_tr;
        p = fftw_plan_dft_2d(dim_sz[1],dim_sz[0], data1, tmp1, FFTW_FORWARD, FFTW_ESTIMATE);
        pinv = fftw_plan_dft_2d(dim_sz[1],dim_sz[0], data1, tmp1, FFTW_BACKWARD, FFTW_ESTIMATE);        
        p_tr = fftw_plan_dft_2d(dim_sz[0],dim_sz[1], data2, tmp2, FFTW_FORWARD, FFTW_ESTIMATE);
        pinv_tr = fftw_plan_dft_2d(dim_sz[0],dim_sz[1], data2, tmp2, FFTW_BACKWARD, FFTW_ESTIMATE);
        double nfac = 1/double(dim_sz[0]*dim_sz[1]);
        
        for (int k = 0 ; k < dim_sz[1];k++)
           for (int j = 0 ; j < dim_sz[0];j++)
           {
                data2[j*dim_sz[1]+k][0] = data1[k*dim_sz[0]+j][0];
                data2[j*dim_sz[1]+k][1] = data1[k*dim_sz[0]+j][1];
           }
        
        fftw_execute_dft(p,data1,tmp1);
        fftw_execute_dft(p_tr,data2,tmp2);
        
        for (int k = 0 ; k < dim_sz[1];k++)
        {
            double ck = (1+cos(2*PI*(double(k)/dim_sz[1])))*0.5 +eps;
            for (int j = 0 ; j < dim_sz[0];j++)
            {                
                double cj = (1+cos(2*PI*(double(j)/dim_sz[0])))*0.5 +eps;
                tmp1[k*dim_sz[0]+j][0] = nfac*(tmp1[k*dim_sz[0]+j][0] * ck) / (ck+cj);        
                tmp1[k*dim_sz[0]+j][1] = nfac*(tmp1[k*dim_sz[0]+j][1] * ck) / (ck+cj);        
                tmp2[j*dim_sz[1]+k][0] = nfac*(tmp2[j*dim_sz[1]+k][0] * cj) / (ck+cj);        
                tmp2[j*dim_sz[1]+k][1] = nfac*(tmp2[j*dim_sz[1]+k][1] * cj) / (ck+cj);        
            }
        }
        
        fftw_execute_dft(pinv,tmp1,data1);
        fftw_execute_dft(pinv_tr,tmp2,data2);
        
        unring_1D(data1,dim_sz[0],dim_sz[1],nsh,minW,maxW);
        unring_1D(data2,dim_sz[1],dim_sz[0],nsh,minW,maxW);
         
  
        fftw_execute_dft(p,data1,tmp1);
        fftw_execute_dft(p_tr,data2,tmp2);

        
        for (int k = 0 ; k < dim_sz[1];k++)
        {
//            double ck = (1+cos(2*PI*(double(k)/dim_sz[1])))*0.5 +eps;
            for (int j = 0 ; j < dim_sz[0];j++)
            {                
  //              double cj = (1+cos(2*PI*(double(j)/dim_sz[0])))*0.5 +eps;
                tmp1[k*dim_sz[0]+j][0] = nfac*(tmp1[k*dim_sz[0]+j][0]  + tmp2[j*dim_sz[1]+k][0] ) ;        
                tmp1[k*dim_sz[0]+j][1] = nfac*(tmp1[k*dim_sz[0]+j][1]  + tmp2[j*dim_sz[1]+k][1] ) ;                                        
     //           tmp1[k*dim_sz[0]+j][0] = nfac*(tmp1[k*dim_sz[0]+j][0]  + tmp2[j*dim_sz[1]+k][0] ) /(ck+cj);        
     //           tmp1[k*dim_sz[0]+j][1] = nfac*(tmp1[k*dim_sz[0]+j][1]  + tmp2[j*dim_sz[1]+k][1] ) /(ck+cj);                                        
//                 tmp1[k*dim_sz[0]+j][0] = nfac*(tmp1[k*dim_sz[0]+j][0]*ck  + tmp2[j*dim_sz[1]+k][0]*cj ) /(ck+cj);        
//                 tmp1[k*dim_sz[0]+j][1] = nfac*(tmp1[k*dim_sz[0]+j][1]*ck  + tmp2[j*dim_sz[1]+k][1]*cj ) /(ck+cj);                                        
            }
        }
        
        fftw_execute_dft(pinv,tmp1,tmp2);
                
        fftw_free(data2);
        fftw_free(tmp1);
}



void print_usage(const string& progname) {
  cout << endl;
  cout << "unring - tool for removal of the Gibbs ringing artefact" << endl;
  cout << "Usage: unring <input> <output> [options]" << endl;
  cout << " Options: -d     slice direction, in the perpendicular plane unringing is performed. Could be either 1,2 or 3  (default 3)" << endl;
  cout << "          -nsh   discretization of subpixel spaceing (default 20)" << endl;
  cout << "          -minW  left border of window used for TV computation (default 1)" << endl;
  cout << "          -maxW  right border of window used for TV computation (default 3)" << endl;
}


template <class T>
int fmrib_main(int argc, char *argv[])
{

  string alongdim = "3";
  int nsh = 20;
  int minW = 1;
  int maxW = 3;
  
  for (int k = 3; k < argc;k+=2)
  {
      string pname = string(argv[k]);
      
      if (k+1 >= argc)
      {
	fprintf(stderr,"missing parameter value\n");
	return -1;
      }
            
      string pval = string(argv[k+1]);                  
      if (pname.compare(string("-d"))==0)
	  alongdim = pval;	
      else if (pname.compare(string("-n"))==0)
	  nsh = atoi(pval.c_str());
      else if (pname.compare(string("-minW"))==0)
	  minW = atoi(pval.c_str());
      else if (pname.compare(string("-maxW"))==0)
	  maxW = atoi(pval.c_str());
      else 
      {
	  fprintf(stderr,"unknown parameter %s \n",pname.c_str());
      }
    
  }
  
  
  if (minW<0 || minW > 10 || minW >= maxW || maxW<0 || maxW > 128)
  {
     fprintf(stderr,"sth wrong with minW/maxW\n");
     return -1 ;
  }
  
  if (nsh< 16 || nsh > 128)
  {
     fprintf(stderr,"sth wrong with nsh\n");
     return -1 ;
  }
  
  
  


  volume4D<T> input_vol,output_vol;
  string input_name=string(argv[1]);
  string output_name=string(argv[2]);
  read_volume4D(input_vol,input_name);
  
  
  output_vol=input_vol;
 
  int dim_sz[4];
  dim_sz[0] = input_vol.xsize();
  dim_sz[1] = input_vol.ysize();
  dim_sz[2] = 1;
  dim_sz[3] = 1 ; 

  



  for (int t=0; t<input_vol.tsize(); t++) {
    fprintf(stderr,"."); fflush(stderr);
    if (alongdim.compare(string("3"))==0)
    {
	dim_sz[0] = input_vol.xsize();
	dim_sz[1] = input_vol.ysize();
	fftw_complex *data_complex =  (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * dim_sz[0]*dim_sz[1]);
	fftw_complex *res_complex  =  (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * dim_sz[0]*dim_sz[1]);

	for (int z=0; z<input_vol.zsize(); z++)       
	{
	  for (int x = 0 ; x < input_vol.xsize();x++)
	      for (int y = 0 ; y < input_vol.ysize();y++)
	      {
		    data_complex[input_vol.xsize()*y+x][0] = double(input_vol.value(x,y,z,t));
		    data_complex[input_vol.xsize()*y+x][1] = 0;
	      }
	    unring_2d(data_complex,res_complex, dim_sz,nsh,minW,maxW);               
	    for (int x = 0 ; x < input_vol.xsize();x++)
		for (int y = 0 ; y < input_vol.ysize();y++)
		  output_vol.value(x,y,z,t) = res_complex[input_vol.xsize()*y+x][0];
	}
	
	fftw_free(data_complex);
	fftw_free(res_complex);
	
    } 
    else if (alongdim.compare(string("1"))==0)
    {
	dim_sz[0] = input_vol.ysize();
	dim_sz[1] = input_vol.zsize();
	fftw_complex *data_complex =  (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * dim_sz[0]*dim_sz[1]);
	fftw_complex *res_complex  =  (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * dim_sz[0]*dim_sz[1]);

	for (int x=0; x<input_vol.xsize(); x++)       
	{
	  for (int z = 0 ; z < input_vol.zsize();z++)
	      for (int y = 0 ; y < input_vol.ysize();y++)
	      {
		    data_complex[input_vol.ysize()*z+y][0] = double(input_vol.value(x,y,z,t));
		    data_complex[input_vol.ysize()*z+y][1] = 0;
	      }
	    unring_2d(data_complex,res_complex, dim_sz,nsh,minW,maxW);               
	    for (int z = 0 ; z < input_vol.zsize();z++)
		for (int y = 0 ; y < input_vol.ysize();y++)
		  output_vol.value(x,y,z,t) = res_complex[input_vol.ysize()*z+y][0];
	}
	fftw_free(data_complex);
	fftw_free(res_complex);
	
     }
    else if (alongdim.compare(string("2"))==0)
    {
	dim_sz[0] = input_vol.xsize();
	dim_sz[1] = input_vol.zsize();
	fftw_complex *data_complex =  (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * dim_sz[0]*dim_sz[1]);
	fftw_complex *res_complex  =  (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * dim_sz[0]*dim_sz[1]);

	for (int y=0; y<input_vol.ysize(); y++)       
	{
	  for (int z = 0 ; z < input_vol.zsize();z++)
	      for (int x = 0 ; x < input_vol.xsize();x++)
	      {
		    data_complex[input_vol.xsize()*z+x][0] = double(input_vol.value(x,y,z,t));
		    data_complex[input_vol.xsize()*z+x][1] = 0;
	      }
	    unring_2d(data_complex,res_complex, dim_sz,nsh,minW,maxW);               
	    for (int z = 0 ; z < input_vol.zsize();z++)
		for (int x = 0 ; x < input_vol.xsize();x++)
		  output_vol.value(x,y,z,t) = res_complex[input_vol.xsize()*z+x][0];
	}
	fftw_free(data_complex);
	fftw_free(res_complex);
	
     }
     else 
     {
	fprintf(stderr,"wrong parameter value for -d\n");
	return -1 ;
     }
  
  }
  
  fprintf(stderr,"\n"); 
   
    
  
  
  
  
  
  save_volume4D(output_vol,output_name);
  return 0;
}


int main(int argc,char *argv[])
{

  Tracer tr("main");

  try {
    string progname=argv[0];
    if (argc <3)
      { 
	print_usage(progname);
	return 1; 
      }
    
    string iname=string(argv[1]);
    return call_fmrib_main(dtype(iname),argc,argv); 
  } catch(std::exception &e) {
    cerr << e.what() << endl;
  }  catch(Exception &e) {
    exit(EXIT_FAILURE);
  } 
}


