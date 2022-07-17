/*
This file is part of ftdock, a program for rigid-body protein-protein docking 
Copyright (C) 1997-2000 Gidon Moont

Biomolecular Modelling Laboratory
Imperial Cancer Research Fund
44 Lincoln's Inn Fields
London WC2A 3PX

+44 (0)20 7269 3348
http://www.bmm.icnet.uk/

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

*/

#include "structures.h"
__global__ void assign_charges_on_GPU(struct Amino_Acid *Residue)
{
  int residue=threadIdx.y;
  int atom=threadIdx.x;


  if((residue>0)&&(atom>0)&&(atom<Residue[residue].size)){

    Residue[residue].Atom[atom].charge = 0.0;
    /* peptide backbone */

    if( strcmp(Residue[residue].Atom[atom].atom_name , " N  " ) == 0 ) {
        if( strcmp(Residue[residue].res_name , "PRO" ) == 0 ) {
          Residue[residue].Atom[atom].charge = -0.10 ;
        } else {
          Residue[residue].Atom[atom].charge =  0.55 ;
          if( residue == 1 )Residue[residue].Atom[atom].charge = 1.00 ;
        }
      }


    if( strcmp( This_Structure.Residue[residue].Atom[atom].atom_name , " O  " ) == 0 ) {
        This_Structure.Residue[residue].Atom[atom].charge = -0.55 ;
        if( residue == This_Structure.length  ) This_Structure.Residue[residue].Atom[atom].charge = -1.00 ;
      }
     /* charged residues */

      if( ( strcmp( This_Structure.Residue[residue].res_name , "ARG" ) == 0 ) && ( strncmp( This_Structure.Residue[residue].Atom[atom].atom_name , " NH" , 3 ) == 0 ) ) This_Structure.Residue[residue].Atom[atom].charge =  0.50 ;
      if( ( strcmp( This_Structure.Residue[residue].res_name , "ASP" ) == 0 ) && ( strncmp( This_Structure.Residue[residue].Atom[atom].atom_name , " OD" , 3 ) == 0 ) ) This_Structure.Residue[residue].Atom[atom].charge = -0.50 ;
      if( ( strcmp( This_Structure.Residue[residue].res_name , "GLU" ) == 0 ) && ( strncmp( This_Structure.Residue[residue].Atom[atom].atom_name , " OE" , 3 ) == 0 ) ) This_Structure.Residue[residue].Atom[atom].charge = -0.50 ;
      if( ( strcmp( This_Structure.Residue[residue].res_name , "LYS" ) == 0 ) && ( strcmp( This_Structure.Residue[residue].Atom[atom].atom_name , " NZ " ) == 0 ) ) This_Structure.Residue[residue].Atom[atom].charge =  1.00 ;

  }
}

void assign_charges( struct Structure This_Structure ) {

/************/

  /* Counters */

  int	residue , atom,a=0 ;

/************/
struct Amino_Acid Residue[This_Structure.length],*d_Residue;
for (int i = 0; i < This_Structure.length; i++)
{
  Residue[i]=This_Structure.Residue[i];
  cudaMalloc(&Residue[i].Atom,This_Structure.Residue[i].size*sizeof(struct Atom));
  cudaMemcpy(Residue[i].Atom,This_Structure.Residue[i].Atom,This_Structure.Residue[i].size*sizeof(struct Atom),cudaMemcpyHostToDevice);
  a=max(a,This_Structure.Residue[i].size);
  
}
cudaMalloc((void**)&d_Residue,This_Structure.length*sizeof(struct Amino_Acid));
cudaMemcpy(d_Residue,Residue,This_Structure.length*sizeof(struct Amino_Acid),cudaMemcpyHostToDevice);

dim3 threadPerBlock(a,This_Structure.length);
assign_charges_on_GPU<<<1,threadPerBlock>>>(d_Residue);
cudaDeviceSynchronize;
cudaMemcpy(This_Structure.Residue,d_Residue,This_Structure.length*sizeof(struct Amino_Acid),cudaMemcpyDeviceToHost);
cudaFree(d_Residue);
free(Residue);
/************/

}



/************************/
__global__ void zero_interaction_grid(cufftReal *grid,int grid_size)
{
    x=threadIdx.x;
    y=threadIdx.y;
    z=threadIdx.z;
    grid[gaddress(x,y,z,grid_size)] = (cufftReal)0;
}
__global__ void field_calculation(float *phi,Amino_Acid *Residue)
{
   int residue=threadIdx.y;
   int atom=threadIdx.x;
   float		distance ;
   float epsilon ;
   if((residue>0)&&(atom>0)&&(atom<Residue[residue].size))
   {
      
            if(Residue[residue].Atom[atom].charge != 0 ) {

              distance = pythagoras( Residue[residue].Atom[atom].coord[1] , Residue[residue].Atom[atom].coord[2] , Residue[residue].Atom[atom].coord[3] , x_centre , y_centre , z_centre ) ;
         
              if( distance < 2.0 ) distance = 2.0 ;

              if( distance >= 2.0 ) {

                if( distance >= 8.0 ) {

                  epsilon = 80 ;

                } else { 

                  if( distance <= 6.0 ) { 

                    epsilon = 4 ;
             
                  } else {

                    epsilon = ( 38 * distance ) - 224 ;

                  }

                }
  
                *phi += (Residue[residue].Atom[atom].charge / ( epsilon * distance ) ) ;

              }

            }

   }
}
__global__ void electric_fieldonGPU(cufftReal *grid,int grid_size,float grid_span,float *phi,dim3 threadPerBlock, Amino_Acid *Residue)
{
    x=threadIdx.x;
    y=threadIdx.y;
    z=threadIdx.z;
    if (y==0&&z==0)
    {
      printf( "." );
    }
    
    x_centre  = gcentre( x , grid_span , grid_size ) ;
    y_centre  = gcentre( y , grid_span , grid_size ) ;
    z_centre  = gcentre( z , grid_span , grid_size ) ;
    *phi=0;
    field_calculation<<<1,threadPerBlock>>>(phi,Residue);
    cudaDeviceSynchronize();
    grid[gaddress(x,y,z,grid_size)] = (cufftReal)*phi ;
}




void electric_field( struct Structure This_Structure , float grid_span , int grid_size , cufftReal *grid ) {

/************/

  /* Counters */

  int	residue , atom ;

  /* Co-ordinates */

  int	x , y , z ;
  float		x_centre , y_centre , z_centre ;//scope for cuda

  /* Variables */

  float		distance ;
  float		*phi , epsilon ;

/************/

dim3 threadsperblock(grid_size,grid_size,grid_size);

cudaMalloc((void**)&phi,sizeof(float));


zero_interaction_grid<<<1,threadsperblock>>>(grid,grid_size);
cudaDeviceSynchronize();
struct Amino_Acid Residue[This_Structure.length],*d_Residue;
int a=0;
for (int i = 0; i < This_Structure.length; i++)
{
  Residue[i]=This_Structure.Residue[i];
  cudaMalloc(&Residue[i].Atom,This_Structure.Residue[i].size*sizeof(struct Atom));
  cudaMemcpy(Residue[i].Atom,This_Structure.Residue[i].Atom,This_Structure.Residue[i].size*sizeof(struct Atom),cudaMemcpyHostToDevice);
  a=max(a,This_Structure.Residue[i].size);
  
}
cudaMalloc((void**)&d_Residue,This_Structure.length*sizeof(struct Amino_Acid));
cudaMemcpy(d_Residue,Residue,This_Structure.length*sizeof(struct Amino_Acid),cudaMemcpyHostToDevice);

  dim3 threadPerBlock1(a,This_Structure.length);


/************/

  setvbuf( stdout , (char *)NULL , _IONBF , 0 ) ;

  printf( "  electric field calculations ( one dot / grid sheet ) " ) ;

  electric_fieldonGPU<<<1,threadPerBlock>>>(grid,grid_size,grid_span,phi,threadPerBlock1,d_Residue);
  cudaDeviceSynchronize();

  printf( "\n" ) ;
  cudaFree(d_Residue);
  free(Residue);
  cudaFree(phi);


/************/

  return ;

}



/************************/



void electric_point_charge( struct Structure This_Structure , float grid_span , int grid_size , cufftReal *grid ) {

/************/

  /* Counters */

  int	residue , atom ;

  /* Co-ordinates */

  int	x , y , z ;
  int	x_low , x_high , y_low , y_high , z_low , z_high ;

  float		a , b , c ;
  float		x_corner , y_corner , z_corner ;
  float		w ;

  /* Variables */

  float		one_span ;

/************/
dim3 threadsperblock(grid_size,grid_size,grid_size);
zero_interaction_grid<<<1,threadsperblock>>>(grid,grid_size);
cudaDeviceSynchronize();


/************/

  one_span = grid_span / (float)grid_size ;

  for( residue = 1 ; residue <= This_Structure.length ; residue ++ ) {
    for( atom = 1 ; atom <= This_Structure.Residue[residue].size ; atom ++ ) {

      if( This_Structure.Residue[residue].Atom[atom].charge != 0 ) {

        x_low = gord( This_Structure.Residue[residue].Atom[atom].coord[1] - ( one_span / 2 ) , grid_span , grid_size ) ;
        y_low = gord( This_Structure.Residue[residue].Atom[atom].coord[2] - ( one_span / 2 ) , grid_span , grid_size ) ;
        z_low = gord( This_Structure.Residue[residue].Atom[atom].coord[3] - ( one_span / 2 ) , grid_span , grid_size ) ;

        x_high = x_low + 1 ;
        y_high = y_low + 1 ;
        z_high = z_low + 1 ;

        a = This_Structure.Residue[residue].Atom[atom].coord[1] - gcentre( x_low , grid_span , grid_size ) - ( one_span / 2 ) ;
        b = This_Structure.Residue[residue].Atom[atom].coord[2] - gcentre( y_low , grid_span , grid_size ) - ( one_span / 2 ) ;
        c = This_Structure.Residue[residue].Atom[atom].coord[3] - gcentre( z_low , grid_span , grid_size ) - ( one_span / 2 ) ;

        for( x = x_low ; x <= x_high  ; x ++ ) {
 
          x_corner = one_span * ( (float)( x - x_high ) + .5 ) ;

          for( y = y_low ; y <= y_high  ; y ++ ) {

            y_corner = one_span * ( (float)( y - y_high ) + .5 ) ;

            for( z = z_low ; z <= z_high  ; z ++ ) {

              z_corner = one_span * ( (float)( z - z_high ) + .5 ) ;

              w = ( ( x_corner + a ) * ( y_corner + b ) * ( z_corner + c ) ) / ( 8.0 * x_corner * y_corner * z_corner ) ;

              grid[gaddress(x,y,z,grid_size)] += (cufftReal)( w * This_Structure.Residue[residue].Atom[atom].charge ) ;

            }
          }
        }

      }

    }
  }

/************/

  return ;

}



/************************/



__global__ void electric_field_zero_core( int grid_size , cufftReal *elec_grid , cufftReal *surface_grid , float internal_value ) {

/************/

  /* Co-ordinates */

  int	x=threadIdx.x,y=threadIdx.y,z=threadIdx.z;

/************/

  
        if( surface_grid[gaddress(x,y,z,grid_size)] == (cufftReal)internal_value ) elec_grid[gaddress(x,y,z,grid_size)] = (cufftReal)0 ;



/************/


}
