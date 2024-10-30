#include"../../include/MEP_include/MEP.h"
QB_tools input;
int MEP_main(int argc,char **argv)
{
		char pattern[1024];
		strcpy(pattern,argv[3]);	
		int step_mep=(int)QB_checkdat(argv[4]);	//points
		printf("%d\n",step_mep);
		if(!strcmp(argv[5],"1"))
		{
			QB_init(&input);
			QB_read_vasp(&input,argv[1]);
			QB_tools input2;
			QB_init(&input2);
			QB_read_vasp(&input2,argv[2]);
			int max_step=(int)QB_checkdat(argv[6]);
			int max_try=(int)QB_checkdat(argv[7]);
			double min_delta=QB_checkdat(argv[8]);
			double step_length=QB_checkdat(argv[9]);
			printf("stageI");
			//get init 
			SPG_tools start;
			start.num=start.ele_n=0;
			QB2SPG(&input,&start);
	
			for(int i=0;i<start.num;i++)
			{
				double in[3]={start.pos[i][0],start.pos[i][1],start.pos[i][2]};
				double out[3];
				QB_f2c(&input,in,out);
				start.pos[i][0]=out[0];
				start.pos[i][1]=out[1];
				start.pos[i][2]=out[2];
			}
			//get end
			SPG_tools end;
			end.num=end.ele_n=0;
			QB2SPG(&input2,&end);
			for(int j=0;j<start.num;j++)
			{
				double in[3]={end.pos[j][0],end.pos[j][1],end.pos[j][2]};
				double out[3];
				QB_f2c(&input,in,out);
				end.pos[j][0]=out[0];
				end.pos[j][1]=out[1];
				end.pos[j][2]=out[2];
			}
			double **m_start;
			double **m_end;
			init_distance_matrix(start.num,&m_start);
			init_distance_matrix(end.num  ,&m_end);
			get_distance_matrix_correction(start,m_start);
			get_distance_matrix_correction(end  ,m_end);
			
			for(int i_mep=0;i_mep<step_mep;i_mep++)
			{
				printf("%d-%d-%d\n",i_mep,step_mep,i_mep<step_mep);
				double delta=(double)(i_mep+1)/(double)(step_mep+1);
				//get linear matrix and linear position
				SPG_tools lin;
				lin.num=lin.ele_n=0;
				QB2SPG(&input,&lin);
		
				for(int j=0;j<lin.num;j++)
					for(int k=0;k<3;k++)
						lin.pos[j][k]=start.pos[j][k]*(1-delta)+end.pos[j][k]*delta;
		
				double **m_lin;
				init_distance_matrix(lin.num,&m_lin);
				for(int j=0;j<lin.num;j++)
					for(int k=0;k<lin.num;k++)
						m_lin[j][k]=m_start[j][k]*(1-delta)+m_end[j][k]*delta;
		
				//init temp postion
				SPG_tools tem;
				tem.num=tem.ele_n=0;
				QB2SPG(&input,&tem);
		
				//get temp matrix
				double ** m_tem;
				init_distance_matrix(lin.num,&m_tem);
		
				double ** m_tem2;
				init_distance_matrix(lin.num,&m_tem2);
 				//begin linesearch
				for(int j=0;j<max_step;j++)
				{
					get_distance_matrix(lin,m_tem);
					double e0=energy_distance_matrix(lin.num,m_lin,m_tem);
					double e1;
					for(int k=0;k<max_try;k++)
					{				
						//calculate force
						for(int l=0;l<tem.num;l++)
						{
							double grad[3]={0,0,0};
							double dis[3],fdis[3];
							for(int m=0;m<tem.num;m++)
								if(l!=m)
								{
									for(int n=0;n<3;n++)
									{
										dis[n]=lin.pos[l][n]-lin.pos[m][n];
									}
									QB_c2f(&input,dis,fdis);
									for(int n=0;n<3;n++)
									{	
										while(fdis[n]> 0.5)fdis[n]-=1;
										while(fdis[n]<-0.5)fdis[n]+=1;
									}
									QB_f2c(&input,fdis,dis);
									double dis_l=sqrt(dis[0]*dis[0]+dis[1]*dis[1]+dis[2]*dis[2]);
									for(int n=0;n<3;n++)
									grad[n]+=2*(m_lin[l][m]-m_tem[l][m])*dis[n]/dis_l*(2*m_lin[l][m]-m_tem[l][m])/pow(m_tem[l][m],5); 
								}
							for(int n=0;n<3;n++)
								tem.pos[l][n]=grad[n];
						}
						for(int l=0;l<tem.num;l++)
						{
							for(int n=0;n<3;n++)
								tem.pos[l][n]=lin.pos[l][n]+tem.pos[l][n]*pow(RATE,k)*step_length;
						}
						get_distance_matrix(tem,m_tem2);
						e1=energydif_distance_matrix(lin.num,m_lin,m_tem2,m_tem);
						//printf("%d-%d-E=%lf\n",j,k,e1);
						//energy decrease
						if(e1<0)
						{
							for(int l=0;l<tem.num;l++)
								for(int m=0;m<3;m++)
									lin.pos[l][m]=tem.pos[l][m];
							break;
						}
					}	
					//step is too short
					if(fabs(e1)<min_delta)
						break;
				}
				free_distance_matrix(lin.num,m_tem);
				free_distance_matrix(lin.num,m_lin);
				QB_tools output;
				QB_init(&output);
				for(int j=0;j<start.num;j++)
				{
					double in[3]={lin.pos[j][0],lin.pos[j][1],lin.pos[j][2]};
					double out[3];
					QB_c2f(&input,in,out);
					lin.pos[j][0]=out[0];
					lin.pos[j][1]=out[1];
					lin.pos[j][2]=out[2];
				}
				SPG2QB(&lin,&output);
				char name[1024];
				sprintf(name,"%s_%d.vasp",pattern,i_mep+1);
				QB_dump_vasp_Direct(&output,name);
				QB_free_atom(&output);
				SPG_free(&lin);
				SPG_free(&tem);
			}
		}

		else if(!strcmp(argv[5],"0"))
		{
			QB_tools input1;
			QB_init(&input1);
			QB_read_file(&input1,argv[1]);
			QB_tools input2;
			QB_init(&input2);
			QB_read_file(&input2,argv[2]);
			QB_tools input3;
			QB_init(&input3);
			QB_read_file(&input3,argv[2]);
			QB_tools output;
			QB_init(&output);
			for(int j=1;j<=step_mep;j++)
			{
      			for(int i=0;i<input1.TotalNumber;i++)	
				{
					if (input1.atom[i].x==input2.atom[i].x && input1.atom[i].y==input2.atom[i].y && input1.atom[i].z==input2.atom[i].z)
					{
						input3.atom[i].x=input1.atom[i].x;
						input3.atom[i].y=input1.atom[i].y;
						input3.atom[i].z=input1.atom[i].z;
					}
 					else if(input1.atom[i].x!=input2.atom[i].x && input1.atom[i].y==input2.atom[i].y && input1.atom[i].z==input2.atom[i].z)
					{
						double d1=input2.atom[i].x-input1.atom[i].x;
						double n1=d1/(step_mep+1);			
						input3.atom[i].x=input1.atom[i].x+j*n1;
						input3.atom[i].y=input1.atom[i].y;
						input3.atom[i].z=input1.atom[i].z;
					}
					else if(input1.atom[i].x==input2.atom[i].x && input1.atom[i].y!=input2.atom[i].y && input1.atom[i].z==input2.atom[i].z)
					{
						double d2=input2.atom[i].y-input1.atom[i].y;
						double n2=d2/(step_mep+1);			
						input3.atom[i].x=input1.atom[i].x;
						input3.atom[i].y=input1.atom[i].y+j*n2;
						input3.atom[i].z=input1.atom[i].z;
					}
					else if(input1.atom[i].x==input2.atom[i].x && input1.atom[i].y==input2.atom[i].y && input1.atom[i].z!=input2.atom[i].z)
					{
						double d3=input2.atom[i].z-input1.atom[i].z;
						double n3=d3/(step_mep+1);			
						input3.atom[i].x=input1.atom[i].x;
						input3.atom[i].y=input1.atom[i].y;
						input3.atom[i].z=input1.atom[i].z+j*n3;
					}
					else if(input1.atom[i].x!=input2.atom[i].x && input1.atom[i].y!=input2.atom[i].y && input1.atom[i].z==input2.atom[i].z)
					{
						double k1=(input1.atom[i].y-input2.atom[i].y)/(input1.atom[i].x-input2.atom[i].x);
						double b1=input1.atom[i].y-k1*input1.atom[i].x;
						double d4=input2.atom[i].x-input1.atom[i].x;
						double n4=d4/(step_mep+1);			
						input3.atom[i].x=input1.atom[i].x+j*n4;
						input3.atom[i].y=k1*input3.atom[i].x+b1;
						input3.atom[i].z=input1.atom[i].z;
					}	
					else if(input1.atom[i].x!=input2.atom[i].x && input1.atom[i].y==input2.atom[i].y && input1.atom[i].z!=input2.atom[i].z)
					{
						double k2=(input1.atom[i].z-input2.atom[i].z)/(input1.atom[i].x-input2.atom[i].x);
						double b2=input1.atom[i].z-k2*input1.atom[i].x;
						double d5=input2.atom[i].x-input1.atom[i].x;
						double n5=d5/(step_mep+1);			
						input3.atom[i].x=input1.atom[i].x+j*n5;
						input3.atom[i].y=input1.atom[i].y;
						input3.atom[i].z=k2*input3.atom[i].x+b2;
					}	
					else if(input1.atom[i].x==input2.atom[i].x && input1.atom[i].y!=input2.atom[i].y && input1.atom[i].z!=input2.atom[i].z)
					{
						double k3=(input1.atom[i].z-input2.atom[i].z)/(input1.atom[i].y-input2.atom[i].y);
						double b3=input1.atom[i].z-k3*input1.atom[i].y;
						double d6=input2.atom[i].y-input1.atom[i].y;
						double n6=d6/(step_mep+1);			
						input3.atom[i].x=input1.atom[i].x;
						input3.atom[i].y=input1.atom[i].y+j*n6;
						input3.atom[i].z=k3*input3.atom[i].y+b3;
					}	
					else		
					{
						double k=(input1.atom[i].y-input2.atom[i].y)/(input1.atom[i].x-input2.atom[i].x);
						double b=input1.atom[i].y-k*input1.atom[i].x;
						double d=input2.atom[i].x-input1.atom[i].x;
						double n=d/(step_mep+1);
						input3.atom[i].x=input1.atom[i].x+j*n;
						input3.atom[i].y=k*input3.atom[i].x+b;
						input3.atom[i].z=(input3.atom[i].x-input1.atom[i].x)/(input2.atom[i].x-input1.atom[i].x)*(input2.atom[i].z-input1.atom[i].z)+input1.atom[i].z;
					}
					char name[1024];
					sprintf(name,"%s_%d.vasp",pattern,j);
					QB_dump_vasp_Direct(&input3,name);
				}
			}
			QB_free_atom(&input1);
			QB_free_atom(&input2);
			QB_free_atom(&input3);
		}
		else
		{
			printf("Error!/n");
		}
}
void init_distance_matrix(int num,double***matrix)
{
	(*matrix)=(double**)malloc(num*sizeof(double*));
	for(int i=0;i<num;i++)
		(*matrix)[i]=(double*)malloc(num*sizeof(double));
}
void get_distance_matrix(SPG_tools spg,double**matrix)
{
	for(int i=0;i<spg.num;i++)
	for(int j=0;j<spg.num;j++)
	{
		double dis[3];
		double fdis[3];
		for(int k=0;k<3;k++)
		{
			dis[k]=spg.pos[i][k]-spg.pos[j][k];
		}
		QB_c2f(&input,dis,fdis);
		for(int k=0;k<3;k++)
		{	
			//while(fdis[k]> 0.5)fdis[k]-=1;
			//while(fdis[k]<-0.5)fdis[k]+=1;
		}
		QB_f2c(&input,fdis,dis);
		matrix[i][j]=sqrt(dis[0]*dis[0]+dis[1]*dis[1]+dis[2]*dis[2]);
	}
}
void get_distance_matrix_correction(SPG_tools spg,double**matrix)
{
	for(int i=0;i<spg.num;i++)
	for(int j=0;j<spg.num;j++)
	{
		double dis[3];
		double fdis[3];
		for(int k=0;k<3;k++)
		{
			dis[k]=spg.pos[i][k]-spg.pos[j][k];
		}
		QB_c2f(&input,dis,fdis);
		for(int k=0;k<3;k++)
		{	
			while(fdis[k]> 0.5)fdis[k]-=1;
			while(fdis[k]<-0.5)fdis[k]+=1;
		}
		QB_f2c(&input,fdis,dis);
		matrix[i][j]=sqrt(dis[0]*dis[0]+dis[1]*dis[1]+dis[2]*dis[2]);
	}
}

void free_distance_matrix(int num,double**matrix)
{
	for(int i=0;i<num;i++)
	{
		free(matrix[i]);
	}
	free(matrix);
}

double energy_distance_matrix(int num,double**m_lin,double**m_cur)
{
	double energy=0;
	for(int i=0;i<num;i++)
	for(int j=0;j<num;j++)
	if(i!=j)
		energy+=(m_lin[i][j]-m_cur[i][j])*(m_lin[i][j]-m_cur[i][j])/pow(m_cur[i][j],4);
	return energy;
}

double energydif_distance_matrix(int num,double**m_lin,double**m_cur,double**m_cur2)
{
	double energy=0;
	for(int i=0;i<num;i++)
	for(int j=0;j<num;j++)
	if(i!=j)
		energy+=(m_lin[i][j]-m_cur[i][j])*(m_lin[i][j]-m_cur[i][j])/pow(m_cur[i][j],4)
		       -(m_lin[i][j]-m_cur2[i][j])*(m_lin[i][j]-m_cur2[i][j])/pow(m_cur2[i][j],4);
	return energy;
}
