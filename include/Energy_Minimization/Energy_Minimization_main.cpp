#include"QB/QB.h"
#include"QSPG/QSPG.h"
#include"Energy_Minimization_main.h"
int Energy_Minimization_main(int argc,char*argv[])
{
	QB_tools QB;
	QB_init(&QB);
	QB_read_file(&QB,argv[1]);
	//QB_add_exdouble(&QB,"v_x");	
	//QB_add_exdouble(&QB,"v_y");	
	//QB_add_exdouble(&QB,"v_z");	
	//int v_sx=QB_slot_get(&QB,"v_x");
	//int v_sy=QB_slot_get(&QB,"v_y");
	//int v_sz=QB_slot_get(&QB,"v_z");
	int max_step=(int)QB_checkdat(argv[3]); //100000
	double scale=QB_checkdat(argv[4]); //0.001
	//double cutoff1=QB_checkdat(argv[5]);
	//double cutoff2=QB_checkdat(argv[6]);
	QB_vector *vel=(QB_vector*)malloc(QB.TotalNumber*sizeof(QB_vector));
	double vecc[3];
	double vecf[3];
	double veca[3];
	double dis;
	double force;
	int active_flag;
	for(int ii=0;ii<max_step;ii++)
	{
		active_flag=0;
		for(int i=0;i<QB.TotalNumber;i++)
			vel[i].x=vel[i].y=vel[i].z=0;
		for(int i=0;i<QB.TotalNumber;i++)	
		for(int j=i+1;j<QB.TotalNumber;j++)
		{
			vecc[0]=QB.atom[j].x-QB.atom[i].x;
			vecc[1]=QB.atom[j].y-QB.atom[i].y;
			vecc[2]=QB.atom[j].z-QB.atom[i].z;
			QB_c2f(&QB,vecc,vecf);
			while(vecf[0]>0.5)vecf[0]-=1;
			while(vecf[1]>0.5)vecf[1]-=1;
			while(vecf[2]>0.5)vecf[2]-=1;
			while(vecf[0]<-0.5)vecf[0]+=1;
			while(vecf[1]<-0.5)vecf[1]+=1;
			while(vecf[2]<-0.5)vecf[2]+=1;
			QB_f2c(&QB,vecf,vecc);
			dis=1.0/sqrt(vecc[0]*vecc[0]+vecc[1]*vecc[1]+vecc[2]*vecc[2]);
			//if(dis<cutoff2)
			{
				force=-scale*dis*dis;
				vel[i].x+=vecc[0]*force;
				vel[i].y+=vecc[1]*force;
				vel[i].z+=vecc[2]*force;
				vel[j].x-=vecc[0]*force;
				vel[j].y-=vecc[1]*force;
				vel[j].z-=vecc[2]*force;
				active_flag=1;
			}	
		}
		for(int i=0;i<QB.TotalNumber;i++)
		{
			QB.atom[i].x+=vel[i].x;
			QB.atom[i].y+=vel[i].y;
			QB.atom[i].z+=vel[i].z;
			QB_wrap(&QB,i);
		}
		if(!active_flag)
			break;
	}
	for(int i=0;i<QB.TotalNumber;i++)
	{ 
		QB_wrap(&QB,i);//wrap an atom back to box;
		//QB_slot_save(&QB,v_sx,i,vel[i].x);
		//QB_slot_save(&QB,v_sy,i,vel[i].y);
		//QB_slot_save(&QB,v_sz,i,vel[i].z);
	}
	QB_dump_vasp_Cartesian(&QB,argv[2]);
	//QB_dump_lmc(&QB,"log.lmc");
	return 0;
}