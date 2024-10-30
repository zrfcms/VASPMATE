#include"./../QB/QB.h"
#include"QSPG.h"

void QSPG_primitive(QB_tools* output,QB_tools* input,double symprec)
{
	SPG_tools spg;
	spg.num=0;
	QB2SPG(input,&spg);
	spg.num=spg_find_primitive(spg.mat,
		      spg.pos,
		       spg.type,
		       spg.num,
		       symprec);
	if(spg.num!=0)
	{
		SPG2QB(&spg,output);
		SPG_free(&spg);
	}
	else if(input!=output)//no primitive cell found
	{
		QB_copy_system(output,input);
	}
}

void QSPG_refined(QB_tools* output,QB_tools* input,double symprec)
{
	SPG_tools spg;
	spg.num=0;
	QB2SPG(input,&spg);
	if(spg.num)
	{
		spg.pos=(double(*)[3])realloc(spg.pos,sizeof(double)*12*spg.num);
		spg.type=(int*)realloc(spg.type,4*sizeof(int)*spg.num);
	}
	spg.num=spg_refine_cell(spg.mat,
		      spg.pos,
		       spg.type,
		       spg.num,
		       symprec);
	if(spg.num!=0)
	{
		SPG2QB(&spg,output);
		SPG_free(&spg);
	}
	else if(input!=output)//no primitive cell found
	{
		QB_copy_system(output,input);
	}
}

void QSPG_get_symmetry(QB_tools*input,QSPGDataset**data,double symprec)
{
	SPG_tools spg;
	spg.num=0;
	QB2SPG(input,&spg);
	if(spg.num!=0)
	{
		spg.pos=(double(*)[3])realloc(spg.pos,sizeof(double)*12*spg.num);
		spg.type=(int*)realloc(spg.type,sizeof(int)*4*spg.num);
	}
	*data=spg_get_dataset(spg.mat,
									 spg.pos,
									 spg.type,
									 spg.num,
									 symprec);
	if(spg.num!=0)
		SPG_free(&spg);
}