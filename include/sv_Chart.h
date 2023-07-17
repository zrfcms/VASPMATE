#pragma once
#ifndef _CHART_
#define _CHART_
#include <stdlib.h>
#include <time.h> // time() - used to seed rand()
#include <cstring>
#include<cstdio>
#include<vector>
#include "Enumerations.H"
using std::vector;
#define FL_NO_SYMBOL_SVMOD 	 			0
#define FL_CIRCLE_SVMOD 	 			1
#define FL_SQUARE_SVMOD 	 			2
#define FL_SQUARE_UP_SVMOD 				3
#define FL_TRIANGLE_UP_SVMOD			4
#define FL_TRIANGLE_DOWN_SVMOD			5
#define FL_TRIANGLE_LEFT_SVMOD			6
#define FL_TRIANGLE_RIGHT_SVMOD 		7
#define FL_CIRCLE_EMPTY_SVMOD 	 		8
#define FL_SQUARE_EMPTY_SVMOD 	 		9
#define FL_SQUARE_UP_EMPTY_SVMOD 		10
#define FL_TRIANGLE_UP_EMPTY_SVMOD		11
#define FL_TRIANGLE_DOWN_EMPTY_SVMOD	12
#define FL_TRIANGLE_LEFT_EMPTY_SVMOD	13
#define FL_TRIANGLE_RIGHT_EMPTY_SVMOD 	14
#define FL_SYMBOL_NUM					15

#define FL_LEGEND_RIGHTTOP_SVMOD		0
#define FL_LEGEND_RIGHTBOTTOM_SVMOD		1
#define FL_LEGEND_LEFTTOP_SVMOD			2
#define FL_LEGEND_LEFTBOTTOM_SVMOD		3
#define FL_LEGEND_MIDDLETOP_SVMOD		4
#define FL_LEGEND_MIDDLEBOTTOM_SVMOD	5

#define FL_TICK_NONE_SVMOD				0
#define FL_TICK_INSIDE_SVMOD			1
#define FL_TICK_OUTSIDE_SVMOD			2

#define FL_LENGEND_NO_BORDER			0
#define FL_LENGEND_BLACK_BORDER			1

#define FL_LENGEND_NO_BORDER			0
#define FL_LENGEND_BLACK_BORDER			1

#define FL_DISPLAY_DECIMAL				0
#define FL_DISPLAY_SCIENTIFIC			1

#define FL_NO_LINE						5

enum {
  FL_SOLID	= 0,		///< line style: <tt>___________</tt>
  FL_DASH	= 1,		///< line style: <tt>_ _ _ _ _ _</tt>
  FL_DOT	= 2,		///< line style: <tt>. . . . . .</tt>
  FL_DASHDOT	= 3,		///< line style: <tt>_ . _ . _ .</tt>
  FL_DASHDOTDOT	= 4,		///< line style: <tt>_ . . _ . .</tt>

  FL_CAP_FLAT	= 0x100,	///< cap style: end is flat
  FL_CAP_ROUND	= 0x200,	///< cap style: end is round
  FL_CAP_SQUARE	= 0x300,	///< cap style: end wraps end point

  FL_JOIN_MITER	= 0x1000,	///< join style: line join extends to a point
  FL_JOIN_ROUND	= 0x2000,	///< join style: line join is rounded
  FL_JOIN_BEVEL	= 0x3000	///< join style: line join is tidied
};

typedef struct
{
	int num;
	float (*x)[4];
	unsigned col;
	int linewidth;
	int linestyle;
	int symbol;
	int symbol_size;
	char name[1024];
	bool name_flag;
	bool show_flag;
	double colbar_min;
	double colbar_max;
	int linemode;
}Fl_Chart_SVMOD_ENTRY;

typedef struct
{
	float x;
	bool yaxis;
	unsigned col;
	int linewidth;
	int linestyle;
	char name[1024];
	bool show_flag;
	bool dash_flag;
	unsigned text_col;
	Fl_Font font;
    Fl_Fontsize size;
}Fl_Chart_SVMOD_LABEL;

class Fl_Chart_SVMOD{
	public:
	int numb;
    Fl_Chart_SVMOD_ENTRY *entries;
	int numbl;
	Fl_Chart_SVMOD_LABEL *labelist;
    double xmin,xmax;
	double ymin,ymax;
	double xspace,yspace;
	double xscale,yscale,hscale,wscale;
	bool   xlabel_show,ylabel_show;
	char   xlabel[1024],ylabel[1024];
	Fl_Color xaxis_color;
	Fl_Color yaxis_color;
	int    axis_width;
	int    legend_mode;
	bool   legend_flag;
	int	   legend_background;
	int    tick_mode;
	
	Fl_Font title_font;
    Fl_Fontsize title_size;
    
	Fl_Font axis_font;
    Fl_Fontsize axis_size;
	
	Fl_Font legend_font;
    Fl_Fontsize legend_size;
	
    Fl_Chart_SVMOD()
	{
		xscale	   = 0.14;
		yscale	   = 0.01;
		wscale	   = 0.85;
		hscale	   = 0.86;
		numb       = 0;
		numbl	   = 0;
		xmin = xmax= 0;
		ymin = ymax= 0;
		tick_mode  = 1;
		
		title_font = FL_TIMES;
		title_size = 22;
    
		axis_font  = FL_TIMES;
		axis_size  = 22;
	
		legend_font= FL_TIMES;
		legend_size= 15;
		legend_mode= FL_LEGEND_RIGHTTOP_SVMOD;
		strcpy(xlabel,"");
		strcpy(ylabel,"");
		legend_flag= 1;
		legend_background = FL_LENGEND_BLACK_BORDER;
		xaxis_color = FL_BLACK;
		yaxis_color = FL_BLACK;
		axis_width = 1;
		xlabel_show = 1;
		ylabel_show = 1;
	}
    ~Fl_Chart_SVMOD();
	void add(int n, vector<double> x, vector<double> y,const char name[1024]="No Label",unsigned col=FL_BLACK,int linestyle=FL_SOLID,int linewidth=2,int symbol=FL_NO_SYMBOL_SVMOD,int symbol_size=3);
	void add(int n, vector<double> x, vector<double> y, vector<double> r,const char name[1024]="No Label",unsigned col=FL_BLACK,int linestyle=FL_SOLID,int linewidth=2,int symbol=FL_NO_SYMBOL_SVMOD,int symbol_size=3);
	
	void add(int n, double* x, double* y, const char name[1024] = "No Label", unsigned col = FL_BLACK, int linestyle = FL_SOLID, int linewidth = 2, int symbol = FL_NO_SYMBOL_SVMOD, int symbol_size = 3);
	void add(int n, double* x, double* y, double* r, const char name[1024] = "No Label", unsigned col = FL_BLACK, int linestyle = FL_SOLID, int linewidth = 2, int symbol = FL_NO_SYMBOL_SVMOD, int symbol_size = 3);

	void clear();
	
	void set_scale(double xs1,double ys1,double xs2,double ys2);
	void set_axis_width(int in_width);
	void set_axis_color(Fl_Color in_x,Fl_Color in_y);
	
	void set_name_flag(int i,bool name_flag);
	void set_show_flag(int i,bool show_flag);
	void set_label_show_flag(int i,bool show_flag);
	
	void add_label(double x,const char name[1024],bool axis=0,Fl_Color text_col=FL_BLACK,Fl_Font font=FL_TIMES,Fl_Fontsize size=15,
				   bool dash_flag=1,unsigned col=FL_BLACK,int linestyle=FL_DASH,int linewidth=1);
	void clear_label();

	void bound(double xmi,double ymi,double xma,double yma);
	void auto_bound(double space);
	
	void titles(const char title_x[1024],const char title_y[1024]);
	void show_titles(bool show_titlex,bool show_titley);
	void set_titles_text(Fl_Font font,Fl_Fontsize size){title_font=font;title_size=size;}
	void space(double xs,double ys);
	
	void read(const char name[1024]);
	void write(const char name[1024]);
};

#endif