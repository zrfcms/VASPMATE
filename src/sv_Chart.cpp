#include"sv_Chart.h"

Fl_Chart_SVMOD::~Fl_Chart_SVMOD()
{
	if (numb > 0)
	{
		for (int i = 0; i < numb; i++)
		{
			free(entries[i].x);
		}
		free(entries);
	}
	if (numbl > 0)free(labelist);
}

void Fl_Chart_SVMOD::clear()
{
	if (numb > 0)
	{
		for (int i = 0; i < numb; i++)
		{
			free(entries[i].x);
		}
		free(entries);
	}
	numb = 0;
}

void Fl_Chart_SVMOD::clear_label()
{
	if (numbl > 0)free(labelist);
	numbl = 0;
}

void Fl_Chart_SVMOD::add(int n, vector<double> x, vector<double> y, const char name[1024], unsigned col, int linestyle, int linewidth, int symbol, int symbol_size)
{
	if (numb == 0)
		entries = (Fl_Chart_SVMOD_ENTRY*)malloc(sizeof(Fl_Chart_SVMOD_ENTRY));
	else
		entries = (Fl_Chart_SVMOD_ENTRY*)realloc(entries, (numb + 1) * sizeof(Fl_Chart_SVMOD_ENTRY));
	entries[numb].x = (float(*)[4])malloc(4 * n * sizeof(float));
	for (int i = 0; i < n; i++)
	{
		entries[numb].x[i][0] = x[i];
		entries[numb].x[i][1] = y[i];
		entries[numb].x[i][2] = 1;
		entries[numb].x[i][3] = 0;
	}
	if (name != NULL)
	{
		strcpy(entries[numb].name, name);
		entries[numb].name_flag = 1;
	}
	else entries[numb].name_flag = 0;
	entries[numb].num = n;
	entries[numb].col = col;
	entries[numb].symbol = symbol;
	entries[numb].symbol_size = symbol_size;
	entries[numb].linestyle = linestyle;
	entries[numb].linewidth = linewidth;
	entries[numb].show_flag = 1;
	entries[numb].linemode = 0;
	numb++;
}

void Fl_Chart_SVMOD::add(int n, double* x, double* y, const char name[1024], unsigned col, int linestyle, int linewidth, int symbol, int symbol_size)
{
	if (numb == 0)
		entries = (Fl_Chart_SVMOD_ENTRY*)malloc(sizeof(Fl_Chart_SVMOD_ENTRY));
	else
		entries = (Fl_Chart_SVMOD_ENTRY*)realloc(entries, (numb + 1) * sizeof(Fl_Chart_SVMOD_ENTRY));
	entries[numb].x = (float(*)[4])malloc(4 * n * sizeof(float));
	for (int i = 0; i < n; i++)
	{
		entries[numb].x[i][0] = x[i];
		entries[numb].x[i][1] = y[i];
		entries[numb].x[i][2] = 1;
		entries[numb].x[i][3] = 0;
	}
	if (name != NULL)
	{
		strcpy(entries[numb].name, name);
		entries[numb].name_flag = 1;
	}
	else entries[numb].name_flag = 0;
	entries[numb].num = n;
	entries[numb].col = col;
	entries[numb].symbol = symbol;
	entries[numb].symbol_size = symbol_size;
	entries[numb].linestyle = linestyle;
	entries[numb].linewidth = linewidth;
	entries[numb].show_flag = 1;
	entries[numb].linemode = 0;
	numb++;
}

void Fl_Chart_SVMOD::add(int n, vector<double> x, vector<double> y, vector<double>r, const char name[1024], unsigned col, int linestyle, int linewidth, int symbol, int symbol_size)
{
	if (numb == 0)
		entries = (Fl_Chart_SVMOD_ENTRY*)malloc(sizeof(Fl_Chart_SVMOD_ENTRY));
	else
		entries = (Fl_Chart_SVMOD_ENTRY*)realloc(entries, (numb + 1) * sizeof(Fl_Chart_SVMOD_ENTRY));
	entries[numb].x = (float(*)[4])malloc(4 * n * sizeof(float));
	for (int i = 0; i < n; i++)
	{
		entries[numb].x[i][0] = x[i];
		entries[numb].x[i][1] = y[i];
		entries[numb].x[i][2] = r[i];
		entries[numb].x[i][3] = 0;
	}
	if (name != NULL)
	{
		strcpy(entries[numb].name, name);
		entries[numb].name_flag = 1;
	}
	else entries[numb].name_flag = 0;
	entries[numb].num = n;
	entries[numb].col = col;
	entries[numb].symbol = symbol;
	entries[numb].symbol_size = symbol_size;
	entries[numb].linestyle = linestyle;
	entries[numb].linewidth = linewidth;
	entries[numb].show_flag = 1;
	entries[numb].linemode = 1;
	numb++;
}

void Fl_Chart_SVMOD::add(int n, double* x, double* y, double* r, const char name[1024], unsigned col, int linestyle, int linewidth, int symbol, int symbol_size)
{
	if (numb == 0)
		entries = (Fl_Chart_SVMOD_ENTRY*)malloc(sizeof(Fl_Chart_SVMOD_ENTRY));
	else
		entries = (Fl_Chart_SVMOD_ENTRY*)realloc(entries, (numb + 1) * sizeof(Fl_Chart_SVMOD_ENTRY));
	entries[numb].x = (float(*)[4])malloc(4 * n * sizeof(float));
	for (int i = 0; i < n; i++)
	{
		entries[numb].x[i][0] = x[i];
		entries[numb].x[i][1] = y[i];
		entries[numb].x[i][2] = r[i];
		entries[numb].x[i][3] = 0;
	}
	if (name != NULL)
	{
		strcpy(entries[numb].name, name);
		entries[numb].name_flag = 1;
	}
	else entries[numb].name_flag = 0;
	entries[numb].num = n;
	entries[numb].col = col;
	entries[numb].symbol = symbol;
	entries[numb].symbol_size = symbol_size;
	entries[numb].linestyle = linestyle;
	entries[numb].linewidth = linewidth;
	entries[numb].show_flag = 1;
	entries[numb].linemode = 1;
	numb++;
}

void Fl_Chart_SVMOD::add_label(double x, const char name[1024], bool yaxis, Fl_Color text_col, Fl_Font font, Fl_Fontsize size,
	bool dash_flag, unsigned col, int linestyle, int linewidth)
{
	if (numbl == 0)
		labelist = (Fl_Chart_SVMOD_LABEL*)malloc(sizeof(Fl_Chart_SVMOD_LABEL));
	else
		labelist = (Fl_Chart_SVMOD_LABEL*)realloc(labelist, (numbl + 1) * sizeof(Fl_Chart_SVMOD_LABEL));
	labelist[numbl].x = x;
	strcpy(labelist[numbl].name, name);
	labelist[numbl].show_flag = 1;
	labelist[numbl].dash_flag = 1;
	labelist[numbl].yaxis = yaxis;
	labelist[numbl].text_col = text_col;
	labelist[numbl].font = font;
	labelist[numbl].size = size;
	labelist[numbl].col = col;
	labelist[numbl].linestyle = linestyle;
	labelist[numbl].linewidth = linewidth;
	numbl++;
}

void Fl_Chart_SVMOD::bound(double xmi, double ymi, double xma, double yma)
{
	xmin = xmi;
	xmax = xma;
	ymin = ymi;
	ymax = yma;
}

void Fl_Chart_SVMOD::auto_bound(double space)
{
	if (numb == 0)return;
	else
	{
		xmin = entries[0].x[0][0];
		xmax = entries[0].x[0][0];
		ymin = entries[0].x[0][1];
		ymax = entries[0].x[0][1];
		for (int i = 0; i < numb; i++)
			for (int j = 0; j < entries[i].num; j++)
			{
				if (xmin > entries[i].x[j][0])xmin = entries[i].x[j][0];
				if (xmax < entries[i].x[j][0])xmax = entries[i].x[j][0];
				if (ymin > entries[i].x[j][1])ymin = entries[i].x[j][1];
				if (ymax < entries[i].x[j][1])ymax = entries[i].x[j][1];
			}
		double spacing = (xmax - xmin) * space;
		xmin -= spacing;
		xmax += spacing;
		spacing = (ymax - ymin) * space;
		ymin -= spacing;
		ymax += spacing;
	}
}

void Fl_Chart_SVMOD::set_scale(double xs1, double ys1, double xs2, double ys2)
{
	xscale = xs1;
	yscale = 1 - ys2;
	wscale = xs2 - xs1;
	hscale = ys2 - ys1;
}

void Fl_Chart_SVMOD::set_axis_width(int in_width)
{
	axis_width = in_width;
}

void Fl_Chart_SVMOD::set_axis_color(Fl_Color in_x, Fl_Color in_y)
{
	xaxis_color = in_x;
	yaxis_color = in_y;
}

void Fl_Chart_SVMOD::space(double xs, double ys)
{
	xspace = xs;
	yspace = ys;
}

void Fl_Chart_SVMOD::titles(const char titles_x[1024], const char titles_y[1024])
{
	strcpy(xlabel, titles_x);
	strcpy(ylabel, titles_y);
}

void Fl_Chart_SVMOD::show_titles(bool show_titlex, bool show_titley)
{
	xlabel_show = show_titlex;
	ylabel_show = show_titley;
}

void Fl_Chart_SVMOD::set_show_flag(int i, bool show_flag)
{
	if (i < 0)return;
	if (i > numb)return;
	entries[i].show_flag = show_flag;
}

void Fl_Chart_SVMOD::set_name_flag(int i, bool name_flag)
{
	if (i < 0)return;
	if (i > numb)return;
	entries[i].name_flag = name_flag;
}

void Fl_Chart_SVMOD::set_label_show_flag(int i, bool show_flag)
{
	if (i < 0)return;
	if (i > numbl)return;
	labelist[i].show_flag = show_flag;
}

void Fl_Chart_SVMOD::read(const char name[1024])
{
	clear();
	clear_label();
	char line[1025];
	char* token;
	char titlex[1024] = "";
	char titley[1024] = "";
	double xmi = 0, ymi = 0, xma = 1, yma = 1;
	double xs = 0, ys = 0;
	bool bound_flag = 0;
	int n;
	double* x, * y, * r;
	unsigned col = 0;
	char line_name[1024] = "";
	int linestyle = 0;
	int linewidth = 2;
	Fl_Color in_xaxis_color = FL_BLACK;
	Fl_Color in_yaxis_color = FL_BLACK;
	int in_axis_width = 0;
	int symbol = 0;
	int symbol_size = 3;
	int legend_flag = 1;
	int show_flag = 1;

	char label_name[1024] = "";
	double label_x = 0;
	bool label_axis = 0;
	int label_show_flag = 1;
	int label_dash_flag = 1;
	unsigned label_col = FL_BLACK;
	unsigned label_textcol = FL_BLACK;
	int label_linestyle = FL_DASH;
	int label_linewidth = 2;
	int label_font = FL_TIMES;
	int label_size = 22;
	FILE* input = fopen(name, "rb");
	if (input == NULL)return;
	while (fgets(line, 1024, input))
	{
		token = strtok(line, " \t\n\r\f");
		if (token != NULL) {
			//main setting	
			if (!strcmp(token, "AERAX:"))
			{
				token = strtok(NULL, " \t\n\r\f");
				xscale = atof(token);
				token = strtok(NULL, " \t\n\r\f");
				wscale = atof(token);
			}
			if (!strcmp(token, "AERAY:"))
			{
				token = strtok(NULL, " \t\n\r\f");
				yscale = atof(token);
				token = strtok(NULL, " \t\n\r\f");
				hscale = atof(token);
			}
			if (!strcmp(token, "TITLEX:"))
			{
				token = strtok(NULL, " \t\n\r\f");
				strcpy(titlex, token);
			}
			if (!strcmp(token, "TITLEY:"))
			{
				token = strtok(NULL, " \t\n\r\f");
				strcpy(titley, token);
			}
			if (!strcmp(token, "TITLE_FONT:"))
			{
				token = strtok(NULL, " \t\n\r\f");
				title_font = atoi(token);
			}
			if (!strcmp(token, "TITLE_SIZE:"))
			{
				token = strtok(NULL, " \t\n\r\f");
				title_size = atoi(token);
			}
			if (!strcmp(token, "SHOW_TITLEX:"))
			{
				token = strtok(NULL, " \t\n\r\f");
				xlabel_show = 1;
			}
			if (!strcmp(token, "SHOW_TITLEY:"))
			{
				token = strtok(NULL, " \t\n\r\f");
				ylabel_show = 1;
			}
			if (!strcmp(token, "MINX:"))
			{
				token = strtok(NULL, " \t\n\r\f");
				xmi = atof(token);
				bound_flag = 1;
			}
			if (!strcmp(token, "MINY:"))
			{
				token = strtok(NULL, " \t\n\r\f");
				ymi = atof(token);
				bound_flag = 1;
			}
			if (!strcmp(token, "MAXX:"))
			{
				token = strtok(NULL, " \t\n\r\f");
				xma = atof(token);
				bound_flag = 1;
			}
			if (!strcmp(token, "MAXY:"))
			{
				token = strtok(NULL, " \t\n\r\f");
				yma = atof(token);
				bound_flag = 1;
			}
			if (!strcmp(token, "AXIS_FONT:"))
			{
				token = strtok(NULL, " \t\n\r\f");
				axis_font = atoi(token);
			}
			if (!strcmp(token, "AXIS_SIZE:"))
			{
				token = strtok(NULL, " \t\n\r\f");
				axis_size = atoi(token);
			}
			if (!strcmp(token, "AXIS_COLOR:"))
			{
				token = strtok(NULL, " \t\n\r\f");
				in_xaxis_color = (Fl_Color)atoi(token);
				token = strtok(NULL, " \t\n\r\f");
				in_yaxis_color = (Fl_Color)atoi(token);
			}
			if (!strcmp(token, "AXIS_WIDTH:"))
			{
				token = strtok(NULL, " \t\n\r\f");
				in_axis_width = atoi(token);
			}
			if (!strcmp(token, "TICKMODE:"))
			{
				token = strtok(NULL, " \t\n\r\f");
				tick_mode = atoi(token);
			}
			if (!strcmp(token, "TICKX:"))
			{
				token = strtok(NULL, " \t\n\r\f");
				xs = atof(token);
			}
			if (!strcmp(token, "TICKY:"))
			{
				token = strtok(NULL, " \t\n\r\f");
				ys = atof(token);
			}
			if (!strcmp(token, "LEGEND_FONT:"))
			{
				token = strtok(NULL, " \t\n\r\f");
				legend_font = atoi(token);
			}
			if (!strcmp(token, "LEGEND_SIZE:"))
			{
				token = strtok(NULL, " \t\n\r\f");
				legend_size = atoi(token);
			}
			if (!strcmp(token, "LEGEND_ALL:"))
			{
				token = strtok(NULL, " \t\n\r\f");
				legend_flag = atoi(token);
			}
			if (!strcmp(token, "LEGEND_POS:"))
			{
				token = strtok(NULL, " \t\n\r\f");
				legend_mode = atoi(token);
			}
			if (!strcmp(token, "LEGEND_MODE:"))
			{
				token = strtok(NULL, " \t\n\r\f");
				legend_background = atoi(token);
			}
			//data setting
			if (!strcmp(token, "NAME:"))
			{
				token = strtok(NULL, " \t\n\r\f");
				strcpy(line_name, token);
			}
			if (!strcmp(token, "COLOR:"))
			{
				token = strtok(NULL, " \t\n\r\f");
				col = (unsigned)atoi(token);
			}
			if (!strcmp(token, "LINESTYLE:"))
			{
				token = strtok(NULL, " \t\n\r\f");
				linestyle = atoi(token);
			}
			if (!strcmp(token, "LINEWIDTH:"))
			{
				token = strtok(NULL, " \t\n\r\f");
				linewidth = atoi(token);
			}
			if (!strcmp(token, "SYMBOL:"))
			{
				token = strtok(NULL, " \t\n\r\f");
				symbol = atoi(token) % FL_SYMBOL_NUM;
			}
			if (!strcmp(token, "SYMBOLSIZE:"))
			{
				token = strtok(NULL, " \t\n\r\f");
				symbol_size = atoi(token);
			}
			if (!strcmp(token, "LEGEND:"))
			{
				token = strtok(NULL, " \t\n\r\f");
				legend_flag = atoi(token);
			}
			if (!strcmp(token, "LINE:"))
			{
				token = strtok(NULL, " \t\n\r\f");
				show_flag = atoi(token);
			}
			if (!strcmp(token, "DATA:"))
			{
				token = strtok(NULL, " \t\n\r\f");
				n = atoi(token);
				if (n > 0)
				{
					x = (double*)malloc(n * sizeof(double));
					y = (double*)malloc(n * sizeof(double));
					for (int i = 0; i < n; i++)
					{
						fgets(line, 1024, input);
						token = strtok(line, " \t\n\r\f");
						x[i] = atof(token);
						token = strtok(NULL, " \t\n\r\f");
						y[i] = atof(token);
					}
					add(n, x, y, line_name, col, linestyle, linewidth, symbol, symbol_size);
					set_name_flag(numb - 1, legend_flag != 0);
					set_show_flag(numb - 1, show_flag != 0);
					free(x);
					free(y);
				}
				sprintf(line_name, "LINE %d", numb);
				col = 0;
				linestyle = 0;
				linewidth = 2;
				symbol = 0;
				symbol_size = 3;
				legend_flag = 1;
				show_flag = 1;
			}

			if (!strcmp(token, "DATA_R:"))
			{
				token = strtok(NULL, " \t\n\r\f");
				n = atoi(token);
				if (n > 0)
				{
					x = (double*)malloc(n * sizeof(double));
					y = (double*)malloc(n * sizeof(double));
					r = (double*)malloc(n * sizeof(double));
					for (int i = 0; i < n; i++)
					{
						fgets(line, 1024, input);
						token = strtok(line, " \t\n\r\f");
						x[i] = atof(token);
						token = strtok(NULL, " \t\n\r\f");
						y[i] = atof(token);
						token = strtok(NULL, " \t\n\r\f");
						r[i] = atof(token);
					}
					add(n, x, y, r, line_name, col, linestyle, linewidth, symbol, symbol_size);
					set_name_flag(numb - 1, legend_flag != 0);
					set_show_flag(numb - 1, show_flag != 0);
					free(x);
					free(y);
					free(r);
				}
				sprintf(line_name, "LINE %d", numb);
				col = 0;
				linestyle = 0;
				linewidth = 2;
				symbol = 0;
				symbol_size = 3;
				legend_flag = 1;
				show_flag = 1;
			}

			//label setting
			if (!strcmp(token, "LABEL:"))
			{
				token = strtok(NULL, " \t\n\r\f");
				label_x = atof(token);
				token = strtok(NULL, " \t\n\r\f");
				strcpy(label_name, token);
				token = strtok(NULL, " \t\n\r\f");
				if (token != NULL)
					label_axis = atoi(token);
				token = strtok(NULL, " \t\n\r\f");
				if (token != NULL)
					label_show_flag = atoi(token);
				token = strtok(NULL, " \t\n\r\f");
				if (token != NULL)
					label_textcol = atoi(token);
				token = strtok(NULL, " \t\n\r\f");
				if (token != NULL)
					label_font = atoi(token);
				token = strtok(NULL, " \t\n\r\f");
				if (token != NULL)
					label_size = atoi(token);
				token = strtok(NULL, " \t\n\r\f");
				if (token != NULL)
					label_dash_flag = atoi(token);
				token = strtok(NULL, " \t\n\r\f");
				if (token != NULL)
					label_col = atoi(token);
				token = strtok(NULL, " \t\n\r\f");
				if (token != NULL)
					label_linestyle = atoi(token);
				token = strtok(NULL, " \t\n\r\f");
				if (token != NULL)
					label_linewidth = atoi(token);


				add_label(label_x, label_name, label_axis, label_textcol, label_font, label_size, label_dash_flag, label_col, label_linestyle, label_linewidth);
				set_label_show_flag(numbl - 1, label_show_flag != 0);
				label_axis = 0;
				label_show_flag = 1;
				label_dash_flag = 1;
				label_col = FL_BLACK;
				label_textcol = FL_BLACK;
				label_linestyle = FL_DASH;
				label_linewidth = 1;
				label_font = FL_TIMES;
				label_size = 22;
			}
		}
	}
	if (bound_flag)
		bound(xmi, ymi, xma, yma);
	else
		auto_bound(0.05);
	titles(titlex, titley);
	space(xs, ys);
	set_axis_color(in_xaxis_color, in_yaxis_color);
	set_axis_width(in_axis_width);
	fclose(input);
}

void Fl_Chart_SVMOD::write(const char name[1024])
{
	FILE* output = fopen(name, "wb");
	if (output == NULL)return;
	fprintf(output, "TITLEX: %s\n", xlabel);
	fprintf(output, "TITLEY: %s\n", ylabel);
	fprintf(output, "TITLE_FONT: %d\n", title_font);
	fprintf(output, "TITLE_SIZE: %d\n", title_size);
	fprintf(output, "AERAX: %lf %lf\n", xscale, wscale);
	fprintf(output, "AERAY: %lf %lf\n", yscale, hscale);
	fprintf(output, "SHOW_TITLEX: %d\n", xlabel_show);
	fprintf(output, "SHOW_TITLEY: %d\n", ylabel_show);
	fprintf(output, "MINX: %lf\nMAXX: %lf\n", xmin, xmax);
	fprintf(output, "MINY: %lf\nMAXY: %lf\n", ymin, ymax);
	fprintf(output, "AXIS_FONT: %d\n", axis_font);
	fprintf(output, "AXIS_SIZE: %d\n", axis_size);
	fprintf(output, "AXIS_COLOR: %d %d\n", xaxis_color, yaxis_color);
	fprintf(output, "AXIS_WIDTH: %d\n", axis_width);
	fprintf(output, "TICKMODE: %d\n", tick_mode);
	fprintf(output, "TICKX: %lf\n", xspace);
	fprintf(output, "TICKY: %lf\n", yspace);
	fprintf(output, "LEGEND_ALL: %d\n", legend_flag);
	fprintf(output, "LEGEND_FONT: %d\n", legend_font);
	fprintf(output, "LEGEND_SIZE: %d\n", legend_size);
	fprintf(output, "LEGEND_POS: %d\n", legend_mode);
	fprintf(output, "LEGEND_MODE: %d\n", legend_background);
	for (int i = 0; i < numb; i++)
	{
		fprintf(output, "NAME: %s\n", entries[i].name);
		fprintf(output, "COLOR: %d\n", entries[i].col);
		fprintf(output, "LINESTYLE: %d\n", entries[i].linestyle);
		fprintf(output, "LINEWIDTH: %d\n", entries[i].linewidth);
		fprintf(output, "SYMBOL: %d\n", entries[i].symbol);
		fprintf(output, "SYMBOLSIZE: %d\n", entries[i].symbol_size);
		fprintf(output, "LEGEND: %d\n", entries[i].name_flag);
		fprintf(output, "LINE: %d\n", entries[i].show_flag);
		if (entries[i].linemode == 0)
		{
			fprintf(output, "DATA: %d\n", entries[i].num);
			for (int j = 0; j < entries[i].num; j++)
				fprintf(output, "%lf\t%lf\n", entries[i].x[j][0], entries[i].x[j][1]);
		}
		else if (entries[i].linemode == 1)
		{
			fprintf(output, "DATA_R: %d\n", entries[i].num);
			for (int j = 0; j < entries[i].num; j++)
				fprintf(output, "%lf\t%lf\t%lf\n", entries[i].x[j][0], entries[i].x[j][1], entries[i].x[j][2]);
		}
	}
	for (int i = 0; i < numbl; i++)
	{
		fprintf(output, "LABEL: %lf %s %d %d %d %d %d %d %d %d %d\n", labelist[i].x, labelist[i].name, labelist[i].yaxis, labelist[i].show_flag,
			labelist[i].text_col, labelist[i].font, labelist[i].size, labelist[i].dash_flag, labelist[i].col, labelist[i].linestyle, labelist[i].linewidth);
	}
	fclose(output);
}
