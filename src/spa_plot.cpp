#include"../include/spa_plot.h"
#include"../include/read_write.h"
using std::string;
using std::min_element;
using std::max_element;
using std::min;
using std::max;

const double min_value = 1000000;
const double max_value = -1000000;
// convert normal data to SPAMD file

vector<Fl_Color> sv_color = { FL_RED, FL_BLUE,FL_GREEN ,FL_YELLOW, FL_MAGENTA, FL_DARK_RED,FL_DARK_BLUE, FL_DARK_GREEN,FL_DARK_YELLOW,FL_DARK_MAGENTA };

string convert_name(const char file[])
{
	string spa_file = "";
	for (int i = 0; i < strlen(file); i++)
	{
		if (file[i] != '.')
			spa_file.push_back(file[i]);
		else
			break;
	}
	spa_file += ".spa";
	return spa_file;
}

void addbenchmarkX(Fl_Chart_SVMOD* chart)
{
	vector<double> y{ 0,0 };
	vector<double> x{ chart->xmin ,chart->xmax };
	chart->add(y.size(), x, y, "NULL", FL_BLACK, FL_DASH, 2, FL_NO_SYMBOL_SVMOD, 3);
	chart->set_name_flag(chart->numb - 1, 0);
}

void addbenchmarkY(Fl_Chart_SVMOD* chart)
{
	vector<double> y{ chart->ymin ,chart->ymax };
	vector<double> x{ 0,0 };
	chart->add(y.size(), x, y, "NULL", FL_BLACK, FL_DASH, 2, FL_NO_SYMBOL_SVMOD, 3);
	chart->set_name_flag(chart->numb - 1, 0);
}

double get_Enthalpy_formation(vector<double> en)
{
	double energy = get_energy();
	FILE* fp = fopen("POSCAR", "r");
	if (fp == NULL)
	{
		printf("POOSCAR IS NOT EXIST!\n");
		return 0;
	}
	POSCAR pos;
	readposcar(fp, pos);
	if (pos.nant[1] != en.size())
	{
		printf("Element type does not match the number of input parameters!\n");
		return 0;
	}
	for (int i = 0; i < pos.nant[1]; i++)
		energy -= pos.typenum[i] * en[i];
	energy /= pos.nant[0];
	return energy;
}

void bound_fix(double& xmin, double& xmax, double& ymin, double& ymax)
{
	double y = ymax - ymin;
	double x = xmax - xmin;
	ymax += y * 0.1;
	ymin -= y * 0.1;
	xmax += x * 0.1;
	xmin -= x * 0.1;
}

void insert_kpath(Fl_Chart_SVMOD* chart)
{
	FILE* fp = fopen("KLABELS", "r");
	if (fp == NULL) return;
	char buf[1024];
	fgets(buf, 1024, fp);
	while (fgets(buf, 1024, fp))
	{
		if (strspn(buf, " \t\n\r") == strlen(buf))
			break;
		char name[1024];
		double x = 0;
		sscanf(buf, "%s%lf", name, &x);
		chart->add_label(x, name);
	}
	fclose(fp);
}

string convert_TDOS(const char file[], int ISPIN)
{
	FILE* fp = fopen(file, "r");
	if (fp == NULL) return "";
	char buf[1024];
	vector<double> Y_up;
	vector<double> Y_dw;
	vector<double> X;
	while (fgets(buf, 1024, fp))
	{
		if (strstr(buf, "#"))
			continue;
		double x = 0, y1 = 0, y2 = 0;
		if (ISPIN == 1)
			sscanf(buf, "%lf%lf", &x, &y1);
		else
			sscanf(buf, "%lf%lf%lf", &x, &y1, &y2);
		X.push_back(x); Y_up.push_back(y1); Y_dw.push_back(y2 * (-1));
	}
	fclose(fp);
	int color_cnt = 0;
	double xmin = *min_element(X.begin(), X.end());
	double ymin = *min_element(Y_up.begin(), Y_up.end());
	double xmax = *max_element(X.begin(), X.end());
	double ymax = *max_element(Y_up.begin(), Y_up.end());
	if (ISPIN == 2)
	{
		ymin = *min_element(Y_dw.begin(), Y_dw.end());
		ymax = max(ymax, fabs(ymin));
		ymin = -ymax;
	}
	bound_fix(xmin, xmax, ymin, ymax);
	double space_x = (int)((xmax - xmin) / 8);
	double space_y = (int)((ymax - ymin) / 5);
	Fl_Chart_SVMOD* chart = new Fl_Chart_SVMOD();
	chart->titles("Energy(ev)", "TDOS");
	chart->space(space_x, space_y);
	chart->bound(xmin, ymin, xmax, ymax);
	if (ISPIN == 1)
		chart->add(X.size(), X, Y_up, "TDOS", sv_color[color_cnt++], FL_SOLID, 2, FL_NO_SYMBOL_SVMOD, 3);
	if (ISPIN == 2)
	{
		chart->add(X.size(), X, Y_up, "TDOS_up", sv_color[color_cnt++], FL_SOLID, 2, FL_NO_SYMBOL_SVMOD, 3);
		chart->add(X.size(), X, Y_dw, "TDOS_dw", sv_color[color_cnt++], FL_SOLID, 2, FL_NO_SYMBOL_SVMOD, 3);
	}
	addbenchmarkY(chart);
	chart->write(convert_name(file).c_str());
	delete chart;
	return convert_name(file);
}

string convert_ADOS(const char file[])
{
	FILE* fp = fopen(file, "r");
	if (fp == NULL) return "";
	char buf[1024];
	string x_label;
	vector<string> y_label;
	vector<double> X;
	vector<vector<double> > Y;
	while (fgets(buf, 1024, fp))
	{
		buf[strlen(buf) - 1] = '\0';
		if (strstr(buf, "#"))
		{
			char* p = strtok(buf, " \t\n\r\f");
			x_label = p;
			x_label.erase(x_label.begin());
			while (p = strtok(NULL, " \t\n\r\f"))
				y_label.push_back(p);
			Y.resize(y_label.size());
			continue;
		}
		char* p = strtok(buf, " \t\n\r\f");
		X.push_back(atof(p));
		int cnt = 0;
		while (p = strtok(NULL, " \t\n\r\f"))
			Y[cnt++].push_back(atof(p));
	}
	fclose(fp);
	Fl_Chart_SVMOD* chart = new Fl_Chart_SVMOD();
	double space_x = 0, space_y = 0, xmin = 0, ymin = min_value, xmax = 0, ymax = max_value;
	xmin = *min_element(X.begin(), X.end());
	xmax = *max_element(X.begin(), X.end());
	int color_cnt = 0;
	for (int i = 0; i < Y.size(); i++)
	{
		ymin = min(ymin, *min_element(Y[i].begin(), Y[i].end()));
		ymax = max(ymax, *max_element(Y[i].begin(), Y[i].end()));
		chart->add(X.size(), X, Y[i], y_label[i].c_str(), sv_color[color_cnt++], FL_SOLID, 2, FL_NO_SYMBOL_SVMOD, 3);
	}
	bound_fix(xmin, xmax, ymin, ymax);
	space_x = (int)((xmax - xmin) / 8);
	space_y = (int)((ymax - ymin) / 5);
	chart->titles(x_label.c_str(), "PDOS");
	chart->space(space_x, space_y);
	chart->bound(xmin, ymin, xmax, ymax);
	addbenchmarkY(chart);
	chart->write(convert_name(file).c_str());
	delete chart;
	return convert_name(file);
}

string convert_ADOS(const char file_up[], const char file_dw[])
{
	vector<string> file = { file_up,file_dw };
	string x_label;
	vector<vector<string> >y_label(2);
	vector<double> X;
	vector<vector<vector<double> > > Y(2);
	for (int i = 0; i < 2; i++)
	{
		FILE* fp = fopen(file[i].c_str(), "r");
		if (fp == NULL)return "";
		char buf[1024];
		while (fgets(buf, 1024, fp))
		{
			buf[strlen(buf) - 1] = '\0';
			if (strstr(buf, "#"))
			{
				char* p = strtok(buf, " \t\n\r\f");
				x_label = p;
				x_label.erase(x_label.begin());
				while (p = strtok(NULL, " \t\n\r\f"))
				{
					string tmp = p;
					if (i == 0)
						tmp += "_up";
					else
						tmp += "_dw";
					y_label[i].push_back(tmp);
				}
				Y[i].resize(y_label[i].size());
				continue;
			}
			char* p = strtok(buf, " \t\n\r\f");
			if (i == 0)X.push_back(atof(p));
			int cnt = 0;
			while (p = strtok(NULL, " \t\n\r\f"))
				Y[i][cnt++].push_back((i == 0) ? atof(p) : -atof(p));
		}
		fclose(fp);
	}
	Fl_Chart_SVMOD* chart = new Fl_Chart_SVMOD();
	double space_x = 0, space_y = 0, xmin = 0, ymin = min_value, xmax = 0, ymax = max_value;
	xmin = *min_element(X.begin(), X.end());
	xmax = *max_element(X.begin(), X.end());
	int color_cnt = 0;
	for (int k = 0; k < 2; k++)
	{
		for (int i = 0; i < Y[k].size(); i++)
		{
			ymin = min(ymin, *min_element(Y[k][i].begin(), Y[k][i].end()));
			ymax = max(ymax, *max_element(Y[k][i].begin(), Y[k][i].end()));
			chart->add(X.size(), X, Y[k][i], y_label[k][i].c_str(), sv_color[color_cnt % Y[k].size()], FL_SOLID, 2, FL_NO_SYMBOL_SVMOD, 3);
			color_cnt++;
		}
	}
	//for (int i = 1; i < chart->numb; i++)
	//	chart->set_show_flag(i, 0);
	ymax = max(ymax, fabs(ymin));
	ymin = -ymax;
	bound_fix(xmin, xmax, ymin, ymax);
	space_x = (int)((xmax - xmin) / 8);
	space_y = (int)((ymax - ymin) / 5);
	chart->titles(x_label.c_str(), "PDOS");
	chart->space(space_x, space_y);
	chart->bound(xmin, ymin, xmax, ymax);
	string spa_name;
	for (int i = 0; i < file[0].size(); i++)
	{
		if (file[0][i] == file[1][i])
			spa_name.push_back(file[0][i]);
		else
			break;
	}
	spa_name.pop_back();
	spa_name += ".spa";
	addbenchmarkY(chart);
	chart->write(spa_name.c_str());
	delete chart;
	return spa_name;
}

string convert_basic_band(const char file[], int ISPIN)
{
	FILE* fp = fopen(file, "r");
	if (fp == NULL) return "";
	char buf[1024];
	fgets(buf, 1024, fp);
	fgets(buf, 1024, fp);
	int nkpt = 0, nband = 0;
	sscanf(buf, "%*[^:]%*[:]%d%d", &nkpt, &nband);
	vector<vector<double> > Y_up(nband, vector<double>(nkpt));
	vector<vector<double> > Y_dw(nband, vector<double>(nkpt));
	vector<vector<double> > X(nband, vector<double>(nkpt));
	int ikpt = -1, iband = -1;
	while (fgets(buf, 1024, fp))
	{
		if (strspn(buf, " \t\n\r") == strlen(buf))
			continue;
		if (strstr(buf, "#"))
		{
			ikpt = -1;
			iband++;
			continue;
		}
		ikpt++;
		double x = 0, y1 = 0, y2 = 0;
		if (ISPIN == 1)
			sscanf(buf, "%lf%lf", &x, &y1);
		else
			sscanf(buf, "%lf%lf%lf", &x, &y1, &y2);
		X[iband][ikpt] = x; Y_up[iband][ikpt] = y1; Y_dw[iband][ikpt] = y2;
	}
	fclose(fp);
	double xmin = *min_element(X[0].begin(), X[0].end());
	double xmax = *max_element(X[0].begin(), X[0].end());
	double ymin = min_value, ymax = max_value;
	for (int i = 0; i < iband; i++)
	{
		ymin = min(ymin, *min_element(Y_up[i].begin(), Y_up[i].end()));
		ymax = max(ymax, *max_element(Y_up[i].begin(), Y_up[i].end()));
	}
	Fl_Chart_SVMOD* chart = new Fl_Chart_SVMOD();
	chart->titles("K-Path", "Energy(ev)");
	bound_fix(xmin, xmax, ymin, ymax);
	double space_x = (int)((xmax - xmin) / 8);
	double space_y = (int)((ymax - ymin) / 5);
	chart->space(space_x, space_y);
	chart->bound(xmin, ymin, xmax, ymax);
	for (int i = 0; i < iband; i++)
	{
		if (ISPIN == 1)
		{
			string name = std::to_string(i);
			chart->add(X[i].size(), X[i], Y_up[i], name.c_str(), sv_color[0], FL_SOLID, 2, FL_NO_SYMBOL_SVMOD, 3);
		}
		if (ISPIN == 2)
		{
			string name_up = "up_" + std::to_string(i);
			string name_dw = "dw_" + std::to_string(i);
			chart->add(X[i].size(), X[i], Y_up[i], name_up.c_str(), sv_color[0], FL_SOLID, 2, FL_NO_SYMBOL_SVMOD, 3);
			chart->add(X[i].size(), X[i], Y_dw[i], name_dw.c_str(), sv_color[1], FL_DASH, 2, FL_NO_SYMBOL_SVMOD, 3);
		}
	}
	for (int i = 0; i < chart->numb; i++)
		chart->set_name_flag(i, 0);
	insert_kpath(chart);
	addbenchmarkX(chart);
	chart->write(convert_name(file).c_str());
	delete chart;
	return convert_name(file);
}

string convert_fat_band(const char file[])
{
	FILE* fp = fopen(file, "r");
	if (fp == NULL) return "";
	char buf[1024];
	fgets(buf, 1024, fp);
	char* p = strtok(buf, " \t\n\r\f");
	string x_label = p;
	string y_label;
	vector<string> z_label;
	x_label.erase(x_label.begin());
	int cnt = 0;
	while (p = strtok(NULL, " \t\n\r\f"))
	{
		if (cnt == 0)
			y_label = p;
		else
			z_label.push_back(p);
		cnt++;
	}
	fgets(buf, 1024, fp);
	int nkpt = 0, nband = 0;
	sscanf(buf, "%*[^:]%*[:]%d%d", &nkpt, &nband);
	int ikpt = -1, iband = -1;
	vector<vector<double> > X(nband, vector<double>(nkpt)), Y(nband, vector<double>(nkpt));
	vector<vector<vector<double> > > Z(z_label.size(), vector<vector<double> >(nband, vector<double>(nkpt)));
	while (fgets(buf, 1024, fp))
	{
		if (strspn(buf, " \t\n\r") == strlen(buf))
			continue;
		if (strstr(buf, "#"))
		{
			ikpt = -1;
			iband++;
			continue;
		}
		cnt = 0;
		ikpt++;
		p = strtok(buf, " \t\n\r\f");
		X[iband][ikpt] = atof(p);
		while (p = strtok(NULL, " \t\n\r\f"))
		{
			if (cnt == 0)
				Y[iband][ikpt] = atof(p);
			else
				Z[cnt - 1][iband][ikpt] = atof(p);
			cnt++;
		}
	}
	fclose(fp);
	Fl_Chart_SVMOD* chart = new Fl_Chart_SVMOD();
	double xmin = *min_element(X[0].begin(), X[0].end());
	double xmax = *max_element(X[0].begin(), X[0].end());
	double ymin = min_value, ymax = max_value;
	for (int i = 0; i < iband; i++)
	{
		ymin = min(ymin, *min_element(Y[i].begin(), Y[i].end()));
		ymax = max(ymax, *max_element(Y[i].begin(), Y[i].end()));
	}
	chart->titles(x_label.c_str(), y_label.c_str());
	bound_fix(xmin, xmax, ymin, ymax);
	double space_x = (int)((xmax - xmin) / 8);
	double space_y = (int)((ymax - ymin) / 5);
	chart->space(space_x, space_y);
	chart->bound(xmin, ymin, xmax, ymax);
	for (int i = 0; i < iband; i++)
	{
		int color_cnt = 0;
		for (int j = 0; j < Z.size(); j++)
			chart->add(X[i].size(), X[i], Y[i], Z[j][i], z_label[j].c_str(), sv_color[color_cnt++], FL_SOLID, 2, FL_CIRCLE_SVMOD, 3);
	}
	insert_kpath(chart);
	for (int i = 0; i < chart->numb; i++)
		chart->set_name_flag(i, 0);
	for (int i = 0; i < chart->numb; i++)
	{
		if ((i + 1) % z_label.size() != 0)
			chart->set_show_flag(i, 0);
	}
	addbenchmarkX(chart);
	chart->write(convert_name(file).c_str());
	delete chart;
	return convert_name(file);
}

void hull_point::spa_convexhull(const char data[])
{
	FILE* fp = fopen(data, "r");
	if (fp == NULL)
	{
		printf("%s IS NOT EXIST!\n", data);
		return;
	}
	vector<double> X, Y;
	vector<double> hullx{ 0,1 }, hully{ 0,0 };
	char buf[1024];
	while (fgets(buf, 1024, fp))
	{
		double x = 0, y = 0;
		sscanf(buf, "%lf%lf", &x, &y);
		if (find(hullx.begin(), hullx.end(), x) == hullx.end())
		{
			hullx.push_back(x);
			hully.push_back(y);
		}
		else
		{
			if (y < hully[find(hullx.begin(), hullx.end(), x) - hullx.begin()])
				hully[find(hullx.begin(), hullx.end(), x) - hullx.begin()] = y;
		}
		X.push_back(x);
		Y.push_back(y);
	}
	fclose(fp);
	vector<hull_point> hull(hullx.size());
	for (int i = 0; i < hullx.size(); i++)
	{
		hull[i].x = hullx[i];
		hull[i].y = hully[i];
	}
	std::sort(hull.begin(), hull.end(), mysort);
	int flag = 1;
	while (flag)
	{
		flag = 0;
		for (auto iter = hull.begin() + 1; iter != hull.end() - 1;)
		{
			if ((((*(iter + 1)).y - (*iter).y) / ((*(iter + 1)).x - (*iter).x)) >= ((*(iter)).y - (*(iter - 1)).y) / ((*(iter)).x - (*(iter - 1)).x))
				iter++;
			else
			{
				iter = hull.erase(iter);
				flag = 1;
			}
		}
	}
	for (int i = 0; i < hull.size(); i++)
	{
		hullx[i] = hull[i].x;
		hully[i] = hull[i].y;
	}
	Fl_Chart_SVMOD* chart = new Fl_Chart_SVMOD();
	chart->titles("Compositon: N fraction", "Formation Enthalpy (eV/atom)");
	double ymax = max(0.0, *max_element(Y.begin(), Y.end()));
	double ymin = min(0.0, *min_element(Y.begin(), Y.end()));
	double delta_y = ymax - ymin;
	chart->space(0.2, 0.5);
	chart->bound(0, ymin - 0.1 * delta_y, 1, ymax + 0.1 * delta_y);
	chart->add(X.size(), X, Y, "NULL", FL_RED, FL_NO_LINE, 2, FL_CIRCLE_EMPTY_SVMOD, 3);
	chart->add(hull.size(), hullx, hully, "NULL", FL_RED, FL_DASH, 2, FL_NO_SYMBOL_SVMOD, 3);
	vector<double> y1 = { 0,0 };
	vector<double> x1 = { 0,1 };
	chart->add(y1.size(), x1, y1, "NULL", FL_BLACK, FL_DASH, 2, FL_NO_SYMBOL_SVMOD, 3);
	for (int i = 0; i < chart->numb; i++)
		chart->set_name_flag(i, 0);
	chart->write("convex_hull.txt");
}