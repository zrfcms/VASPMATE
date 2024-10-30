#include"../../include/VASPMATE_include/tools.h"

string GetInfoINCAR(char key[])
{
	FILE* fp = fopen("INCAR", "r");
	string value;
	if (fp == NULL)
	{
		perror("INCAR IS NOT EXIST!\n");
		exit(1);
	}
	char buf[1024];
	int iat = 0;
	while (fgets(buf, 1024, fp) != NULL)
	{
		if (strstr(buf, key) != NULL)
		{
			for (int i = 0; i < strlen(buf); i++)
			{
				if (buf[i] == '#')
					break;
				if (buf[i] == '=')
					iat = 1;
				if (iat)
				{
					if (buf[i] != '=' && buf[i] != ' ' && buf[i] != '\n' && buf[i] != '\r' && buf[i] != '\t')
						value.push_back(buf[i]);
				}
			}
		}
	}
	fclose(fp);
	return value;
}

string GetInfoOUTCAR(char key[])
{
	FILE* fp = fopen("OUTCAR", "r");
	string value;
	if (fp == NULL)
	{
		perror("OUTCAR IS NOT EXIST!\n");
		exit(1);
	}
	char buf[1024];
	int iat = 0;
	while (fgets(buf, 1024, fp) != NULL)
	{
		if (strstr(buf, key) != NULL)
		{
			for (int i = 0; i < strlen(buf); i++)
			{
				if (buf[i] == '#')
					break;
				if (buf[i] == '=')
					iat = 1;
				if (iat)
				{
					if (buf[i] != '=' && buf[i] != ' ' && buf[i] != '\n' && buf[i] != '\r' && buf[i] != '\t')
						value.push_back(buf[i]);
				}
			}
		}
	}
	fclose(fp);
	return value;
}

vector<string> GetInfoINCAR(char key[], char vector_type[])
{
	FILE* fp = fopen("INCAR", "r");
	vector<string> value_key; //return
	string value_line;
	if (fp == NULL)
	{
		perror("INCAR IS NOT EXIST!\n");
		exit(1);
	}
	char buf[1024];
	int iat = 0;
	while (fgets(buf, 1024, fp) != NULL)
	{
		if (strstr(buf, key) != NULL)
		{
			for (int i = 0; i < strlen(buf); i++)
			{
				if (buf[i] == '#')
					break;
				if (buf[i] == '=')
					iat = 1;
				if (iat)
				{
					if (buf[i] != '=')
						value_line.push_back(buf[i]);
				}
			}
		}
	}
	fclose(fp);
	size_t pos = 0;
    if (pos != string::npos) 
	{
        string part = value_line.substr(pos);
        size_t start = part.find_first_not_of(" ");
        if (start != string::npos)
            part = part.substr(start);
        size_t begin = 0;
        while (begin < part.length()) 
		{
            while (begin < part.length() && std::isspace(part[begin])) 
                begin++;
            if (begin >= part.length()) 
				break;
            size_t end = begin;
            while (end < part.length() && !std::isspace(part[end]))
                end++;
            string value_temp = part.substr(begin, end - begin);
            size_t starPos = value_temp.find('*');
            if (starPos != string::npos) 
			{
                int count = stoi(value_temp.substr(0, starPos));
                string value = value_temp.substr(starPos + 1);
                for (int i = 0; i < count; i++)
                    value_key.push_back(value);
            } 
			else
                value_key.push_back(value_temp);
            begin = end;
        }
    }
	return value_key;
}


int IsNum(string symbol)
{
	for (int i = 0; i < symbol.length(); i++)
	{
		int tmp = (int)symbol[i];
		if (tmp >= 48 && tmp <= 57)
		{
			continue;
		}
		else
		{
			return false;
		}
	}
	return true;
}

vector<double> AddVec(vector<double> v1, vector<double> v2)
{
	vector<double> ret;
	if (v1.size() != v2.size())
		return ret;
	ret.resize(v1.size());
	for (int i = 0; i < v1.size(); i++)
		ret[i] = v1[i] + v2[i];
	return ret;
}

vector<double> operator-(vector<double> v1, double value)
{
	vector<double> ret;
	ret.resize(v1.size());
	for (int i = 0; i < v1.size(); i++)
		ret[i] = v1[i] - value;
	return ret;
}

vector<double> operator+(vector<double> v1, vector<double> v2)
{
	vector<double> ret;
	if (v1.size() != v2.size())
		return ret;
	ret.resize(v1.size());
	for (int i = 0; i < v1.size(); i++)
		ret[i] = v1[i] + v2[i];
	return ret;
}

vector<double> operator*(vector<double> v1, vector<double> v2)
{
	vector<double> ret;
	if (v1.size() != v2.size())
		return ret;
	ret.resize(v1.size());
	for (int i = 0; i < v1.size(); i++)
		ret[i] = v1[i] * v2[i];
	return ret;
}

vector<double> Intergral(vector<double> v1, vector<double> v2)
{
	vector<double> v3;
	if (v1.size() != v2.size())
		return v3;
	v3.push_back(0);
	for (int i = 1; i < v1.size(); i++)
		v3.push_back((v1[i] - v1[i - 1]) * (v2[i] + v2[i - 1]) / 2 + v3[i - 1]);
	return v3;
}

double SumVec(vector<double> v)
{
	double sum = 0;
	for (int i = 0; i < v.size(); i++)
		sum += v[i];
	return sum;
}

double volume(double a[3][3])
{
	double odd = 0;
	double even = 0;
	for (int i = 0; i < 3; i++)
		odd += a[0][i] * a[1][(i + 1) % 3] * a[2][(i + 2) % 3];
	for (int i = 0; i < 3; i++)
		even += a[0][i] * a[1][(i + 2) % 3] * a[2][(i + 1) % 3];
	return fabs(odd - even);
}

double* getangle( double vec[3][3])
{
	double* angle = new double[3];
	angle[0] = (180 / Pi) * acos(( vec[1][0] * vec[2][0] + vec[1][1] * vec[2][1] + vec[1][2] * vec[2][2])/
									sqrt( vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2])/
									sqrt( vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]));
	angle[1] = (180 / Pi) * acos(( vec[0][0] * vec[2][0] + vec[0][1] * vec[2][1] + vec[0][2] * vec[2][2])/
									sqrt( vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2])/
									sqrt( vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]));
	angle[2] = (180 / Pi) * acos(( vec[0][0] * vec[1][0] + vec[0][1] * vec[1][1] + vec[0][2] * vec[1][2])/
									sqrt( vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2])/
									sqrt( vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]));
	return angle;
}

double* getcevec( double vec[3][3]) //cell vector
{
	double* cevec = new double[3];
	cevec[0] = sqrt( vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]); //a
	cevec[1] = sqrt( vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]); //b
	cevec[2] = sqrt( vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]); //c
	return cevec;
}

void RecMat(double POSI[3][3], double ret[3][3])
{
	double v = volume(POSI);
	if (v == 0)
		v == 1;
	ret[0][0] = (POSI[1][1] * POSI[2][2] - POSI[1][2] * POSI[2][1]) * (2 * Pi) / v;
	ret[0][1] = -(POSI[1][0] * POSI[2][2] - POSI[1][2] * POSI[2][0]) * (2 * Pi) / v;
	ret[0][2] = (POSI[1][0] * POSI[2][1] - POSI[1][1] * POSI[2][0]) * (2 * Pi) / v;
	ret[1][0] = -(POSI[0][1] * POSI[2][2] - POSI[0][2] * POSI[2][1]) * (2 * Pi) / v;
	ret[1][1] = (POSI[0][0] * POSI[2][2] - POSI[0][2] * POSI[2][0]) * (2 * Pi) / v;
	ret[1][2] = -(POSI[0][0] * POSI[2][1] - POSI[0][1] * POSI[2][0]) * (2 * Pi) / v;
	ret[2][0] = (POSI[0][1] * POSI[1][2] - POSI[1][1] * POSI[0][2]) * (2 * Pi) / v;
	ret[2][1] = -(POSI[0][0] * POSI[1][2] - POSI[1][0] * POSI[0][2]) * (2 * Pi) / v;
	ret[2][2] = (POSI[1][1] * POSI[0][0] - POSI[1][0] * POSI[0][1]) * (2 * Pi) / v;
}

void copy_vec3d(double vec1[3][3], double vec2[3][3])
{
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			vec2[i][j] = vec1[i][j];
}

void multi_mat(double mat1[3][3], double mat2[3][3], double multi_mat[3][3])
{
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			multi_mat[i][j] = mat1[i][0] * mat2[0][j] + mat1[i][1] * mat2[1][j] + mat1[i][2] * mat2[2][j];
		}
	}
}

bool Isequal(double a, double b, double SYMPREC)
{
	return fabs(a - b) > SYMPREC ? 0 : 1;
}

vector<int> ExtractNumbersFromString(char a[], int len)
{
	int i, j, count = 0, wei[20], times = 0;
	vector<int> num(len);
	bool ctoi = 0, befctoi = 0;
	for (i = 0; i < len + 1; i++)
	{
		if (a[i] >= '0' && a[i] <= '9')
		{
			ctoi = 1;
		}
		else
		{
			ctoi = 0;
		}
		if (befctoi == 0 && ctoi == 1)
		{
			wei[count] = a[i] - '0';
			befctoi = 1;
			count++;
		}
		else if (befctoi == 1 && ctoi == 1)
		{
			wei[count] = a[i] - '0';
			count++;
		}
		else if (befctoi == 1 && ctoi == 0)
		{
			for (j = 0; j < count; j++)
			{
				num[times] += wei[j] * pow(10, count - j - 1);
			}
			times++;
			befctoi = 0;
			count = 0;
		}
	}
	return num;
}

double maxvalue(vector<double> v)
{
	double ret = v[0];
	for (int i = 1; i < v.size(); i++)
		ret = max(ret, v[i]);
	return ret;
}

double minvalue(vector<double> v)
{
	double ret = v[0];
	for (int i = 1; i < v.size(); i++)
		ret = min(ret, v[i]);
	return ret;
}

double maxvalue(vector<vector<double> >v)
{
	double ret = v[0][0];
	for (int i = 0; i < v.size(); i++)
		for (int j = 0; j < v[i].size(); j++)
			ret = max(ret, v[i][j]);
	return ret;
}

double minvalue(vector<vector<double> >v)
{
	double ret = v[0][0];
	for (int i = 0; i < v.size(); i++)
		for (int j = 0; j < v[i].size(); j++)
			ret = min(ret, v[i][j]);
	return ret;
}

int zdgys(int a, int b) {
	for (int i = a > b ? a : b; i >= 1; i--) {
		if (a % i == 0 && b % i == 0) return i;
	}
}

string TranDecimalTofraction(string x) {
	int len = x.length(), y, z,
		zs = 0, xs = 0;
	for (int i = 0; i < len; i++) {
		if (x[i] == '.') {
			z = i - 1;
			y = i + 1;
			break;
		}
	}
	for (int i = 0; i <= z; i++) {
		zs += pow(10, z - i) * (x[i] - '0');
	}
	for (int i = y; i < len; i++) {
		xs += pow(10, len - i - 1) * (x[i] - '0');
	}
	int b = pow(10, len - y), a = xs + zs * b, gys = zdgys(a, b);
	string ans = to_string(a / gys) + "/" + to_string(b / gys);
	return ans;
}

string RemovePostfixNum(string elem)
{
	int length = elem.size();
	for (int i = 0; i < elem.length(); i++)
	{
		int tmp = (int)elem[i];
		if (tmp >= 48 && tmp <= 57)
		{
			length = i;
			break;
		}
	}
	string ret(elem.begin(), elem.begin() + length);
	return ret;
}

vector<int> closest::getclosest(vector<pair<int, int> > v, int target)
{
	globaldelta = target;
	backtrack(v, 0, 0, target, target);
	return result;
}

bool mysort(pair<int, int> a, pair<int, int> b)
{
	return a.second < b.second;
}
void closest::backtrack(vector<pair<int, int> > v, int start, int sum, int delta, int target)
{
	set<int> record;
	if (delta < globaldelta) {
		globaldelta = delta;
		result = path;
	}
	if (delta == 0)
		return;
	int limit = INT_MAX;
	for (int i = start; i < v.size(); i++) {
		if (v[i].second >= limit)
			continue;
		if (find(record.begin(), record.end(), v[i].second) != record.end())
			continue;
		record.insert(v[i].second);
		sum += v[i].second;
		int newdelta = abs(sum - target);
		if (newdelta > delta) {
			limit = v[i].second;
			sum -= v[i].second;
			continue;
		}
		path.push_back(i);
		backtrack(v, i + 1, sum, newdelta, target);
		sum -= v[i].second;
		path.pop_back();
	}
}

map<char, int> out = { {'#', 0}, {'(', 8}, {'^', 6}, {'*', 4}, {'/', 4}, {'%', 4},
	{'+', 2},{'-', 2},{')', 1} };
map<char, int> inside = { {'#', 0}, {'(', 1}, {'^', 7}, {'*', 5}, {'/', 5}, {'%', 5},
	{'+', 3},{'-', 3},{')', 8} };


vector<Node> midToBack(string mid) {
	vector<Node> back;
	if (mid.size() == 0) return back;

	string s;
	for (char c : mid) {
		if (c != ' ')
			s.push_back(c);
	}

	if (s[0] == '-') {
		s.insert(0, "0");
	}
	for (int i = 0; i < s.size() - 1; i++) {
		if (s[i] == '(' && s[i + 1] == '-')
			s.insert(i + 1, "0");
	}
	mid = s;

	stack<char> op;
	op.push('#');
	bool flag = false;
	for (int i = 0; i < mid.size(); i++) {
		if (mid[i] >= '0' && mid[i] <= '9') {
			if (flag == false) {
				flag = true;
				int tmp = mid[i] - '0';
				back.push_back(Node(tmp, false));
			}
			else {
				back.back().num *= 10;
				back.back().num += mid[i] - '0';
			}
		}
		else if (out[mid[i]] > inside[op.top()]) {
			flag = false;
			op.push(mid[i]);
		}
		else if (out[mid[i]] < inside[op.top()]) {
			flag = false;
			back.push_back(Node(op.top(), true));
			op.pop();
			i--;
		}
		else {
			flag = false;
			op.pop();
		}
	}
	while (op.top() != '#') {
		back.push_back(Node(op.top(), true));
		op.pop();
	}
	return back;
}

double getValue(vector<Node> back) {
	if (back.size() == 0) return 0;
	stack<double> nums;
	for (auto e : back) {
		if (!e.isop) {
			nums.push(e.num);
		}
		else {
			double b = nums.top();
			nums.pop();
			double a = nums.top();
			nums.pop();
			double res;
			switch (e.op)
			{
			case '+':
				res = a + b; break;
			case '-':
				res = a - b; break;
			case '*':
				res = a * b; break;
			case '/':
				res = a / b; break;
			case '^':
				res = pow(a, b); break;
			case '%':
				res = int(a) % int(b); break;
			default:
				break;
			}
			nums.push(res);
		}
	}
	return nums.top();
}

void copy_file(const  char* target, const  char* source)
{
	FILE* fp1 = fopen(source, "r");
	if (fp1 == NULL)
	{
		printf("%s IS NOT EXIST!\n", source);
		return;
	}
	FILE* fp2 = fopen(target, "w");
	char buff[1024];
	while (fgets(buff, 1024, fp1) != NULL)
		fputs(buff, fp2);
	fclose(fp1);
	fclose(fp2);
}

int len(char buf[])
{
	int flag = 0;
	int cnt = 0;
	for (int i = 0; i < strlen(buf); i++)
	{
		if (buf[i] == '\n' || buf[i] == '\r')
			break;
		if (!flag && (buf[i] == ' ' || buf[i] == '\t'))
			continue;
		if (!flag && buf[i] != ' ' && buf[i] != '\t')
		{
			cnt++;
			flag = 1;
			continue;
		}
		if (flag && (buf[i] == ' ' || buf[i] == '\t'))
		{
			flag = 0;
			continue;
		}
	}
	return cnt;
}

void getlength(double vec[3][3], double length[3])
{
	length[0] = length[1] = length[2] = 0;
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			length[i] += vec[i][j] * vec[i][j];
		}
		length[i] = sqrt(length[i]);
	}
}

void brinv(double vec[3][3], double inv_vec[3][3])
{
	double v = volume(vec);
	if (v == 0)
		v == 1;
	inv_vec[0][0] = (vec[1][1] * vec[2][2] - vec[1][2] * vec[2][1]) / v;
	inv_vec[0][1] = -(vec[1][0] * vec[2][2] - vec[1][2] * vec[2][0]) / v;
	inv_vec[0][2] = (vec[1][0] * vec[2][1] - vec[1][1] * vec[2][0]) / v;
	inv_vec[1][0] = -(vec[0][1] * vec[2][2] - vec[0][2] * vec[2][1]) / v;
	inv_vec[1][1] = (vec[0][0] * vec[2][2] - vec[0][2] * vec[2][0]) / v;
	inv_vec[1][2] = -(vec[0][0] * vec[2][1] - vec[0][1] * vec[2][0]) / v;
	inv_vec[2][0] = (vec[0][1] * vec[1][2] - vec[1][1] * vec[0][2]) / v;
	inv_vec[2][1] = -(vec[0][0] * vec[1][2] - vec[1][0] * vec[0][2]) / v;
	inv_vec[2][2] = (vec[1][1] * vec[0][0] - vec[1][0] * vec[0][1]) / v;
}

void transpose_matrix(vector<vector<double> > v, vector<vector<double> >& tran_v)
{
	for (int i = 0; i < tran_v.size(); i++)
	{
		for (int j = 0; j < tran_v[i].size(); j++)
		{
			tran_v[i][j] = v[j][i];
		}
	}
}

vector<vector<double>> mutil(vector<vector<double>> m1, vector<vector<double>> m2)
{
	int m = m1.size();
	int n = m1[0].size();
	int p = m2[0].size();
	vector<vector<double>> arr;
	vector<double> temparay;
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < p; j++) {
			double sum = 0;
			for (int k = 0; k < n; k++) {
				sum += m1[i][k] * m2[k][j];
			}
			temparay.push_back(sum);
		}
		arr.push_back(temparay);
		temparay.erase(temparay.begin(), temparay.end());
	}
	return arr;
}

void cross(double a[3], double b[3], double c[3])
{
	c[0] = a[1] * b[2] - a[2] * b[1];
	c[1] = a[2] * b[0] - a[0] * b[2];
	c[2] = a[0] * b[1] - a[1] * b[0];
}

void ITP_french(int init_size[3], vector<double> init_pot, int out_size[3], vector<double> &out_pot, int scale)
{
	if (init_pot.size() != init_size[0] * init_size[1] * init_size[2])
		return;
	out_pot.resize(scale * scale * scale * init_size[0] * init_size[1] * init_size[2]);
	double step1[4][4];
	double step2[4];
	double par1[3][4];
	double par0[3];
	double iscale = 1.0 / (double)(scale * scale * scale);
	for (int i = 0; i < 3; i++)
		out_size[i] = init_size[i] * scale;

	for (double ii = 0; ii < scale; ii++)
		for (double jj = 0; jj < scale; jj++)
			for (double kk = 0; kk < scale; kk++)
			{
				par0[0] = ii;
				par0[1] = jj;
				par0[2] = kk;
				for (int i = 0; i < 3; i++)
				{
					par1[i][0] = -0.5 * par0[i] * (scale - par0[i]) * (scale - par0[i]) * iscale;
					par1[i][1] = ((scale - par0[i]) * (scale - par0[i]) + 3 * par0[i] * (scale - par0[i]) + 0.5 * par0[i] * par0[i]) * (scale - par0[i]) * iscale;
					par1[i][2] = (0.5 * (scale - par0[i]) * (scale - par0[i]) + 3 * par0[i] * (scale - par0[i]) + par0[i] * par0[i]) * par0[i] * iscale;
					par1[i][3] = -0.5 * par0[i] * par0[i] * (scale - par0[i]) * iscale;
				}
				for (int i = 0; i < init_size[0]; i++)
					for (int j = 0; j < init_size[1]; j++)
						for (int k = 0; k < init_size[2]; k++)
						{
							for (int o = 0; o < 4; o++)
							{
								step2[o] = 0;
								for (int p = 0; p < 4; p++)
									step1[o][p] = 0;
							}
							for (int o = 0; o < 4; o++)
								for (int p = 0; p < 4; p++)
									for (int q = 0; q < 4; q++)
										step1[o][p] += par1[2][q] * init_pot[((k - 1 + q + init_size[2]) % init_size[2]) * init_size[1] * init_size[0]
										+ ((j - 1 + o + init_size[1]) % init_size[1]) * init_size[0]
										+ (i - 1 + p + init_size[0]) % init_size[0]];
							for (int o = 0; o < 4; o++)
								for (int p = 0; p < 4; p++)
									step2[p] += par1[1][o] * step1[o][p];
							for (int p = 0; p < 4; p++)
								out_pot[(scale * k + (int)kk) * out_size[0] * out_size[1] + (scale * j + (int)jj) * out_size[0] + (scale * i + (int)ii)] += par1[0][p] * step2[p];
						}
			}
}

Eigen::VectorXd FitterLeastSquareMethod(vector<double>& X, vector<double>& Y, uint8_t orders)
{
	// abnormal input verification
	if (X.size() < 2 || Y.size() < 2 || X.size() != Y.size() || orders < 1)
		return Eigen::VectorXd{};

	// map sample data from STL vector to eigen vector
	Eigen::Map<Eigen::VectorXd> sampleX(X.data(), X.size());
	Eigen::Map<Eigen::VectorXd> sampleY(Y.data(), Y.size());

	Eigen::MatrixXd mtxVandermonde(X.size(), orders + 1);  // Vandermonde matrix of X-axis coordinate vector of sample data
	Eigen::VectorXd colVandermonde = sampleX;              // Vandermonde column

	// construct Vandermonde matrix column by column
	for (size_t i = 0; i < orders + 1; ++i)
	{
		if (0 == i)
		{
			mtxVandermonde.col(0) = Eigen::VectorXd::Constant(X.size(), 1, 1);
			continue;
		}
		if (1 == i)
		{
			mtxVandermonde.col(1) = colVandermonde;
			continue;
		}
		colVandermonde = colVandermonde.array() * sampleX.array();
		mtxVandermonde.col(i) = colVandermonde;
	}

	// calculate coefficients vector of fitted polynomial
	Eigen::VectorXd result = (mtxVandermonde.transpose() * mtxVandermonde).inverse() * (mtxVandermonde.transpose()) * sampleY;

	return result;
}

int check_para(int argc, char* argv[], const char* par_name, int par_numb ,int *param_position)
{
	*param_position = -1;
	for (int i = 0; i < argc; i++) 
	{
        if (!strcmp(argv[i], par_name)) 
		{
            *param_position = i; // Update the position tracker
            if (i >= (argc - par_numb)) 
                return 0; //Too few arguments, please enter the correct number of parameters!
			else 
                return 1; //Includes '-par'&'-parameter'
        }
    }
    return 2; // '-par' not found
}

int check_para(int argc, char* argv[], const char* par_name_1, const char* par_name_2, int par_numb ,int *param_position)
{
	*param_position = -1;
	for (int i = 0; i < argc; i++) 
	{
        if (!strcmp(argv[i], par_name_1) || !strcmp(argv[i], par_name_2))
		{
            *param_position = i; // Update the position tracker
            if (i >= (argc - par_numb)) 
                return 0; //Too few arguments, please enter the correct number of parameters!
			else 
                return 1; //Includes '-par'
        }
    }
    return 2; // '-par' not found
}

int check_mulpara(int argc, char* argv[], const vector<const char*> par_names, int par_numb, int *param_position)
{
    for (int i = 0; i < argc; i++) 
	{
        for (const auto& par_name : par_names) 
		{
            if (!strcmp(argv[i], par_name)) 
			{
                *param_position = i; // Update the position tracker
                if (i >= (argc - par_numb))
                    return 0; // Too few arguments, please enter the correct number of parameters!
                else
                    return 1; // Includes '-par'&'-parameter'
            }
        }
    } 
    return 2; // No '-par' found
}

int check_number(int argc, char* argv[], int *param_position)
{
	*param_position = -1;
	for (int i = 0; i < argc; i++) 
	{      
		bool is_number = true;
		string str = string(argv[i]);
		size_t start = (str[0] == '-' || str[0] == '+') ? 1 : 0;
        for (size_t j = start; j < str.size(); j++)
        {
            if (!isdigit(str[j]))
            {
                is_number = false;
                break;
            }
        }
        if (is_number)
        {
            *param_position = i; 
        }
    }
    if(*param_position == -1)
		return 2; // none number has been found
    else
		return 1; // number has been found
}

//Convert the contents of the file 'table' into a 'map<string, double>'
std::map<std::string, double> createMap(const char table[]) 
{
	std::map<string, double> myMap;
	FILE* fpt = fopen(table, "r");
	if (fpt == nullptr)
	{
		printf("%s IS NOT EXIST!\n", table);
		return myMap;
	}
	else
	{
		char buf[1024];
		while (fgets(buf, 1024, fpt))
		{
			char key[10];
			double doub;
			if (sscanf(buf, "%s%lf", key, &doub) != 2)
			{
				remove(std::begin(buf), std::end(buf), '\n');
				printf("[%s] in [%s] file may be wrong format!\n", buf, table);
				continue;
			}
			if (!myMap.count(string(key)))
			{
				myMap[string(key)] = doub;
			}
			else
			{
				printf("VASPMATE has detected that you have entered two identical elem in the table.\n");
				continue;
			}
		}
	}
	fclose(fpt);
	return myMap;
}

void runVASPMATE(const string& command, const string& logFile) 
{
    string fullCommand = command + " >> " + logFile;
    system(fullCommand.c_str());
}

void Output_to_terminal(const char file[])
{
	FILE* fp = fopen(file, "r");
	if (fp == NULL)
	{
		printf("%s IS NOT EXIST!\n", file);
		return;
	}
	char buf[1024];
	while (fgets(buf, 1024, fp) != NULL)
		cout << buf;
	fclose(fp);
	cout << endl;
}

bool VM_fileexist(const char file[])
{
    ifstream file_(file, ios::binary | ios::ate);
    if (!file_.is_open())
	{
        return false;
    }
	streamsize size = file_.tellg();
    file_.close();
    return size > 0;
}

bool VM_isnumb(const char numb[])
{
    if (!numb || *numb == '\0') 
		return false;
    char* end = nullptr;
    errno = 0;
    double val = std::strtod(numb, &end);
    if (errno == 0 && end != numb && *end == '\0')
        return true;
    return false;
}