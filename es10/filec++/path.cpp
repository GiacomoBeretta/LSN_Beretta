#include "path.hpp"

using namespace std;

Path ::Path()
{
  n_cities_ = 0;
  vector<Position> path;
  vector<int> indexes;
  path_ = path;
  indexes_ = indexes;
  norm_power_ = 1;
  sphere_ = false;
  path_length_ = 0;
  temp_ = 1;
  accepted_ = 0;
  attempted_ = 0;
}

Path ::Path(const int n_cities)
{
  n_cities_ = n_cities;
  vector<Position> path(n_cities_ + 1);
  vector<int> indexes(n_cities_ + 1);
  for (unsigned int i = 0; i < n_cities_ + 1; i++)
  {
    indexes[i] = i;
    path[i].setX(0);
    path[i].setY(0);
  }
  indexes[n_cities_] = 0;
  path_ = path;
  norm_power_ = 1;
  sphere_ = false;
  path_length_ = 0;
  indexes_ = indexes;
  temp_ = 1;
  accepted_ = 0;
  attempted_ = 0;
}

Path::Path(const int n_cities, const int norm_power) : Path(n_cities)
{
  norm_power_ = norm_power;
}

Path::Path(const int n_cities, const int norm_power, const double temp) : Path(n_cities, norm_power)
{
  temp_ = temp;
}

Path ::Path(const vector<vector<double>> &v, const int norm_power, const double temp)
{
  n_cities_ = v.size();
  vector<Position> path(n_cities_ + 1);
  vector<int> indexes(n_cities_ + 1);
  for (unsigned int i = 0; i < n_cities_; i++)
  {
    Position pos(v[i][0], v[i][1]);
    path[i] = pos;
    indexes[i] = i;
  }
  path[n_cities_] = path[0];
  indexes[n_cities_] = 0;
  path_ = path;
  indexes_ = indexes;
  norm_power_ = norm_power;
  sphere_ = false;
  computePathLength();
  temp_ = temp;
  accepted_ = 0;
  attempted_ = 0;
}

Path ::Path(const std::vector<Position> &path, const int norm_power, const double temp)
{
  n_cities_ = path.size();
  path_ = path;
  path_.push_back(path[0]);
  norm_power_ = norm_power;
  sphere_ = false;
  computePathLength();
  vector<int> indexes(n_cities_ + 1);
  for (unsigned int i = 0; i < n_cities_; i++)
  {
    indexes[i] = i;
  }
  indexes[n_cities_] = 0;
  indexes_ = indexes;
  temp_ = temp;
  accepted_ = 0;
  attempted_ = 0;
}

Path ::Path(const Path &path)
{
  n_cities_ = path.getNcities();
  path_ = path.getPath();
  norm_power_ = path.getNormPower();
  sphere_ = path.getSphere();
  path_length_ = path.getPathLength();
  indexes_ = path.getIndexes();
  temp_ = path.getTemp();
  accepted_ = path.getAccepted();
  attempted_ = path.getAttempted();
}

void Path::setXcity(const int i, const double x) { path_[i].setX(x); }

void Path::setYcity(const int i, const double y) { path_[i].setX(y); }

void Path ::setCity(const int i, const Position &p) { path_[i] = p; }

void Path ::setPath(const std::vector<Position> &path) { path_ = path; }

void Path::setNormPower(const int norm_power) { norm_power_ = norm_power; }

void Path::setSphere(const bool sphere) { sphere_ = sphere; }

void Path::setPathLength(const double L) { path_length_ = L; }

void Path::setIndex(const int i, const int index) { indexes_[i] = index; }

void Path ::setIndexes(const vector<int> &indexes) { indexes_ = indexes; }

void Path::setTemp(const double temp) { temp_ = temp; }

void Path::setAccepted(const int accepted) { accepted_ = accepted; }

void Path::setAttempted(const int attempted) { attempted_ = attempted; }

void Path ::randomOrder(Random &rand)
{
  vector<int> indexes(n_cities_ + 1);
  for (unsigned int i = 0; i < n_cities_ + 1; i++)
  {
    indexes[i] = 0;
  }
  unsigned int j = 1;
  int r;
  while (j < n_cities_) // n_cities_=34, j=1,...,33
  {
    r = rand.Integer(1, n_cities_ - 1); // r=1,...,33
    if (indexes[r] == 0)
    {
      indexes[r] = j;
      j++;
    }
  }
  indexes_ = indexes;
  vector<Position> path_copy = path_;
  for (unsigned int i = 1; i < n_cities_; i++)
  {
    path_[i] = path_copy[indexes[i]];
  }
  computePathLength();
}

void Path ::setCitiesRandomlyOnCircumference(Random &rand)
{
  double r;
  for (unsigned int i = 0; i < n_cities_; i++)
  {
    r = rand.Rannyu(0, 2 * pi);
    path_[i].setX(cos(r)); // x on the unit circle
    path_[i].setY(sin(r)); // y on the unit circle
  }
  path_[n_cities_] = path_[0];
  computePathLength();
}

void Path ::setCitiesRandomlyInSquare(Random &rand)
{
  for (unsigned int i = 0; i < n_cities_; i++)
  {
    path_[i].setX(rand.Rannyu()); // x,y in the square (0,1)  (1,1)
    path_[i].setY(rand.Rannyu()); //                  (0,0)  (1,0)
  }
  path_[n_cities_] = path_[0];
  computePathLength();
}

void Path::setAmericanCapitals()
{
  sphere_ = true;
  ifstream input("American_capitals.dat");
  string line;
  getline(input, line); // first line doesn't contain data
  string word[3];
  double phi, theta;
  for (unsigned int i = 0; i < n_cities_; i++)
  {
    getline(input, line);
    istringstream iss(line);
    while (iss >> word[2])
    {
      word[0] = word[1];
      word[1] = word[2];
    }
    phi = stod(word[0]);   // -180°<phi<180° longitude
    theta = stod(word[1]); // latitude

    phi += 180.0;         // 0°<phi<360°
    theta = 90.0 - theta; // so that theta is 0° at the north pole and 90° at the equator

    phi = phi * pi / 180.0; // phi and theta in radians
    theta = theta * pi / 180.0;
    path_[i].setX(phi);
    path_[i].setY(theta);
  }
  input.close();
  path_[n_cities_] = path_[0];
  computePathLength();
}

void Path ::swapCities(const int a, const int b)
{
  Position temp_pos = path_[a];
  int temp_index = indexes_[a];
  path_[a] = path_[b];
  indexes_[a] = indexes_[b];
  path_[b] = temp_pos;
  indexes_[b] = temp_index;
}

void Path ::permutateCities(const int a, const unsigned int n, const int b)
{
  int c, d;
  for (unsigned int i = 0; i < n; i++)
  {
    c = boundaryCondition(a + i, n_cities_);
    d = boundaryCondition(b + i, n_cities_);
    swapCities(c, d);
  }
}

void Path ::reverseCities(const int a, const unsigned int n)
{
  int b, c;
  for (unsigned int i = 0; i < n / 2; i++)
  {
    b = boundaryCondition(a + i, n_cities_);
    c = boundaryCondition(a + n - 1 - i, n_cities_);
    swapCities(b, c);
  }
}

void Path ::checkPath() const
{
  // check by index
  if (n_cities_ != indexes_.size() - 1)
  {
    cout << "error: the indexes vector hasn't the right dimension" << endl;
    exit(-1);
  }

  if (indexes_[0] != indexes_[n_cities_])
  {
    cout << "error: the first city's index doesn't match the last city's"
         << endl;
    exit(-1);
  }

  for (unsigned int i = 0; i < n_cities_ - 1; i++)
  {
    for (unsigned int j = i + 1; j < n_cities_; j++)
    {
      if (indexes_[i] == indexes_[j])
      {
        cout << "error: there are two cities with same index, the " << i
             << " city and the " << j << " city" << endl;
        exit(-1);
      }
    }
  }
  // check by position
  if (path_[0] != path_[n_cities_])
  {
    cout << "error: the first city doesn't match the last city" << endl;
    exit(-1);
  }
  for (unsigned int i = 0; i < n_cities_ - 1; i++)
  {
    for (unsigned int j = i + 1; j < n_cities_; j++)
    {
      if (path_[i] == path_[j])
      {
        cout << "error: there are two cities with same position, the " << i
             << " city and the " << j << " city" << endl;
        exit(-1);
      }
    }
  }
}

void Path ::computePathLength()
{
  double d;
  double sum = 0;
  for (unsigned int i = 0; i < n_cities_; i++)
  {
    if (sphere_ == false)
    {
      d = path_[i].distance(path_[i + 1]);
    }
    else
    {
      d = path_[i].distanceOnSphere(path_[i + 1]);
    }
    sum += pow(d, norm_power_);
  }
  path_length_ = pow(sum, 1 / (double)norm_power_);
}

double Path::probability()
{
  computePathLength();
  double length = getPathLength();
  return exp(-length / temp_);
}

void Path::metropolis(Random &rand)
{
  Path new_path;
  int a, b, d, n;
  double p_old, p_new;
  for (int i = 0; i < 3; i++)
  {
    new_path = *this;
    switch (i)
    {
    case 0:
      a = rand.Integer(1, n_cities_ - 1);
      b = rand.Integer(1, n_cities_ - 1);
      new_path.swapCities(a, b);
      break;
    case 1:
      a = rand.Integer(1, n_cities_ - 1);
      b = rand.Integer(1, n_cities_ - 1);
      d = boundaryDistance(a, b, n_cities_);
      if (d != 0)
      {
        n = rand.Integer(1, d);
        new_path.permutateCities(a, n, b);
      }
      break;
    case 2:
      a = rand.Integer(1, n_cities_ - 1);
      n = rand.Integer(1, n_cities_ - 1);
      new_path.reverseCities(a, n);
      break;
    }
    p_old = probability();
    p_new = new_path.probability();
    if (p_old < p_new)
    {
      *this = new_path;
      accepted_++;
    }
    else if (rand.Rannyu() < p_new / p_old)
    {
      *this = new_path;
      accepted_++;
    }
    attempted_++;
  }
}

Path &Path ::operator=(const Path &path)
{
  if (this != &path)
  {
    n_cities_ = path.getNcities();
    path_ = path.getPath();
    norm_power_ = path.getNormPower();
    sphere_ = path.getSphere();
    path_length_ = path.getPathLength();
    indexes_ = path.getIndexes();
    temp_ = path.getTemp();
    accepted_ = path.getAccepted();
    attempted_ = path.getAttempted();
  }
  return *this;
}

Position &Path ::operator[](const int i) { return path_[i]; }

const Position &Path ::operator[](const int i) const { return path_[i]; }

void swap(Path &a, Path &b)
{
  Path temp = a;
  a = b;
  b = temp;
}

void Path ::printPath() const
{
  for (unsigned int i = 0; i < path_.size(); i++)
  {
    cout << "i = " << i << ", x = " << getXcity(i) << ", y = " << getYcity(i)
         << endl;
  }
}

void Path ::printPathFile(ofstream &out) const
{
  for (unsigned int i = 0; i < path_.size(); i++)
  {
    out << i << " " << getXcity(i) << " " << getYcity(i) << endl;
  }
}

void Path ::printIndexes() const
{
  for (unsigned int i = 0; i < path_.size(); i++)
  {
    cout << "i = " << i << ", index = " << indexes_[i] << endl;
  }
}

/*
void Path :: shiftCities(int a, int n, int shift)
{
        int n_cities = getNcities();
        int c, d;

        for(int i=0; i<n; i++)
        {
                if( i+a < n_cities )
                {
                        if( i+a+shift < n_cities )
                        {
                                swapCities(i+a, i+a+shift);
                        }else
                        {
                                c = (i+a+shift) % n_cities;
                                swapCities(i+a, i+a+shift-n_cities*c+1);
                        }
                }else
                {
                        c = (i+a+shift) % n_cities;
                        d = (i+a)%n_cities;
                        swapCities(i+a-n_cities*d+1, i+a+shift-n_cities*c+1);
                }
        }
}*/

int boundaryCondition(const int a, const int n_cities)
{
  int b = a % n_cities;
  if (b == 0)
    return 1;
  else
    return b;
}

int boundaryDistance(const int a, const int b, const int n_cities)
{
  int d1 = abs(a - b);
  int d2 = 0;
  if (a < b)
  {
    d2 = n_cities - b + a - 2; //(n_cities - 1) - b + a - 1;
  }
  else
  {
    d2 = n_cities - a + b - 2;
  }
  return min(d1, d2);
}
