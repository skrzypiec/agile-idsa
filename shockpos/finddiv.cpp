#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <cstdlib>

using namespace std;


double div(double r1, double r2, double r3, double v1, double v2, double v3)
{
	//return abs(((v3-v1)/(r3-r1)));
	return abs(1/(r2*r2) * (2*r2*v2 + r2*r2*(v3-v1)/(r3-r1)));
}

int main(int argc, char** argv) {
  
  vector<vector<double> > allData;
  
  ifstream fin(argv[1]);
  string line;
  while (getline(fin, line)) {      // for each line
    vector<double> lineData;           // create a new row
    double val;
    istringstream lineStream(line); 
    while (lineStream >> val) {          // for each value in line
      lineData.push_back(val);           // add to the current row
    }
    allData.push_back(lineData);         // add row to allData
  }

///////////find maximal value of v

  double vel = allData[0][1];
  double rad = allData[0][0];

  for(int i = 0; i < allData.size(); i++)
  {
	 if(abs(allData[i][1]) > vel)
		{
 	    	vel = abs(allData[i][1]);
  			rad = allData[i][0];
  		}
  }
//  cout << min << "\n" << rad << endl;

//////////////////////////////////////

double max1 = 0;
double rad1 = 0;
double vel1 = 0;

/////////find maximal divergence value of v
for(int i = 2; i < allData.size()-1; i++)
  {
	 if(div(allData[i-1][0], allData[i][0], allData[i+1][0], allData[i-1][1], allData[i][1], allData[i+1][1]) > max1)
		{
			max1 = div(allData[i-1][0], allData[i][0], allData[i+1][0], allData[i-1][1], allData[i][1], allData[i+1][1]);
 	    	vel1 = allData[i][1];
  			rad1 = allData[i][0];
  		}
  }
//////////


    string linia;
    fstream plik;
 
    plik.open(argv[1], ios::out | ios::trunc);
    if(plik.good() == true)
    {
		plik << rad << "\t" << vel << "\t" << rad1 << "\t" << vel1 << endl;
        plik.close();
    }


return 0;
}
