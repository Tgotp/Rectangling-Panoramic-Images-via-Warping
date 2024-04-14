#include <opencv2/core/core.hpp> 
#include <opencv2/highgui/highgui.hpp> 
#include <opencv2/opencv.hpp>
using namespace std;
using namespace cv;

vector<vector<pair<int,int> > > Mesh_Placement(Mat img,Mat image,vector<vector<pair<int,int> > > U,bool is_show=1);