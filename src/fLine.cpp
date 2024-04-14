#include "fLine.h"
double* fline::line_img(Mat img,int *n_line,bool is_show)
{
    int n = img.rows,m = img.cols;
    double *input_img = solve_img(img);
    double *out = lsd(n_line,input_img,m,n);
    if(is_show)
    {
        Mat line_img = Mat::zeros(n,m,CV_8UC1);
        line_img = ~line_img;
        int num = * n_line;
        for(int i = 0;i < num;++ i)
        {
            // cout << out[i*7]<< ' ' << out[i*7+1] << ' ' << out[i*7+2] << ' ' << out[i*7+3] << endl;
            line(line_img,Point(out[i*7],out[i*7+1]),Point(out[i*7+2],out[i*7+3]),0,1);
        }
        namedWindow("line",WINDOW_FREERATIO);
        imshow("line",line_img);
        waitKey(0);
    }
    return out;
}

bool onSegment(Point p, Point q, Point r) {
    if (q.x <= max(p.x, r.x) && q.x >= min(p.x, r.x) &&
        q.y <= max(p.y, r.y) && q.y >= min(p.y, r.y))
        return true;
    return false;
}

int orientation(Point p, Point q, Point r) {
    double val = (q.y - p.y) * (r.x - q.x) - (q.x - p.x) * (r.y - q.y);
    if (val == 0) return 0; 
    return (val > 0) ? 1 : 2; 
}

bool doIntersect(Line l1, Line l2) {
    int o1 = orientation(l1.x, l1.y, l2.x);
    int o2 = orientation(l1.x, l1.y, l2.y);
    int o3 = orientation(l2.x, l2.y, l1.x);
    int o4 = orientation(l2.x, l2.y, l1.y);

    if (o1 != o2 && o3 != o4)
        return true;

    if (o1 == 0 && onSegment(l1.x, l2.x, l1.y)) return true;
    if (o2 == 0 && onSegment(l1.x, l2.y, l1.y)) return true;
    if (o3 == 0 && onSegment(l2.x, l1.x, l2.y)) return true;
    if (o4 == 0 && onSegment(l2.x, l1.y, l2.y)) return true;

    return false;
}
Point lineIntersection(Line l1, Line l2) {
    // a1x + b1y = c1
    double a1 = l1.y.y - l1.x.y;
    double b1 = l1.x.x - l1.y.x;
    double c1 = a1*(l1.x.x) + b1*(l1.x.y);

    // a2x + b2y = c2
    double a2 = l2.y.y - l2.x.y;
    double b2 = l2.x.x - l2.y.x;
    double c2 = a2*(l2.x.x) + b2*(l2.x.y);

    double determinant = a1*b2 - a2*b1;

    int x = (b2*c1 - b1*c2) / determinant;
    int y = (a1*c2 - a2*c1) / determinant;
    return Point(x, y);
}

Mat lineImg;
Line getLine(Line x,vector<Point> p)
{
    for(int i = 0;i <= 3;++ i)
    {
        Line y = Line(p[i],p[(i+1)%4]);
        
        if(doIntersect(x, y))
        {
            // cout << "i:" << i << endl;
            // cout << "x: "<< x.x.x << " " << x.x.y << " " << x.y.x << " "<< x.y.y << endl;
            // cout << "y: "<< y.x.x << " " << y.x.y << " " << y.y.x << " "<< y.y.y << endl;
            // cout << "get line intersection: " << Line(lineIntersection(x, y),x.x).x.x << " " << Line(lineIntersection(x, y),x.x).x.y << " "<< 
            //     Line(lineIntersection(x, y),x.x).y.x << " " << Line(lineIntersection(x, y),x.x).y.y << endl; 
            
            // line(lineImg,Point(x.x.y,x.x.x),Point(x.y.y,x.y.x),Scalar(0,0,255),1);
            // line(lineImg,Point(y.x.y,y.x.x),Point(y.y.y,y.y.x),Scalar(0,255,0),1);
            // line(lineImg,Point(lineIntersection(x, y).y,lineIntersection(x, y).x),Point(x.x.y,x.x.x),Scalar(255,0,0),1);
            // namedWindow("lineImg",WINDOW_FREERATIO);
            // imshow("lineImg",lineImg);
            // waitKey(0);
            Point L = lineIntersection(x, y);
            if(L.x == INT_MIN && L.y == INT_MIN)
                return x;
            return Line(L,x.x);
        }
    }
    cout << "error exit" << endl;
    exit(0);
}

vector<vector<vector<Line > > > fline::init_line(vector<vector<Point > > mesh,double *line,int n)
{
    // cout << "init line" << endl;
    // lineImg = Mat::zeros(291,640,CV_8UC3);
    vector<vector<vector<Line > > > Mesh_Line;
    vector<vector<Line > > Mesh_line;
    vector<Line > mesh_line;
    vector<Point > contours;
    for(int i = 0;i < 20; ++ i)
    {
        Mesh_line.clear();
        for(int j = 0;j < 20;++ j)
        {
            mesh_line.clear();
            Point lt = Point(mesh[i][j].x,mesh[i][j].y);
            Point rt = Point(mesh[i][j+1].x,mesh[i][j+1].y);
            Point lb = Point(mesh[i+1][j].x,mesh[i+1][j].y);
            Point rb = Point(mesh[i+1][j+1].x,mesh[i+1][j+1].y);
            contours.clear();
            contours.push_back(lt);
            contours.push_back(rt);
            contours.push_back(rb);
            contours.push_back(lb);
            for(int k = 0;k < n;++ k)
            {
                Point st = Point(line[k*7+1],line[k*7]);
                Point ed = Point(line[k*7+3],line[k*7+2]);
                int z1 = pointPolygonTest(contours,st,false);
                int z2 = pointPolygonTest(contours,ed,false);

                if(z1 >= 0 && z2 >= 0)
                    mesh_line.push_back(Line(st,ed,line[k*7+5]));
                else if(z1 >= 0)
                    mesh_line.push_back(getLine(Line(st,ed,line[k*7+5]),contours));
                else if(z2 >= 0)
                    mesh_line.push_back(getLine(Line(ed,st,line[k*7+5]),contours));
                
            }
            Mesh_line.push_back(mesh_line);
        }
        Mesh_Line.push_back(Mesh_line);
    }
    // cout << "init line end" << endl;
    return Mesh_Line;
}

void fline::check_line(Mat img,vector<vector<vector<Line > > >mesh_line)
{
    cout << "check line" << endl;
    Mat fimg;
    img.copyTo(fimg);
    for(int i = 0;i < 20;++ i)
        for(int j = 0;j < 20;++ j)
            for(auto k : mesh_line[i][j])
                // line(fimg,k.x,k.y,Scalar(0,255,0),1);
                line(fimg,Point(k.x.y,k.x.x),Point(k.y.y,k.y.x),Scalar(0,255,0),1);
        
    namedWindow("mesh line",WINDOW_FREERATIO);
    imshow("mesh line",fimg);
    waitKey(0);
}

double fline::get_gray(Mat&img,int i,int j) { return img.at<Vec3b>(i,j)[0] * 0.299 + img.at<Vec3b>(i,j)[1] * 0.587 + img.at<Vec3b>(i,j)[2] * 0.114; }

double* fline::solve_img(Mat img)
{
    int n = img.rows,m = img.cols;
    double*ans = new double[n * m];
    for(int j = 0;j < m; ++ j)
        for(int i = 0;i < n; ++ i)
        {
            // cout << n << ' ' << m <<' '<< i << ' ' << j <<' ' <<i +j*n <<endl;
            ans[i*m+j] = get_gray(img,i,j);
        }
    return ans;
}


double* fline::Line_rotate_count(int *num,vector<vector<vector<Line > > > mesh_line,int number)
{
    double *bins;
    bins = new double[number];
    for(int i = 0;i < 20;++ i)
        for(int j = 0;j < 20;++ j)
            for(auto k:mesh_line[i][j])
            {
                int pos = ceil(k.theta * 50);
                bins[pos] += k.theta * PI;
                num[pos] ++;
            }
    for(int i = 0;i < number; ++ i)
        bins[i] /= num[i];
    return bins;
}